# Parallel script for regridding the forcing data to desired inpu
# OVERVIEW:
import sys
from mpi4py import MPI
import netCDF4 as nc
import ast
import numpy as np
import os
import rasterio
import subprocess
import datetime

# MPI4PY Stuff
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# -------------- #
# LOAD ARGUMENTS #
# -------------- #
# load in the arguments from wrf_preprocessing.py 
argv = sys.argv[1:]
w_input = argv[0]. # location of the PCF files; I use '/stor/soteria/hydro/shared/data/PCF/1hr/daily/'
w_dir = argv[1] # a workspace directory; blank directory you don't care about i.e. ../workspace/
forcing_dir = argv[2] # Directory for the output forcing files; I use FORCING/
run_dir = argv[3] # Run Directory for WRF (in my setup, this is often ../wrf_hydro_run/)
geogrid = run_dir+argv[4] # Geogrid file i.e. /path/to/the/file/geo_em.d01.nc
gdal_cmd = argv[5] # see comment below for details
filelist =ast.literal_eval(argv[6]) # see below for details

### OVERALL ###
# the run should look like this:
# mpiexec -n N python wrf_regrid.py w_input w_dir forcing_dir run_dir geogrid gdal_cmd filelist
# pay close attention to use of '' and "" below, as having these wrong can mess up the run

### GDAL_CMD ###
# gdal_cmd looks like this: 
# gdal_cmd = '\"gdalwarp -s_srs '+src_proj+' -t_srs \''+proj+'\' -te '+ex_string+' -tr '+str(dx)+' '+str(dy)+' -r bilinear \"'

# src_proj doesn't need to change from my code and looks like this: src_proj = "'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'"

# proj will be dependent on your specific projection used to run wrf; for me this is lcc and looks like this:
# proj = "+proj=lcc +a=6370000 +b=6370000 +lon_0="+str(stand_lon)+" +lat_1="+str(truelat1)+" +lat_2="+str(truelat2)
# you will need to determine stand_lon and truelate1 and truelat2 for your setup if using lcc; otherwise you
# will need to figure out the proj4 statement for how you are running the code. The infromation for the proj4
# statement should be in namelist.wps; I use the same naming conventions (truelat1,stand_lon etc.) as the namelist.wps file

# ex_string is coded inefficiently in my script, but is... 
# extent = [0,0,0,0]
# extent[0]=np.round(np.mean(x_val))-dx*x/2 # get left extent for x
# extent[2]=np.round(np.mean(x_val))+dx*x/2 # get right extent for x
# extent[1]=np.round(np.mean(y_val))-dy*y/2 # get lower extent for y
# extent[3]=np.round(np.mean(y_val))+dy*y/2 # get upper extent for y
# ex_string = str(extent[0])+' '+str(extent[1])+' '+str(extent[2])+' '+str(extent[3])

# dx and dy are probably in meters depends on the projection

### FILELIST ###
# this is a list of the files that you will generate forcing data for. 
# each file in the list is formatted as yyyymmdd.nc (i.e. 20170612.nc)
# Here is the code I use to generate filelist:

'''
filelist = []
dt = timedelta(1)
num_days=(end_dt-start_dt)/dt #start and end datetime for simulation
num_days = round(num_days+1)
n_files= num_days*24
geogrid = dom_dir+'geo_em.d01.nc'

for i in range(num_days):
	date = start_dt+dt*i
	if date.month > 9:
		month = date.month
	else:
		month = '0'+str(date.month)
	if date.day > 9:
		day = date.day
	else:
		day = '0'+str(date.day)
	filelist.append(str(date.year)+str(month)+str(day)+'.nc')
'''
# and this filelist gets passed to wrf_regrid.py as ... str(filelist) ...



# --------------------- #
# INITIAL SETUP STUFFS  #
# --------------------- #

# pull info from geogrid
fp_geo = nc.Dataset(geogrid,'r')
corner_lats = fp_geo.corner_lats[0:4]
corner_lons = fp_geo.corner_lons[0:4]
dx = fp_geo.DX
dy = fp_geo.DY

# Set up approximate extent map
fp_1=nc.Dataset(w_input+filelist[0],'r')
lats = fp_1['lat'][:]
lons = fp_1['lon'][:]
dlat = lats[1]-lats[0]
dlon = lons[1]-lons[0]
min_lon = np.min(corner_lons)-2*dlon
min_lat = np.min(corner_lats)-2*dlat
max_lon = np.max(corner_lons)+2*dlon
max_lat = np.max(corner_lats)+2*dlat
lat_m = (lats>=min_lat) & (lats<=max_lat) # mask for lat extent
lon_m = (lons>min_lon) & (lons<=max_lon)  # mask for lon extent
trim_lats = lats[lat_m]
trim_lons = lons[lon_m]
fill_value=fp_1['swdown']._FillValue
fp_1.close()


# --------------------------------- #
# REGRID THE FILES NETCDF to NETCDF #
# --------------------------------- #
# NOTE: this is not an optimized method to accomplish the regridding
# figuring out how to regrid with CDO or with GDAL directly from netCDF
# would be optimal, but unfortunately I had issues with this and getting
# it to work is a higher priority than optimization

# Create subdirectory for this parallel filegroup
out_dir = 'out'+str(rank)+'/'
subprocess.run('mkdir '+run_dir+w_dir+out_dir,shell=True)

# Some variable lists
var_names={'T2D':'tair','Q2D':'spfh','U2D':'wind','V2D':'wind','PSFC':'psurf',
           'RAINRATE':'precip','SWDOWN':'swdown','LWDOWN':'lwdown'}
var_max={'T2D':335,'Q2D':.05,'U2D':75,'V2D':75,'PSFC':110000,
           'RAINRATE':10,'SWDOWN':2000,'LWDOWN':600}

# Delcare proj4 and Create transform (west, north, xsize, ysize)
proj4_src = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

print('Rank '+str(rank)+' regridding...',flush=True)
for file in filelist[rank::size]:
	for var in var_names.keys():
		src_file =w_input+file
		src_name = "NETCDF:\""+src_file+"\":"+var_names[var]
		new_name = run_dir+w_dir+out_dir+file[0:8]+"_"+var+".nc"
		regrid_cmd = gdal_cmd+' '+src_name+' '+new_name+' >/dev/null'
		subprocess.run(regrid_cmd,shell=True)
	data = nc.Dataset(new_name,'r')['Band20'][:]
	for t in range(24):
		name=file[0:8]+format(t,'02d')+'.LDASIN_DOMAIN1'
		fp_out=nc.Dataset(run_dir+forcing_dir+name,'w')
		fp_out.createDimension('south_north',data.shape[0])
		fp_out.createDimension('west_east',data.shape[1])
		fp_out.createDimension('Time',1)
		fp_out.createDimension('DateStrLen',20)

		# read in data
		for var in var_names.keys():
			fp_out.createVariable(var,'f4',('Time','south_north','west_east'),\
                                              fill_value=fill_value)
			var_file =  run_dir+w_dir+out_dir+file[0:8]+"_"+var+".nc"
			data = nc.Dataset(var_file,'r')['Band'+str(t+1)][:]
			msk=(data>=var_max[var])
			smmsk=np.sum(msk)
			data[msk]=np.median(data)
			if(smmsk>0):
				print(str(smmsk)+' BAD VALUES FILLED FOR '+str(var)\
                    +' ON '+file[0:8]+format(t,'02d'))
			fp_out[var][0,:,:]=data[:]
			fp_out[var].missing_value=fill_value
		# load in other data
		fp_out.createVariable('Times','c',('DateStrLen'))
		time_str=file[0:4]+'-'+file[4:6]+'-'+file[6:8]+'_'+format(t,'02d')+':00:00 '
		fp_out['Times'][:]=time_str[:]
		fp_out.createVariable('valid_time','i4',('Time'))
		dt = datetime.datetime(int(file[0:4]),int(file[4:6]),int(file[6:8]),t)
		fp_out['valid_time'][0]=dt.replace(tzinfo=datetime.timezone.utc).timestamp()
		fp_out.createVariable('lat','f4',('south_north','west_east'))
		fp_out.createVariable('lon','f4',('south_north','west_east'))
		fp_out['lat'].FieldType=104
		fp_out['lat'].units="degrees latitude"
		fp_out['lat'].sr_x=1
		fp_out['lat'].sr_y=1
		fp_out['lat'].stagger='M'
		fp_out['lat'].description = 'Latitude on mass grid'
		fp_out['lon'].FieldType=104
		fp_out['lon'].units="degrees longitude"
		fp_out['lon'].sr_x=1
		fp_out['lon'].sr_y=1
		fp_out['lon'].stagger='M'
		fp_out['lon'].description = 'Longitude on mass grid'

		# load in lat and long
		fp_out['lat'][:]=fp_geo['XLAT_M'][0,:,:]
		fp_out['lon'][:]=fp_geo['XLONG_M'][0,:,:]
		fp_out.close()
	subprocess.run('rm '+run_dir+w_dir+out_dir+file[0:8]+"_*",shell=True)

			
	
