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
argv = sys.argv[1:]
w_input = argv[0]
w_dir = argv[1]
forcing_dir = argv[2]
run_dir = argv[3]
geogrid = run_dir+argv[4]
gdal_cmd = argv[5]
filelist =ast.literal_eval(argv[6])

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

# Delcare proj4 and Create transform (west, north, xsize, ysize)
proj4_src = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'


print('Rank '+str(rank)+' regridding...',flush=True)
for file in filelist[rank::size]:
	for var in var_names.keys():
		src_file =w_input+file
		src_name = "NETCDF:\""+src_file+"\":"+var_names[var]
		new_name = run_dir+w_dir+out_dir+file[0:8]+"_"+var+".nc"
		regrid_cmd = gdal_cmd+'-srcnodata '+str(fill_value)+\
                     ' '+src_name+' '+new_name+' >/dev/null'
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
			fp_out[var][0,:,:]=nc.Dataset(var_file,'r')['Band'+str(t+1)][:]
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

			
	
