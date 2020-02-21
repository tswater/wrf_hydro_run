# A regridding script for the weather data in the cluster
import netCDF4 as nc
import numpy as np
import subprocess
import rasterio
import datetime
from datetime import timedelta
import os
import sys

start_date = [2011,8,26,0] # format: [year,month,day,hour]
end_date = [2011,9,2,0] # format: [year,month,day,hour]
geogrid = '../../DOMAIN/geo_em.d01.nc'
w_input = '/stor/soteria/hydro/shared/data/PCF/1hr/daily/' # weather input location
clean = True # if True, will remove intermediary nc files and only keep needed output
# -------------------------------------------------------------------




# -------------------------------------------------------------------------- #
# -------------------------- HELPER FUNCTION ------------------------------- #
# A helper function to display progress in the evaluation of the files       #
# -------------------------------------------------------------------------- #
def count_help(c,n):
	if np.floor(c/n*10)>np.floor((c-1)/n*10):
		print(np.round(c/n*100),end="",flush=True)
	elif np.floor(c/n*50)>np.floor((c-1)/n*50):
		print('.',end='',flush=True)

# -------------------------------------------------------------------





# -------------------------------------------------------------------------- #
# ------------------------------ SETUP  ------------------------------------ #
# Set everything up and read in arguments from main script                   #
# -------------------------------------------------------------------------- #
print('Initial Regrid Setup...',end='',flush=True)

# read in arguments
argv = sys.argv[1:]
start_date = argv[0]
end_date = argv[1]
w_input = argv[2]
proj = argv[3] # proj statement
w_dir = argv[4]
dom_dir = argv[5]
forcing_dir = argv[6]
log_dir = argv[7]
run_dir = argv[8]

# Create the working directory
subprocess.run('mkdir '+run_dir+w_dir+'out',shell=True)

# declare some names
var_names = {'tair':'TMP_110_HTGL','spfh':'SPF_H_110_HTGL',
    'wind':['U_GRD_110_HTGL','V_GRD_110_HTGL'],'psurf':'PRES_110_SFC',
    'precip':'A_PCP_110_SFC_acc1h','swdown':'DSWRF_110_SFC','lwdown':'DLWRF_110_SFC'}

# set up datetimes
start_dt = datetime.datetime.strptime(start_date,'%Y-%m-%d_%X')
end_dt = datetime.datetime.strptime(end_date,'%Y-%m-%d_%X')

# load corner lats and lons
fp_geo = nc.Dataset(run_dir+dom_dir+'geo_em.d01.nc','r')
corner_lats = fp_geo.corner_lats[0:4]
corner_lons = fp_geo.corner_lons[0:4]

# Set up approximate extent map
fp_1=nc.Dataset(w_input+start_date[0:4]+start_date[5:7]+start_date[8:10]+'.nc','r')
lats = fp_1['lat'][:]
lons = fp_1['lon'][:]
min_lon = np.floor(np.min(corner_lons))
min_lat = np.floor(np.min(corner_lats))
max_lon = np.ceil(np.max(corner_lons))
max_lat = np.ceil(np.max(corner_lats))
lat_m = (lats>=min_lat) & (lats<=max_lat) # mask for lat extent
lon_m = (lons>min_lon) & (lons<=max_lon)  # mask for lon extent

# Create a filelist using datetime
filelist = []
dt = datetime.timedelta(1)
num_days=(end_dt-start_dt)/dt
num_days = round(num_days+1)
n_files= num_days*24

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


# Setup counting scheme
counter = 0 # simple counter to keep track of process pace
print('COMPLETE',flush=True)



# -------------------------------------------------------------------------- #
# ----------------------- REFORMAT NETCDF FILES ---------------------------- #
# Reformat the netCDF files to get ready for regridding                      #
# -------------------------------------------------------------------------- #
print('Reformating Input:',end='')

# Reformat to match input for ncl scripts, trim to extent, and split by hour
for file in filelist:
	fp = nc.Dataset(w_input+file,'r')

	# extract day, month, year. this is dependent on input filename
	year = file[0:4]
	month = file[4:6]
	day = file[6:8]

	# Loop through the time dimension in the nc file
	for t in range(24):
		counter = counter + 1
		count_help(counter,n_files)
		if t > 9:
			hour = str(t)
		else: 
			hour = '0'+str(t)
		
		fpout=nc.Dataset(run_dir+w_dir+'out/NLDAS_FORA0125_H.A'+year+month+day+'.'+hour+'00.002.nc','w')
		fpout.createDimension('lat_110',len(lats[lat_m]))
		fpout.createDimension('lon_110',len(lons[lon_m]))
		fpout.createVariable('lat_110','f4',('lat_110'))
		fpout.createVariable('lon_110','f4',('lon_110'))
		fpout['lat_110'][:]=lats[lat_m]
		fpout['lon_110'][:]=lons[lon_m]
		for var in var_names.keys():
			data = fp[var][t,lat_m,lon_m]
			if var == 'wind':
				data = np.sqrt((data**2)/2)
				name = var_names[var][0]
				fpout.createVariable(name,'f4',('lat_110','lon_110'))
				fpout[name][:]=data[:]
				name = var_names[var][1]
			else:
				name = var_names[var]
			fpout.createVariable(name,'f4',('lat_110','lon_110'))
			fpout[name][:]=data[:]
		fpout.close()
	fp.close()
print('')
print('Reformating Complete')




# -------------------------------------------------------------------------- #
# ------------------------- REGRID THE FILES ------------------------------- #
# Regrid the incomming forcing data. This step will take a while             #
# -------------------------------------------------------------------------- #

# Some filenames
grid_file = run_dir+w_dir+'out/NLDAS_FORA0125_H.A'+start_date[0:4]+start_date[5:7]+start_date[8:10]+'.0100.002.nc'
geogrid = run_dir+dom_dir+'geo_em.d01.nc'

# Edit the NCL scripts
fp_ncl = open(run_dir+w_dir+"NLDAS2WRFHydro_regrid.ncl",'r')
ncl_write = ''
for line in fp_ncl:
	if "./output_files/" in line:
		ncl_write += '  outdir  = \"'+run_dir+forcing_dir+'\"'+'\n'
	elif "./out/" in line:
		ncl_write += '  dirm    = \"'+run_dir+w_dir+'out/'+'\"'+'\n'
	elif "./NLDAS2WRFHydro_weight_bilinear.nc" in line:
		ncl_write += '  wgtFileName_bilinear = \"'+run_dir+w_dir+'NLDAS2WRFHydro_weight_bilinear.nc\"\n'
	else:
		ncl_write += line
fp_ncl.close()
fp_out=open(run_dir+w_dir+"temp.ncl",'w')
fp_out.write(ncl_write)
subprocess.run('rm '+run_dir+w_dir+"NLDAS2WRFHydro_regrid.ncl",shell=True)
subprocess.run('mv '+run_dir+w_dir+"temp.ncl "+run_dir+w_dir+"NLDAS2WRFHydro_regrid.ncl", shell=True)
fp_out.close()

# Write Commands
print(os.getcwd())
cmd1 = "ncl 'interp_opt=\"bilinear\"' 'srcGridName=\""+grid_file+"\"' 'dstGridName=\""+geogrid+"\"' "+run_dir+w_dir+"NLDAS2WRFHydro_generate_weights.ncl >"+run_dir+log_dir+"regrid_log1.txt"
cmd2 = "ncl 'srcFileName=\"NLDAS_FORA0125_H.*\"' 'dstGridName=\""+geogrid+"\"' "+run_dir+w_dir+"NLDAS2WRFHydro_regrid.ncl "+'>'+run_dir+log_dir+"regrid_log2.txt"
print(cmd2)
# Run Commands
print('Generating Regrid Weights...',end='',flush=True)
subprocess.run(cmd1,shell=True)
print('COMPLETE',flush=True)
print('Regridding...',end='',flush=True)
subprocess.run(cmd2,shell=True)
print('COMPLETE',flush=True)

# Clean up
clean_cmd=['rm -r out']
if clean:
	for cmd in clean_cmd:
		subprocess.run(cmd,shell=True)




