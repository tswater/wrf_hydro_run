# A regridding script for the weather data in the cluster
# NOTE almost any major change to this script requires similar or other
# changes to the parallel_regrid.py, especially if there is a change in the 
# input netCDF forcing files, since both scripts read them.
import netCDF4 as nc
import numpy as np
import subprocess
import rasterio
import datetime
from datetime import timedelta
import os
import sys

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
ncores = int(argv[3])
proj = argv[4] # proj statement
w_dir = argv[5]
dom_dir = argv[6]
forcing_dir = argv[7]
log_dir = argv[8]
run_dir = argv[9]

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
fp_1.close()

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


print('COMPLETE',flush=True)



# -------------------------------------------------------------------------- #
# ----------------------- REFORMAT NETCDF FILES ---------------------------- #
# Reformat one netCDF file to use as template for weight generation/regrid   #
# -------------------------------------------------------------------------- #
file = filelist[0]
fp = nc.Dataset(w_input+file,'r')
year = file[0:4]
month = file[4:6]
day = file[6:8]
hour = '12'
out_name=run_dir+w_dir+'out/NLDAS_FORA0125_H.A'
out_name=out_name+year+month+day+'.'+hour+'00.002.nc'
fpout=nc.Dataset(out_name,'w')
fpout.createDimension('lat_110',len(lats[lat_m]))
fpout.createDimension('lon_110',len(lons[lon_m]))
fpout.createVariable('lat_110','f4',('lat_110'))
fpout.createVariable('lon_110','f4',('lon_110'))
fpout['lat_110'][:]=lats[lat_m]
fpout['lon_110'][:]=lons[lon_m]
for var in var_names.keys():
	data = fp[var][12,lat_m,lon_m]
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



# -------------------------------------------------------------------------- #
# ------------------------- REGRID THE FILES ------------------------------- #
# Regrid the incomming forcing data. This step will take a while             #
# -------------------------------------------------------------------------- #

# Some filenames
grid_file = out_name
geogrid = run_dir+dom_dir+'geo_em.d01.nc'


# Regrid Weights
cmd1 = "ncl 'interp_opt=\"bilinear\"' 'srcGridName=\""+grid_file+\
       "\"' 'dstGridName=\""+geogrid+"\"' "+run_dir+w_dir+\
       "NLDAS2WRFHydro_generate_weights.ncl >"+\
       run_dir+log_dir+"regrid_weights_log.txt"

# Run Commands
print('Generating Regrid Weights...',end='',flush=True)
subprocess.run(cmd1,shell=True)
print('COMPLETE',flush=True)
print('Regridding...',end='',flush=True)

# Run Parallel Regrid
print(os.getcwd(),flush=True)
runargs=w_input+' '+w_dir+' '+forcing_dir+' '+log_dir+' '+run_dir+' '+geogrid\
     +' '+'"'+str(filelist)+'"'
mpi_cmd = 'mpiexec -n '+str(ncores)+' python parallel_regrid.py '+runargs
subprocess.run(mpi_cmd,shell=True)
print('COMPLETE',flush=True)



