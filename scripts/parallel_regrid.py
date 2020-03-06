# -------------------------------------------------------------------------- #
# Script meant to be run parallel from wrf_weather_regrid.py
import glob, sys
from mpi4py import MPI
import ast
import netCDF4 as nc
import numpy as np
import os
import subprocess
import mpi4py.rc
#import parUtils as par
#mpi4py.rc.finalize = False
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# move to directory, create a list of subfolders using glob
# dir = directories[rank]

# WORKFLOW 
# Give each one a portion of the initial filelist
# Put each portion (reformated) into its own directory
# Place a copy of the ncl script unique for each rank 
#	edit the out part and the name of the script
# Run each script for each subdir
# Yay done


# -------------- #
# LOAD ARGUMENTS #
# -------------- #
argv = sys.argv[1:]
w_input = argv[0]
w_dir = argv[1]
forcing_dir = argv[2]
log_dir = argv[3]
run_dir = argv[4]
geogrid = argv[5]
filelist = ast.literal_eval(argv[6])

# --------------------- #
# INITIAL SETUP STUFFS  #
# --------------------- #
if (rank == 0):
	sys.exit()

fp_geo = nc.Dataset(geogrid,'r')
corner_lats = fp_geo.corner_lats[0:4]
corner_lons = fp_geo.corner_lons[0:4]
fp_geo.close()

# Set up approximate extent map
fp_1=nc.Dataset(w_input+filelist[0],'r')
lats = fp_1['lat'][:]
lons = fp_1['lon'][:]
min_lon = np.floor(np.min(corner_lons))
min_lat = np.floor(np.min(corner_lats))
max_lon = np.ceil(np.max(corner_lons))
max_lat = np.ceil(np.max(corner_lats))
lat_m = (lats>=min_lat) & (lats<=max_lat) # mask for lat extent
lon_m = (lons>min_lon) & (lons<=max_lon)  # mask for lon extent
fp_1.close()

# declare some names
var_names = {'tair':'TMP_110_HTGL','spfh':'SPF_H_110_HTGL',
    'wind':['U_GRD_110_HTGL','V_GRD_110_HTGL'],'psurf':'PRES_110_SFC',
    'precip':'A_PCP_110_SFC_acc1h','swdown':'DSWRF_110_SFC','lwdown':'DLWRF_110_SFC'}



# --------------------- #
# REFORMAT NETCDF FILES #
# --------------------- #
out_dir = 'out'+str(rank)+'/'
subprocess.run('mkdir '+run_dir+w_dir+out_dir,shell=True)
n = len(filelist)
for file in filelist[rank-1::size-1]:
	fp = nc.Dataset(w_input+file,'r') 
	
	# extract day, month year
	year = file[0:4]
	month = file[4:6]
	day = file[6:8]
	
	for t in range(24):
		if t > 9:
			hour = str(t)
		else:
			hour = '0'+str(t)
		out_name=run_dir+w_dir+out_dir+'NLDAS_FORA0125_H.A'
		out_name=out_name+year+month+day+'.'+hour+'00.002.nc'
		fpout=nc.Dataset(out_name,'w')
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

# ---------------- #
# REGRID THE FILES #
# ---------------- #
# Rewrite the ncl file to match the specific filelist
regrid_ncl = run_dir+w_dir+"NLDAS2WRFHydro_regrid"+str(rank)+".ncl"
fp_ncl = open(regrid_ncl,'r')
ncl_write = ''
for line in fp_ncl:
	if "./output_files/" in line:
		ncl_write += '  outdir  = \"'+run_dir+forcing_dir+'\"'+'\n'
	elif "./out/" in line:
		ncl_write += '  dirm    = \"'+run_dir+w_dir+out_dir+'\"'+'\n'
	elif "./NLDAS2WRFHydro_weight_bilinear.nc" in line:
		ncl_write += '  wgtFileName_bilinear = \"'+run_dir+w_dir+'NLDAS2WRFHydro_weight_bilinear.nc\"\n'
	else:
		ncl_write += line
fp_ncl.close()

# Write to the file
fp_out=open(run_dir+w_dir+"temp"+str(rank)+".ncl",'w')
fp_out.write(ncl_write)
subprocess.run('rm '+regrid_ncl,shell=True)
subprocess.run('mv '+run_dir+w_dir+"temp"+str(rank)+".ncl "+regrid_ncl, shell=True)
fp_out.close()

# Regrid
print('Regrid for rank: '+str(rank),flush=True)
cmd = "ncl 'srcFileName=\"NLDAS_FORA0125_H.*\"' 'dstGridName=\""+geogrid+"\"' "+ \
       regrid_ncl+' >'+run_dir+log_dir+"regrid_log"+str(rank)+".txt"
print(cmd,flush=True)
subprocess.run(cmd,shell=True)
print(str(rank)+': complete',flush=True)
