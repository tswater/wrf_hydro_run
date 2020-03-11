# A script to pre process inputs for wrf_hydro standalone

import os
import subprocess
import numpy as np
import netCDF4 as nc
import rasterio
from datetime import datetime
from datetime import timedelta
import matplotlib.pyplot as plt
# -------------------------------------------------------------------------- #
# ----------------------------- PARAMETERS --------------------------------- #
# this is the only section of the code that the user should have to edit for #
# most basic runs. check the readme.txt file for more details on changes to  #
# other sections of the script.                                              #
# -------------------------------------------------------------------------- #

# Description of Domain
lat_center     = 41.471 # latitude of the center of the domain
lon_center     = -73.74365 # longitude of the center of the domain
e_we           = 16 # number of gridcells in west_east direction
e_sn           = 17 # number of gridcells in north_south direction
dx             = 1000 # gridcell size for mass grid
dy             = 1000 # gridcell size for mass grid
grid_ratio     = 4 # ratio between mass and routing gridcell size
start_date     = '2011-08-26_00:00:00' # format 'yyyy-mm-dd_hh:mm:ss'
end_date       = '2011-08-27_23:00:00' # format 'yyyy-mm-dd_hh:mm:ss'

# Projection information (always lambert conformal conic)
truelat1       = 30.0 # first standard parallel
truelat2       = 60.0 # second standard parallel
stand_lon      = -97.0 # standard_longitude

# File locations
dem_loc        = '/stor/tyche/hydro/private/nc153/nwc/GAEA/data/NEDtiled/NED.vrt' # location of dem.vrt or dem.tif
wps_geo_loc    = '../../software/WPS_GEOG/'# location of WRF Preprocessing Script geological data
wps_loc        = '../../software/WPS/'

# Other Script Parameters 
n_cores        = 4 # number to use in parallelized processes MINIMUM OF 2
clean          = True # set to False to keep the working directory


# -------------------------------------------------------------------------- #
# --------------------------- GENERAL SETUP -------------------------------- #
# General setup of the script and declaration of some constants              #
# -------------------------------------------------------------------------- #
# Lets get started
print('WRF Hydro Preprocessing Script by Tyler Waterman, Duke University')
print('Setup...',end='',flush=True)

# Setup Directories
tbl_dir = 'tables/'
namelist_dir = 'namelists/'
scripts_dir = 'scripts/'
w_dir = 'workspace/'
dom_dir = w_dir+'DOMAIN/'
log_dir = 'logs/'
run_dir = os.getcwd()+'/'
subprocess.run('mkdir '+w_dir,shell=True)
subprocess.run('mkdir '+dom_dir,shell=True)

# Make datetimes
start_dt = datetime.strptime(start_date,'%Y-%m-%d_%X')
end_dt = datetime.strptime(end_date,'%Y-%m-%d_%X')
dt = timedelta(hours=1)
sim_time = (end_dt-start_dt)/dt # simulation time in hours

# Declare proj4 statement
proj = "+proj=lcc +a=6370000 +b=6370000 +lon_0="+str(stand_lon)+" +lat_1="+\
       str(truelat1)+" +lat_2="+str(truelat2)

print('Complete',flush=True)

# Set mass grid x and y
x = e_we - 1
y = e_sn - 1



# -------------------------------------------------------------------------- #
# ---------------------------- DOMAIN SETUP -------------------------------- #
# Edit the namelist.wps file, and then run geogrid.exe in the WPS folder to  #
# create a geo_emd01.nc file which is then placed in DOMAIN. From the        #
# geogrid file, extract the extent of the mass grid                          #
# -------------------------------------------------------------------------- #

print('Editing namelist.wps...',end='',flush=True)

# edit the namelist.wps
wps_write =''
fp_wps = open(namelist_dir+'template/namelist.wps','r')
for line in fp_wps:
	if 'start_date' in line:
		wps_write += ' start_date = \''+start_date+'\','+'\n'
	elif 'end_date' in line:
		wps_write += ' end_date   = \''+end_date+'\','+'\n'
	elif 'e_we' in line:
		wps_write += ' e_we              ='+str(e_we)+','+'\n'
	elif 'e_sn' in line:
		wps_write += ' e_sn              ='+str(e_sn)+','+'\n'
	elif 'dx' in line:
		wps_write += 'dx = '+str(dx)+','+'\n'
	elif 'dy' in line:
		wps_write += 'dy = '+str(dy)+','+'\n'
	elif 'ref_lat' in line:
		wps_write += ' ref_lat   = '+str(lat_center)+','+'\n'
	elif 'ref_lon' in line:
		wps_write += ' ref_lon   = '+str(lon_center)+','+'\n'
	elif 'truelat1' in line:
		wps_write += ' truelat1  = '+str(truelat1)+','+'\n'
	elif 'truelat2' in line:
		wps_write += ' truelat2  = '+str(truelat2)+','+'\n'
	elif 'stand_lon' in line:
		wps_write += ' stand_lon = '+str(stand_lon)+','+'\n'
	elif 'geog_data_path' in line:
		wps_write += ' geog_data_path = \''+wps_geo_loc+'\''+'\n'
	else:
		wps_write += line
fp_wps.close()
fp_out = open(namelist_dir+'namelist.wps','w')
fp_out.write(wps_write)
fp_out.close()
print(' Complete',flush=True)

# Geogrid generation
print('Geogrid generation...',end='',flush=True)
os.chdir(wps_loc)
subprocess.run('./geogrid.exe >'+run_dir+w_dir+'log.txt',shell=True)
os.chdir(run_dir)
subprocess.run('mv '+wps_loc+'geo_em.d01.nc '+run_dir+dom_dir,shell=True)
print(' Complete',flush=True)
print('Skipping geogrid generation',flush=True)

# Extract Extent from Geogrid
print('Extracting extent from geogrid...',end='',flush=True)
fp_geo = nc.Dataset(run_dir+dom_dir+'geo_em.d01.nc','r')
os.chdir(w_dir)

# Get corner values in lat lon and then reproject
corner_lats = fp_geo.corner_lats[0:4]
corner_lons = fp_geo.corner_lons[0:4]
fp_geo.close()

fp_co_ll = open('lat_lon.txt','w')
for i in range(4):
	fp_co_ll.write(str(corner_lons[i])+" "+str(corner_lats[i])+'\n')
fp_co_ll.close()
cmd = "gdaltransform -output_xy -s_srs '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' -t_srs '"+proj+"' < lat_lon.txt > xy_out.txt"
subprocess.run(cmd,shell=True)

fp_xy=open('xy_out.txt','r')
x_val=np.ones((4,))
y_val=np.ones((4,))
for i in range(4):
	xy = fp_xy.readline().split()
	x_val[i]=xy[0]
	y_val[i]=xy[1]
fp_xy.close()

# Calculate extent from corner values
extent=[0,0,0,0]
extent[0]=np.round(np.mean(x_val))-dx*x/2 # get left extent for x
extent[2]=np.round(np.mean(x_val))+dx*x/2 # get right extent for x
extent[1]=np.round(np.mean(y_val))-dy*y/2 # get lower extent for y
extent[3]=np.round(np.mean(y_val))+dy*y/2 # get upper extent for y

# Clean up text files
subprocess.run('rm xy_out.txt',shell=True)
subprocess.run('rm lat_lon.txt',shell=True)
os.chdir(run_dir)
print(' Complete',flush=True)
print('FINISHED DOMAIN PREPROCESSING')



# -------------------------------------------------------------------------- #
# ----------------------- GIS PROCESSING/ROUTING --------------------------- #
# Extract the appropriate domain from a DEM, then using TauDEM create grids  #
# with routing information and other hydrology information. Compile these    #
# grids into Fulldom_hires.nc
# -------------------------------------------------------------------------- #
ex_string = str(extent[0])+' '+str(extent[1])+' '+str(extent[2])+' '+str(extent[3])
print('Trimming DEM...',end='',flush=True)
os.chdir(w_dir)
cmd = "gdalwarp -t_srs '"+proj+"' -te "+ex_string+" -tr "+str(dx/grid_ratio)+ \
              " "+str(dy/grid_ratio)+" "+dem_loc+" dem.tif >" \
              +run_dir+w_dir+'demlog.txt'
subprocess.run(cmd,shell=True)
print('Complete',flush=True)
os.chdir(run_dir) 
dem=run_dir+w_dir+'dem.tif'

thresh=[5,10,15,25,50,75,100,150,250,500]

t_file = w_dir+'junk/wrfT.tif' # file location for topography tif
s_file = w_dir+'junk/wrfS.tif'
f_file = w_dir+'junk/wrfF.tif'
subprocess.run('mkdir '+w_dir+'junk',shell=True)

# Run beginning steps
n_cores = str(n_cores)

cmd1 = 'mpiexec -n '+n_cores+' pitremove '+dem+' >/dev/null'
cmd_rename = 'mv '+dem.split('.')[0]+'fel.tif '+t_file+' >/dev/null'
cmd2 = 'mpiexec -n '+n_cores+' d8flowdir -p '+f_file+' -fel '+t_file+' -sd8 '+w_dir+'junk/junk.tif'+' >/dev/null'
cmd3 = 'mpiexec -n '+n_cores+' aread8 -p '+f_file+' -ad8 '+w_dir+'junk/areaad8.tif'+' >/dev/null'

subprocess.run(cmd1,shell=True)
subprocess.run(cmd_rename,shell=True)
subprocess.run(cmd2,shell=True)
subprocess.run(cmd3,shell=True)
print('INITIAL SETUP SUCCESSFUL')

# Retrieve the array and clear
fp = rasterio.open(w_dir+'junk/areaad8.tif')
a = fp.read(1)
fp.close()
subprocess.run('rm -r '+w_dir+'junk/',shell=True)

# Iterate through thresholds
i=0
for t in thresh:
        plt.subplot(2,5,i+1)
        i = i+1
        a_p = np.copy(a)
        a_p[a<t]=0
        a_p[a>=t]=1
        plt.imshow(a_p)
        plt.title(str(t))

plt.show()

if clean:
	subprocess.run('rm -r '+w_dir,shell=True)


