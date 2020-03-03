# A script to pre process inputs for wrf_hydro standalone

import os
import subprocess
import numpy as np
import netCDF4 as nc
import rasterio
from datetime import datetime
from datetime import timedelta
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
restart_loc    = '' # location of LSM restart files if any
restart_loch   = '' # location of hydro restart files if any
dem_loc        = '/stor/tyche/hydro/private/nc153/nwc/GAEA/data/NEDtiled/NED.vrt' # location of dem.vrt or dem.tif
forcing_loc    = '/stor/soteria/hydro/shared/data/PCF/1hr/daily/' # location of forcing files
wps_geo_loc    = '../../software/WPS_GEOG/'# location of WRF Preprocessing Script geological data
wps_loc        = '../../software/WPS-4.1/'
dem_loc        = '/home/tsw35/soteria/software/WRF_hydro_standalone/tst_ny/dem/dem.tif'

# Other Script Parameters 
thresh         = 50 # channel determination threshold; see readme.txt
n_cores        = 3 # number to use in parallelized processes MINIMUM OF 2
clean          = False # set to False to keep the working directory
restart        = False # set to True to look for restart files
setup_domain   = True # set to False to skip geogrid generation
setup_wrfinput = True # set to False to skip creation of wrfinput
setup_routing  = True # set to False to skip generation of fulldom_hires.nc
setup_forcing  = True # set to False to skip generation of forcing files
setup_hydro    = True # set to False to skip autofill of hydro.namelist
setup_hrldas   = True # set to False to skip autofill of namelist.hrldas


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
dom_dir = 'DOMAIN/'
forc_dir= 'FORCING/'
log_dir = 'logs/'
run_dir = os.getcwd()+'/'
subprocess.run('mkdir '+w_dir,shell=True)
subprocess.run('mkdir '+forc_dir,shell=True)
subprocess.run('mkdir '+dom_dir,shell=True)
subprocess.run('mkdir '+log_dir,shell=True)

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
if setup_domain:
	print('Geogrid generation...',end='',flush=True)
	os.chdir(wps_loc)
	subprocess.run('./geogrid.exe >'+run_dir+log_dir+'geogrid_log.txt',shell=True)
	os.chdir(run_dir)
	subprocess.run('mv '+wps_loc+'geo_em.d01.nc '+dom_dir,shell=True)
	print(' Complete',flush=True)
else:
	print('Skipping geogrid generation',flush=True)

# Extract Extent from Geogrid
print('Extracting extent from geogrid...',end='',flush=True)
fp_geo = nc.Dataset(dom_dir+'geo_em.d01.nc','r')
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
# ------------------------ INITIAL INPUT SETUP ----------------------------- #
# Create the initial input for the model Yes this is necessary with restart  #
# -------------------------------------------------------------------------- #
if setup_wrfinput:
	print('Creating wrfinput...', end='',flush=True)
	subprocess.run('mv '+scripts_dir+'create_wrfinput.R '+w_dir,shell=True)
	os.chdir(w_dir)
	subprocess.run('./create_wrfinput.R --geogrid='+run_dir+dom_dir+'geo_em.d01.nc >'+run_dir+log_dir+'wrfinput_log.txt',shell=True)
	os.chdir(run_dir)
	subprocess.run('mv '+w_dir+'create_wrfinput.R '+scripts_dir,shell=True)
	subprocess.run('mv '+w_dir+'wrfinput_d01.nc '+dom_dir,shell=True)
	subprocess.run('rm '+w_dir+'*',shell=True)
	print(' Complete',flush=True)



# -------------------------------------------------------------------------- #
# ----------------------- GIS PROCESSING/ROUTING --------------------------- #
# Extract the appropriate domain from a DEM, then using TauDEM create grids  #
# with routing information and other hydrology information. Compile these    #
# grids into Fulldom_hires.nc
# -------------------------------------------------------------------------- #
ex_string = str(extent[0])+' '+str(extent[1])+' '+str(extent[2])+' '+str(extent[3])
if setup_routing:
	print('Creating DEM...',end='',flush=True)
	os.chdir(w_dir)
	cmd = "gdalwarp -t_srs '"+proj+"' -te "+ex_string+" -tr "+str(dx/grid_ratio)+ \
              " "+str(dy/grid_ratio)+" "+dem_loc+" dem.tif >" \
              +run_dir+log_dir+'demlog.txt'
	subprocess.run(cmd,shell=True)
	print('Complete',flush=True)
	print('BEGIN GIS PROCESSING SCRIPT')
	gis_cmd = 'python '+run_dir+scripts_dir+'wrf_gisprocess.py '+str(x)+' '+str(y)+' '+str(dx)+' '+ \
	          str(dy)+' '+str(grid_ratio)+' \''+proj+'\' dem.tif '+str(thresh)+' '+\
	          str(n_cores)+' '+dom_dir+' '+w_dir+' '+run_dir+' '+ex_string \
		  +' '+str(truelat1)+' '+str(truelat2)+' '+str(stand_lon)
	subprocess.run(gis_cmd,shell=True)
	os.chdir(run_dir) 

# -------------------------------------------------------------------------- #
# ------------------------ SETUP FORCING FILES ----------------------------- #
# Regrid the incomming forcing data. This step will take a while
# -------------------------------------------------------------------------- #
if setup_forcing:
	for i in range(n_cores):
		subprocess.run('cp '+scripts_dir+'NLDAS2WRFHydro_regrid.ncl '\
		     +w_dir+'NLDAS2WRFHydro_regrid'+str(i)+'.ncl ',shell=True)
	subprocess.run('cp '+scripts_dir+'parallel_regrid.py '+w_dir,shell=True)
	subprocess.run('cp '+scripts_dir+'NLDAS2WRFHydro_generate_weights.ncl '\
	     +w_dir,shell=True)
	os.chdir(w_dir)
	forcing_cmd = 'python ' +run_dir+scripts_dir+'wrf_weather_regrid.py '\
             +start_date+' '+end_date+' '+forcing_loc+' '+str(n_cores)+' \''+\
             proj+'\' '+w_dir+' '+dom_dir+' '+forc_dir+' '+log_dir+' '+run_dir
	subprocess.run(forcing_cmd,shell=True)
	os.chdir(run_dir)
	if clean:
		subprocess.run('rm '+w_dir+'NLDAS2WRFHydro_*',shell=True)


# -------------------------------------------------------------------------- #
# --------------------------- SETUP NAMELISTS ------------------------------ #
# Autofills what is possible for the namelist.hrldas and hydro.namelist file #
# -------------------------------------------------------------------------- #
if setup_hydro:
	print('Autofilling hydro.namelist...',end='',flush=True)
	hydro_write=''
	fp_hydro = open(namelist_dir+'template/hydro.namelist','r')
	for line in fp_hydro:
		if 'RESTART_FILE' in line: # restart filename
			if restart:
				hydro_write+= 'RESTART_FILE  = \''+\
				              restart_loch+'\''+'\n'
			else:
				hydro_write+=line
		elif 'DXRT' in line: # DXRT grid spacing routing
			hydro_write+= 'DXRT = '+str(dx)+'\n'
		elif 'AGGFACTRT' in line: # AGGFACTRT grid ratio
			hydro_write+= 'AGGFACTRT = '+str(grid_ratio)+'\n'
		else:
			hydro_write+= line
	fp_hydro.close()
	fp_out = open(namelist_dir+'hydro.namelist','w')
	fp_out.write(hydro_write)
	fp_out.close()
	print('COMPLETE',flush=True)
if setup_hrldas:
	print('Autofilling namelist.hrldas...',end='',flush=True)
	hrldas_write=''
	fp_hrldas = open(namelist_dir+'template/namelist.hrldas','r')
	for line in fp_hrldas:
		if 'START_YEAR' in line: # start time
			hrldas_write+=' START_YEAR  = '+str(start_dt.year)+'\n'
		elif 'START_MONTH' in line:
			hrldas_write+=' START_MONTH = '+start_date[5:7]+'\n'
		elif 'START_DAY' in line:
			hrldas_write+=' START_DAY   = '+start_date[8:10]+'\n'
		elif 'START_HOUR' in line:
			hrldas_write+=' START_HOUR  = '+start_date[11:13]+'\n'
		elif 'START_MIN' in line:
			hrldas_write+=' START_MIN   = '+start_date[14:16]+'\n'
		elif ('RESTART_FILENAME_REQUESTED' in line) and (restart):
			hrldas_write+=' RESTART_FILENAME_REQUESTED = \"'+\
			              restart_loc+'\"'+'\n'
		elif 'KHOUR' in line: # runlength
			hrldas_write+=' KHOUR = '+str(int(round(sim_time)))+'\n'
		else:
			hrldas_write+= line
	fp_hrldas.close()
	fp_out = open(namelist_dir+'namelist.hrldas','w')
	fp_out.write(hrldas_write)
	fp_out.close()
	print('COMPLETE',flush=True)



# -------------------------------------------------------------------------- #
# --------------------------- FINAL SETUP ---------------------------------- #
# Creates 3D soil properties from a table, final cleaning, add to path 
# -------------------------------------------------------------------------- #
print('Creating 3D Soil Properties...',end='',flush=True)
soil_write =''
fp_soil = open(scripts_dir+'create_SoilProperties.R','r')
for line in fp_soil:
	if 'geoFile <- ' in line:
		soil_write += 'geoFile <- \"'+run_dir+dom_dir+'geo_em.d01.nc\"'+'\n'
	elif 'soilParamFile <- ' in line:
		soil_write += 'soilParamFile <- \"'+run_dir+tbl_dir+'SOILPARM.TBL\"'+'\n'
	elif 'mpParamFile <- 'in line:
		soil_write += 'mpParamFile <- \"'+run_dir+tbl_dir+'MPTABLE.TBL\"'+'\n'
	elif 'genParamFile <- 'in line:
		soil_write += 'genParamFile <- \"'+run_dir+tbl_dir+'GENPARM.TBL\"'+'\n'
	elif 'hydParamFile <- 'in line:
		soil_write += 'hydParamFile <- \"'+run_dir+tbl_dir+'HYDRO.TBL\"'+'\n'
	elif 'slpropFile <- 'in line:
		soil_write += 'slpropFile <- \"'+run_dir+dom_dir+'soil_properties.nc\"'+'\n'
	elif 'hyd2dFile <-' in line:
		soil_write += 'hyd2dFile <- \"'+run_dir+'hydro2dtbl.nc\"'+'\n'
	else:
		soil_write+= line
fp_soil.close()
fp_out = open(w_dir+'create_SoilProperties.R','w')
fp_out.write(soil_write)
fp_out.close()

# Actually run the file
os.chdir(w_dir)
subprocess.run('chmod +x create_SoilProperties.R',shell=True)
subprocess.run('Rscript create_SoilProperties.R >'+run_dir+log_dir+'soilproplog.txt',shell=True)
os.chdir(run_dir)
print('COMPLETE',flush=True)
subprocess.run('mkdir OUT',shell=True)
subprocess.run('export PATH=$PATH:'+run_dir,shell=True)
if clean:
	subprocess.run('rm -r '+w_dir,shell=True)

print('----- WRF PREPROCESSING COMPLETE -----')

