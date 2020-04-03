# A function to extract important variables from WRF output, chunk, and consolidate by month
# goal: monthly netcdf files with only the variables below, chunked as 
# Each (monthly) netCDF4 file should have dimensions of approximately 720,y_dim,x_dim

# DETERMINING CHUNKSIZE
# assuming the user wants relatively equal i/o for 2D and 1D (time), then a
# okay default chunk shape is time_dim/N**2,y_dim/N,x_dim/N 
# for some N (assuming x and y relatively equal sizes)
# The resulting block size should be ????
# But apparently less than 4096 for most desktop systems

# VARIABLES WE WANT
# SNEQ : Snow Water Equivalent
# SOIL_T : Soil Temperature
# TRAD : Surface Temperature
# ALBEDO : albedo
# LH : latent heat to atm
# HFX : sensible heat to atm
# GRDFLX: heat flux to soil
# FIRA : net LW to atm
# FSA : total absorbed SW
# EMISS : emissivity
# RAINRATE : precip
# LWFORC : LW forcing
# SWFORC : SW forcing

from os import listdir
from datetime import datetime
from datetime import timedelta
import netCDF4 as nc
import subprocess

# ----------------------------- #
# CONSTANTS and USER PARAMETERS #
# ----------------------------- #
# Define Variables to Keep from Output
var_list = ['LH','HFX','TRAD','GRDFLX','FSA','SOIL_M','CM','CH','SWFORC','LWFORC','RAINRATE']
# soil moisture needs to be handled specially, (time, y, soillayer, x) 

# Other Paramters
chunk_shape = (7,12,12)
remove_old  = False # mark true to remove original output
in_dir      = 'OUT/'
out_dir     = 'COMPRESSED_OUTPUT/'

# ---------------------- #
# POST PROCESSING BASICS #
# ---------------------- #
# Extract a proper file list with only LSM output
files = listdir(in_dir)
LSM_files = []
for file in files:
	if 'LDASOUT' in file:
		LSM_files.append(file)
LSM_files.sort()

# Extract start and end date
start_date = LSM_files[0][0:10]
end_date = LSM_files[-1][0:10]

# Extract some basic information
xsize=nc.Dataset(in_dir+file,'r').dimensions['x'].size
ysize=nc.Dataset(in_dir+file,'r').dimensions['y'].size

# Initialize some variables for the loop
date = '200010'
fp_out = ''
t = 0

# Create out_dir
subprocess.run('mkdir '+out_dir,shell=True)

# Loop through all the files
for file in LSM_files:
	fp_in = nc.Dataset(in_dir+file,'r')
	if date != file[0:6]:
		# Close old file and Initialize new output file...
		try:
			fp_out.close()
		except Exception:
			pass
		print(date[0:4]+'-'+date[4:6]+' Processed',flush=True)
		date = file[0:6]
		fp_out= nc.Dataset(out_dir+'LSM_OUT_'+date+'.nc','w')
		t = 0
		# Determine Size of time dimension
		tsize = 0
		for file in LSM_files:
			if date in file:
				tsize += 1
		
		# Add dimensions
		fp_out.createDimension('time',tsize)
		fp_out.createDimension('y',ysize)
		fp_out.createDimension('x',xsize)
		fp_out.createDimension('reference_time',1)

		# Create Variables
		for var in var_list:
			fp_out.createVariable(var,'f4',('time','y','x'),chunksizes=chunk_shape,fill_value=-9999)
			
			# Load in attributes
			for attr in fp_in[var].ncattrs():
				if attr == '_FillValue':
					continue
				else:
					exec_str='fp_out[\''+var+'\'].'+attr+'=fp_in[\''+var+'\'].'+attr
					exec(exec_str)
			if var == 'SOIL_M':
				fp_out['SOIL_M'].descrip='Weighted Average over 2m'
				
		# Manage time independent variables (x,y,crs,time,reference_time)
		# x
		var = 'x'
		fp_out.createVariable('x','f4',('x'))
		for attr in fp_in[var].ncattrs():
			exec_str='fp_out[\''+var+'\'].'+attr+'=fp_in[\''+var+'\'].'+attr
			exec(exec_str)
		fp_out[var][:]=fp_in[var][:]
		
		# y 
		var = 'y'
		fp_out.createVariable('y','f4',('y'))
		for attr in fp_in[var].ncattrs():
			exec_str='fp_out[\''+var+'\'].'+attr+'=fp_in[\''+var+'\'].'+attr
			exec(exec_str)
		fp_out[var][:]=fp_in[var][:]

		# crs
		var = 'crs'
		fp_out.createVariable('crs','c')
		for attr in fp_in[var].ncattrs():
			exec_str='fp_out[\''+var+'\'].'+attr+'=fp_in[\''+var+'\'].'+attr
			exec(exec_str)
		
		# time
		var = 'time'
		fp_out.createVariable('time','i2',('time'))
		for attr in fp_in[var].ncattrs():
			exec_str='fp_out[\''+var+'\'].'+attr+'=fp_in[\''+var+'\'].'+attr
			exec(exec_str)
		
		# reference_time
		var = 'reference_time'
		fp_out.createVariable('reference_time','i2',('reference_time'))
		for attr in fp_in[var].ncattrs():
			if 'valid' in attr:
				pass
			else:
				exec_str='fp_out[\''+var+'\'].'+attr+'=fp_in[\''+var+'\'].'+attr
				exec(exec_str)
		fp_out[var][0]=fp_in[var][0]

	# Load in the variables for
	for var in var_list:
		if var == 'SOIL_M':
			data = fp_in[var][0,:,0,:]*.1/2+fp_in[var][0,:,1,:]*.3/2+ \
                               fp_in[var][0,:,2,:]*.6/2+fp_in[var][0,:,3,:]*1/2
			fp_out[var][t,:,:]=data
		else:
			fp_out[var][t,:,:]=fp_in[var][0,:,:]
	# Increase t
	t = t+1


	
