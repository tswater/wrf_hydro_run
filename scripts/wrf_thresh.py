# Script to help (arbitrarily) determine the appropriate threshold for stream channel determination
# all this actuall does is set up an areaad8.tif for the second script to work with
import subprocess
import numpy as np
import rasterio
import matplotlib.pyplot as plt

# THIS IS THE PARAMETER SECTION AND THE ONLY SECTION YOU SHOULD NEED TO CHANGE
# ----------------------------------------------------------------------------
n_cores = '4' #number of cores for tauDEM 
dem = 'dem.tif' # dem location; DEM is assumed to be trimed to extent
# ----------------------------------------------------------------------------
# END OF PARAMETER SECTION




# MAIN PROGRAM
# ----------------------------------------------------------------------------
thresh=[5,10,15,25,50,75,100,150,250,500]

t_file = 'junk/wrfT.tif' # file location for topography tif
s_file = 'junk/wrfS.tif'
f_file = 'junk/wrfF.tif'
subprocess.run('mkdir junk',shell=True)

# Run beginning steps
cmd1 = 'mpiexec -n '+n_cores+' pitremove '+dem+' >/dev/null'
cmd_rename = 'mv '+dem.split('.')[0]+'fel.tif '+t_file+' >/dev/null'
cmd2 = 'mpiexec -n '+n_cores+' d8flowdir -p '+f_file+' -fel '+t_file+' -sd8 junk/junk.tif'+' >/dev/null'
cmd3 = 'mpiexec -n '+n_cores+' aread8 -p '+f_file+' -ad8 junk/areaad8.tif'+' >/dev/null'

subprocess.run(cmd1,shell=True)
subprocess.run(cmd_rename,shell=True)
subprocess.run(cmd2,shell=True)
subprocess.run(cmd3,shell=True)
print('INITIAL SETUP SUCCESSFUL')

# Retrieve the array and clear
fp = rasterio.open('junk/areaad8.tif')
a = fp.read(1)
fp.close()
subprocess.run('rm -r junk/',shell=True)

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
