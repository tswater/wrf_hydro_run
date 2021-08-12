# Creates the Fulldomhires.nc
# tif files :named CHANNELGRID, FLOWDIRECTION (d8 flow dir from tauDEM), STREAMORDER (from tauDEM streamnet function), 
# and TOPOGRAPHY (pit filled elevation grid) are required input, as is a param.txt textfile which contain parameters
# necessary to fill out the metadata of Fulldomhires.nc
# NOTE: this script assumes a Lambert Conformal Conic Projection wwith spherical earth radius 6370000
import netCDF4 as nc
import numpy as np
import subprocess
import rasterio
import sys
import os

# -------------------------------------------------------------------------- #
# ------------------------------ GET INPUT --------------------------------- #
# Recieve input from outside. Order of args must match order of appearance   #
# in the wrf_preprocess.py script. Also some general setup
# -------------------------------------------------------------------------- #
print('Initial GIS configuration...',end='',flush=True)

# Read in argv
argv = sys.argv[1:] 
xp = int(argv[0]) # mass grid physical x/east_west comp y/col/j
yp = int(argv[1]) # mass grid physical y/north_south comp x/row/i
dxp = int(argv[2]) # mass grid dx
dyp = int(argv[3]) # mass grid dy
grid_ratio = int(argv[4]) # ratio between mass grid and routing grid
proj = argv[5] # proj statement
dem = argv[6] # dem location
thresh = int(argv[7])
n_cores = argv[8]
fh_loc = argv[9] # DOMAIN directory
w_dir = argv[10]
run_dir = argv[11]
extent = [0,0,0,0] #
extent[0]=float(argv[12])
extent[1]=float(argv[13])
extent[2]=float(argv[14])
extent[3]=float(argv[15])
truelat1 =float(argv[16])
truelat2 =float(argv[17])
stand_lon=float(argv[18])


# change working directory
os.chdir(run_dir+w_dir)
fh_loc = run_dir+fh_loc

# Load in Parameters from geogrid file
x = xp*grid_ratio # routing grid physical x
y = yp*grid_ratio # routing grid physical y
dx = dxp/grid_ratio
dy = dyp/grid_ratio

# Set up NC file for Fulldom_hires
fp_out=nc.Dataset(fh_loc+"Fulldom_hires.nc",'w')
fp_out.createDimension('y',size=y)
fp_out.createDimension('x',size=x)

# Set up NC file for Metadata
fp_met = nc.Dataset(fh_loc+'GEOGRID_LDASOUT_Spatial_Metadata.nc','w')
fp_met.createDimension('y',size=yp)
fp_met.createDimension('x',size=xp)

# declare some filename
t_file = 'TOPOGRAPHY.tif' # file location for topography tif
c_file = 'CHANNELGRID.tif'
s_file = 'STREAMORDER.tif'
f_file = 'FLOWDIRECTION.tif'

print('COMPLETE',flush=True)

# -------------------------------------------------------------------------- #
# ------------------------------- TAU DEM ---------------------------------- #
# Do TauDEM stuff
# -------------------------------------------------------------------------- #

print('Create Routing Files...',end='',flush=True)

# declare TauDEM commands
cmd1 = 'mpiexec -n '+n_cores+' pitremove '+dem
cmd_rename = 'mv '+dem.split('.')[0]+'fel.tif '+t_file
cmd2 = 'mpiexec -n '+n_cores+' d8flowdir -p '+f_file+' -fel '+t_file+' -sd8 junk.tif'
cmd3 = 'mpiexec -n '+n_cores+' aread8 -p '+f_file+' -ad8 areaad8.tif'

# cmd4 = 'mpiexec -n '+n_cores+' threshold -ssa areaad8.tif -src '+c_file+' -thresh '+thresh
cmd5 = 'mpiexec -n 4 streamnet -fel '+t_file+' -p '+f_file+' -ad8 areaad8.tif -src '+c_file+' -ord '+s_file+' -tree junk.dat -coord junk2.dat -net junk3.shp -w junk4.tif'

# run TauDEM stuff
subprocess.run(cmd1,shell=True)
subprocess.run(cmd_rename,shell=True)
subprocess.run(cmd2,shell=True)
subprocess.run(cmd3,shell=True)
# subprocess.run(cmd4,shell=True) #we had some issues with threshold

# apply threshold and save tif
fp_ad8 = rasterio.open('areaad8.tif','r')
a = fp_ad8.read(1)
a[a<int(thresh)]=0
a[a>=int(thresh)]=1

with rasterio.Env():
    profile = fp_ad8.profile
    profile.update(
        count=1,
        compress='lzw')
    with rasterio.open(c_file, 'w', **profile) as dst:
        dst.write(a, 1)

subprocess.run(cmd5,shell=True)

# clean
subprocess.run('rm junk*',shell=True)
subprocess.run('rm areaad8.tif',shell=True)

print('COMPLETE',flush=True)

# -------------------------------------------------------------------------- #
# ------------------------- MAKE FULLDOM_HIRES ----------------------------- #
# Make the fulldom_hires.nc file from the inputs. Also make latlon grid      #
# -------------------------------------------------------------------------- #

print('Making Fulldom_hires.nc...',end='',flush=True)

# X Y
x_lin = np.linspace(extent[0]+dx/2,extent[2]-dx/2,x)
y_lin = np.linspace(extent[3]-dy/2,extent[1]+dy/2,y)
y_g, x_g = np.meshgrid(y_lin,x_lin,indexing='ij')

fp_out.createVariable('x','d',dimensions=('x'))
fp_out['x'].standard_name = 'projection_x_coordinate'
fp_out['x'].long_name = 'x coordinate of projection'
fp_out['x'].units = 'm'
fp_out['x']._CoordinateAxisType = 'GeoX'
fp_out['x'].resolution=float(dx)
fp_out['x'][:] = x_lin[:]

fp_out.createVariable('y','d',dimensions=('y'))
fp_out['y'].standard_name = 'projection_y_coordinate'
fp_out['y'].long_name = 'y coordinate of projection'
fp_out['y'].units = 'm'
fp_out['y']._CoordinateAxisType = 'GeoY'
fp_out['y'].resolution=float(dy)
fp_out['y'][:] = y_lin[:]

# crs
fp_out.createVariable('crs','c')
fp_out['crs'].transform_name="lambert_conformal_conic"
fp_out['crs'].grid_mapping_name="lambert_conformal_conic"
fp_out['crs'].long_name="CRS definition"
fp_out['crs'].longitude_of_prime_meridian = 0.0
fp_out['crs']._CoordinateAxes = "y x"
fp_out['crs']._CoordinateTransformType = "Projection"
fp_out['crs'].standard_parallel = [float(truelat1),float(truelat2)]
fp_out['crs'].longitude_of_central_meridian=float(stand_lon)
fp_out['crs'].latitude_of_projection_origin=0.0
fp_out['crs'].false_easting=0.0
fp_out['crs'].false_northing=0.0
fp_out['crs'].earth_radius=6370000.0
fp_out['crs'].semi_major_axis=6370000.0
fp_out['crs'].inverse_flattening=0.0

# CHANNELGRID
fp_out.createVariable('CHANNELGRID','i',dimensions=('y','x'))
fp_out['CHANNELGRID'].grid_mapping = 'crs'
cg = rasterio.open(c_file).read(1)
cg[cg==0]=-9999
cg[cg==1]=0
print(cg.shape)
print(str(y)+','+str(x))
fp_out['CHANNELGRID'][:]=cg[:]


# FLOWDIRECTION
# NOTE: tauDEM does not specify a nodata value; wrf requires nodata 255, potential source of error
fp_out.createVariable('FLOWDIRECTION','short',dimensions=('y','x'))
fp_out['FLOWDIRECTION'].grid_mapping='crs'
fd = rasterio.open(f_file).read(1)
fd[(fd>9)|(fd<0)]=255
fp_out['FLOWDIRECTION'][:]=fd[:]

# TOPOGRAPHY
fp_out.createVariable('TOPOGRAPHY','f',dimensions=('y','x'))
fp_out['TOPOGRAPHY'].grid_mapping='crs'
fp_out['TOPOGRAPHY'][:]=rasterio.open(t_file).read(1)[:]

# STREAMORDER
fp_out.createVariable('STREAMORDER','b',dimensions=('y','x'))
fp_out['STREAMORDER'].grid_mapping='crs'
fp_out['STREAMORDER'][:]=rasterio.open(s_file).read(1)

# LAT LON PREPROCESSING
f_xy = open('xy_coords.txt','w')
for i in range(y):
	for j in range(x):
		f_xy.write(str(x_g[i,j])+" "+str(y_g[i,j])+'\n')
f_xy.close()
cmd = "gdaltransform -output_xy -s_srs '"+proj+"' -t_srs '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' < xy_coords.txt > lat_lon.txt"
subprocess.run(cmd,shell=True)

lon_grid = np.zeros((y,x))
lat_grid = np.zeros((y,x))
f_ll = open('lat_lon.txt','r')
for i in range(y):
	for j in range(x):
		ll = f_ll.readline().split()
		lon_grid[i,j] = ll[0]
		lat_grid[i,j] = ll[1]
f_ll.close()

# Remove Intermediary Files
subprocess.run('rm xy_coords.txt',shell=True)
subprocess.run('rm lat_lon.txt',shell=True)

# LATITUDE
fp_out.createVariable('LATITUDE','f',dimensions=('y','x'))
fp_out['LATITUDE'].long_name = "latitude coordinate"
fp_out['LATITUDE'].units = "degrees_north"
fp_out['LATITUDE'].grid_mapping="crs"
fp_out['LATITUDE']._CoordinateAxisType="Lat"
fp_out['LATITUDE']._CoordinateSystems="crs"
fp_out['LATITUDE'][:]=lat_grid[:]

# LONGITUDE
fp_out.createVariable('LONGITUDE','f',dimensions=('y','x'))
fp_out['LONGITUDE'].long_name = "longitude coordinate"
fp_out['LONGITUDE'].units = "degrees_east"
fp_out['LONGITUDE'].grid_mapping="crs"
fp_out['LONGITUDE']._CoordinateAxisType="Lon"
fp_out['LONGITUDE']._CoordinateSystems="crs"
fp_out['LONGITUDE'][:]=lon_grid[:]

# LKSATFAC
fp_out.createVariable('LKSATFAC','f',dimensions=('y','x'))
fp_out['LKSATFAC'].grid_mapping="crs"
fp_out['LKSATFAC'][:]=np.ones((y,x))*1000

# EMPTIES -9999: frxst_pts, basn_msk, LAKEGRID,
fp_out.createVariable('frxst_pts','i',dimensions=('y','x'))
fp_out['frxst_pts'][:]=np.ones((y,x))*-9999
fp_out['frxst_pts'].grid_mapping='crs'

fp_out.createVariable('basn_msk','i',dimensions=('y','x'))
fp_out['basn_msk'][:]=np.ones((y,x))*-9999
fp_out['basn_msk'].grid_mapping='crs'

fp_out.createVariable('LAKEGRID','i',dimensions=('y','x'))
fp_out['LAKEGRID'][:]=np.ones((y,x))*-9999
fp_out['LAKEGRID'].grid_mapping='crs'

# EMPTIES 1: OVROUGHRTFAC, RETDEPRTFAC
fp_out.createVariable('OVROUGHRTFAC','f',dimensions=('y','x'))
fp_out['OVROUGHRTFAC'][:]=np.ones((y,x))
fp_out['OVROUGHRTFAC'].grid_mapping='crs'

fp_out.createVariable('RETDEPRTFAC','f',dimensions=('y','x'))
fp_out['RETDEPRTFAC'][:]=np.ones((y,x))
fp_out['RETDEPRTFAC'].grid_mapping = "crs"

# ADD GLOBAL ATTRIBUTES
fp_out.Conventions = "CF-1.5"
fp_out.GDAL_DataType = "Generic"
fp_out.proj4 = proj

# CLEAN
fp_out.close()
subprocess.run('rm '+t_file,shell=True)
subprocess.run('rm '+s_file,shell=True)
subprocess.run('rm '+f_file,shell=True)
subprocess.run('rm '+c_file,shell=True)

print('COMPLETE',flush=True)

# -------------------------------------------------------------------------- #
# ------------------------ MAKE SPATIAL METADATA --------------------------- #
# Make the spatial metadata file from the inputs.                            #
# -------------------------------------------------------------------------- #

print('Creating Spatial Metadata file...',end='',flush=True)

# X and Y
x_lin = np.linspace(extent[0]+dxp/2,extent[2]-dxp/2,xp)
y_lin = np.linspace(extent[3]-dyp/2,extent[1]+dyp/2,yp)
fp_met.createVariable('x','d',dimensions=('x'))
fp_met['x'].standard_name = 'projection_x_coordinate'
fp_met['x'].long_name = 'x coordinate of projection'
fp_met['x'].units = 'm'
fp_met['x']._CoordinateAxisType = 'GeoX'
fp_met['x'].resolution=float(dxp)
fp_met['x'][:] = x_lin[:]

fp_met.createVariable('y','d',dimensions=('y'))
fp_met['y'].standard_name = 'projection_y_coordinate'
fp_met['y'].long_name = 'y coordinate of projection'
fp_met['y'].units = 'm'
fp_met['y']._CoordinateAxisType = 'GeoY'
fp_met['y'].resolution=float(dyp)
fp_met['y'][:] = y_lin[:]

# crs
fp_met.createVariable('crs','c')
fp_met['crs'].transform_name="lambert_conformal_conic"
fp_met['crs'].grid_mapping_name="lambert_conformal_conic"
fp_met['crs'].long_name="CRS definition"
fp_met['crs'].longitude_of_prime_meridian = 0.0
fp_met['crs']._CoordinateAxes = "y x"
fp_met['crs']._CoordinateTransformType = "Projection"
fp_met['crs'].standard_parallel = [float(truelat1),float(truelat2)]
fp_met['crs'].longitude_of_central_meridian=float(stand_lon)
fp_met['crs'].latitude_of_projection_origin=0.0
fp_met['crs'].false_easting=0.0
fp_met['crs'].false_northing=0.0
fp_met['crs'].earth_radius=6370000.0
fp_met['crs'].semi_major_axis=6370000.0
fp_met['crs'].inverse_flattening=0.0

# ADD GLOBAL ATTRIBUTES
fp_met.Conventions = "CF-1.5"
fp_met.GDAL_DataType = "Generic"
fp_met.proj4 = proj
fp_met.close()

print('COMPLETE',flush=True)
