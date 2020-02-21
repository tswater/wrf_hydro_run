RUN WRF HYDRO!

REQUIRED: 
- conda environment with the following packages installed:
	ncl
	r-ncdf4
	netcdf4
	nco
	r-optparse
	r-stringr
	r-plyr
	rasterio
	gdal
- Compiled WRF and WRF hydro
- Compiled WPS
- Compiled TauDEM (and gdal library added to LD_LIBRARY_PATH)
- A DEM
- Forcing data



First, fill in the text on the wrf_preprocess.py file in the top section. No
other changes should be necessary in this file. 

Run wrf_preprocess.py

Check the output in the DOMAIN and FORCING folders. In DOMAIN you should have
4 netcdf files, metadata, fulldom_hires, wrfinput and geo_emd01.nc

In the Forcing folder, there should be many files

Add the working folder to the path. (export PATH=$PATH:[run directory])

Run run_wrf.py [n]  where n is the number of processes requested

All output should be in the OUT folder created by running run_wrf.py

To clean the output AND preprocessing from the folder, run scripts/clean.py
