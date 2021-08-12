RUN WRF HYDRO!

REQUIRED: 
- conda environment with the following packages installed (this is innacurrate, ask me for updated conda environment):
	ncl
	r-ncdf4
	netcdf4
	nco
	r-optparse
	r-stringr
	r-plyr
	rasterio
	gdal
- Compiled WRF hydro
- Compiled WPS
- Compiled TauDEM (and the gdal library added to LD_LIBRARY_PATH)
- A DEM from NED
- Forcing data
- Optional: Hiresolution data from NLCD and Polaris Soil types
	+ note: the polaris soil type input will have to be 
	  specially created to match expectations of the script
	  
NOTE: script only allows for lambert conformal conic projection, and a
significant rewrite will likely be required to get it to work for other
projections, although it does allow selection of standard latitude and
longitude

First, move wrf_hydro.exe from the wrf hydro run folder to this location

Fill in the text on the wrf_preprocess.py file in the top section. No
other changes should be necessary in this file. 

Run wrf_preprocess.py

Check the output in the DOMAIN and FORCING folders. In DOMAIN you should have
4 netcdf files, metadata, fulldom_hires, wrfinput and geo_emd01.nc

Specifically check the FLOW DIRECTION variable in fulldom_hires

In the Forcing folder, there should be many files. Run ncdump [filename]
to verify the forcing data was written properly.

Add the working folder to the path. (export PATH=$PATH:[run directory])

You may wish to make changes to hydro.namelist and namelist.hrldas. Specifically,
the defaults assume that you you are running without subsurface or surface
routing and without any restart files. Look at the docummentation for these
namelists online.

Run mpiexec -n [num_processes] wrf_hydro.exe

Alternatively, run wrf_hydro_d.exe for the diagnostic wrf hydro

To Postprocess the data (improve chunking for i/o, remove unecessary variables)
run (int the main folder) scripts/wrf_postprocess.py

To clean the output AND preprocessing from the folder, run scripts/clean.py
Note that you may need to change the run directory listed in clean.py

KNOWN ISSUES:
	- wrf_regrid.py sometimes fails to properly regrid forcing files
	  no known solution, other than manually editing the reprojected
	  files to remove the abnormal gridcells. Problem will cause 
	  wrf_hydro to crash part way through on or immediately after
	  the date of the improperly regrided file.
	- TauDEM has no reported "nodata" value for the flow direction
	  output (and corresponding variable in Fulldom_hires.nc) but 
	  wrf does have a nodata value. No known problems have occured,
	  but may become an issue

