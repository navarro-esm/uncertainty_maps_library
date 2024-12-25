This folder contains three programs: two written in NCL and one in R.
The workflow for these programs is as follows:

Part 01 (whittaker_read.ncl)
	This program creates 2 intermediate files.The script is customized to run
	with different models and scenarios.
		
	INPUT: netcdf files for precipitation and temperature (e.g.
	pr_amon_CESM2_historical_r1i1p1_185001-201412.nc
	ts_amon_CESM2_historical_r1i1p1_185001-201412.nc).
	OUTPUT: intermediate files (pr.csv and ts.csv) to be to be used in part02.

Part 02 (whittaker.R)
	This program takes the output from part01 and performs the classifiying operation.
	Required libraries: devtools, rgeos, plotbiomes (https://github.com/valentinitnelav/plotbiomes)

	INPUT FILES: pr.csv, ts.csv
	OUTPUT FILE: whittaker.csv


Part 03 (whittaker_plot.ncl)
	This program prepares the end user files by taking the output from part02.
	The script is customized to run with different models and scenarios.
   
	INPUT FILES: whittaker.csv
	OUTPUT VARIABLES: whittaker (all), APP, AT (only netcdf file).
	OUTPUT FILES: netcdf file (x1), bil file (x1), hdr file (x1), plot (x1).


The gdal tools are used to create GeoTIFF files from the output netcdf files.
	gdal_translate -of Gtiff -a_srs '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' NETCDF:"whittaker_class_historical_"$MODELNAME".nc":whittaker  "whittaker_class_historical_"$MODELNAME".tif"
