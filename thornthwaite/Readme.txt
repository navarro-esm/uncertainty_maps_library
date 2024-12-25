This folder contains one program written in NCL.
The workflow is as follows:

Part 01 (thornthwaite.ncl)
	This program takes the netcdf files and performs the classification
	according to Feddema (2005). It then outputs the results in various formats.
	The script is customized to run with different models and scenarios.

	INPUT FILES: netcdf files for precipitation and temperature (e.g.
 	pr_amon_CESM2_historical_r1i1p1_185001-201412.nc
	ts_amon_CESM2_historical_r1i1p1_185001-201412.nc).
	OUTPUT VARIABLES: thorn (all), MI, EP (only netcdf file).
	OUTPUT FILES: netcdf file (x1), bil file (x1), hdr file (x1), plot (x1).

The gdal tools are used to create GeoTIFF files from the output netcdf files.
	gdal_translate -of Gtiff -a_srs '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' NETCDF:"thornthwaite_class_historical_"$MODELNAME".nc":thorn  "thornthwaite_class_historical_"$MODELNAME".tif"
