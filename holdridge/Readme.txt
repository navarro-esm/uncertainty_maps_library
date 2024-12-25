This folder contains one program written in NCL.
The workflow is as follows:

Part 01 (hlz.ncl)
	This program takes the netcdf files and performs the HLZ classification
	according to Sisneros et al. (2011) and Monserud & Leemans (1992).
	It then outputs the results in various formats.
	The script is customized to run with different models and scenarios.

	INPUT FILES: netcdf files for precipitation and temperature (e.g.
	pr_amon_CESM2_historical_r1i1p1_185001-201412.nc
	ts_amon_CESM2_historical_r1i1p1_185001-201412.nc).
	OUTPUT VARIABLES: HLZ (all), APP, ABT, PER (only netcdf file).
	OUTPUT FILES: netcdf file (x1), bil file (x1), hdr file (x1), plot (x1).

The gdal tools are used to create GeoTIFF files from the output netcdf files.
	gdal_translate -of Gtiff -a_srs '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' NETCDF:"koppen_class_historical_"$MODELNAME".nc":koppen "koppen_class_historical_"$MODELNAME".tif"
