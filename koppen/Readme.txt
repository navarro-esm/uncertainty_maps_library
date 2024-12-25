This folder contains three programs two written in C and one written in NCL.
The workflow is as follows:

Part 01 (bil.out & classifier.out)
	The programs are written in C and need to be compiled. We have provided
	a makefile. Please, make sure that the LAPACK and NETCDF libraries are installed first.
	Further instructions can be found in the Readme-classifier.txt.

	OUTPUT FILES: intermediate bil files with koppen climate types (e.g.:
	koppen_class_CESM2.bil).

Part 02 (koppen_read-plot.ncl)
	This program takes the intermediate bil files from part01 and then outputs
	the results in various formats. The script is customized to run with
	different models and scenarios.

	INPUT FILES:intermediate bil files for individual models (e.g. koppen_class_CESM2.bil).
	OUTPUT VARIABLES: koppen (all).
	OUTPUT FILES: netcdf file (x1), bil file (x1), hdr file (x1), plot (x1).

The gdal tools are used to create GeoTIFF files from the output netcdf files.
	gdal_translate -of Gtiff -a_srs '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' NETCDF:"holdridge_class_historical_"$MODELNAME".nc":HLZ "holdridge_class_historical_"$MODELNAME".tif"
