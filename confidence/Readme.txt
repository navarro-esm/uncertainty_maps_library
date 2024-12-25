This folder contains one program written in NCL.
The workflow is as follows:

Part 01 (confidence_intervals.ncl)
	This program performs the uncertainty analysis of the top-10 CMIP6 models by
	calculating two variables: confidence and modvar.
	The program takes bil files of reference data and top-10 models.
	The script is customized to run with different classification schemes and scenarios.
	Before running, a TXT file should be prepared,which include the filename list
	of top 10 models. e.g. "Z-per80.txt"

  INPUT FILES: Bil files of top-10 models (x10). Bil file of reference data.
	OUTPUT VARIABLES: confidence, modvar.
	OUTPUT FILES: netcdf file (x1), bil files (x2), hdr files (x2), plot (x1).

The gdal tools are used to create GeoTIFF files from the output netcdf files.
	gdal_translate -of Gtiff -a_srs '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' NETCDF:$CLASS"_consensus_"$CASE".nc":confidence  $CLASS"_confidence_"$CASE".tif"
	gdal_translate -of Gtiff -a_srs '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' NETCDF:$CLASS"_consensus_"$CASE".nc":modvar  $CLASS"_modvar_"$CASE".tif"
