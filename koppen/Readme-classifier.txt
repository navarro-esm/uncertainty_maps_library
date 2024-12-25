# Koppen / K-Means classifier

This program is intended for obtaining Koppen and K-Means climate classifications from the processing of CMIP6-based files.

## Description

Input NetCDF files from CMPI6 project are used in this program.
A CMIP6 file follows the next naming convention:
*ts_Amon_MPI-ESM-LR_historical_r1i1p1_185001-201412.nc*.
Where the *ts* attribute means is a temperature file,
the *Amon* attribute means is a montly data file,
the *MPI-ESM-LR* attribute indicates the model which generated this data,
the *historical* attribute means the experiment,
the *r1i1p1* attribute means the ensemble and the *185001-201412* attribute means the range of dates provided in the file.

The input files are converted to an intermediate **.bil** format based on a sequence of float point values.
A *.bil* file can contain precipitation or temperature data with these conditions:

* Only one month of data.
* Only data from one model.
* A grid with the same dimensions as the grid provided by the model is stored.

Also, the files are stored in folders, one folder for every model.
For example, the file *CESM1-CAM5/ts_197005.bil* will contain the monthly temperatures of May 1970 for the CESM1-CAM5 model.

The intermediate *.bil* files can be used for computing a Koppen and K-Means classification from a range of dates, if the right files are present.
If one data file is missing in this range, the model is not classified.
The first step is to iterate the models in the working folder and compute the averages and standard deviations, and then perform a Principal Component Analysis.
Next, the models are again iterated to compute the required grid of averages for a Koppen classification (such as annual precipitation or the coldest temperature in the year).
Afterwards the Koppen, Extended Koppen and K-Means classification are computed in that order.

With the previously computed values, the Koppen and an Extended Koppen (Trewartha) classification is computed and stored in an output file in the output folder.
For example, the file *koppen/koppen_class_ACCESS1-0.bil* will contain the Koppen classes for the model ACCESS1-0 and the file *koppen_ex/koppen_ex_class_ACCESS1-0.bil* will contain the Extended Koppen classes.

After computing the Koppen classes, their centroids are obtained, standardized and used for initializing the K-Means algorithm.
Then the values of the grids are standardized.
The algorithm iterates until the classes converge, which means that the classification is the same as in the previous iteration.
The distances are Euclidean-based and can be used the distance between the climate parameters or the distance between Principal Components.
In the last case a set of weights obtained from the eigenvalues are used for computing the distance.

## Usage

This software is composed by 2 executables:

* **bil.out**: Used for obtaining the intermediate **.bil** files.
* **classifier.out**: Used for obtaining the Koppen and K-Means classification from **.bil** files.

First, a set of input files should be given to the **bil.out** executable through a namelist with the path of every NetCDF file:
*bin/bil.out namelist.input /bil_folder*
The *bil_folder* will be created with a folder inside for every model containing its intermediate *.bil* files.

To classify the data, the classifier needs the *bil_folder* in order to work:
*bin/classifier.out /bil_folder 1960 2004*

The previous command will use the data in the *bil_folder* to perform the classification of the models from 1960 to 2004.
The following folders will be created:

* *bil_folder/koppen*: With the output classification from Koppen.
* *bil_folder/koppen_ex*: With the output classification from Extended Koppen.
* *bil_folder/kmeans*: With the output classification from K-Means algorithm.
* *bil_folder/kmeans_dist*: With the output distances from K-Means algorithm.

One additional parameter can be passed to the classifier:
*bin/classifier.out /bil_folder 2010 2100 /init_bil_folder*
In this case, the centroids from another *init_bil_folder* will be used to initialize the K-Means algorithm.

## Compilation

1. Install LAPACK library in your system if it is not present.
2. Modify the desired parameters in the *Makefile*.
3. Execute make in the classifier source folder.
