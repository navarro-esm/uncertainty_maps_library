/*
 * =====================================================================================
 *
 *       Filename:  bil.c
 *
 *    Description:  Generate .bil intermediate files for later classification.
 *
 *        Version:  1.0
 *        Created:  27/07/17 11:49:08
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Raul Moreno,
 *   Organization:  University of Castilla-La Mancha
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libgen.h>
#include <assert.h>
#include <netcdf.h>

#include <sys/stat.h>

#include "log.h"
#include "struct.h"
#include "utils.h"

// Log file
FILE *__flog;

#define VAR_PR 1
#define VAR_TS 2

#define MIN_VAL -9999.f   // Minimum valid data value
#define MAX_VAL 9999.f // Maximum valid data value

// Helper functions
int _file_metadata(char *fpath, file_t *d);
int _set_missing(float *data, const size_t nlat, const size_t nlon, const int variable);
int _generate_bil(const char *input_file, file_t *metadata, const char *output_dir, int variable);

/* *
 * Iterate namelist files and generate monthly .bil files in the output folder.
 * */
int main(int argc, char* argv[])
{
    const char *namelist;
    FILE *fnamelist;

    // File metadata
    char input_file[4096];
    char var_name[16];
    file_t d;

    // Output dir
    char *output_dir;
    char model_dir[4096];

    // Logging
    char log_file[4096];

    // Check arguments
    if(argc < 3)
    {
        LOG_ERROR("Please provide a namelist and output folder.\n");
        return(-1);
    }

    // Get namelist from arguments
    namelist = basename(argv[1]);
    LOG_INFO("Using '%s' namelist.\n", namelist);

    // Get namelist from arguments
    output_dir = argv[2];
    LOG_INFO("Using output folder:%s.\n", output_dir);

    // Open namelist
    fnamelist = fopen(namelist, "r");
    if(!fnamelist)
    {
        LOG_ERROR("Cannot open '%s'\n", namelist);
        return(-1);
    }

    // Create output folder
    // Error check should be implemented
    mkdir(output_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    // Iterate files
    while(fgets(input_file, 4096, fnamelist))
    {
        // Remove newline
        input_file[strcspn(input_file, "\r\n")] = '\0';

        // Get date range from metadata
        LOG_INFO("Processing file '%s'...\n", input_file);
        _file_metadata(input_file, &d); // Read metadata from NetCDF
        if(d.variable == VAR_PR) strcpy(var_name, "VAR_PR");
        if(d.variable == VAR_TS) strcpy(var_name, "VAR_TS");
        LOG_INFO("%s -- %s\n", d.model, var_name);
        LOG_INFO("Date: %d/%02d -> %d/%02d\n", d.year_ini, d.mon_ini, d.year_end, d.mon_end);

        // Create directory for this model
        snprintf(model_dir, 4096, "%s/%s/", output_dir, d.model);
        mkdir(model_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        // Set log
        snprintf(log_file, 4096, "%s/%s.bil.log", output_dir, d.model);

        // Generate .bil files
        if(d.variable == VAR_PR) _generate_bil(input_file, &d, model_dir, VAR_PR); // Precipitation
        if(d.variable == VAR_TS) _generate_bil(input_file, &d, model_dir, VAR_TS); // Temperature
    }

    // Close files
    fclose(fnamelist);

    return(0);
}

/* *
 * Reads metadata from NetCDF file.
 * */
int _file_metadata(char *fpath, file_t *d)
{
    const char *bname;
    char *fname;
    char *ptr_tok;

    // Auxiliar
    int n;

    // Flag for CRU files
    int cru_file = 0;

    // Set output to zero
    memset(d, 0, sizeof(file_t));

    // Get file name
    bname = (char *)basename(fpath);
    fname = (char *)malloc(strlen(bname) + 2);
    strcpy(fname, bname);

    // Split "_"
    ptr_tok = strtok(fname, "_");

    // First element is the variable
    if(strcmp(ptr_tok, "pr") == 0)
        d->variable = VAR_PR;
    else if(strcmp(ptr_tok, "ts") == 0)
        d->variable = VAR_TS;
    else if(strcmp(ptr_tok, "cru") == 0)
        cru_file = 1;

    // CMIP file case
    if(cru_file == 0)
    {
        // We are suppossing the date is in position 5 in the filename
        // Ex: ts_Amon_MPI-ESM-LR_decadal1990_r10i1p1_199101-200012.nc
        n = 0;
        while(ptr_tok != NULL && n < 5)
        {
            // Get model name
            if(n==2) strncpy(d->model, ptr_tok, 31);

            // Next token
            ptr_tok = strtok(NULL, "_");
            n++;
        }

        // Parse date
        if(ptr_tok)
            sscanf(ptr_tok, "%4d%2d-%4d%2d.nc", &d->year_ini, &d->mon_ini, &d->year_end, &d->mon_end);
    }
    // CRU file case
    else
    {
        LOG_INFO("CRU file detected.\n");

        // Next token
        ptr_tok = strtok(NULL, "_");

        // Split "."
        ptr_tok = strtok(ptr_tok, ".");

        // Set model name to CRU
        strncpy(d->model, "CRU", 31);

        // Set months
        d->mon_ini = 1;
        d->mon_end = 12;

        // We are suppossing the date is in position 3 and 4 in the filename
        // Ex: cru_ts4.00.1971.1980.tmp.dat.nc
        n = 0;
        char *ptr;
        while(ptr_tok != NULL && n < 5)
        {
            // Get initial date
            if(n==2) d->year_ini = strtol(ptr_tok, &ptr, 10);

            // Get final date
            if(n==3) d->year_end = strtol(ptr_tok, &ptr, 10);

            // Get variable
            if(n==4)
            {
                if(strcmp(ptr_tok, "pre") == 0)
                    d->variable = VAR_PR;
                else if(strcmp(ptr_tok, "tmp") == 0)
                    d->variable = VAR_TS;
                else
                    LOG_ERROR("Error: Invalid variable name '%s' on file '%s'\n", ptr_tok, fname);
            }

            // Next token
            ptr_tok = strtok(NULL, ".");
            n++;
        }

    }

    // Free memory
    free(fname);

    return 0;
}

/* *
 * Checks for valid float-point values in data and set missing values.
 * */
int _set_missing(float *data, const size_t nlat, const size_t nlon, const int variable)
{
    // Check validity of float values
    int idx = 0;
    int non_zero_found = 0;
    for(int lat = 0; lat < nlat; lat++)
    for(int lon = 0; lon < nlon; lon++)
    {
        switch(fpclassify(data[idx]))
        {
            // Invalid float
            case FP_INFINITE:
                LOG_ERROR("Error in FP: Infinite float number found when retrieving data, pos (%d,%d)\n", lat, lon);
                return(-1);
                break;
            case FP_NAN:
                LOG_ERROR("Error in FP: NaN float number found when retrieving data, pos (%d,%d)\n", lat, lon);
                return(-1);
                break;
            case FP_SUBNORMAL:
                LOG_ERROR("Error in FP: Underflow float number found when retrieving data, pos (%d,%d)\n", lat, lon);
                return(-1);
                break;
            // Valid float
            default:
                // Check bounds
                if(data[idx] < MIN_VAL || data[idx] > MAX_VAL)
                    data[idx] = MISSING_VAL; // Set as missing value
                else if(data[idx] != 0.f) non_zero_found = 1;
                break;
        }

        // Increase index
        idx++;
    }

    // All the data is zero?
    if(!non_zero_found)
    {
        LOG_ERROR("Error: all data is zero!\n");
        return(-1); // Something is wrong if everything is zero
    }

    return(0);
}

/* *
 * Read the contents of a NetCDF input_file and write every month data to a separate .bil file.
 * */
int _generate_bil(const char *input_file, file_t *metadata, const char *output_dir, int variable)
{
    assert(input_file);
    assert(metadata);
    assert(output_dir);

    // NetCDF IDs
    int ncid;
    int lat_id, lon_id;
    int time_id;
    char str_var[16];
    char output_str_var[16];

    // Dimensions
    size_t nlat, nlon;
    size_t ntime;
    int ndims;
    size_t start[] = {0,0,0};
    size_t dims[3];

    // Output file
    FILE *fout;
    char output_file[4096];

    // Error code for NC
    int errcode;

    // Error code for function
    int retcode = 0;

    // Select variable
    if(variable == VAR_TS) strcpy(str_var, "ts");
    else if(variable == VAR_PR) strcpy(str_var, "pr");
    else strcpy(str_var, "VAR_UNDEF");

    // Select output variable name
    if(variable == VAR_TS) strcpy(output_str_var, "ts");
    else if(variable == VAR_PR) strcpy(output_str_var, "pr");
    else strcpy(output_str_var, "VAR_UNDEF");

    // Select variable if CRU file
    if(strcmp(metadata->model, "CRU") == 0)
    {
        if(variable == VAR_TS) strcpy(str_var, "tmp");
        else if(variable == VAR_PR) strcpy(str_var, "pre");
        else strcpy(str_var, "VAR_UNDEF");
    }

    // Open NC file
    LOG_INFO("Opening file '%s'...\n", input_file);
    errcode = nc_open(input_file, NC_NOWRITE, &ncid);
    if(errcode)
    {
        LOG_ERROR("Error opening the file '%s': %s\n", input_file, nc_strerror(errcode));
        return(-1);
    }

    // Read dimensions
    errcode = nc_inq(ncid, &ndims, NULL, NULL, NULL);
    if(errcode)
    {
        LOG_ERROR("Error inquiring the file '%s': %s\n", input_file, nc_strerror(errcode));
        return(-1);
    }
    LOG_INFO("Dimensions: %d\n", ndims);

    // Read time
    errcode = nc_inq_dimid(ncid, "time", &time_id);
    if(errcode)
    {
        LOG_ERROR("Error retrieving time ID in '%s': %s\n", input_file, nc_strerror(errcode));
        return(-1);
    }
    errcode = nc_inq_dim(ncid, time_id, NULL, &ntime);
    if(errcode)
    {
        LOG_ERROR("Error retrieving latitude dimension in '%s': %s\n", input_file, nc_strerror(errcode));
        return(-1);
    }
    LOG_INFO("Time dimension: %lu\n", ntime);

    // Read coordinates
    errcode = nc_inq_dimid(ncid, "lat", &lat_id);
    if(errcode)
    {
        LOG_ERROR("Error retrieving latitude ID in '%s': %s\n", input_file, nc_strerror(errcode));
        return(-1);
    }
    errcode = nc_inq_dim(ncid, lat_id, NULL, &nlat);
    if(errcode)
    {
        LOG_ERROR("Error retrieving latitude dimension in '%s': %s\n", input_file, nc_strerror(errcode));
        return(-1);
    }
    LOG_INFO("Latitude dimension:  %lu\n", nlat);
    errcode = nc_inq_dimid(ncid, "lon", &lon_id);
    if(errcode)
    {
        LOG_ERROR("Error retrieving longitude ID in '%s': %s\n", input_file, nc_strerror(errcode));
        return(-1);
    }
    errcode = nc_inq_dim(ncid, lon_id, NULL, &nlon);
    if(errcode)
    {
        LOG_ERROR("Error retrieving longitude dimension in '%s': %s\n", input_file, nc_strerror(errcode));
        return(-1);
    }
    LOG_INFO("Longitude dimension: %lu\n", nlon);

    // Allocate arrays
    float *f_data = (float *)malloc(ntime * nlat * nlon * sizeof(float));
    assert(f_data);
    float *f_time = (float *)malloc(ntime * sizeof(float));
    assert(f_time);
    float *f_lat = (float *)malloc(nlat * sizeof(float));
    assert(f_lat);
    float *f_lon = (float *)malloc(nlon * sizeof(float));
    assert(f_lon);

    // Setup dimensions for time
    dims[0] = ntime;

    // Read time
    int var_time_id;
    errcode = nc_inq_varid(ncid, "time", &var_time_id);
    if(errcode)
    {
        LOG_ERROR("Error retrieving %s ID in '%s': %s\n", "time", input_file, nc_strerror(errcode));
        free(f_data);
        free(f_time);
        free(f_lat);
        free(f_lon);
        return(-1);
    }
    errcode = nc_get_vara_float(ncid, var_time_id, start, dims, f_time);
    if(errcode)
    {
        LOG_ERROR("Error retrieving %s data in '%s': %s\n", "time", input_file, nc_strerror(errcode));
        free(f_data);
        free(f_time);
        free(f_lat);
        free(f_lon);
        return(-1);
    }

    // Setup dimensions for latitude
    dims[0] = nlat;

    // Read latitude
    int var_lat_id;
    errcode = nc_inq_varid(ncid, "lat", &var_lat_id);
    if(errcode)
    {
        LOG_ERROR("Error retrieving %s ID in '%s': %s\n", "lat", input_file, nc_strerror(errcode));
        free(f_data);
        free(f_time);
        free(f_lat);
        free(f_lon);
        return(-1);
    }
    errcode = nc_get_vara_float(ncid, var_lat_id, start, dims, f_lat);
    if(errcode)
    {
        LOG_ERROR("Error retrieving %s data in '%s': %s\n", "lat", input_file, nc_strerror(errcode));
        free(f_data);
        free(f_time);
        free(f_lat);
        free(f_lon);
        return(-1);
    }

    // Setup dimensions for longitude
    dims[0] = nlon;

    // Read longitude
    int var_lon_id;
    errcode = nc_inq_varid(ncid, "lon", &var_lon_id);
    if(errcode)
    {
        LOG_ERROR("Error retrieving %s ID in '%s': %s\n", "lon", input_file, nc_strerror(errcode));
        free(f_data);
        free(f_time);
        free(f_lat);
        free(f_lon);
        return(-1);
    }
    errcode = nc_get_vara_float(ncid, var_lon_id, start, dims, f_lon);
    if(errcode)
    {
        LOG_ERROR("Error retrieving %s data in '%s': %s\n", "lon", input_file, nc_strerror(errcode));
        free(f_data);
        free(f_time);
        free(f_lat);
        free(f_lon);
        return(-1);
    }

    // Setup dimensions
    dims[0] = ntime;
    dims[1] = nlat;
    dims[2] = nlon;

    // Read data
    int var_id;
    errcode = nc_inq_varid(ncid, str_var, &var_id);
    if(errcode)
    {
        LOG_ERROR("Error retrieving %s ID in '%s': %s\n", str_var, input_file, nc_strerror(errcode));
        free(f_data);
        free(f_time);
        free(f_lat);
        free(f_lon);
        return(-1);
    }
    errcode = nc_get_vara_float(ncid, var_id, start, dims, f_data);
    if(errcode)
    {
        LOG_ERROR("Error retrieving %s data in '%s': %s\n", str_var, input_file, nc_strerror(errcode));
        free(f_data);
        free(f_time);
        free(f_lat);
        free(f_lon);
        return(-1);
    }

    // Write header
    snprintf(output_file, 4096, "%s/header.txt", output_dir);
    fout = fopen(output_file, "w");
    fprintf(fout, "%lu %lu", nlat, nlon);
    fclose(fout);

    // Write coordinates binary file
    snprintf(output_file, 4096, "%s/coordinates.bil", output_dir);
    fout = fopen(output_file, "w");
    fwrite(f_lat, sizeof(float), nlat, fout);
    fwrite(f_lon, sizeof(float), nlon, fout);
    fclose(fout);

    // Iterate months
    int m;
    int current_year = metadata->year_ini;
    int current_mon = metadata->mon_ini;
    int prev_year = metadata->year_ini;
    int prev_mon = metadata->mon_ini;
    for(m = 0; m < ntime; m++)
    {
        // If CMIP model convert units to Celsius and mm/month
        if(strcmp(metadata->model, "CRU") != 0)
        {
            // Convert temperatures to Celsius (Not CRU case)
            if(variable == VAR_TS)
            {
                for(int i = 0; i < nlat * nlon; i++)
                    f_data[m * nlat * nlon + i] = celsius(f_data[m * nlat * nlon + i]);
            }
            else if (variable == VAR_PR)
            {
                // Convert precipitation to mm/month (Not CRU case)
                // Get conversion constant for precipitation
                float sec_in_month = 3600.f * 24.f * 31.f;
                for (int i = 0; i < nlat * nlon; i++)
                    f_data[m * nlat * nlon + i] *= sec_in_month; // 1 kg m-2 s-1 == 1 mm/s ===> mm / month
            }
            else
            {
                LOG_ERROR("Error: Unknown variable while converting units.");
            }
        }

        // Check data validity before generate the .bil and set missing values
        if(_set_missing(&f_data[m * nlat * nlon], nlat, nlon, variable))
        {
            // Something is wrong in this month
            // Skip this corrupt .bil
            LOG_ERROR("Failed to read data from '%s', date: %04d/%02d\n", input_file, current_year, current_mon);
        }
        else
        {
            if(f_time[m] <= 0.f)
            {
                LOG_ERROR("Warning: Invalid time value (%.1f) from '%s', date: %04d/%02d\n", f_time[m], input_file, current_year, current_mon);
            }

            // Write data to .bil file
            snprintf(output_file, 4096, "%s/%s_%04d%02d.bil", output_dir, output_str_var, current_year, current_mon);
            fout = fopen(output_file, "w");
            fwrite(&f_data[m * nlat * nlon], sizeof(float), nlat * nlon, fout);
            fclose(fout);
        }

        // Update previous date
        prev_year = current_year;
        prev_mon = current_mon;

        // Update current date
        current_mon++;
        if(current_mon > 12) { current_mon = 1; current_year++; }
    }

    // Check that the end date match the metadata
    if(prev_year != metadata->year_end || prev_mon != metadata->mon_end)
    {
        // WARNING, something is wrong.
        // The last m should match the last date in the metadata.
        LOG_ERROR("WARNING: The last date in the data does not match the last date in the filename: '%s'\n", input_file);
        retcode = -1;
    }

    // Free arrays
    free(f_data);
    free(f_lat);
    free(f_lon);

    // Close NC file
    LOG_INFO("Closing file '%s'...\n", input_file);
    nc_close(ncid);

    return(retcode);
}
