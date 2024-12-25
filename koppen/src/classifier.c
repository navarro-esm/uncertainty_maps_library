/*
 * =====================================================================================
 *
 *       Filename:  classifier.c
 *
 *    Description:  Classify climate.
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <libgen.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>

#include <sys/stat.h>
#include <dirent.h>

#include "log.h"
#include "struct.h"

#include "koppen.h"
#include "kmeans.h"
#include "pca.h"

/* *
 * K-Means general configuration.
 * */
#define K 28
#define KMEANS_MAX_RAND_IT 50
#define KMEANS_INIT_STRATEGY KMEANS_INIT
//#define KMEANS_INIT_STRATEGY KMEANS_DISTRIBUTED
//#define KMEANS_INIT_STRATEGY KMEANS_RAND_UNIFORM

/* *
 * Helper functions.
 * */
int _model_metadata(const char *model_path, model_t *metadata, int year_ini, int mon_ini);
int _model_get_stats(const char *bil_dir, const char *model_path, const model_t *model_meta, const int initial_year, const int final_year, stats_t *model_table);
int _load_data(const char *data_file, size_t n, float *data);
int _get_models_list(const char *dir, const int initial_year, char ***models, int *n_models);
int _load_kmeans_init_file(const char *filename, stats_t **out_k_centroids, int *out_k);

/* *
 * Process Koppen and K-Means classification of a given model.
 * bil_dir: Directory of the .bil files.
 * model_path: Relative path of the model from bil_dir.
 * model_meta: Metadata of the model, previously loaded.
 * k_classes: How many classes for K-Means classification.
 * k_strategy: Initialization strategy for K-Means (KMEANS_INIT, KMEANS_RAND_UNIFORM, KMEANS_DISTRIBUTED. See kmeans.h.).
 * k_init_centroids: Initial centroids in case of KMEANS_INIT strategy.
 * st_avg: Average values (used for standardization).
 * st_dev: Standard deviation values (used for standardization).
 * weights: Weights used in K-Means distance computation.
 * eigenvectors: Used for computing Principal Components.
 * model_table: Grid with weather values.
 * */
int model_process(const char *bil_dir, const char *model_path, const model_t *model_meta, const int k_classes, const int k_strategy, const stats_t *k_init_centroids, const stats_t *st_avg, const stats_t *st_dev, const double *weights, const double *eigenvectors, stats_t *model_table)
{
    // Error codes
    int errcode;

    // Coordinates
    float *coords;

    // Auxiliar arrays
    char data_file[4096];

    // K-Means
    stats_t *k_centroids = NULL;
    stats_t *koppen_centroids = NULL;
    int *k_count = NULL;
    float k_energy;
    int koppen_k = KOPPEN_EX_N;

    // Output data dir
    char output_dir[4096];
    char output_bil[4096];

    // Output files
    FILE *fheader;
    FILE *fclass;

    // Allocate k_centroids and k_count
    k_centroids = (stats_t *)malloc(k_classes * sizeof(stats_t));
    k_count = (int *)malloc(k_classes * sizeof(int));

    // Allocate coordinates
    coords = (float *)malloc((model_meta->nlat + model_meta->nlon) * sizeof(float));
    if(!coords)
    {
        LOG_ERROR("Coordinate array allocation failed.\n");
        return(-1);
    }

    // Load latitude/longitude data
    snprintf(data_file, 4096, "%s/coordinates.bil", model_path);
    if(_load_data(data_file, model_meta->nlat + model_meta->nlon, coords))
    {
        LOG_ERROR("Failed to read data from '%s'.\n", data_file);
        free(coords);
        return(-1);
    }

    /* *****************************************
     * KOPPEN
     * *****************************************/

    // Compute Koppen
    errcode = koppen(model_table, model_meta->nlat, model_meta->nlon);
    if(errcode == 0)
    {
        // Create output folder
        snprintf(output_dir, 4096, "%s/koppen", bil_dir);
        mkdir(output_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        // Write output header
        snprintf(output_bil, 4096, "%s/koppen/%s.txt", bil_dir, model_meta->name);
        fheader = fopen(output_bil, "w");
        fprintf(fheader, "%lu %lu", model_meta->nlat, model_meta->nlon);
        fclose(fheader);

        // Write output file
        snprintf(output_bil, 4096, "%s/koppen/koppen_class_%s.bil", bil_dir, model_meta->name);
        fclass = fopen(output_bil, "w");
        fwrite(coords, sizeof(float), model_meta->nlat + model_meta->nlon, fclass);
        for(int i = 0; i < model_meta->nlat * model_meta->nlon; i++)
        {
            float class = (float)model_table[i].class;
            fwrite(&class, sizeof(float), 1, fclass);
        }
        fclose(fclass);

        LOG_INFO("%s - Computed Koppen.\n", model_meta->name);
    }
    else
    {
        LOG_ERROR("%s - Error when computing Koppen.\n", model_meta->name);
        return(-1);
    }

    /* *****************************************
     * EXTENDED KOPPEN (Trewartha)
     * *****************************************/

    // Compute Koppen - Trewartha
    errcode = koppen_extended(model_table, model_meta->nlat, model_meta->nlon);
    if(errcode == 0)
    {
        // Create output folder
        snprintf(output_dir, 4096, "%s/koppen_ex", bil_dir);
        mkdir(output_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        // Write output header
        snprintf(output_bil, 4096, "%s/koppen_ex/%s.txt", bil_dir, model_meta->name);
        fheader = fopen(output_bil, "w");
        fprintf(fheader, "%lu %lu", model_meta->nlat, model_meta->nlon);
        fclose(fheader);

        // Write output file with classes
        snprintf(output_bil, 4096, "%s/koppen_ex/koppen_ex_class_%s.bil", bil_dir, model_meta->name);
        fclass = fopen(output_bil, "w");
        fwrite(coords, sizeof(float), model_meta->nlat + model_meta->nlon, fclass);
        for(int i = 0; i < model_meta->nlat * model_meta->nlon; i++)
        {
            float class = (float)model_table[i].class;
            fwrite(&class, sizeof(float), 1, fclass);
        }
        fclose(fclass);

        // Compute centroids
        {
            // Allocate
            koppen_centroids = (stats_t *) malloc(koppen_k * sizeof(stats_t));
            stats_t *__centroids = (stats_t *) malloc(koppen_k * sizeof(stats_t));
            int *__count = (int *) malloc(koppen_k * sizeof(int));

            // Initialize to zero
            memset(koppen_centroids, 0, koppen_k * sizeof(stats_t));
            memset(__count, 0, koppen_k * sizeof(int));

            // Compute centroids using the extended Koppen classification
            errcode = kmeans_centroids(model_table, model_meta->nlat * model_meta->nlon, koppen_k, koppen_centroids, __count);

            // Standardize and compute PC
            standardize_stats(koppen_centroids, k_classes, st_avg, st_dev, __centroids);
            compute_PC(__centroids, k_classes, eigenvectors);

            // Destandardize again
            destandardize_stats(__centroids, k_classes, st_avg, st_dev, koppen_centroids);

            // Save centroids to file
            snprintf(output_bil, 4096, "%s/koppen_ex/%s.centroids", bil_dir, model_meta->name);
            LOG_INFO("Saving Koppen centroids to: %s\n", output_bil);
            errcode = kmeans_centroids_to_file(koppen_centroids, koppen_k, output_bil);

            // Deallocate
            free(__centroids);
            free(__count);
        }

        LOG_INFO("%s - Computed extended Koppen.\n", model_meta->name);
    }
    else
    {
        LOG_ERROR("%s - Error when computing extended Koppen.\n", model_meta->name);
        return(-1);
    }

    /* *****************************************
     * K-MEANS (and distances)
     * *****************************************/

    // Auxiliar variables
    stats_t *aux_k_centroids;
    stats_t *aux_k_init_centroids;
    int *aux_k_count;
    stats_t kmin;
    stats_t kmax;

    // Allocate
    aux_k_init_centroids = (stats_t *)malloc(k_classes * sizeof(stats_t));
    aux_k_centroids = (stats_t *)malloc(k_classes * sizeof(stats_t));
    aux_k_count = (int *)malloc(k_classes * sizeof(stats_t));

    // Standardize data
    // This is a neccessary step for K-means to work.
    stats_t *std_table;
    std_table = (stats_t *)malloc(model_meta->nlat * model_meta->nlon * sizeof(stats_t));
    memset(std_table, 0, model_meta->nlat * model_meta->nlon * sizeof(stats_t));
    standardize_stats(model_table, model_meta->nlat * model_meta->nlon, st_avg, st_dev, std_table);
    compute_PC(std_table, model_meta->nlat * model_meta->nlon, eigenvectors);

    // Standardize initial centroids
    stats_t *std_init_centroids = NULL;
    if(k_init_centroids != NULL)
    {
        std_init_centroids = (stats_t *)malloc(k_classes * sizeof(stats_t));
        standardize_stats(k_init_centroids, k_classes, st_avg, st_dev, std_init_centroids);
        compute_PC(std_init_centroids, k_classes, eigenvectors);
    }
    else
    {
        // Use standardized Koppen centroids
        std_init_centroids = (stats_t *)malloc(k_classes * sizeof(stats_t));
        standardize_stats(koppen_centroids, k_classes, st_avg, st_dev, std_init_centroids);
        compute_PC(std_init_centroids, k_classes, eigenvectors);
        LOG_INFO("Using Koppen centroids for %s\n", model_meta->name);
    }

    // Compute min and max values, needed for distributed or random initialization
    kinit_min_max(std_table, model_meta->nlat * model_meta->nlon, &kmin, &kmax);

    // Compute K-means
    // The loop is for the case of random initialization
    float best_k_energy = FLT_MAX;
    for(int kit = 0; kit < KMEANS_MAX_RAND_IT; kit++)
    {
        // Init centroids
        memset(aux_k_centroids, 0, k_classes * sizeof(stats_t));
        memset(aux_k_count, 0, k_classes * sizeof(int));

        switch(k_strategy)
        {
            case KMEANS_DISTRIBUTED:
#ifdef _DEBUG
                LOG_INFO("KMEANS_DISTRIBUTED\n");
                LOG_INFO("%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
                        "Temp", "Tmin", "Tmax", "Prec", "Pmin", "Pmax", "Psummer", "PSmin", "PSmax", "Pwinter", "PWmin", "PWmax", "Pannual");
                LOG_INFO("Min: %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
                        kmin.t_avg,
                        kmin.t_min,
                        kmin.t_max,
                        kmin.r_avg,
                        kmin.r_min,
                        kmin.r_max,
                        kmin.r_s,
                        kmin.r_smin,
                        kmin.r_smax,
                        kmin.r_w,
                        kmin.r_wmin,
                        kmin.r_wmax,
                        kmin.r_annual
                        );
                LOG_INFO("Max: %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
                        kmax.t_avg,
                        kmax.t_min,
                        kmax.t_max,
                        kmax.r_avg,
                        kmax.r_min,
                        kmax.r_max,
                        kmax.r_s,
                        kmax.r_smin,
                        kmax.r_smax,
                        kmax.r_w,
                        kmax.r_wmin,
                        kmax.r_wmax,
                        kmax.r_annual
                        );
#endif
                // Only do one iteration
                kit = KMEANS_MAX_RAND_IT;
                // Generate initial centroids using distributed strategy
                kinit_distributed(k_classes, &kmin, &kmax, aux_k_init_centroids);
                // Update PC
                compute_PC(aux_k_init_centroids, k_classes, eigenvectors);
                // Compute kmeans
                k_energy = kmeans(std_table, k_classes,
                        model_meta->nlat, model_meta->nlon,
                        aux_k_centroids, aux_k_count,
                        aux_k_init_centroids,
                        eigenvectors, weights);
                break;

            case KMEANS_RAND_UNIFORM:
#ifdef _DEBUG
                LOG_INFO("KMEANS_RAND_UNIFORM\n");
#endif
                // Generate initial using random strategy
                kinit_random(k_classes, &kmin, &kmax, aux_k_init_centroids);
                compute_PC(aux_k_init_centroids, k_classes, eigenvectors);
                // Compute kmeans
                k_energy = kmeans(std_table, k_classes,
                        model_meta->nlat, model_meta->nlon,
                        aux_k_centroids, aux_k_count,
                        aux_k_init_centroids,
                        eigenvectors, weights);
                break;

            case KMEANS_INIT:
#ifdef _DEBUG
                LOG_INFO("KMEANS_INIT\n");

                // Allocate auxiliar centroids and get them
                stats_t *__k_centroids = (stats_t *)malloc(k_classes * sizeof(stats_t));
                destandardize_stats(std_init_centroids, k_classes, st_avg, st_dev, __k_centroids);

                LOG_INFO("Initial centroids:\n");
                LOG_INFO("%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
                        "Class", "Temp", "Tmin", "Tmax", "Prec", "Pmin", "Pmax", "Psummer", "PSmin", "PSmax", "Pwinter", "PWmin", "PWmax", "Pannual");
                for(int i = 0; i < k_classes; i++)
                {
                    LOG_INFO("%12d %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
                            __k_centroids[i].class,
                            __k_centroids[i].t_avg,
                            __k_centroids[i].t_min,
                            __k_centroids[i].t_max,
                            __k_centroids[i].r_avg,
                            __k_centroids[i].r_min,
                            __k_centroids[i].r_max,
                            __k_centroids[i].r_s,
                            __k_centroids[i].r_smin,
                            __k_centroids[i].r_smax,
                            __k_centroids[i].r_w,
                            __k_centroids[i].r_wmin,
                            __k_centroids[i].r_wmax,
                            __k_centroids[i].r_annual
                          );
                }
                free(__k_centroids);
#endif
                // Only do one iteration
                kit = KMEANS_MAX_RAND_IT;
                // Compute kmeans (Koppen centroids are used)
                k_energy = kmeans(std_table, k_classes,
                        model_meta->nlat, model_meta->nlon,
                        aux_k_centroids, aux_k_count,
                        std_init_centroids,
                        eigenvectors, weights);
                break;
            default:
                LOG_ERROR("Unknown Kmeans strategy: %d\n", k_strategy);
                return(-1);
        }

        // Check result
        // Error if energy (error) is negative
        if(k_energy >= 0.f)
        {
#ifdef _DEBUG
            LOG_INFO("K-means last energy: %8.4f\n", k_energy);
#endif
            // Better?
            if(k_energy < best_k_energy)
            {
                // Update best energy
                best_k_energy = k_energy;

                // Get destandardized centroids and copy count array
                destandardize_stats(aux_k_centroids, k_classes, st_avg, st_dev, k_centroids);
                memcpy(k_count, aux_k_count, k_classes * sizeof(int));

#ifdef TRANSFORM_KMEANS_TO_KOPPEN
                // Transform classes to Koppen
                {
                    koppen_extended(k_centroids, k_classes, 1);

                    // Copy classes to auxiliar centroids
                    for(int i = 0; i < k_classes; i++)
                        aux_k_centroids[i].class = k_centroids[i].class;

                    // Apply Koppen classes to table again
                    kclassify(model_meta->nlat * model_meta->nlon, aux_k_centroids, k_classes, weights, std_table);
                }
#endif

#ifdef _DEBUG
                LOG_INFO("K-means best energy: %8.4f\n", k_energy);
                LOG_INFO("%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
                        "Class", "Temp", "Tmin", "Tmax", "Prec", "Pmin", "Pmax", "Psummer", "PSmin", "PSmax", "Pwinter", "PWmin", "PWmax", "Pannual");
                for(int i = 0; i < k_classes; i++)
                {
                    LOG_INFO("%12d %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
                            k_centroids[i].class,
                            k_centroids[i].t_avg,
                            k_centroids[i].t_min,
                            k_centroids[i].t_max,
                            k_centroids[i].r_avg,
                            k_centroids[i].r_min,
                            k_centroids[i].r_max,
                            k_centroids[i].r_s,
                            k_centroids[i].r_smin,
                            k_centroids[i].r_smax,
                            k_centroids[i].r_w,
                            k_centroids[i].r_wmin,
                            k_centroids[i].r_wmax,
                            k_centroids[i].r_annual
                          );
                }
#endif

                // Create output folder
                snprintf(output_dir, 4096, "%s/kmeans", bil_dir);
                mkdir(output_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

                // Same for distances
                snprintf(output_dir, 4096, "%s/kmeans_dist", bil_dir);
                mkdir(output_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

                // Write output header
                snprintf(output_bil, 4096, "%s/kmeans/%s.txt", bil_dir, model_meta->name);
                fheader = fopen(output_bil, "w");
                fprintf(fheader, "%lu %lu\n", model_meta->nlat, model_meta->nlon); // Dimensions
                fclose(fheader);

                // Write output centroids
                snprintf(output_bil, 4096, "%s/kmeans/%s.centroids", bil_dir, model_meta->name);
                LOG_INFO("Saving Kmeans centroids to: %s\n", output_bil);
                errcode = kmeans_centroids_to_file(k_centroids, k_classes, output_bil);

                // Write output class file
                snprintf(output_bil, 4096, "%s/kmeans/kmeans_class_%s.bil", bil_dir, model_meta->name);
                fclass = fopen(output_bil, "w");
                fwrite(coords, sizeof(float), model_meta->nlat + model_meta->nlon, fclass);
                for(int i = 0; i < model_meta->nlat * model_meta->nlon; i++)
                {
                    float class = (float)std_table[i].class;
                    fwrite(&class, sizeof(float), 1, fclass);
                }
                fclose(fclass);
                snprintf(output_bil, 4096, "%s/kmeans/kmeans_class_%s.class", bil_dir, model_meta->name);
                fclass = fopen(output_bil, "w");
                for(int i = 0; i < model_meta->nlat; i++)
                {
                    for(int j = 0; j < model_meta->nlon; j++)
                    {
                        fprintf(fclass, "%02d ", std_table[i * model_meta->nlon + j].class);
                    }
                    fprintf(fclass, "\n");
                }
                fclose(fclass);

                // Write output header for distances dimensions
                snprintf(output_bil, 4096, "%s/kmeans_dist/%s.txt", bil_dir, model_meta->name);
                fheader = fopen(output_bil, "w");
                fprintf(fheader, "%lu %lu\n", model_meta->nlat, model_meta->nlon); // Dimensions
                fclose(fheader);

                // Get maximum distance
                float __max_dist = 0.f;
                for(int i = 0; i < model_meta->nlat * model_meta->nlon; i++)
                    if(__max_dist < std_table[i].class_dist) __max_dist = std_table[i].class_dist;

                // Bound distances [0,1]
                for(int i = 0; i < model_meta->nlat * model_meta->nlon; i++)
                    std_table[i].class_dist /= __max_dist;

                // Write output distance to file
                snprintf(output_bil, 4096, "%s/kmeans_dist/kmeans_dist_%s.bil", bil_dir, model_meta->name);
                fclass = fopen(output_bil, "w");
                fwrite(coords, sizeof(float), model_meta->nlat + model_meta->nlon, fclass);
                for(int i = 0; i < model_meta->nlat * model_meta->nlon; i++)
                    fwrite(&std_table[i].class_dist, sizeof(float), 1, fclass);
                fclose(fclass);
            }
        }
        else
        {
            LOG_ERROR("%s - Error when computing K-means.\n", model_meta->name);
            return(-1);
        }

    }

    // Free memory
    free(std_table);
    if(std_init_centroids != NULL) free(std_init_centroids);
    free(aux_k_init_centroids);
    free(aux_k_centroids);
    free(aux_k_count);
    free(k_centroids);
    free(k_count);

    return(0);
}


/* *
 * Iterate monthly .bil files and generate Koppen and K-means classification.
 * */
int main(int argc, char* argv[])
{
    // Intermediate directory
    const char *bil_dir;
    struct stat sb;
    DIR *bil_fd;

    // K-means
    unsigned int k_classes = K;
    stats_t *k_init_centroids = NULL;
    char *k_init_folder = NULL;

    // Averages and covariance
    double *avg = NULL;
    double *dev = NULL;
    double *cov = NULL;
    double *weights = NULL;
    stats_t *all_table = NULL;
    size_t n_all_table = 0;
    stats_t st_avg;
    stats_t st_dev;

    // PC
    double *eigenvalues;
    double *eigenvectors;

    // Models directories
    char **model_dir;
    int n_models;

    // Error code
    int errcode;

    // Years
    int initial_year;
    int final_year;

    // Model metadata
    model_t model_meta;

    // Init random
    srand(time(NULL));

    // Check arguments
    if(argc < 4)
    {
        LOG_ERROR("Please provide a .bil directory and start/end year.\n");
        return(-1);
    }

    // Get directory from arguments
    bil_dir = argv[1];
    LOG_INFO("Using data in '%s'.\n", bil_dir);

    // Check if it is a valid directory
    if(stat(bil_dir, &sb) != 0 || !S_ISDIR(sb.st_mode))
    {
        LOG_ERROR("Invalid directory!\n");
        return(-1);
    }

    // Info
#ifdef _OPENMP
    LOG_INFO("Using OpenMP.\n");
#endif
#ifdef KMEANS_USE_PC
    LOG_INFO("Using Principal Components for K-Means distances.\n");
#endif

    // Get years
    initial_year = strtol(argv[2], NULL, 10);
    final_year = strtol(argv[3], NULL, 10);
    LOG_INFO("From %d to %d.\n", initial_year, final_year);

    // List all the models and their date range
    _get_models_list(bil_dir, initial_year, &model_dir, &n_models);

    // Check if argument 4 (RCP centroids) is present
    if(argc > 4)
    {
        // RCP folder
        // Use centroids from files
        k_init_folder = argv[5];
    }

    // Open dir
    if((bil_fd = opendir(bil_dir)) == NULL)
    {
        LOG_ERROR("Error accessing %s\n", bil_dir);
        return(-1);
    }

#ifdef _OPENMP
#pragma omp parallel \
    firstprivate(model_dir, n_models, initial_year, final_year, k_classes, k_init_centroids) \
    private(errcode, model_meta) \
    shared(bil_fd, bil_dir, k_init_folder, eigenvalues, eigenvectors, weights, all_table, n_all_table, avg, dev, cov, st_avg, st_dev, stdout, stderr) \
    default(none)
#endif
    {

#ifdef _OPENMP
#pragma omp master
#endif
        {
            /* *
             * Compute Principal Components
             * */

            // Allocate averages, covariance matrix and eigenvalues
            avg = (double *)malloc(STATS_T_N * sizeof(stats_t));
            dev = (double *)malloc(STATS_T_N * sizeof(stats_t));
            cov = (double *)malloc(STATS_T_N * STATS_T_N * sizeof(stats_t));
            eigenvalues = (double *)malloc(STATS_T_N * sizeof(stats_t));
            eigenvectors = (double *)malloc(STATS_T_N * STATS_T_N * sizeof(stats_t));
            weights = (double *)malloc(STATS_T_N * sizeof(stats_t));

        } // OMP master

            // Iterate models
#ifdef _OPENMP
#pragma omp barrier
#pragma omp for schedule(dynamic)
#endif
            for(int m = 0; m < n_models; m++)
            {
                // Get model metadata
                if(_model_metadata(model_dir[m], &model_meta, initial_year, 1) != 0)
                {
                    LOG_ERROR("Model: %s, INVALID DATA\n", model_meta.name);
                }
                else
                {
                    //if(strcmp(model_meta.name, "CRU") == 0)
                    {
                    // Check if the model is in range
                    if(initial_year < model_meta.year_ini || final_year > model_meta.year_end)
                    {
                        // Not enough data of this model
                        LOG_INFO("'%s' ---> Not enough data for this range.\n", model_meta.name);
                    }
                    else
                    {
                        LOG_INFO("Loading data ---> %s\n", model_meta.name);

                        // Allocate model grid
                        stats_t *model_table = (stats_t *)malloc(model_meta.nlat * model_meta.nlon * sizeof(stats_t));
                        memset(model_table, 0, model_meta.nlat * model_meta.nlon * sizeof(stats_t));

                        // Get model table
                        errcode = _model_get_stats(bil_dir, model_dir[m], &model_meta, initial_year, final_year, model_table);

                        // Append values
#ifdef _OPENMP
#pragma omp critical
#endif
                        {
                            all_table = realloc(all_table, (n_all_table + (model_meta.nlat * model_meta.nlon)) * sizeof(stats_t));
                            memcpy(&all_table[n_all_table], model_table, model_meta.nlat * model_meta.nlon * sizeof(stats_t));
                            n_all_table += model_meta.nlat * model_meta.nlon;
                        }

                        // Free memory
                        free(model_table);
                    }
                    }
                }
            }

#ifdef _OPENMP
#pragma omp master
#endif
            {

            // Allocate values
            double *raw_values = (double *) malloc(n_all_table * STATS_T_N * sizeof(double));
            size_t n_rows;

            // Transform models to array
            errcode = stats_to_array(all_table, n_all_table, raw_values, &n_rows);

            LOG_INFO("Getting averages and standard deviation...\n");
            // Compute standardized values
            {
                // Compute averages
                errcode = compute_avg(raw_values, n_rows, STATS_T_N, avg);

                // Compute standard deviations
                errcode = compute_dev(raw_values, n_rows, STATS_T_N, avg, dev);

                // Set st_avg
                st_avg.t_avg    = avg[0 ];
                st_avg.t_min    = avg[1 ];
                st_avg.t_max    = avg[2 ];
                st_avg.r_avg    = avg[3 ];
                st_avg.r_min    = avg[4 ];
                st_avg.r_max    = avg[5 ];
                st_avg.r_s      = avg[6 ];
                st_avg.r_smin   = avg[7 ];
                st_avg.r_smax   = avg[8 ];
                st_avg.r_w      = avg[9 ];
                st_avg.r_wmin   = avg[10];
                st_avg.r_wmax   = avg[11];
                st_avg.r_annual = avg[12];

                // Set st_dev
                st_dev.t_avg    = dev[0 ];
                st_dev.t_min    = dev[1 ];
                st_dev.t_max    = dev[2 ];
                st_dev.r_avg    = dev[3 ];
                st_dev.r_min    = dev[4 ];
                st_dev.r_max    = dev[5 ];
                st_dev.r_s      = dev[6 ];
                st_dev.r_smin   = dev[7 ];
                st_dev.r_smax   = dev[8 ];
                st_dev.r_w      = dev[9 ];
                st_dev.r_wmin   = dev[10];
                st_dev.r_wmax   = dev[11];
                st_dev.r_annual = dev[12];
            }

            // Summary (avg)
            LOG_INFO("Averaging done (RAW, %lu count):\n", n_rows);
            for(int i = 0; i < STATS_T_N; i++)
                LOG_INFO("%12.2lf ", avg[i]);
            LOG_INFO("\n");
            LOG_INFO("Standard deviations (RAW):\n");
            for(int i = 0; i < STATS_T_N; i++)
                LOG_INFO("%12.2lf ", dev[i]);
            LOG_INFO("\n");

            // Allocate standardized table
            stats_t *std_table = malloc(n_all_table * sizeof(stats_t));
            double *std_values = (double *) malloc(n_all_table * STATS_T_N * sizeof(double));

            LOG_INFO("Computing eigenvalues and eigenvectors...\n");
            // Compute Principal Components over standarized values
            {
                // Standardize table
                standardize_stats(all_table, n_all_table, &st_avg, &st_dev, std_table);

                // Obtain standardized values array
                errcode = stats_to_array(std_table, n_all_table, std_values, &n_rows);

                // Compute averages
                errcode = compute_avg(std_values, n_rows, STATS_T_N, avg);

                // Compute covariance
                errcode = compute_cov(std_values, n_rows, STATS_T_N, avg, cov);

                // Compute eigen
                errcode = compute_eigen(cov, STATS_T_N, eigenvalues, eigenvectors);
            }

            // Summary (avg)
            LOG_INFO("Averaging done (STANDARDIZED, %lu count):\n", n_rows);
            for(int i = 0; i < STATS_T_N; i++)
                LOG_INFO("%12.2lf ", avg[i]);
            LOG_INFO("\n");
            LOG_INFO("Standard deviations (STANDARDIZED):\n");
            for(int i = 0; i < STATS_T_N; i++)
                LOG_INFO("%12.2lf ", sqrt(cov[(i*STATS_T_N) + i]));
            LOG_INFO("\n");

            LOG_INFO("Covariance done (%lu count):\n", n_rows);
            for(int i = 0; i < STATS_T_N; i++)
            {
                for(int j = 0; j < STATS_T_N; j++)
                    LOG_INFO("%12.2lf ", cov[(i*STATS_T_N) + j]);
                LOG_INFO("\n");
            }

            // Summary (eigen)
            LOG_INFO("Eigenvalues:\n");
            float __total = 0.f;
            float __acum = 0.f;
            for(int i = 0; i < STATS_T_N; i++)
                __total += eigenvalues[i];
            for(int i = 0; i < STATS_T_N; i++)
            {
                LOG_INFO("%12.2f ", eigenvalues[i]);
            }
            LOG_INFO("\nEigenvalues %%:\n");
            for(int i = 0; i < STATS_T_N; i++)
            {
                LOG_INFO("%6.2f ", eigenvalues[i] / __total);
            }
            LOG_INFO("\nEigenvalues accumulated %%:\n");
            for(int i = 0; i < STATS_T_N; i++)
            {
                __acum += eigenvalues[i] / __total;
                LOG_INFO("%6.2f ", __acum);
            }
            LOG_INFO("\nEigenvectors:\n");
            for(int i = 0; i < STATS_T_N; i++)
            {
                for(int j = 0; j < STATS_T_N; j++)
                    LOG_INFO("%6.2f ", eigenvectors[i * STATS_T_N + j]);
                LOG_INFO("\n");
            }

            // Set weights
            for(int i = 0; i < STATS_T_N; i++)
                weights[i] = eigenvalues[i] / __total;

            // Free memory
            free(all_table);
            free(std_table);
            free(raw_values);
            free(std_values);
            free(avg);
            free(dev);
            free(cov);
            avg = NULL;
            dev = NULL;
            cov = NULL;
        } // OMP Master

    /* *
     * Process models with K-means and Koppen
     * */
#ifdef _OPENMP
#pragma omp barrier
#pragma omp for schedule(dynamic)
#endif
    for(int m = 0; m < n_models; m++)
    {
        // Get model metadata
        if(_model_metadata(model_dir[m], &model_meta, initial_year, 1) != 0)
        {
            LOG_ERROR("Model: %s, INVALID DATA\n", model_meta.name);
        }
        else
        {
#ifdef _PROCESS_ONLY
            if(strcmp(model_meta.name, _PROCESS_ONLY) == 0)
#endif
            {
                // Check if the model is in range
                if(initial_year < model_meta.year_ini || final_year > model_meta.year_end)
                {
                    // Not enough data of this model
                    LOG_INFO("Not enough data for processing '%s'.\n", model_meta.name);
                }
                else
                {
                    // Allocate model grid
                    stats_t *model_table = (stats_t *)malloc(model_meta.nlat * model_meta.nlon * sizeof(stats_t));
                    memset(model_table, 0, model_meta.nlat * model_meta.nlon * sizeof(stats_t));

                    // Get model table
                    errcode = _model_get_stats(bil_dir, model_dir[m], &model_meta, initial_year, final_year, model_table);

                    // Set Kmeans strategy to default
                    int k_strategy = KMEANS_INIT_STRATEGY;

                    // Selected initial centroids
                    stats_t *sel_initial_centroids = NULL;

                    // Try to load centroids file for this model
                    if(k_init_folder != NULL)
                    {
                        int readed_k;
                        char k_init_file[4096];
                        snprintf(k_init_file, 4096, "%s/kmeans/%s.centroids", k_init_folder, model_meta.name);
                        // Try to read the file
                        if(_load_kmeans_init_file(k_init_file, &k_init_centroids, &readed_k) == 0)
                        {
                            // Change strategy
                            k_strategy = KMEANS_INIT;
                            LOG_INFO("Kmeans init file (%d classes): %s\n", k_classes, k_init_file);

                            // Compute PC of readed centroids
                            if(k_init_centroids != NULL)
                            {
                                compute_PC(k_init_centroids, k_classes, eigenvectors);
                                sel_initial_centroids = k_init_centroids;
                            }
                        }
                    }

                    // Use or not the weights for PC?
                    double *__weights = NULL;
#ifdef KMEANS_USE_PC
                    __weights = weights;
#endif

                    // Process model
                    errcode = model_process(bil_dir, model_dir[m], &model_meta, k_classes, k_strategy, sel_initial_centroids, &st_avg, &st_dev, __weights, eigenvectors, model_table);
                    if(errcode == 0)
                    {
                        LOG_INFO("Model %s finished.\n", model_meta.name);
                    }
                    else
                    {
                        LOG_ERROR("Model %s FAILED!.\n", model_meta.name);
                    }

                    // Deallocate
                    free(model_table);
                    if(k_init_centroids != NULL) free(k_init_centroids);
                }
            }
        }
    }

    } // OMP End

    // Close dir
    closedir(bil_fd);

    return(0);
}

/* *
 * Read model metadata from file.
 * The range of dates will be obtained from year_ini and mon_ini and will be contigous.
 * If some month is missing, the final date is set to the previous month.
 * */
int _model_metadata(const char *model_path, model_t *metadata, int year_ini, int mon_ini)
{
    assert(model_path);
    assert(metadata);

    // Basename
    const char *bname;

    // Header file
    FILE *fheader;

    // File name
    char fname[4096];
    struct stat sb;

    // Set output to zero
    memset(metadata, 0, sizeof(model_t));

    // Get file name
    strncpy(fname, model_path, 4096);
    bname = (char *)basename(fname);
    strncpy(metadata->name, bname, 128);

    // Read header
    snprintf(fname, 4096, "%s/header.txt", model_path);
    fheader = fopen(fname, "r");
    if(!fheader)
    {
        LOG_ERROR("Failed to open '%s'.\n", fname);
        return(-1);
    }

    // Read dimensions
    fscanf(fheader, "%lu %lu", &metadata->nlat, &metadata->nlon);

    // Close file
    fclose(fheader);

    // Check contiguous range from year_ini
    int current_year = year_ini;
    int current_mon = mon_ini;
    int last_year = current_year;
    int last_mon = current_mon;
    while(1)
    {
        // Check continuity
        // Precipitation file exists for this month?
        snprintf(fname, 4096, "%s/pr_%04d%02d.bil", model_path, current_year, current_mon);
        if(stat(fname, &sb) != 0 || !S_ISREG(sb.st_mode))
        {
            // Continuity is broken here
            metadata->year_ini = year_ini;
            metadata->mon_ini = 1;
            metadata->year_end = last_year;
            metadata->mon_end = last_mon;
            break;
        }
        // Temperature file exists for this month?
        snprintf(fname, 4096, "%s/ts_%04d%02d.bil", model_path, current_year, current_mon);
        if(stat(fname, &sb) != 0 || !S_ISREG(sb.st_mode))
        {
            // Continuity is broken here
            metadata->year_ini = year_ini;
            metadata->mon_ini = 1;
            metadata->year_end = last_year;
            metadata->mon_end = last_mon;
            break;
        }

        // Update last date
        last_mon = current_mon;
        last_year = current_year;

        // Update current date
        current_mon++;
        if(current_mon > 12) { current_mon = 1; current_year++; }
    }

    return(0);
}

/* *
 * Gets the grid with the montly data.
 * Only data between initial year and final year is going to be collected.
 * */
int _model_get_stats(const char *bil_dir, const char *model_path, const model_t *model_meta, const int initial_year, const int final_year, stats_t *model_table)
{
    // Current model dir
    DIR *model_fd;

    // Auxiliar arrays
    char data_file[4096];

    // Global table
    // Dimensions (lat, lon)
    // All the monthly data is combined here and will be copied to the output model_table
    stats_t *annual_table;

    // Monthly tables
    // Dimensions (12, lat, lon)
    // Every first dimension is the average for the corresponding month
    // Only the averages of precipitation and temperature is used
    stats_t *mon_table;

    // Open model dir
    if((model_fd = opendir(model_path)) == NULL)
    {
        LOG_ERROR("Error accessing %s\n", model_path);
        return(-1);
    }

    // Allocate annual table
    annual_table = (stats_t *)malloc(model_meta->nlat * model_meta->nlon * sizeof(stats_t));
    if(!annual_table)
    {
        LOG_ERROR("Stats data allocation failed.\n");
        return(-1);
    }
    memset(annual_table, 0, model_meta->nlat * model_meta->nlon * sizeof(stats_t));

    // Allocate monthly tables
    mon_table = (stats_t *)malloc(12 * model_meta->nlat * model_meta->nlon * sizeof(stats_t));
    if(!mon_table)
    {
        LOG_ERROR("Stats data allocation failed.\n");
        return(-1);
    }
    memset(mon_table, 0, 12 * model_meta->nlat * model_meta->nlon * sizeof(stats_t));

    // Initialize tables
    for(int i = 0; i < 12 * model_meta->nlat * model_meta->nlon; i++)
    {
        // Averages
        mon_table[i].t_avg = 0.f;
        mon_table[i].r_avg = 0.f;
        mon_table[i].n = 0;
    }

    // Annual table min/max
    for(int i = 0; i < model_meta->nlat * model_meta->nlon; i++)
    {
        annual_table[i].t_avg    = 0.f;
        annual_table[i].r_avg    = 0.f;
        annual_table[i].r_s      = 0.f;
        annual_table[i].r_w      = 0.f;
        annual_table[i].r_annual = 0.f;

        annual_table[i].t_max    = -9999.f;
        annual_table[i].r_max    = -9999.f;
        annual_table[i].r_smax   = -9999.f;
        annual_table[i].r_wmax   = -9999.f;

        annual_table[i].t_min    =  9999.f;
        annual_table[i].r_min    =  9999.f;
        annual_table[i].r_smin   =  9999.f;
        annual_table[i].r_wmin   =  9999.f;

        annual_table[i].n = 0;
        annual_table[i].PC1 = 0.f;
        annual_table[i].PC2 = 0.f;
        annual_table[i].PC3 = 0.f;
        annual_table[i].PC4 = 0.f;
        annual_table[i].PC5 = 0.f;
        annual_table[i].PC6 = 0.f;
        annual_table[i].PC7 = 0.f;
        annual_table[i].PC8 = 0.f;
        annual_table[i].PC9 = 0.f;
        annual_table[i].PC10 = 0.f;
        annual_table[i].PC11 = 0.f;
        annual_table[i].PC12 = 0.f;
        annual_table[i].PC13 = 0.f;
    }

    // Allocate auxiliar grid
    float *pr = (float *)malloc(model_meta->nlat * model_meta->nlon * sizeof(float));
    float *ts = (float *)malloc(model_meta->nlat * model_meta->nlon * sizeof(float));
    if(!ts || !pr)
    {
        LOG_ERROR("Intermediate grid allocation failed.\n");
        return(-1);
    }

    // We need to accumulate montly data for every different month from temperature and precipitation
    // Iterate months
    int model_failed = 0;
    for(int m = 0; m < 12; m++)
    {
        // Iterate years
        for(int y = initial_year; y <= final_year; y++)
        {
            // Load precipitation data
            snprintf(data_file, 4096, "%s/pr_%04d%02d.bil", model_path, y, m+1);
            if(_load_data(data_file, model_meta->nlat * model_meta->nlon, pr))
            {
                LOG_ERROR("Failed to read data from '%s'.\n", data_file);
                model_failed = 1;
                break;
            }

            // Load temperature data
            snprintf(data_file, 4096, "%s/ts_%04d%02d.bil", model_path, y, m+1);
            if(_load_data(data_file, model_meta->nlat * model_meta->nlon, ts))
            {
                LOG_ERROR("Failed to read data from '%s'.\n", data_file);
                model_failed = 1;
                break;
            }

            // Get sums and count for later averaging
            stats_t *current_mon = &mon_table[m * model_meta->nlat * model_meta->nlon]; // Get the grid/table for the current month
            for(int i = 0; i < model_meta->nlat * model_meta->nlon; i++)
            {
                // Check for missing values and skip if neccessary
                if(ts[i] <= MISSING_VAL || pr[i] <= MISSING_VAL) continue;

                // Add precipitation and temperature (averaged later)
                current_mon[i].t_avg += ts[i];
                current_mon[i].r_avg += pr[i];

                // Increase counter
                current_mon[i].n++;
            }
        }

        // Do not continue if this model failed.
        if(model_failed) break;
    }

    // Check for errors
    if(model_failed)
    {
        LOG_ERROR("Error when processing '%s'\n", model_meta->name);

        // Free memory
        free(annual_table);
        free(mon_table);
        free(pr);
        free(ts);

        // Close model dir
        closedir(model_fd);

        return(-1);
    }

    /* *************************************************
     * COMBINE MONTHLY DATA INTO ONE TABLE
     * *************************************************/

    // Iterate months and combine
    for(int m = 0; m < 12; m++)
    {
        // Select month
        stats_t *current_mon = &mon_table[m * model_meta->nlat * model_meta->nlon];

        // Iterate grid
        for(int lat = 0; lat < model_meta->nlat; lat++)
        for(int lon = 0; lon < model_meta->nlon; lon++)
        {
            // Index in the grid
            int i = lat * model_meta->nlon + lon;

            // Missing data?
            if(current_mon[i].n == 0 || current_mon[i].t_avg <= MISSING_VAL) continue;

            // Increase n (so we can know if the stats are valid)
            annual_table[i].n += current_mon[i].n;

            // Now we can compute montly averages
            current_mon[i].t_avg /= current_mon[i].n;
            current_mon[i].r_avg /= current_mon[i].n;

            // Compute partial averages for the final table (12 months)
            annual_table[i].t_avg += current_mon[i].t_avg / 12.f;
            annual_table[i].r_avg += current_mon[i].r_avg / 12.f;

            // Annual precipitation accumulation
            annual_table[i].r_annual += current_mon[i].r_avg;

            // Warmest/coldest month temperature in annual table
            if(annual_table[i].t_max < current_mon[i].t_avg) annual_table[i].t_max = current_mon[i].t_avg; // Warmest average
            if(annual_table[i].t_min > current_mon[i].t_avg) annual_table[i].t_min = current_mon[i].t_avg; // Coldest average

            // Driest/wettest month precipitation in annual table
            if(annual_table[i].r_max < current_mon[i].r_avg) annual_table[i].r_max = current_mon[i].r_avg; // Wettest average
            if(annual_table[i].r_min > current_mon[i].r_avg) annual_table[i].r_min = current_mon[i].r_avg; // Driest average

            // North hemisphere case
            if(lat < (model_meta->nlat / 2))
            {
                // Summer (Jun - Aug)
                if(m >= 5 && m <= 7)
                {
                    annual_table[i].r_s += current_mon[i].r_avg; // Summer precipitation
                    if(annual_table[i].r_smax < current_mon[i].r_avg) annual_table[i].r_smax = current_mon[i].r_avg; // Wettest summer precipitation
                    if(annual_table[i].r_smin > current_mon[i].r_avg) annual_table[i].r_smin = current_mon[i].r_avg; // Driest summer precipitation
                }
                // Winter (Dec - Feb)
                if(m == 11 || m == 0 || m == 1)
                {
                    annual_table[i].r_w += current_mon[i].r_avg; // Winter precipitation
                    if(annual_table[i].r_wmax < current_mon[i].r_avg) annual_table[i].r_wmax = current_mon[i].r_avg; // Wettest winter precipitation
                    if(annual_table[i].r_wmin > current_mon[i].r_avg) annual_table[i].r_wmin = current_mon[i].r_avg; // Driest winter precipitation
                }
            }
            // South hemisphere case
            else
            {
                // Winter (Jun - Aug)
                if(m >= 5 && m <= 7)
                {
                    annual_table[i].r_w += current_mon[i].r_avg; // Winter precipitation
                    if(annual_table[i].r_wmax < current_mon[i].r_avg) annual_table[i].r_wmax = current_mon[i].r_avg; // Wettest winter precipitation
                    if(annual_table[i].r_wmin > current_mon[i].r_avg) annual_table[i].r_wmin = current_mon[i].r_avg; // Driest winter precipitation
                }
                // Summer (Dec - Feb)
                if(m == 11 || m == 0 || m == 1)
                {
                    annual_table[i].r_s += current_mon[i].r_avg; // Summer precipitation
                    if(annual_table[i].r_smax < current_mon[i].r_avg) annual_table[i].r_smax = current_mon[i].r_avg; // Wettest summer precipitation
                    if(annual_table[i].r_smin > current_mon[i].r_avg) annual_table[i].r_smin = current_mon[i].r_avg; // Driest summer precipitation
                }
            }
        }
    }

    // Set missing values
    // Iterate grid
    for(int i = 0; i < model_meta->nlat * model_meta->nlon; i++)
    {
        int set_missing = 0;

        // No data?
        if(annual_table[i].n == 0) set_missing = 1;
        // Check precipitation
        // Precipitation cant be negative
        else if(
                annual_table[i].r_avg    < 0.f ||
                annual_table[i].r_min    < 0.f ||
                annual_table[i].r_max    < 0.f ||
                annual_table[i].r_s      < 0.f ||
                annual_table[i].r_w      < 0.f ||
                annual_table[i].r_annual < 0.f
                )
        {
            LOG_ERROR("Detected negative precipitation in model %s\n", model_meta->name);
            LOG_ERROR("%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
                    annual_table[i].t_avg,
                    annual_table[i].t_min,
                    annual_table[i].t_max,
                    annual_table[i].r_avg,
                    annual_table[i].r_min,
                    annual_table[i].r_max,
                    annual_table[i].r_s,
                    annual_table[i].r_smin,
                    annual_table[i].r_smax,
                    annual_table[i].r_w,
                    annual_table[i].r_wmin,
                    annual_table[i].r_wmax,
                    annual_table[i].r_annual,
                    annual_table[i].PC1,
                    annual_table[i].PC2,
                    annual_table[i].PC3,
                    annual_table[i].PC4,
                    annual_table[i].PC5,
                    annual_table[i].PC6,
                    annual_table[i].PC7,
                    annual_table[i].PC8,
                    annual_table[i].PC9,
                    annual_table[i].PC10,
                    annual_table[i].PC11,
                    annual_table[i].PC12,
                    annual_table[i].PC13
                   );
            set_missing = 1;
        }

        // Set the data to missing?
        if(set_missing == 1)
        {
            annual_table[i].t_avg    = MISSING_VAL;
            annual_table[i].t_min    = MISSING_VAL;
            annual_table[i].t_max    = MISSING_VAL;
            annual_table[i].r_avg    = MISSING_VAL;
            annual_table[i].r_min    = MISSING_VAL;
            annual_table[i].r_max    = MISSING_VAL;
            annual_table[i].r_s      = MISSING_VAL;
            annual_table[i].r_smin   = MISSING_VAL;
            annual_table[i].r_smax   = MISSING_VAL;
            annual_table[i].r_w      = MISSING_VAL;
            annual_table[i].r_wmin   = MISSING_VAL;
            annual_table[i].r_wmax   = MISSING_VAL;
            annual_table[i].r_annual = MISSING_VAL;
            annual_table[i].PC1      = MISSING_VAL;
            annual_table[i].PC2      = MISSING_VAL;
            annual_table[i].PC3      = MISSING_VAL;
            annual_table[i].PC4      = MISSING_VAL;
            annual_table[i].PC5      = MISSING_VAL;
            annual_table[i].PC6      = MISSING_VAL;
            annual_table[i].PC7      = MISSING_VAL;
            annual_table[i].PC8      = MISSING_VAL;
            annual_table[i].PC9      = MISSING_VAL;
            annual_table[i].PC10     = MISSING_VAL;
            annual_table[i].PC11     = MISSING_VAL;
            annual_table[i].PC12     = MISSING_VAL;
            annual_table[i].PC13     = MISSING_VAL;
        }
    }

    // Copy output
    memcpy(model_table, annual_table, model_meta->nlat * model_meta->nlon * sizeof(stats_t));

    // Free memory
    free(annual_table);
    free(mon_table);
    free(pr);
    free(ts);

    // Close model dir
    closedir(model_fd);

    return(0);
}

/* *
 * Loads a float-based array from file.
 * */
int _load_data(const char *data_file, size_t n, float *data)
{
    assert(data_file);
    assert(n > 0);
    assert(data);

    // Data file descriptor
    FILE *fdata;

    // Number of elements readed
    size_t readed;

    // Try to open file
    fdata = fopen(data_file, "r");
    if(fdata == NULL) return(-1);

    // Read data
    readed = fread(data, sizeof(float), n, fdata);

    // Check success
    if(readed != n)
    {
        LOG_ERROR("Failed to read file '%s', expected %u elements, readed %u\n", data_file, (unsigned int)n, (unsigned int)readed);
        fclose(fdata);
        return(-1);
    }

    // Close file
    fclose(fdata);

    return(0);
}

/* *
 * Obtain the list of model folders from a directory.
 * */
int _get_models_list(const char *dir, const int initial_year, char ***models, int *n_models)
{
    assert(dir);
    assert(models);

    DIR *bil_fd;
    struct stat sb;
    struct dirent *dp;
    char model_dir[4096];

    char **models_dir;
    size_t max_models_dir;
    size_t n_models_dir;

    // Set output
    *models = NULL;
    *n_models = 0;

    // Open dir
    if((bil_fd = opendir(dir)) == NULL)
    {
        LOG_ERROR("Error accessing %s\n", dir);
        return(-1);
    }

    // First allocation
    max_models_dir = 10;
    models_dir = malloc(max_models_dir * sizeof(char *));
    n_models_dir = 0;

    // Iterate models
    while((dp = readdir(bil_fd)) != NULL)
    {
        // Skip current dir and parent
        if(strcmp(dp->d_name, ".") == 0 || strcmp(dp->d_name, "..") == 0) continue;

        // Skip koppen dir
        if(strcmp(dp->d_name, "koppen") == 0) continue;

        // Skip koppen dir
        if(strcmp(dp->d_name, "koppen_ex") == 0) continue;

        // Skip kmeans dir
        if(strcmp(dp->d_name, "kmeans") == 0) continue;

        // Skip kmeans dir
        if(strcmp(dp->d_name, "kmeans_dist") == 0) continue;

        // Get path and check
        snprintf(model_dir, 4096, "%s/%s", dir, dp->d_name);
        if(stat(model_dir, &sb) != 0)
        {
            LOG_ERROR("Error accessing %s\n", model_dir);
            continue;
        }
        // Is directory?
        if(!S_ISDIR(sb.st_mode)) continue;

        // Get model metadata
        model_t model_meta;
        if(_model_metadata(model_dir, &model_meta, initial_year, 1) != 0)
        {
            LOG_ERROR("Model: %s, INVALID DATA\n", model_meta.name);
            continue;
        }

        LOG_INFO("Model: %15s -- (%03lu,%03lu) -- %04d/%02d -> %04d/%02d\n",
                model_meta.name,
                model_meta.nlat, model_meta.nlon,
                model_meta.year_ini, model_meta.mon_ini,
                model_meta.year_end, model_meta.mon_end
                );

        // Move slots needed?
        if(n_models_dir >= max_models_dir)
        {
            models_dir = (char **)realloc(models_dir, (max_models_dir + 10) * sizeof(char *));
            max_models_dir += 10;
        }

        // Save model to list
        models_dir[n_models_dir] = (char *)malloc((strlen(model_dir) + 1) * sizeof(char));
        strcpy(models_dir[n_models_dir], model_dir);
        n_models_dir++;
    }

    // Close dir
    closedir(bil_fd);

    // Set output
    *models = models_dir;
    *n_models = n_models_dir;

    return(0);
}

/* *
 * Load centroids from a file (.centroids or not).
 * */
int _load_kmeans_init_file(const char *filename, stats_t **out_k_centroids, int *out_k)
{
    assert(filename);
    assert(out_k_centroids);
    assert(out_k);

    FILE *fd;
    stats_t *k_centroids;

    // Try open file (if failed it will be NULL anyways)
    fd = fopen(filename, "r");
    if(fd)
    {
        int n = 0;
        // Count file lines
        for(char c = getc(fd); c != EOF; c = getc(fd))
            if(c == '\n') n++;

        // Rewind file
        rewind(fd);

        // Allocate
        k_centroids = (stats_t *)malloc(n * sizeof(stats_t));

        // Iterate file and initialize centroids
        char aux_line[4096];
        int i = 0;
        while(fgets(aux_line, 4096, fd))
        {
            int nelem = sscanf(aux_line, "%f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                    &k_centroids[i].t_avg,
                    &k_centroids[i].t_min,
                    &k_centroids[i].t_max,
                    &k_centroids[i].r_avg,
                    &k_centroids[i].r_min,
                    &k_centroids[i].r_max,
                    &k_centroids[i].r_s,
                    &k_centroids[i].r_smin,
                    &k_centroids[i].r_smax,
                    &k_centroids[i].r_w,
                    &k_centroids[i].r_wmin,
                    &k_centroids[i].r_wmax,
                    &k_centroids[i].r_annual
                    );

            if(nelem < STATS_T_N)
            {
                LOG_ERROR("Skipping line in file '%s'.\n", filename);
                continue;
            }

            // Set class
            k_centroids[i].class = i;

            // Increase index
            i++;
        }

        if(i > 0)
        {
            // Output
            *out_k = i;
            *out_k_centroids = k_centroids;
        }
        else
        {
            free(k_centroids);
            return(-1);
        }
    }
    else
    {
        LOG_ERROR("Failed to open Kmeans init file '%s'\n", filename);
        return(-1);
    }

    return(0);
}
