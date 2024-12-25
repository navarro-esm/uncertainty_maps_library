/*
 * =====================================================================================
 *
 *       Filename:  kmeans.c
 *
 *    Description:  Compute K-means classification.
 *
 *        Version:  1.0
 *        Created:  07/09/17 10:15:35
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Raul Moreno Galdon
 *   Organization:  University of Castilla-La Mancha
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <float.h>
#include <limits.h>

#include "kmeans.h"
#include "pca.h"
#include "log.h"
#include "utils.h"
#include "struct.h"

/* *
 * Compute K-means classes from a sample grid.
 * Returns the total energy (accumulated distances) of the new classification.
 *  grid: Sample grid (dimensions nlat * nlon)
 *  k: Number of output classes
 *  nlat: Latitude dimension of the grid
 *  nlon: Longitude dimension of the grid
 *  k_centroids: Output centroids for every class (dimension k)
 *  k_count: Output number of occurrences for every class (dimension k)
 *  initial_centroids: Used with KMEANS_FROM_FILE and can be NULL otherwise.
 *  koppen_centroids: Used for transformation.
 *  koppen_centroids_n: Number of centroids in koppen_centroids.
 *  eigenv: Eigenvectors to compute PC after updating classes.
 *  weights: Weigths for computing the distances
 * */
float kmeans(stats_t *grid, const int k, const size_t nlat, const size_t nlon, stats_t *k_centroids, int *k_count, const stats_t *initial_centroids, const double *eigenv, const double *weights)
{
    assert(grid);
    assert(nlat > 0);
    assert(nlon > 0);
    assert(k > 0);
    assert(eigenv);

    // Stats means
    stats_t *kst;

    // Iterations
    int it;
    int max_it = 1000;

    // Count
    int *n_k;

    // Distances
    float prev_dist;
    float curr_dist;

    // Return
    float ret = 0.f;

    // Use PC?
    int use_PC = 0;
#ifdef KMEANS_USE_PC
    use_PC = 1;
#endif

    // Allocate memory for every k
    kst = (stats_t *)malloc(k * sizeof(stats_t));
    if(!kst)
    {
        LOG_ERROR("Failed to allocate stats in kmeans.\n");
        return(-1.f);
    }
    memset(kst, 0, k * sizeof(stats_t)); // Set all to zero
    n_k = (int *)malloc(k * sizeof(int));
    if(!n_k)
    {
        LOG_ERROR("Failed to allocate n_k in kmeans.\n");
        return(-1.f);
    }
    memset(n_k, 0, k * sizeof(int)); // Set all to zero

    // Copy initial centroids to the working copy of the centroids
    memcpy(kst, initial_centroids, k * sizeof(stats_t));

    // Iterate
    it = 0;
    do
    {
        /**
         * Classify
         */
        int b_changed = 0;
        for(int pos = 0; pos < nlat * nlon; pos++)
        {
            // Save previous class
            int old_class = grid[pos].class;

            // Set class as missing at first
            grid[pos].class = CLASS_MASK;

            // Skip missing
            if(grid[pos].t_avg <= MISSING_VAL || grid[pos].r_avg <= MISSING_VAL) continue;

            // Iterate clusters and find nearest distance
            prev_dist = FLT_MAX;
            for(int c = 0; c < k; c++)
            {
                // Get distance to this cluster (no weights)
                // Use Principal Components distance
                curr_dist = kdistance(&grid[pos], &kst[c], weights, use_PC);
                if(curr_dist < prev_dist)
                {
                    // Assign class to this position and update distance
                    grid[pos].class = kst[c].class;
                    grid[pos].class_dist = curr_dist;
                    prev_dist = curr_dist;
                }
            }

            // Check convergence (changes in the classification)
            if(grid[pos].class != old_class) b_changed = 1;

            // Check class
            if(grid[pos].class == CLASS_MASK)
            {
                LOG_INFO("Dist: %f\n", prev_dist);
                LOG_INFO("WARNING: %08d - T:%8.4f, %8.4f, %8.4f, P:%8.4f, %8.4f, %8.4f, %8.4f, %8.4f - %1d\n", grid[pos].class,
                        grid[pos].t_avg,
                        grid[pos].t_min,
                        grid[pos].t_max,
                        grid[pos].r_avg,
                        grid[pos].r_min,
                        grid[pos].r_max,
                        grid[pos].r_s,
                        grid[pos].r_w,
                        k
                        );
                LOG_ERROR("Failed to classify, unknown error in kmeans.\n");
                ret = -1.f;
                break;
            }
        }

        /* *
         * Check for convergence
         * */
        if(b_changed == 0)
        {
            // The classification did not change
            LOG_INFO("K-means converged after %d iterations\n", it);
            break;
        }

        /**
         * Compute new centroids
         */
        kmeans_centroids(grid, nlat * nlon, k, kst, n_k);

        // Update PC
        compute_PC(kst, k, eigenv);

        // Increate iteration
        it++;
    } while(it < max_it);

    /**
     * Classify again with new classes
     */
    ret = kclassify(nlat * nlon, kst, k, weights, grid);

    // Output centroids
    if(k_centroids) memcpy(k_centroids, kst, k * sizeof(stats_t));
    if(k_count) memcpy(k_count, n_k, k * sizeof(int));

    // Divide energy
    ret /= (float)(nlat * nlon);

    // Free memory
    free(kst);
    free(n_k);

    return(ret);
}

/* *
 * Classify values depending on the centroids.
 * The class of the grid will be changed and the distances will be stored too.
 * Returns the total energy (accumulated distances) of the new classification.
 *  grid_n: Number of elements in the grid.
 *  k_centroids: Class centroids to use on the classification.
 *  n_k_centroids: Number of classes.
 *  n_k_centroids: Number of classes.
 *  weights: Weights to use in distances computation (see kdistance function).
 *  grid: Input grid, output attributes are class and class_dist.
 * */
float kclassify(const size_t grid_n, const stats_t *k_centroids, const int n_k_centroids, const double *weights, stats_t *grid)
{
    assert(grid_n > 0);
    assert(k_centroids);
    assert(n_k_centroids > 0);
    assert(grid);

    // Distance values
    float prev_dist;
    float curr_dist;

    // Returned energy
    float ret = 0.f;

    // Use PC?
    int use_PC = 0;
#ifdef KMEANS_USE_PC
    use_PC = 1;
#endif

    // Iterate stats and classify them
    for(int pos = 0; pos < grid_n; pos++)
    {
        // Set class as missing at first
        grid[pos].class = CLASS_MASK;
        grid[pos].class_dist = MISSING_VAL;

        // Skip missing
        if(grid[pos].t_avg <= MISSING_VAL || grid[pos].r_avg <= MISSING_VAL) continue;

        // Iterate clusters and find nearest distance
        prev_dist = FLT_MAX;
        for(int c = 0; c < n_k_centroids; c++)
        {
            // Skip if class is not valid/missing
            if(k_centroids[c].t_avg <= MISSING_VAL) continue;

            // Get distance to this cluster
            curr_dist = kdistance(&grid[pos], &k_centroids[c], weights, use_PC);
            if(curr_dist < prev_dist)
            {
                // Set class
                int class = c;

                // Assign centroid class to this position and update distance
                grid[pos].class = k_centroids[c].class;
                grid[pos].class_dist = curr_dist;
                prev_dist = curr_dist;
            }
        }

        // Add energy
        ret += grid[pos].class_dist;
    }

    return(ret);
}

/* *
 * Compute the Euclidean distance between a sample and a class.
 * Returns the Euclidean distance.
 *  sample: Input sample.
 *  class: Reference class centroid.
 *  weights: Weights for the different variables, useful when using PC.
 *  use_PC: Use PC attributes to compute the distances the distance.
 * */
float kdistance(const stats_t *sample, const stats_t *class, const double *weights, const int use_PC)
{
    assert(sample);
    assert(class);

    // Euclidean distance
    float dist = 0.f;

    // Apply weights
    if(use_PC == 1)
    {
        if(weights != NULL)
        {
            dist += weights[0]  * pow(class->PC1  - sample->PC1 , 2);
            dist += weights[1]  * pow(class->PC2  - sample->PC2 , 2);
            dist += weights[2]  * pow(class->PC3  - sample->PC3 , 2);
            dist += weights[3]  * pow(class->PC4  - sample->PC4 , 2);
            dist += weights[4]  * pow(class->PC5  - sample->PC5 , 2);
            dist += weights[5]  * pow(class->PC6  - sample->PC6 , 2);
            dist += weights[6]  * pow(class->PC7  - sample->PC7 , 2);
            dist += weights[7]  * pow(class->PC8  - sample->PC8 , 2);
            dist += weights[8]  * pow(class->PC9  - sample->PC9 , 2);
            dist += weights[9]  * pow(class->PC10 - sample->PC10, 2);
            dist += weights[10] * pow(class->PC11 - sample->PC11, 2);
            dist += weights[11] * pow(class->PC12 - sample->PC12, 2);
            dist += weights[12] * pow(class->PC13 - sample->PC13, 2);
        }
        else
        {
            dist += pow(class->PC1  - sample->PC1 , 2);
            dist += pow(class->PC2  - sample->PC2 , 2);
            dist += pow(class->PC3  - sample->PC3 , 2);
            dist += pow(class->PC4  - sample->PC4 , 2);
            dist += pow(class->PC5  - sample->PC5 , 2);
            dist += pow(class->PC6  - sample->PC6 , 2);
            dist += pow(class->PC7  - sample->PC7 , 2);
            dist += pow(class->PC8  - sample->PC8 , 2);
            dist += pow(class->PC9  - sample->PC9 , 2);
            dist += pow(class->PC10 - sample->PC10, 2);
            dist += pow(class->PC11 - sample->PC11, 2);
            dist += pow(class->PC12 - sample->PC12, 2);
            dist += pow(class->PC13 - sample->PC13, 2);
        }
    }
    else
    {
        if(weights != NULL)
        {
            // Compute distance
            dist += weights[0]  * pow(class->t_avg    - sample->t_avg   , 2);
            dist += weights[1]  * pow(class->t_min    - sample->t_min   , 2);
            dist += weights[2]  * pow(class->t_max    - sample->t_max   , 2);
            dist += weights[3]  * pow(class->r_avg    - sample->r_avg   , 2);
            dist += weights[4]  * pow(class->r_min    - sample->r_min   , 2);
            dist += weights[5]  * pow(class->r_max    - sample->r_max   , 2);
            dist += weights[6]  * pow(class->r_s      - sample->r_s     , 2);
            dist += weights[7]  * pow(class->r_smin   - sample->r_smin  , 2);
            dist += weights[8]  * pow(class->r_smax   - sample->r_smax  , 2);
            dist += weights[9]  * pow(class->r_w      - sample->r_w     , 2);
            dist += weights[10] * pow(class->r_wmin   - sample->r_wmin  , 2);
            dist += weights[11] * pow(class->r_wmax   - sample->r_wmax  , 2);
            dist += weights[12] * pow(class->r_annual - sample->r_annual, 2);
        }
        else
        {
            // Compute distance
            dist += pow(class->t_avg    - sample->t_avg   , 2);
            dist += pow(class->t_min    - sample->t_min   , 2);
            dist += pow(class->t_max    - sample->t_max   , 2);
            dist += pow(class->r_avg    - sample->r_avg   , 2);
            dist += pow(class->r_min    - sample->r_min   , 2);
            dist += pow(class->r_max    - sample->r_max   , 2);
            dist += pow(class->r_s      - sample->r_s     , 2);
            dist += pow(class->r_smin   - sample->r_smin  , 2);
            dist += pow(class->r_smax   - sample->r_smax  , 2);
            dist += pow(class->r_w      - sample->r_w     , 2);
            dist += pow(class->r_wmin   - sample->r_wmin  , 2);
            dist += pow(class->r_wmax   - sample->r_wmax  , 2);
            dist += pow(class->r_annual - sample->r_annual, 2);
        }
    }

    return sqrt(dist);
}

/* *
 * Compute centroids from an already classified grid.
 * It uses the class attribute of the stats_t structure.
 * The number of ocurrences are also computed.
 *  grid: Grid data to process (1D).
 *  n_grid: Number of elements in the grid.
 *  k: Number of output classes.
 *  k_centroids: Output centroids, must have dimension k.
 *  n_k: Output centroids coincidences, must have dimension k.
 * */
int kmeans_centroids(const stats_t *grid, const size_t n_grid, const int k, stats_t *k_centroids, int *n_k)
{
    assert(grid);
    assert(k > 0);
    assert(n_grid > 0);
    assert(k_centroids);
    assert(n_k);

    // Count clusters
    memset(n_k, 0, k * sizeof(int));
    for(int pos = 0; pos < n_grid; pos++)
    {
        unsigned char class = grid[pos].class;
        if(class >= k)
        {
            if(class != CLASS_MASK) LOG_ERROR("Error: Unknown class found in K-Means %d\n", class);
        }
        else
        {
            n_k[class]++;
        }
    }

    // Prepare to compute new means and set everything to zero
    memset(k_centroids, 0, k * sizeof(stats_t));

    // Set class number
    // Set to missing value the classes with zero occurences
    for(int c = 0; c < k; c++)
    {
        // Set class number
        k_centroids[c].class = c;

        // Set missing values if no ocurrences
        if(n_k[c] <= 0)
        {
            // This class is missing
            k_centroids[c].t_avg    = MISSING_VAL;
            k_centroids[c].t_min    = MISSING_VAL;
            k_centroids[c].t_max    = MISSING_VAL;
            k_centroids[c].r_avg    = MISSING_VAL;
            k_centroids[c].r_min    = MISSING_VAL;
            k_centroids[c].r_max    = MISSING_VAL;
            k_centroids[c].r_s      = MISSING_VAL;
            k_centroids[c].r_smin   = MISSING_VAL;
            k_centroids[c].r_smax   = MISSING_VAL;
            k_centroids[c].r_w      = MISSING_VAL;
            k_centroids[c].r_wmin   = MISSING_VAL;
            k_centroids[c].r_wmax   = MISSING_VAL;
            k_centroids[c].r_annual = MISSING_VAL;
            k_centroids[c].PC1      = MISSING_VAL;
            k_centroids[c].PC2      = MISSING_VAL;
            k_centroids[c].PC3      = MISSING_VAL;
            k_centroids[c].PC4      = MISSING_VAL;
            k_centroids[c].PC5      = MISSING_VAL;
            k_centroids[c].PC6      = MISSING_VAL;
            k_centroids[c].PC7      = MISSING_VAL;
            k_centroids[c].PC8      = MISSING_VAL;
            k_centroids[c].PC9      = MISSING_VAL;
            k_centroids[c].PC10     = MISSING_VAL;
            k_centroids[c].PC11     = MISSING_VAL;
            k_centroids[c].PC12     = MISSING_VAL;
            k_centroids[c].PC13     = MISSING_VAL;
        }
    }

    // Accumulate new values to update the means
    for(int pos = 0; pos < n_grid; pos++)
    {
        unsigned char class = grid[pos].class; // Get cluster

        // Skip if missing
        if(class >= CLASS_MASK || class >= k) continue;

        // Skip if n = 0
        if(n_k[class] <= 0) continue;

        // Update means
        k_centroids[class].t_avg    += grid[pos].t_avg    / (float)(n_k[class]);
        k_centroids[class].t_min    += grid[pos].t_min    / (float)(n_k[class]);
        k_centroids[class].t_max    += grid[pos].t_max    / (float)(n_k[class]);
        k_centroids[class].r_avg    += grid[pos].r_avg    / (float)(n_k[class]);
        k_centroids[class].r_min    += grid[pos].r_min    / (float)(n_k[class]);
        k_centroids[class].r_max    += grid[pos].r_max    / (float)(n_k[class]);
        k_centroids[class].r_s      += grid[pos].r_s      / (float)(n_k[class]);
        k_centroids[class].r_smin   += grid[pos].r_smin   / (float)(n_k[class]);
        k_centroids[class].r_smax   += grid[pos].r_smax   / (float)(n_k[class]);
        k_centroids[class].r_w      += grid[pos].r_w      / (float)(n_k[class]);
        k_centroids[class].r_wmin   += grid[pos].r_wmin   / (float)(n_k[class]);
        k_centroids[class].r_wmax   += grid[pos].r_wmax   / (float)(n_k[class]);
        k_centroids[class].r_annual += grid[pos].r_annual / (float)(n_k[class]);
    }

    return(0);
}

/* *
 * Save centroids data to file.
 *  k_centroids: Centroids array to save.
 *  k: Number of centroids.
 *  filename: Path to the output file.
 * */
int kmeans_centroids_to_file(const stats_t *k_centroids, const int k, const char *filename)
{
    assert(k_centroids);
    assert(k);
    assert(filename);

    FILE *fd;

    // Open file
    fd = fopen(filename, "w");
    if(fd == NULL)
    {
        LOG_ERROR("Error opening file '%s'\n", filename);
        return(-1);
    }

    // Print header
    fprintf(fd, "%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
            "Temp", "Tmin", "Tmax", "Prec", "Pmin", "Pmax", "Psummer", "PSmin", "PSmax", "Pwinter", "PWmin", "PWmax", "Pannual", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13");

    // Iterate classes and print to file
    for(int i = 0; i < k; i++)
    {
        fprintf(fd, "%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
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
                k_centroids[i].r_annual,
                k_centroids[i].PC1,
                k_centroids[i].PC2,
                k_centroids[i].PC3,
                k_centroids[i].PC4,
                k_centroids[i].PC5,
                k_centroids[i].PC6,
                k_centroids[i].PC7,
                k_centroids[i].PC8,
                k_centroids[i].PC9,
                k_centroids[i].PC10,
                k_centroids[i].PC11,
                k_centroids[i].PC12,
                k_centroids[i].PC13
               );
    }

    // Close file
    fclose(fd);

    return(0);
}

/* *
 * Get minimum and maximum stats.
 *  samples: Input samples data.
 *  n_samples: number of samples.
 *  min: output minimum values for the data (optional).
 *  max: output maximum values for the data (optional).
 * */
int kinit_min_max(const stats_t *samples, const int n_samples, stats_t *min, stats_t *max)
{
    assert(samples);
    assert(n_samples > 0);

    // Initialize min
    if(min != NULL)
    {
        min->t_avg    = 9999.f;
        min->t_min    = 9999.f;
        min->t_max    = 9999.f;
        min->r_avg    = 9999.f;
        min->r_min    = 9999.f;
        min->r_max    = 9999.f;
        min->r_s      = 9999.f;
        min->r_smin   = 9999.f;
        min->r_smax   = 9999.f;
        min->r_w      = 9999.f;
        min->r_wmin   = 9999.f;
        min->r_wmax   = 9999.f;
        min->r_annual = 9999.f;
        min->PC1      = 9999.f;
        min->PC2      = 9999.f;
        min->PC3      = 9999.f;
        min->PC4      = 9999.f;
        min->PC5      = 9999.f;
        min->PC6      = 9999.f;
        min->PC7      = 9999.f;
        min->PC8      = 9999.f;
        min->PC9      = 9999.f;
        min->PC10     = 9999.f;
        min->PC11     = 9999.f;
        min->PC12     = 9999.f;
        min->PC13     = 9999.f;
    }

    // Initialize max
    if(max != NULL)
    {
        max->t_avg    = -9999.f;
        max->t_min    = -9999.f;
        max->t_max    = -9999.f;
        max->r_avg    = -9999.f;
        max->r_min    = -9999.f;
        max->r_max    = -9999.f;
        max->r_s      = -9999.f;
        max->r_smax   = -9999.f;
        max->r_smin   = -9999.f;
        max->r_w      = -9999.f;
        max->r_wmax   = -9999.f;
        max->r_wmin   = -9999.f;
        max->r_annual = -9999.f;
        max->PC1      = -9999.f;
        max->PC2      = -9999.f;
        max->PC3      = -9999.f;
        max->PC4      = -9999.f;
        max->PC5      = -9999.f;
        max->PC6      = -9999.f;
        max->PC7      = -9999.f;
        max->PC8      = -9999.f;
        max->PC9      = -9999.f;
        max->PC10     = -9999.f;
        max->PC11     = -9999.f;
        max->PC12     = -9999.f;
        max->PC13     = -9999.f;
    }

    // Get min stats
    if(min != NULL)
    {
        for(int i = 0; i < n_samples; i++)
        {
            // Skip missing
            if(samples[i].t_avg <= MISSING_VAL) continue;

            // Min
            min->t_avg    = MIN(samples[i].t_avg     , min->t_avg);
            min->t_min    = MIN(samples[i].t_min     , min->t_min);
            min->t_max    = MIN(samples[i].t_max     , min->t_max);
            min->r_avg    = MIN(samples[i].r_avg     , min->r_avg);
            min->r_min    = MIN(samples[i].r_min     , min->r_min);
            min->r_max    = MIN(samples[i].r_max     , min->r_max);
            min->r_s      = MIN(samples[i].r_s       , min->r_s);
            min->r_smax   = MIN(samples[i].r_smax    , min->r_smax);
            min->r_smin   = MIN(samples[i].r_smin    , min->r_smin);
            min->r_w      = MIN(samples[i].r_w       , min->r_w);
            min->r_wmax   = MIN(samples[i].r_wmax    , min->r_wmax);
            min->r_wmin   = MIN(samples[i].r_wmin    , min->r_wmin);
            min->r_annual = MIN(samples[i].r_annual  , min->r_annual);
            min->PC1      = MIN(samples[i].PC1       , min->PC1);
            min->PC2      = MIN(samples[i].PC2       , min->PC2);
            min->PC3      = MIN(samples[i].PC3       , min->PC3);
            min->PC4      = MIN(samples[i].PC4       , min->PC4);
            min->PC5      = MIN(samples[i].PC5       , min->PC5);
            min->PC6      = MIN(samples[i].PC6       , min->PC6);
            min->PC7      = MIN(samples[i].PC7       , min->PC7);
            min->PC8      = MIN(samples[i].PC8       , min->PC8);
            min->PC9      = MIN(samples[i].PC9       , min->PC9);
            min->PC10     = MIN(samples[i].PC10      , min->PC10);
            min->PC11     = MIN(samples[i].PC11      , min->PC11);
            min->PC12     = MIN(samples[i].PC12      , min->PC12);
            min->PC13     = MIN(samples[i].PC13      , min->PC13);
        }
    }

    // Get max stats
    if(max != NULL)
    {
        for(int i = 0; i < n_samples; i++)
        {
            // Max
            max->t_avg    = MAX(samples[i].t_avg   , max->t_avg);
            max->t_min    = MAX(samples[i].t_min   , max->t_min);
            max->t_max    = MAX(samples[i].t_max   , max->t_max);
            max->r_avg    = MAX(samples[i].r_avg   , max->r_avg);
            max->r_min    = MAX(samples[i].r_min   , max->r_min);
            max->r_max    = MAX(samples[i].r_max   , max->r_max);
            max->r_s      = MAX(samples[i].r_s     , max->r_s);
            max->r_smax   = MAX(samples[i].r_smax  , max->r_smax);
            max->r_smin   = MAX(samples[i].r_smin  , max->r_smin);
            max->r_w      = MAX(samples[i].r_w     , max->r_w);
            max->r_wmax   = MAX(samples[i].r_wmax  , max->r_wmax);
            max->r_wmin   = MAX(samples[i].r_wmin  , max->r_wmin);
            max->r_annual = MAX(samples[i].r_annual, max->r_annual);
            max->PC1      = MAX(samples[i].PC1     , max->PC1);
            max->PC2      = MAX(samples[i].PC2     , max->PC2);
            max->PC3      = MAX(samples[i].PC3     , max->PC3);
            max->PC4      = MAX(samples[i].PC4     , max->PC4);
            max->PC5      = MAX(samples[i].PC5     , max->PC5);
            max->PC6      = MAX(samples[i].PC6     , max->PC6);
            max->PC7      = MAX(samples[i].PC7     , max->PC7);
            max->PC8      = MAX(samples[i].PC8     , max->PC8);
            max->PC9      = MAX(samples[i].PC9     , max->PC9);
            max->PC10     = MAX(samples[i].PC10    , max->PC10);
            max->PC11     = MAX(samples[i].PC11    , max->PC11);
            max->PC12     = MAX(samples[i].PC12    , max->PC12);
            max->PC13     = MAX(samples[i].PC13    , max->PC13);
        }
    }

    return(0);
}

/* *
 * Initialize K-means using a uniform distribution between maximum and minimum values.
 *  k: number of classes.
 *  min: minimum values for the class.
 *  max: maximum values for the class.
 *  kclasses: output classes, must have dimension k.
 * */
int kinit_distributed(const int k, const stats_t *min, const stats_t *max, stats_t *kclasses)
{
    assert(k > 0);
    assert(min);
    assert(max);
    assert(kclasses);

    // Compute initial averages
    for(int i = 0; i < k; i++)
    {
        kclasses[i].t_avg    = ((max->t_avg    - min->t_avg)    / (float)k * (float)i) + min->t_avg;
        kclasses[i].t_min    = ((max->t_min    - min->t_min)    / (float)k * (float)i) + min->t_min;
        kclasses[i].t_max    = ((max->t_max    - min->t_max)    / (float)k * (float)i) + min->t_max;
        kclasses[i].r_avg    = ((max->r_avg    - min->r_avg)    / (float)k * (float)i) + min->r_avg;
        kclasses[i].r_min    = ((max->r_min    - min->r_min)    / (float)k * (float)i) + min->r_min;
        kclasses[i].r_max    = ((max->r_max    - min->r_max)    / (float)k * (float)i) + min->r_max;
        kclasses[i].r_s      = ((max->r_s      - min->r_s)      / (float)k * (float)i) + min->r_s  ;
        kclasses[i].r_smin   = ((max->r_smin   - min->r_smin)   / (float)k * (float)i) + min->r_smin;
        kclasses[i].r_smax   = ((max->r_smax   - min->r_smax)   / (float)k * (float)i) + min->r_smax;
        kclasses[i].r_w      = ((max->r_w      - min->r_w)      / (float)k * (float)i) + min->r_w  ;
        kclasses[i].r_wmin   = ((max->r_wmin   - min->r_wmin)   / (float)k * (float)i) + min->r_wmin;
        kclasses[i].r_wmax   = ((max->r_wmax   - min->r_wmax)   / (float)k * (float)i) + min->r_wmax;
        kclasses[i].r_annual = ((max->r_annual - min->r_annual) / (float)k * (float)i) + min->r_annual;
    }

    return(0);
}

// Return random number between 0 and 1
static inline float random_uniform()
{
    float ret;

    ret = (float)(rand() % 10000);
    ret /= 10000.f;

    return(ret);
}

/* *
 * Initialize K-means using a random (uniform) distribution between maximum and minimum values.
 *  k: number of classes.
 *  min: minimum values for the class.
 *  max: maximum values for the class.
 *  kclasses: output classes, must have dimension k.
 * */
int kinit_random(const int k, const stats_t *min, const stats_t *max, stats_t *kclasses)
{
    assert(k > 0);
    assert(min);
    assert(max);
    assert(kclasses);

    // Compute initial averages
    for(int i = 0; i < k; i++)
    {
        kclasses[i].t_avg    = (random_uniform() * (max->t_avg    - min->t_avg))    + min->t_avg;
        kclasses[i].t_min    = (random_uniform() * (max->t_min    - min->t_min))    + min->t_min;
        kclasses[i].t_max    = (random_uniform() * (max->t_max    - min->t_max))    + min->t_max;
        kclasses[i].r_avg    = (random_uniform() * (max->r_avg    - min->r_avg))    + min->r_avg;
        kclasses[i].r_min    = (random_uniform() * (max->r_min    - min->r_min))    + min->r_min;
        kclasses[i].r_max    = (random_uniform() * (max->r_max    - min->r_max))    + min->r_max;
        kclasses[i].r_s      = (random_uniform() * (max->r_s      - min->r_s)  )    + min->r_s  ;
        kclasses[i].r_smin   = (random_uniform() * (max->r_smin   - min->r_smin))   + min->r_smin;
        kclasses[i].r_smax   = (random_uniform() * (max->r_smax   - min->r_smax))   + min->r_smax;
        kclasses[i].r_w      = (random_uniform() * (max->r_w      - min->r_w)  )    + min->r_w  ;
        kclasses[i].r_wmin   = (random_uniform() * (max->r_wmin   - min->r_wmin))   + min->r_wmin;
        kclasses[i].r_wmax   = (random_uniform() * (max->r_wmax   - min->r_wmax))   + min->r_wmax;
        kclasses[i].r_annual = (random_uniform() * (max->r_annual - min->r_annual)) + min->r_annual;
    }

    return(0);
}
