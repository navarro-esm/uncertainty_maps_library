/*
 * =====================================================================================
 *
 *       Filename:  kmeans.h
 *
 *    Description:  Compute K-means classification.
 *
 *        Version:  1.0
 *        Created:  07/09/17 10:36:01
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Raul Moreno Galdon
 *   Organization:  University of Castilla-La Mancha
 *
 * =====================================================================================
 */

#ifndef KMEANS_H
#define KMEANS_H

#include "struct.h"

#define CLASS_MASK 99 // Missing class

/* *
 * Initialization strategies.
 * */
#define KMEANS_DISTRIBUTED 0
#define KMEANS_RAND_UNIFORM 1
#define KMEANS_INIT 2

/* *
 * Compute K-means classes from a sample grid.
 * Returns the total energy (accumulated distances) of the new classification.
 *  st: Sample grid (dimensions nlat * nlon)
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
 *  */
float kmeans(stats_t *grid, const int k, const size_t nlat, const size_t nlon, stats_t *k_centroids, int *k_count, const stats_t *initial_centroids, const double *eigenv, const double *weights);

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
float kclassify(const size_t grid_n, const stats_t *k_centroids, const int n_k_centroids, const double *weights, stats_t *grid);

/* *
 * Compute the Euclidean distance between a sample and a class.
 * Returns the Euclidean distance.
 *  sample: Input sample.
 *  class: Reference class centroid.
 *  weights: Weights for the different variables, useful when using PC.
 *  use_PC: Use PC attributes to compute the distances the distance.
 * */
float kdistance(const stats_t *sample, const stats_t *class, const double *weights, const int use_PC);

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
int kmeans_centroids(const stats_t *grid, const size_t n_grid, const int k, stats_t *k_centroids, int *n_k);

/* *
 * Save centroids data to file.
 *  k_centroids: Centroids array to save.
 *  k: Number of centroids.
 *  filename: Path to the output file.
 * */
int kmeans_centroids_to_file(const stats_t *k_centroids, const int k, const char *filename);

/* *****************
 *
 * INITIALIZATION FUNCTIONS.
 *
 * *****************/

/* *
 * Get minimum and maximum stats.
 *  samples: Input samples data.
 *  n_samples: number of samples.
 *  min: output minimum values for the data (optional).
 *  max: output maximum values for the data (optional).
 * */
int kinit_min_max(const stats_t *st, const int n_st, stats_t *min, stats_t *max);

/* *
 * Initialize K-means using a uniform distribution between maximum and minimum values.
 *  k: number of classes.
 *  min: minimum values for the class.
 *  max: maximum values for the class.
 *  kclasses: output classes, must have dimension k.
 * */
int kinit_distributed(const int k, const stats_t *min, const stats_t *max, stats_t *st);

/* *
 * Initialize K-means using a random (uniform) distribution between maximum and minimum values.
 *  k: number of classes.
 *  min: minimum values for the class.
 *  max: maximum values for the class.
 *  kclasses: output classes, must have dimension k.
 * */
int kinit_random(const int k, const stats_t *min, const stats_t *max, stats_t *st);



#endif /* KMEANS_H */
