/*
 * =====================================================================================
 *
 *       Filename:  pca.h
 *
 *    Description:  Functions for using Principal Component Analysis
 *
 *        Version:  1.0
 *        Created:  02/10/17 11:56:57
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Raul Moreno Galdon
 *   Organization:  University of Castilla-La Mancha
 *
 * =====================================================================================
 */

#ifndef PCA_H
#define PCA_H

#include <stdlib.h>

#include "struct.h"

/* *
 * Convert an array of stats_t into an contiguous array of double elements.
 *  data: Array of data to convert.
 *  n_data: Number of stats_t elements in data.
 *  stats_array: Output array with dimensions (count, 13) with at least count > n_data allocated.
 *  count: Output records in the stats_array.
 * */
int stats_to_array(const stats_t *data, const size_t n_data, double *stats_array, size_t *count);

/* *
 * Compute the averages array of the attributes in data.
 *  data: Input (n_data, n_values) array to compute.
 *  n_data: Dimension of the input data.
 *  n_values: Number of attributes per record.
 *  avg: Output averages with dimension n_values.
 * */
int compute_avg(const double *data, const size_t n_data, const size_t n_values, double *avg);

/* *
 * Compute the standard deviation array of the attributes in data.
 *  data: Input (n_data, n_values) array to compute.
 *  n_data: Dimension of the input data.
 *  n_values: Number of attributes per record.
 *  dev: Output standard deviations with dimension n_values.
 * */
int compute_dev(const double *data, const size_t n_data, const size_t n_values, const double *avg, double *dev);

/* *
 * Compute the matrix of covariances.
 *  data: Input (n_data, n_values) array to compute.
 *  n_data: Dimension of the input data.
 *  n_values: Number of attributes per record.
 *  avg: Averages array (n_values) for computing the covariances.
 *  cov: Output (n_values, n_values) covariances.
 * */
int compute_cov(const double *data, const size_t n_data, const size_t n_values, const double *avg, double *cov);

/* *
 * Compute the eigenvectors and eigenvalues of a square matrix using LAPACK library.
 *  cov: Covariance matrix (order, order).
 *  order: Leading dimension of the covariance matrix.
 *  eigenvalues: Output eigenvalues (order), ordered by size.
 *  eigenvectors: Output eigenvectors (order, order), ordered by eigenvalues.
 * */
int compute_eigen(const double *cov, const size_t order, double *eigenvalues, double *eigenvectors);

/* *
 * Standardize the input values / table.
 *  table: Records of stats_t to standardize (n_stats).
 *  n_stats: Number of elements of table.
 *  avg: Averages used for the standardization.
 *  dev: Standard deviations used for the standardization.
 *  std_table: Output table with standardize values (n_stats), preallocated with at least n_stats records.
 * */
void standardize_stats(const stats_t *table, const size_t n_stats, const stats_t *avg, const stats_t *dev, stats_t *std_table);

/* *
 * Revert the standardization of the input values / table.
 *  std_table: Records of stats_t to revert (n_stats).
 *  n_stats: Number of elements of table.
 *  avg: Averages used for the standardization.
 *  dev: Standard deviations used for the standardization.
 *  table: Output table with reverted values (n_stats), preallocated with at least n_stats records.
 * */
void destandardize_stats(const stats_t *std_table, const size_t n_stats, const stats_t *avg, const stats_t *dev, stats_t *table);

/* *
 * Compute the 13 Principal Components of the input data using the provided eigenvectors.
 * TODO: The PC should be stored and iterated as an array.
 *  data: Input data (n_data).
 *  n_data: Number of elements in data.
 *  eigenv: Eigenvectors for computing the PC.
 * */
void compute_PC(stats_t *st, const int n_st, const double *eigenv);

#endif
