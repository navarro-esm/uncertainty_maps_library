/*
 * =====================================================================================
 *
 *       Filename:  pca.c
 *
 *    Description:  Functions for using Principal Component Analysis
 *
 *        Version:  1.0
 *        Created:  02/10/17 11:59:05
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Raul Moreno Galdon
 *   Organization:  University of Castilla-La Mancha
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "pca.h"
#include "log.h"

/* *
 * Convert an array of stats_t into an contiguous array of double elements.
 *  data: Array of data to convert.
 *  n_data: Number of stats_t elements in data.
 *  stats_array: Output array with dimensions (count, 13) with at least count > n_data allocated.
 *  count: Output records in the stats_array.
 * */
int stats_to_array(const stats_t *data, const size_t n_data, double *stats_array, size_t *count)
{
    assert(data);
    assert(n_data > 0);
    assert(stats_array);
    assert(count);

    // Iterate values
    int aux_count = 0;
    for(int i = 0; i < n_data; i++)
    {
        // Skip if data is not valid/missing
        if(data[i].t_avg <= MISSING_VAL) continue;

        // Get index for this data
        int idx = aux_count * STATS_T_N;

        // Save
        stats_array[idx + 0 ] = data[i].t_avg;
        stats_array[idx + 1 ] = data[i].t_min;
        stats_array[idx + 2 ] = data[i].t_max;
        stats_array[idx + 3 ] = data[i].r_avg;
        stats_array[idx + 4 ] = data[i].r_min;
        stats_array[idx + 5 ] = data[i].r_max;
        stats_array[idx + 6 ] = data[i].r_s;
        stats_array[idx + 7 ] = data[i].r_smin;
        stats_array[idx + 8 ] = data[i].r_smax;
        stats_array[idx + 9 ] = data[i].r_w;
        stats_array[idx + 10] = data[i].r_wmin;
        stats_array[idx + 11] = data[i].r_wmax;
        stats_array[idx + 12] = data[i].r_annual;

        // Increase n
        aux_count++;
    }

    // Set output
    *count = aux_count;

    return(0);
}

/* *
 * Compute the averages array of the attributes in data.
 *  data: Input (n_data, n_values) array to compute.
 *  n_data: Dimension of the input data.
 *  n_values: Number of attributes per record.
 *  avg: Output averages with dimension n_values.
 * */
int compute_avg(const double *data, const size_t n_data, const size_t n_values, double *avg)
{
    assert(data);
    assert(n_data > 0);
    assert(n_values > 0);
    assert(avg);

    // Set all to zero
    memset(avg, 0, n_values * sizeof(double));

    // Iterate records
    for(int i = 0; i < n_data; i++)
    {
        // Get index
        int idx = i * n_values;

        // Iterate values and accumulate average
        for(int j = 0; j < n_values; j++)
        {
            avg[j] += data[idx + j] / (double)n_data;
        }
    }

    return(0);
}

/* *
 * Compute the standard deviation array of the attributes in data.
 *  data: Input (n_data, n_values) array to compute.
 *  n_data: Dimension of the input data.
 *  n_values: Number of attributes per record.
 *  dev: Output standard deviations with dimension n_values.
 * */
int compute_dev(const double *data, const size_t n_data, const size_t n_values, const double *avg, double *dev)
{
    assert(data);
    assert(n_data > 0);
    assert(n_values > 0);
    assert(avg);
    assert(dev);

    // Set all to zero
    memset(dev, 0, n_values * sizeof(double));

    // Iterate records
    for(int i = 0; i < n_data; i++)
    {
        // Get index
        int idx = i * n_values;

        // Iterate values and accumulate deviation
        for(int j = 0; j < n_values; j++)
        {
            dev[j] += ((data[idx + j] - avg[j]) * (data[idx + j] - avg[j])) / n_data;
        }
    }

    // Square root
    for(int i = 0; i < n_values; i++)
        dev[i] = sqrt(dev[i]);

    return(0);
}

/* *
 * Compute the matrix of covariances.
 *  data: Input (n_data, n_values) array to compute.
 *  n_data: Dimension of the input data.
 *  n_values: Number of attributes per record.
 *  avg: Averages array (n_values) for computing the covariances.
 *  cov: Output (n_values, n_values) covariances.
 * */
int compute_cov(const double *data, const size_t n_data, const size_t n_values, const double *avg, double *cov)
{
    assert(data);
    assert(n_data > 0);
    assert(n_values > 0);
    assert(avg);
    assert(cov);

    // Set all to zero
    memset(cov, 0, n_values * n_values * sizeof(double));

    // Iterate records
    for(int i = 0; i < n_data; i++)
    {
        int idx = i * n_values;

        // Compute covariances for each value
        for(int v1 = 0; v1 < n_values; v1++)
        for(int v2 = 0; v2 < n_values; v2++)
        {
            cov[v1 * n_values + v2] += ((data[idx + v1] - avg[v1]) * (data[idx + v2] - avg[v2])) / n_data;
        }
    }

    return(0);
}

/* *
 * Compute the transpose of a square matrix.
 * */
static inline void matrix_transpose(const double *x, const size_t order, double *t)
{
	assert(x);
	assert(order > 1);
	assert(t);

	// Transpose matrix
	for(int i = 0; i < order; i++)
	for(int j = 0; j < order; j++)
	{
		t[j * order + i] = x[i * order + j];
	}
}


// LAPACK functions
void dsyev_(char *JOBZ, char *UPLO, int *N, double *A, int *LDA, double *E_VAL, double *WORK, int *LWORK, int* INFO);

/* *
 * Compute the eigenvectors and eigenvalues of a square matrix using LAPACK library.
 *  cov: Covariance matrix (order, order).
 *  order: Leading dimension of the covariance matrix.
 *  eigenvalues: Output eigenvalues (order), ordered by size.
 *  eigenvectors: Output eigenvectors (order, order), ordered by eigenvalues.
 * */
int compute_eigen(const double *cov, const size_t order, double *eigenvalues, double *eigenvectors)
{
    int n = order;
	int lwork = 32 * order;
	double work[32 * order];
	int info;
	char jobz = 'V'; // Eigen values and vectors
	char uplo = 'U'; // Does not matter, x is simmetrical

    // Sorting
    int swap;
    double *aux_eigenvector;

	// Transpose matrix (FORTRAN)
	matrix_transpose(cov, order, eigenvectors);

	// Call LAPACK routine
	dsyev_(&jobz, &uplo, &n, eigenvectors, &n, eigenvalues, work, &lwork, &info);
    if(info)
    {
        LOG_ERROR("Failed to compute DSYEV: %d", info);
        return(-1);
    }

    // Allocations
    aux_eigenvector = (double *)malloc(order * sizeof(double));

    // Sort eigenvalues and eigenvectors
	do
	{
		swap = 0;
		for(int i = 1; i < order; i++)
		{
			// Check greater
			if(eigenvalues[i-1] < eigenvalues[i])
			{
				// Swap eigen value
				double aux_evalue = eigenvalues[i];
				eigenvalues[i] = eigenvalues[i-1];
				eigenvalues[i-1] = aux_evalue;

				// Swap eigen vector
				memcpy(aux_eigenvector, &eigenvectors[i * order], order * sizeof(double));
				memcpy(&eigenvectors[i * order], &eigenvectors[(i-1) * order], order * sizeof(double));
				memcpy(&eigenvectors[(i-1) * order], aux_eigenvector, order * sizeof(double));

				swap = 1;
			}
		}
	} while(swap == 1);

    // Free
    free(aux_eigenvector);

    return(0);
}

/* *
 * Standardize the input values / table.
 *  table: Records of stats_t to standardize (n_stats).
 *  n_stats: Number of elements of table.
 *  avg: Averages used for the standardization.
 *  dev: Standard deviations used for the standardization.
 *  std_table: Output table with standardize values (n_stats), preallocated with at least n_stats records.
 * */
void standardize_stats(const stats_t *table, const size_t n_stats, const stats_t *avg, const stats_t *dev, stats_t *std_table)
{
    assert(table);
    assert(avg);
    assert(dev);
    assert(n_stats > 0);
    assert(std_table);

    for(int i = 0; i < n_stats; i++)
    {
        // Copy values
        memcpy(&std_table[i], &table[i], sizeof(stats_t));

        // Skip missing values
        if(table[i].t_avg <= MISSING_VAL) continue;

        // Standardize table
        std_table[i].t_avg    = (table[i].t_avg   - avg->t_avg   ) / dev->t_avg   ;
        std_table[i].t_min    = (table[i].t_min   - avg->t_min   ) / dev->t_min   ;
        std_table[i].t_max    = (table[i].t_max   - avg->t_max   ) / dev->t_max   ;
        std_table[i].r_avg    = (table[i].r_avg   - avg->r_avg   ) / dev->r_avg   ;
        std_table[i].r_min    = (table[i].r_min   - avg->r_min   ) / dev->r_min   ;
        std_table[i].r_max    = (table[i].r_max   - avg->r_max   ) / dev->r_max   ;
        std_table[i].r_s      = (table[i].r_s     - avg->r_s     ) / dev->r_s     ;
        std_table[i].r_smin   = (table[i].r_smin  - avg->r_smin  ) / dev->r_smin  ;
        std_table[i].r_smax   = (table[i].r_smax  - avg->r_smax  ) / dev->r_smax  ;
        std_table[i].r_w      = (table[i].r_w     - avg->r_w     ) / dev->r_w     ;
        std_table[i].r_wmin   = (table[i].r_wmin  - avg->r_wmin  ) / dev->r_wmin  ;
        std_table[i].r_wmax   = (table[i].r_wmax  - avg->r_wmax  ) / dev->r_wmax  ;
        std_table[i].r_annual = (table[i].r_annual- avg->r_annual) / dev->r_annual;
    }
}

/* *
 * Revert the standardization of the input values / table.
 *  std_table: Records of stats_t to revert (n_stats).
 *  n_stats: Number of elements of table.
 *  avg: Averages used for the standardization.
 *  dev: Standard deviations used for the standardization.
 *  table: Output table with reverted values (n_stats), preallocated with at least n_stats records.
 * */
void destandardize_stats(const stats_t *std_table, const size_t n_stats, const stats_t *avg, const stats_t *dev, stats_t *table)
{
    assert(table);
    assert(avg);
    assert(dev);
    assert(n_stats > 0);
    assert(std_table);

    for(int i = 0; i < n_stats; i++)
    {
        // Copy values
        memcpy(&table[i], &std_table[i], sizeof(stats_t));

        // Skip missing values
        if(std_table[i].t_avg <= MISSING_VAL) continue;

        // Standardize table
        table[i].t_avg    = (std_table[i].t_avg    * dev->t_avg   ) + avg->t_avg   ;
        table[i].t_min    = (std_table[i].t_min    * dev->t_min   ) + avg->t_min   ;
        table[i].t_max    = (std_table[i].t_max    * dev->t_max   ) + avg->t_max   ;
        table[i].r_avg    = (std_table[i].r_avg    * dev->r_avg   ) + avg->r_avg   ;
        table[i].r_min    = (std_table[i].r_min    * dev->r_min   ) + avg->r_min   ;
        table[i].r_max    = (std_table[i].r_max    * dev->r_max   ) + avg->r_max   ;
        table[i].r_s      = (std_table[i].r_s      * dev->r_s     ) + avg->r_s     ;
        table[i].r_smin   = (std_table[i].r_smin   * dev->r_smin  ) + avg->r_smin  ;
        table[i].r_smax   = (std_table[i].r_smax   * dev->r_smax  ) + avg->r_smax  ;
        table[i].r_w      = (std_table[i].r_w      * dev->r_w     ) + avg->r_w     ;
        table[i].r_wmin   = (std_table[i].r_wmin   * dev->r_wmin  ) + avg->r_wmin  ;
        table[i].r_wmax   = (std_table[i].r_wmax   * dev->r_wmax  ) + avg->r_wmax  ;
        table[i].r_annual = (std_table[i].r_annual * dev->r_annual) + avg->r_annual;
    }

}

/* *
 * Compute the 13 Principal Components of the input data using the provided eigenvectors.
 * TODO: The PC should be stored and iterated as an array.
 *  data: Input data (n_data).
 *  n_data: Number of elements in data.
 *  eigenv: Eigenvectors for computing the PC.
 * */
void compute_PC(stats_t *data, const int n_data, const double *eigenv)
{
    assert(data);
    assert(n_data > 0);
    assert(eigenv);

    for(int i = 0; i < n_data; i++)
    {
        // Missing?
        if(data[i].t_avg <= MISSING_VAL)
        {
            data[i].PC1 = MISSING_VAL;
            data[i].PC2 = MISSING_VAL;
            data[i].PC3 = MISSING_VAL;
            data[i].PC4 = MISSING_VAL;
            data[i].PC5 = MISSING_VAL;
            data[i].PC6 = MISSING_VAL;
            data[i].PC7 = MISSING_VAL;
            data[i].PC8 = MISSING_VAL;
            data[i].PC9 = MISSING_VAL;
            data[i].PC10 = MISSING_VAL;
            data[i].PC11 = MISSING_VAL;
            data[i].PC12 = MISSING_VAL;
            data[i].PC13 = MISSING_VAL;
            continue;
        }

        // First component
        data[i].PC1 =
            (data[i].t_avg    * eigenv[0 ]) +
            (data[i].t_min    * eigenv[1 ]) +
            (data[i].t_max    * eigenv[2 ]) +
            (data[i].r_avg    * eigenv[3 ]) +
            (data[i].r_min    * eigenv[4 ]) +
            (data[i].r_max    * eigenv[5 ]) +
            (data[i].r_s      * eigenv[6 ]) +
            (data[i].r_smin   * eigenv[7 ]) +
            (data[i].r_smax   * eigenv[8 ]) +
            (data[i].r_w      * eigenv[9 ]) +
            (data[i].r_wmin   * eigenv[10]) +
            (data[i].r_wmax   * eigenv[11]) +
            (data[i].r_annual * eigenv[12]);

        // Second component
        data[i].PC2 =
            (data[i].t_avg    * eigenv[STATS_T_N + 0 ]) +
            (data[i].t_min    * eigenv[STATS_T_N + 1 ]) +
            (data[i].t_max    * eigenv[STATS_T_N + 2 ]) +
            (data[i].r_avg    * eigenv[STATS_T_N + 3 ]) +
            (data[i].r_min    * eigenv[STATS_T_N + 4 ]) +
            (data[i].r_max    * eigenv[STATS_T_N + 5 ]) +
            (data[i].r_s      * eigenv[STATS_T_N + 6 ]) +
            (data[i].r_smin   * eigenv[STATS_T_N + 7 ]) +
            (data[i].r_smax   * eigenv[STATS_T_N + 8 ]) +
            (data[i].r_w      * eigenv[STATS_T_N + 9 ]) +
            (data[i].r_wmin   * eigenv[STATS_T_N + 10]) +
            (data[i].r_wmax   * eigenv[STATS_T_N + 11]) +
            (data[i].r_annual * eigenv[STATS_T_N + 12]);

        // Third component
        data[i].PC3 =
            (data[i].t_avg    * eigenv[(STATS_T_N * 2) + 0 ]) +
            (data[i].t_min    * eigenv[(STATS_T_N * 2) + 1 ]) +
            (data[i].t_max    * eigenv[(STATS_T_N * 2) + 2 ]) +
            (data[i].r_avg    * eigenv[(STATS_T_N * 2) + 3 ]) +
            (data[i].r_min    * eigenv[(STATS_T_N * 2) + 4 ]) +
            (data[i].r_max    * eigenv[(STATS_T_N * 2) + 5 ]) +
            (data[i].r_s      * eigenv[(STATS_T_N * 2) + 6 ]) +
            (data[i].r_smin   * eigenv[(STATS_T_N * 2) + 7 ]) +
            (data[i].r_smax   * eigenv[(STATS_T_N * 2) + 8 ]) +
            (data[i].r_w      * eigenv[(STATS_T_N * 2) + 9 ]) +
            (data[i].r_wmin   * eigenv[(STATS_T_N * 2) + 10]) +
            (data[i].r_wmax   * eigenv[(STATS_T_N * 2) + 11]) +
            (data[i].r_annual * eigenv[(STATS_T_N * 2) + 12]);

        // Forth component
        data[i].PC4 =
            (data[i].t_avg    * eigenv[(STATS_T_N * 3) + 0 ]) +
            (data[i].t_min    * eigenv[(STATS_T_N * 3) + 1 ]) +
            (data[i].t_max    * eigenv[(STATS_T_N * 3) + 2 ]) +
            (data[i].r_avg    * eigenv[(STATS_T_N * 3) + 3 ]) +
            (data[i].r_min    * eigenv[(STATS_T_N * 3) + 4 ]) +
            (data[i].r_max    * eigenv[(STATS_T_N * 3) + 5 ]) +
            (data[i].r_s      * eigenv[(STATS_T_N * 3) + 6 ]) +
            (data[i].r_smin   * eigenv[(STATS_T_N * 3) + 7 ]) +
            (data[i].r_smax   * eigenv[(STATS_T_N * 3) + 8 ]) +
            (data[i].r_w      * eigenv[(STATS_T_N * 3) + 9 ]) +
            (data[i].r_wmin   * eigenv[(STATS_T_N * 3) + 10]) +
            (data[i].r_wmax   * eigenv[(STATS_T_N * 3) + 11]) +
            (data[i].r_annual * eigenv[(STATS_T_N * 3) + 12]);

        // Fifth component
        data[i].PC5 =
            (data[i].t_avg    * eigenv[(STATS_T_N * 4) + 0 ]) +
            (data[i].t_min    * eigenv[(STATS_T_N * 4) + 1 ]) +
            (data[i].t_max    * eigenv[(STATS_T_N * 4) + 2 ]) +
            (data[i].r_avg    * eigenv[(STATS_T_N * 4) + 3 ]) +
            (data[i].r_min    * eigenv[(STATS_T_N * 4) + 4 ]) +
            (data[i].r_max    * eigenv[(STATS_T_N * 4) + 5 ]) +
            (data[i].r_s      * eigenv[(STATS_T_N * 4) + 6 ]) +
            (data[i].r_smin   * eigenv[(STATS_T_N * 4) + 7 ]) +
            (data[i].r_smax   * eigenv[(STATS_T_N * 4) + 8 ]) +
            (data[i].r_w      * eigenv[(STATS_T_N * 4) + 9 ]) +
            (data[i].r_wmin   * eigenv[(STATS_T_N * 4) + 10]) +
            (data[i].r_wmax   * eigenv[(STATS_T_N * 4) + 11]) +
            (data[i].r_annual * eigenv[(STATS_T_N * 4) + 12]);

        // Sixth component
        data[i].PC6 =
            (data[i].t_avg    * eigenv[(STATS_T_N * 5) + 0 ]) +
            (data[i].t_min    * eigenv[(STATS_T_N * 5) + 1 ]) +
            (data[i].t_max    * eigenv[(STATS_T_N * 5) + 2 ]) +
            (data[i].r_avg    * eigenv[(STATS_T_N * 5) + 3 ]) +
            (data[i].r_min    * eigenv[(STATS_T_N * 5) + 4 ]) +
            (data[i].r_max    * eigenv[(STATS_T_N * 5) + 5 ]) +
            (data[i].r_s      * eigenv[(STATS_T_N * 5) + 6 ]) +
            (data[i].r_smin   * eigenv[(STATS_T_N * 5) + 7 ]) +
            (data[i].r_smax   * eigenv[(STATS_T_N * 5) + 8 ]) +
            (data[i].r_w      * eigenv[(STATS_T_N * 5) + 9 ]) +
            (data[i].r_wmin   * eigenv[(STATS_T_N * 5) + 10]) +
            (data[i].r_wmax   * eigenv[(STATS_T_N * 5) + 11]) +
            (data[i].r_annual * eigenv[(STATS_T_N * 5) + 12]);

        // Seventh component
        data[i].PC7 =
            (data[i].t_avg    * eigenv[(STATS_T_N * 6) + 0 ]) +
            (data[i].t_min    * eigenv[(STATS_T_N * 6) + 1 ]) +
            (data[i].t_max    * eigenv[(STATS_T_N * 6) + 2 ]) +
            (data[i].r_avg    * eigenv[(STATS_T_N * 6) + 3 ]) +
            (data[i].r_min    * eigenv[(STATS_T_N * 6) + 4 ]) +
            (data[i].r_max    * eigenv[(STATS_T_N * 6) + 5 ]) +
            (data[i].r_s      * eigenv[(STATS_T_N * 6) + 6 ]) +
            (data[i].r_smin   * eigenv[(STATS_T_N * 6) + 7 ]) +
            (data[i].r_smax   * eigenv[(STATS_T_N * 6) + 8 ]) +
            (data[i].r_w      * eigenv[(STATS_T_N * 6) + 9 ]) +
            (data[i].r_wmin   * eigenv[(STATS_T_N * 6) + 10]) +
            (data[i].r_wmax   * eigenv[(STATS_T_N * 6) + 11]) +
            (data[i].r_annual * eigenv[(STATS_T_N * 6) + 12]);

        // Eighth component
        data[i].PC8 =
            (data[i].t_avg    * eigenv[(STATS_T_N * 7) + 0 ]) +
            (data[i].t_min    * eigenv[(STATS_T_N * 7) + 1 ]) +
            (data[i].t_max    * eigenv[(STATS_T_N * 7) + 2 ]) +
            (data[i].r_avg    * eigenv[(STATS_T_N * 7) + 3 ]) +
            (data[i].r_min    * eigenv[(STATS_T_N * 7) + 4 ]) +
            (data[i].r_max    * eigenv[(STATS_T_N * 7) + 5 ]) +
            (data[i].r_s      * eigenv[(STATS_T_N * 7) + 6 ]) +
            (data[i].r_smin   * eigenv[(STATS_T_N * 7) + 7 ]) +
            (data[i].r_smax   * eigenv[(STATS_T_N * 7) + 8 ]) +
            (data[i].r_w      * eigenv[(STATS_T_N * 7) + 9 ]) +
            (data[i].r_wmin   * eigenv[(STATS_T_N * 7) + 10]) +
            (data[i].r_wmax   * eigenv[(STATS_T_N * 7) + 11]) +
            (data[i].r_annual * eigenv[(STATS_T_N * 7) + 12]);

        // Ninth component
        data[i].PC9 =
            (data[i].t_avg    * eigenv[(STATS_T_N * 8) + 0 ]) +
            (data[i].t_min    * eigenv[(STATS_T_N * 8) + 1 ]) +
            (data[i].t_max    * eigenv[(STATS_T_N * 8) + 2 ]) +
            (data[i].r_avg    * eigenv[(STATS_T_N * 8) + 3 ]) +
            (data[i].r_min    * eigenv[(STATS_T_N * 8) + 4 ]) +
            (data[i].r_max    * eigenv[(STATS_T_N * 8) + 5 ]) +
            (data[i].r_s      * eigenv[(STATS_T_N * 8) + 6 ]) +
            (data[i].r_smin   * eigenv[(STATS_T_N * 8) + 7 ]) +
            (data[i].r_smax   * eigenv[(STATS_T_N * 8) + 8 ]) +
            (data[i].r_w      * eigenv[(STATS_T_N * 8) + 9 ]) +
            (data[i].r_wmin   * eigenv[(STATS_T_N * 8) + 10]) +
            (data[i].r_wmax   * eigenv[(STATS_T_N * 8) + 11]) +
            (data[i].r_annual * eigenv[(STATS_T_N * 8) + 12]);

        // Tenth component
        data[i].PC10 =
            (data[i].t_avg    * eigenv[(STATS_T_N * 9) + 0 ]) +
            (data[i].t_min    * eigenv[(STATS_T_N * 9) + 1 ]) +
            (data[i].t_max    * eigenv[(STATS_T_N * 9) + 2 ]) +
            (data[i].r_avg    * eigenv[(STATS_T_N * 9) + 3 ]) +
            (data[i].r_min    * eigenv[(STATS_T_N * 9) + 4 ]) +
            (data[i].r_max    * eigenv[(STATS_T_N * 9) + 5 ]) +
            (data[i].r_s      * eigenv[(STATS_T_N * 9) + 6 ]) +
            (data[i].r_smin   * eigenv[(STATS_T_N * 9) + 7 ]) +
            (data[i].r_smax   * eigenv[(STATS_T_N * 9) + 8 ]) +
            (data[i].r_w      * eigenv[(STATS_T_N * 9) + 9 ]) +
            (data[i].r_wmin   * eigenv[(STATS_T_N * 9) + 10]) +
            (data[i].r_wmax   * eigenv[(STATS_T_N * 9) + 11]) +
            (data[i].r_annual * eigenv[(STATS_T_N * 9) + 12]);

        // Eleventh component
        data[i].PC11 =
            (data[i].t_avg    * eigenv[(STATS_T_N * 10) + 0 ]) +
            (data[i].t_min    * eigenv[(STATS_T_N * 10) + 1 ]) +
            (data[i].t_max    * eigenv[(STATS_T_N * 10) + 2 ]) +
            (data[i].r_avg    * eigenv[(STATS_T_N * 10) + 3 ]) +
            (data[i].r_min    * eigenv[(STATS_T_N * 10) + 4 ]) +
            (data[i].r_max    * eigenv[(STATS_T_N * 10) + 5 ]) +
            (data[i].r_s      * eigenv[(STATS_T_N * 10) + 6 ]) +
            (data[i].r_smin   * eigenv[(STATS_T_N * 10) + 7 ]) +
            (data[i].r_smax   * eigenv[(STATS_T_N * 10) + 8 ]) +
            (data[i].r_w      * eigenv[(STATS_T_N * 10) + 9 ]) +
            (data[i].r_wmin   * eigenv[(STATS_T_N * 10) + 10]) +
            (data[i].r_wmax   * eigenv[(STATS_T_N * 10) + 11]) +
            (data[i].r_annual * eigenv[(STATS_T_N * 10) + 12]);

        // Twelfth component
        data[i].PC12 =
            (data[i].t_avg    * eigenv[(STATS_T_N * 11) + 0 ]) +
            (data[i].t_min    * eigenv[(STATS_T_N * 11) + 1 ]) +
            (data[i].t_max    * eigenv[(STATS_T_N * 11) + 2 ]) +
            (data[i].r_avg    * eigenv[(STATS_T_N * 11) + 3 ]) +
            (data[i].r_min    * eigenv[(STATS_T_N * 11) + 4 ]) +
            (data[i].r_max    * eigenv[(STATS_T_N * 11) + 5 ]) +
            (data[i].r_s      * eigenv[(STATS_T_N * 11) + 6 ]) +
            (data[i].r_smin   * eigenv[(STATS_T_N * 11) + 7 ]) +
            (data[i].r_smax   * eigenv[(STATS_T_N * 11) + 8 ]) +
            (data[i].r_w      * eigenv[(STATS_T_N * 11) + 9 ]) +
            (data[i].r_wmin   * eigenv[(STATS_T_N * 11) + 10]) +
            (data[i].r_wmax   * eigenv[(STATS_T_N * 11) + 11]) +
            (data[i].r_annual * eigenv[(STATS_T_N * 11) + 12]);

        // Thirdteen component
        data[i].PC13 =
            (data[i].t_avg    * eigenv[(STATS_T_N * 12) + 0 ]) +
            (data[i].t_min    * eigenv[(STATS_T_N * 12) + 1 ]) +
            (data[i].t_max    * eigenv[(STATS_T_N * 12) + 2 ]) +
            (data[i].r_avg    * eigenv[(STATS_T_N * 12) + 3 ]) +
            (data[i].r_min    * eigenv[(STATS_T_N * 12) + 4 ]) +
            (data[i].r_max    * eigenv[(STATS_T_N * 12) + 5 ]) +
            (data[i].r_s      * eigenv[(STATS_T_N * 12) + 6 ]) +
            (data[i].r_smin   * eigenv[(STATS_T_N * 12) + 7 ]) +
            (data[i].r_smax   * eigenv[(STATS_T_N * 12) + 8 ]) +
            (data[i].r_w      * eigenv[(STATS_T_N * 12) + 9 ]) +
            (data[i].r_wmin   * eigenv[(STATS_T_N * 12) + 10]) +
            (data[i].r_wmax   * eigenv[(STATS_T_N * 12) + 11]) +
            (data[i].r_annual * eigenv[(STATS_T_N * 12) + 12]);
    }
}
