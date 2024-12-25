/*
 * =====================================================================================
 *
 *       Filename:  koppen.h
 *
 *    Description:  Compute Koppen classification.
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

#ifndef KOPPEN_H
#define KOPPEN_H

/* *
 * Koppen classes as in Lohman.
 * */
#define CLASS_AF 0
#define CLASS_AW 1
#define CLASS_BS 2
#define CLASS_BW 3
#define CLASS_CS 4
#define CLASS_CW 5
#define CLASS_CF 6
#define CLASS_DW 7
#define CLASS_DF 8
#define CLASS_ET 9
#define CLASS_EF 10

/* *
 * Koppen classes as in Trewartha.
 * Koppen Extended.
 * */

// Tropical
#define CLASS_EX_Af 0
#define CLASS_EX_Aw 1
#define CLASS_EX_Am 2
#define CLASS_EX_As 3
#define CLASS_EX_Ai 4
#define CLASS_EX_Ag 5
// Dry
#define CLASS_EX_BWh 6
#define CLASS_EX_BWk 7
#define CLASS_EX_BWs 8
#define CLASS_EX_BWw 9
#define CLASS_EX_BWn 10
#define CLASS_EX_BSh 11
#define CLASS_EX_BSk 12
#define CLASS_EX_BSs 13
#define CLASS_EX_BSw 14
#define CLASS_EX_BSn 15
// Warm temperate
#define CLASS_EX_Cfa 16
#define CLASS_EX_Cfb 17
#define CLASS_EX_Cwa 18
#define CLASS_EX_Cwb 19
#define CLASS_EX_Csa 20
#define CLASS_EX_Csb 21
// Cold-snowy forest
#define CLASS_EX_Df 22
#define CLASS_EX_Dfd 23
#define CLASS_EX_Dw 24
#define CLASS_EX_Dwd 25
// Polar
#define CLASS_EX_ET 26
#define CLASS_EX_EF 27

// Number of classes
#define KOPPEN_EX_N (CLASS_EX_EF + 1)

/* *
 * Missing value for classes.
 * */
#define CLASS_MASK 99

/* *
 * Compute Koppen classes from a sample grid, as in Lohman.
 *  grid: Sample grid, dimensions (nlat, nlon)
 *  nlat: Latitude dimension of the grid
 *  nlon: Longitude dimension of the grid
 * */
int koppen(stats_t *grid, size_t nlat, size_t nlon);

/* *
 * Compute Extended Koppen classes from a sample grid, as in Trewartha.
 *  grid: Sample grid, dimensions (nlat, nlon)
 *  nlat: Latitude dimension of the grid
 *  nlon: Longitude dimension of the grid
 * */
int koppen_extended(stats_t *grid, size_t nlat, size_t nlon);

#endif /* KOPPEN_H */
