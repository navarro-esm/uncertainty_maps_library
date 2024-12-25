/*
 * =====================================================================================
 *
 *       Filename:  utils.h
 *
 *    Description:  Auxiliar functions.
 *
 *        Version:  1.0
 *        Created:  07/09/17 10:21:59
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Raul Moreno Galdon
 *   Organization:  University of Castilla-La Mancha
 *
 * =====================================================================================
 */

#ifndef UTILS_H
#define UTILS_H

/* *
 * Get min/max of two numbers.
 * */
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

/* *
 * Converts from Kelvin to Celsius.
 * */
static inline float celsius(float kelvin)
{
    return kelvin - 273.15f;
}

/* *
 * Converts from Celsius to Kelvin.
 * */
static inline float kelvin(float celsius)
{
    return celsius + 273.15f;
}

/* *
 * Converts from mm to inches.
 * */
static inline float mm_to_inches(float mm)
{
    return mm / 25.4f;
}

/* *
 * Convert 3D indexes of a matrix to a 1D index.
 * */
static inline int matrix_3D_to_1D(const int x, const int y, const int z, const size_t ldb, const size_t ldc)
{
    return((x * ldb * ldc) + (y * ldc) + z);
}

#endif /* UTILS_H */
