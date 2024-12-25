/*
 * =====================================================================================
 *
 *       Filename:  struct.h
 *
 *    Description:  Data structures used for classification.
 *
 *        Version:  1.0
 *        Created:  07/09/17 10:17:21
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Raul Moreno Galdon
 *   Organization:  University of Castilla-La Mancha
 *
 * =====================================================================================
 */

#ifndef STRUCT_H
#define STRUCT_H

/* *
 * Missing value
 * */
#define MISSING_VAL -9999.0

/* *
 * NetCDF file metadata.
 * */
typedef struct file {
    char model[32];
    int variable;
    int year_ini;
    int year_end;
    int mon_ini;
    int mon_end;
} file_t;

/* *
 * Model metadata
 * */
typedef struct model {
    char name[128];
    size_t nlat, nlon;
    int year_ini;
    int mon_ini;
    int year_end;
    int mon_end;
} model_t;

/* *
 * Statistics needed for classification
 * */
typedef struct stats {
    float t_avg;    // Average temperature
    float t_min;    // Coldest temperature
    float t_max;    // Warmest temperature

    float r_avg;    // Average precipitation
    float r_min;    // Dryiest precipitation
    float r_max;    // Wettest precipitation

    float r_s;      // Summer precipitation
    float r_smin;   // Driest summer month
    float r_smax;   // Wettest summer month

    float r_w;      // Winter precipitation
    float r_wmin;   // Driest winter month
    float r_wmax;   // Wettest winter month

    float r_annual; // Annual precipitation

    // Auxiliar
    int n;

    // Classification
    unsigned char class;
    float class_dist;

    // PC
    // TODO: Use array instead of variables
    float PC1;
    float PC2;
    float PC3;
    float PC4;
    float PC5;
    float PC6;
    float PC7;
    float PC8;
    float PC9;
    float PC10;
    float PC11;
    float PC12;
    float PC13;
} stats_t;

/* *
 * Number of attributes in stats_t
 * */
#define STATS_T_N 13

#endif /* STRUCT_H */
