/*
 * =====================================================================================
 *
 *       Filename:  log.h
 *
 *    Description:  Helper functions for logging.
 *
 *        Version:  1.0
 *        Created:  28/07/17 09:57:40
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Raul Moreno,
 *   Organization:  University of Castilla-La Mancha
 *
 * =====================================================================================
 */

#ifndef LOG_H_
#define LOG_H_

#include <stdio.h>
#include <stdarg.h>

static inline void LOG_INFO(char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
#ifdef _OPENMP
#pragma omp critical
#endif
    vprintf(fmt, args);
    va_end(args);
}

static inline void LOG_ERROR(char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
#ifdef _OPENMP
#pragma omp critical
#endif
    vfprintf(stderr, fmt, args);
    va_end(args);
}

#endif
