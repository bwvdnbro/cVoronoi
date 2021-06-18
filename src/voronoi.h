//
// Created by yuyttenh on 08/06/2021.
//

#ifndef CVORONOI_VORONOI_H

#include <stdlib.h>
#include <stdio.h>

#define voronoi_error(s, ...) \
  fprintf(stderr, s, ##__VA_ARGS__); \
  abort();

/*! @brief Activate runtime assertions. */
#define VORONOI_DO_ASSERTIONS

/**
 *@brief Evaluate the given condition and abort if it evaluates to true.
 *
 * This macro is similar to the standard assert() macro.
 * This macro is only defined when VORONOI_DO_ASSERTIONS is active.
 */
#ifdef VORONOI_DO_ASSERTIONS
#define voronoi_assert(condition)                                    \
  if (!(condition)) {                                                 \
    fprintf(stderr, "%s:%s():%i: Condition failed: " #condition "\n", \
            __FILE__, __FUNCTION__, __LINE__);                        \
    abort();                                                          \
  }
#else
#define delaunay_assert(condition)
#endif

#include "dimensionality.h"

#if defined(DIMENSIONALITY_2D)
#include "voronoi2d.h"
#else
#include "voronoi3d.h"
#endif

#define CVORONOI_VORONOI_H

#endif  // CVORONOI_VORONOI_H
