//
// Created by yuyttenh on 08/06/2021.
//

/**
 * @file voronoi.h
 *
 * @brief includes the correct Voronoi implementation based on the
 * dimensionality.
 */

#ifndef CVORONOI_VORONOI_H

#include <stdlib.h>
#include <stdio.h>
#include "delaunay.h"
#include "geometry.h"
#include "dimensionality.h"

#define voronoi_error(s, ...) \
  fprintf(stderr, s, ##__VA_ARGS__); \
  abort();

/*! @brief Store the edges of faces (so that the actual Voronoi grid can be
 *  reconstructed). */
#define VORONOI_STORE_CONNECTIONS

/*! @brief Store information about the number of faces per cell. */
#define VORONOI_STORE_CELL_STATS

/*! @brief Store cell generators. */
#define VORONOI_STORE_GENERATORS

/*! @brief Activate runtime assertions. */

#define VORONOI_DO_ASSERTIONS

#define VORONOI_CHECKS

/**
 *@brief Evaluate the given condition and abort if it evaluates to true.
 *
 * This macro is similar to the standard assert() macro.
 * This macro is only defined when VORONOI_DO_ASSERTIONS is active.
 */
#ifdef VORONOI_DO_ASSERTIONS
#define voronoi_assert(condition)                                     \
  if (!(condition)) {                                                 \
    fprintf(stderr, "%s:%s():%i: Condition failed: " #condition "\n", \
            __FILE__, __FUNCTION__, __LINE__);                        \
    abort();                                                          \
  }
#else
#define voronoi_assert(condition)
#endif

#if defined(DIMENSIONALITY_2D)
#include "voronoi2d.h"
#else
#include "voronoi3d.h"
#endif

#define CVORONOI_VORONOI_H

#endif  // CVORONOI_VORONOI_H
