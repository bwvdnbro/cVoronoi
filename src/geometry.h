//
// Created by yuyttenh on 09/06/2021.
//

/**
 * @file geometry.h
 *
 * @brief Includes the correct geometry implementation based on the
 * dimensionality.
 */

#ifndef CVORONOI_GEOMETRY_H
#define CVORONOI_GEOMETRY_H

#include "geometry2d.h"
#include "dimensionality.h"
#if defined(DIMENSIONALITY_3D)
#include "geometry3d.h"
#endif

#endif  // CVORONOI_GEOMETRY_H
