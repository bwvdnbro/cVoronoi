//
// Created by yuyttenh on 09/06/2021.
//

#ifndef CVORONOI_GEOMETRY_H
#define CVORONOI_GEOMETRY_H

#include "dimensionality.h"

#if defined(DIMENSIONALITY_2D)
#include "geometry2d.h"
#else
#include "geometry3d.h"
#endif

#endif  // CVORONOI_GEOMETRY_H
