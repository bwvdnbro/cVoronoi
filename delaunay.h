//
// Created by yuyttenh on 08/06/2021.
//

#ifndef CVORONOI_DELAUNAY_H
#define CVORONOI_DELAUNAY_H

#include "dimensionality.h"

#if defined(DIMENSIONALITY_2D)
#include "delaunay2d.h"
#elif defined(DIMENSIONALITY_3D)
#include "delaunay3d.h"
#else
#error "Invalid or undefined dimensionality"
#endif

#endif  // CVORONOI_DELAUNAY_H
