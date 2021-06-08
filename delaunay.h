//
// Created by yuyttenh on 08/06/2021.
//

#ifndef CVORONOI_DELAUNAY_H
#define CVORONOI_DELAUNAY_H

#include "dimensionality.h"

#if defined(DIMENSIONALITY_2D)
#include "delaunay2d.h"
#else
#include "delaunay3d.h"
#endif

#endif  // CVORONOI_DELAUNAY_H
