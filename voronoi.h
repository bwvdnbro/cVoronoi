//
// Created by yuyttenh on 08/06/2021.
//

#ifndef CVORONOI_VORONOI_H

#include "dimensionality.h"

#if defined(DIMENSIONALITY_2D)
#include "voronoi2d.h"
#else
#include "voronoi3d.h"
#endif

#define CVORONOI_VORONOI_H

#endif  // CVORONOI_VORONOI_H
