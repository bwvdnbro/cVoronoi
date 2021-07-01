//
// Created by yuyttenh on 08/06/2021.
//
/**
 * @file dimensionality.h
 *
 * @brief Sets the dimensionality of the code
 */

#ifndef CVORONOI_DIMENSIONALITY_H

/* Define dimensionality of code */
//#define DIMENSIONALITY_2D
#define DIMENSIONALITY_3D

#if !defined(DIMENSIONALITY_2D) && !defined(DIMENSIONALITY_3D)
#error "Invalid or undefined dimensionality"
#endif

#define CVORONOI_DIMENSIONALITY_H

#endif  // CVORONOI_DIMENSIONALITY_H
