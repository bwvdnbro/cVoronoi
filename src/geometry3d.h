//
// Created by yuyttenh on 09/06/2021.
//

#ifndef CVORONOI_GEOMETRY3D_H
#define CVORONOI_GEOMETRY3D_H

#include <gmp.h>

/**
 * @brief Auxiliary variables used by the arbirary exact tests. Since allocating
 * and deallocating these variables poses a significant overhead, they are best
 * reused.
 */
struct geometry {
  /*! @brief Arbitrary exact vertex coordinates */
  mpz_t aix, aiy, bix, biy, cix, ciy, dix, diy, eix, eiy;

  /*! @brief Temporary variables used to store relative vertex coordinates. */
  mpz_t s1x, s1y, s2x, s2y, s3x, s3y;

  /*! @brief Temporary variables used to store intermediate results. */
  mpz_t tmp1, tmp2;

  /*! @brief Temporary variable used to store final exact results, before their
   *  sign is evaluated and returned as a finite precision integer. */
  mpz_t result;
};

/**
 * @brief Initialize the geometry object.
 *
 * This allocates and initialises the auxiliary arbitrary precision variables.
 *
 * @param g Geometry object.
 */
inline static void geometry_init(struct geometry* restrict g) {
  mpz_inits(g->aix, g->aiy, g->bix, g->biy, g->cix, g->ciy, g->dix, g->diy,
            g->eix, g->eiy, g->s1x, g->s1y, g->s2x, g->s2y, g->s3x, g->s3y,
            g->tmp1, g->tmp2, g->result, NULL);
}

/**
 * @brief Deallocate all memory occupied by the geometry object.
 *
 * @param g Geometry object.
 */
inline static void geometry_destroy(struct geometry* restrict g) {
  mpz_clears(g->aix, g->aiy, g->bix, g->biy, g->cix, g->ciy, g->dix, g->diy,
             g->eix, g->eiy, g->s1x, g->s1y, g->s2x, g->s2y, g->s3x, g->s3y,
             g->tmp1, g->tmp2, g->result, NULL);
}

inline static double geometry_orient() {
  // TODO
  return -1.;
}

inline static int geometry_orient_exact() {
  // TODO
  return -1.;
}

inline static double geometry_in_sphere() {
  // TODO
  return -1.;
}

inline static int geometry_in_sphere_exact() {
  // TODO
  return -1;
}


#endif  // CVORONOI_GEOMETRY3D_H
