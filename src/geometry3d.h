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
  mpz_t aix, aiy, aiz, bix, biy, biz, cix, ciy, ciz, dix, diy, diz, eix, eiy,
      eiz;

  /*! @brief Temporary variables used to store relative vertex coordinates. */
  mpz_t s1x, s1y, s1z, s2x, s2y, s2z, s3x, s3y, s3z, s4x, s4y, s4z;

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
  mpz_inits(g->aix, g->aiy, g->aiz, g->bix, g->biy, g->biz, g->cix, g->ciy,
            g->ciz, g->dix, g->diy, g->diz, g->eix, g->eiy, g->eiz, g->s1x,
            g->s1y, g->s1z, g->s2x, g->s2y, g->s2z, g->s3x, g->s3y, g->s3z,
            g->s4x, g->s4y, g->s4z, g->tmp1, g->tmp2, g->result, NULL);
}

/**
 * @brief Deallocate all memory occupied by the geometry object.
 *
 * @param g Geometry object.
 */
inline static void geometry_destroy(struct geometry* restrict g) {
  mpz_clears(g->aix, g->aiy, g->aiz, g->bix, g->biy, g->biz, g->cix, g->ciy,
             g->ciz, g->dix, g->diy, g->diz, g->eix, g->eiy, g->eiz, g->s1x,
             g->s1y, g->s1z, g->s2x, g->s2y, g->s2z, g->s3x, g->s3y, g->s3z,
             g->s4x, g->s4y, g->s4z, g->tmp1, g->tmp2, g->result, NULL);
}

inline static double geometry_orient() {
  // TODO
  return -1.;
}

/**
 * @brief Test the orientation of the tetrahedron that has the four given
 * points as vertices.
 *
 * The test returns a positive result if the fourth vertex is below the plane
 * through the three other vertices, with above the direction from which the
 * three points are ordered counterclockwise.
 *
 * E.g. if the four points are (0, 0, 0), (0, 0, 1), (0, 1, 0), and (1, 0, 0),
 * then this function returns 1.
 *
 * If the four points are exactly coplanar, then this function returns 0.
 *
 * @param a First vertex.
 * @param b Second vertex.
 * @param c Third vertex.
 * @param d Fourth vertex.
 * @return -1, 0, or 1, depending on the orientation of the tetrahedron.
 */
inline static int geometry_orient_exact(
    struct geometry* g, const unsigned long ax, const unsigned long ay,
    const unsigned long az, const unsigned long bx, const unsigned long by,
    const unsigned long bz, const unsigned long cx, const unsigned long cy,
    const unsigned long cz, const unsigned long dx, const unsigned long dy,
    const unsigned long dz) {

  /* store the input coordinates into the temporary large integer variables */
  mpz_set_ui(g->aix, ax);
  mpz_set_ui(g->aiy, ay);
  mpz_set_ui(g->aiz, az);

  mpz_set_ui(g->bix, bx);
  mpz_set_ui(g->biy, by);
  mpz_set_ui(g->biz, bz);

  mpz_set_ui(g->cix, cx);
  mpz_set_ui(g->ciy, cy);
  mpz_set_ui(g->ciz, cz);

  mpz_set_ui(g->dix, dx);
  mpz_set_ui(g->diy, dy);
  mpz_set_ui(g->diz, dz);

  /* compute large integer relative coordinates */
  mpz_sub(g->s1x, g->aix, g->dix);
  mpz_sub(g->s1y, g->aiy, g->diy);
  mpz_sub(g->s1z, g->aiz, g->diz);

  mpz_sub(g->s2x, g->bix, g->dix);
  mpz_sub(g->s2y, g->biy, g->diy);
  mpz_sub(g->s2z, g->biz, g->diz);

  mpz_sub(g->s3x, g->cix, g->dix);
  mpz_sub(g->s3y, g->ciy, g->diy);
  mpz_sub(g->s3z, g->ciz, g->diz);

  /* Compute the result in 3 steps */
  mpz_set_ui(g->result, 0);

  mpz_mul(g->tmp1, g->s2x, g->s3y);
  mpz_submul(g->tmp1, g->s3x, g->s2y);
  mpz_addmul(g->result, g->s1z, g->tmp1);

  mpz_mul(g->tmp1, g->s3x, g->s1y);
  mpz_submul(g->tmp1, g->s1x, g->s3y);
  mpz_addmul(g->result, g->s2z, g->tmp1);

  mpz_mul(g->tmp1, g->s1x, g->s2y);
  mpz_submul(g->tmp1, g->s2x, g->s1y);
  mpz_addmul(g->result, g->s3z, g->tmp1);

  return mpz_sgn(g->result);
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
