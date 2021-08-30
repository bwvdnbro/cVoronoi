/*******************************************************************************
 * This file is part of cVoronoi.
 * Copyright (c) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/**
 * @file geometry2d.h
 *
 * @brief Arbitrary exact and non-exact geometrical tests.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#ifndef SWIFT_GEOMETRY_H
#define SWIFT_GEOMETRY_H

#include <gmp.h>
#include <math.h>

#include "utils.h"

/**
 * @brief Auxiliary variables used by the arbirary exact tests. Since allocating
 * and deallocating these variables poses a significant overhead, they are best
 * reused.
 */
struct geometry2d {
  /*! @brief Arbitrary exact vertex coordinates */
  mpz_t aix, aiy, bix, biy, cix, ciy, dix, diy;

  /*! @brief Temporary variables used to store relative vertex coordinates. */
  mpz_t s1x, s1y, s2x, s2y, s3x, s3y;

  /*! @brief Temporary variables used to store intermediate results. */
  mpz_t tmp1, tmp2;

  /*! @brief Temporary variable used to store final exact results, before their
   *  sign is evaluated and returned as a finite precision integer. */
  mpz_t result;
};

/**
 * @brief Initialize the geometry2d object.
 *
 * This allocates and initialises the auxiliary arbitrary precision variables.
 *
 * @param g Geometry object.
 */
inline static void geometry2d_init(struct geometry2d *restrict g) {
  mpz_inits(g->aix, g->aiy, g->bix, g->biy, g->cix, g->ciy, g->dix, g->diy,
            g->s1x, g->s1y, g->s2x, g->s2y, g->s3x, g->s3y, g->tmp1, g->tmp2,
            g->result, NULL);
}

/**
 * @brief Deallocate all memory occupied by the geometry2d object.
 *
 * @param g Geometry object.
 */
inline static void geometry2d_destroy(struct geometry2d *restrict g) {
  mpz_clears(g->aix, g->aiy, g->bix, g->biy, g->cix, g->ciy, g->dix, g->diy,
             g->s1x, g->s1y, g->s2x, g->s2y, g->s3x, g->s3y, g->tmp1, g->tmp2,
             g->result, NULL);
}

/**
 * @brief Non-exact 2D orientation test.
 *
 * This test returns a positive value if the points a, b and c are ordered
 * counterclockwise, a negative value if they are ordered clockwise, and 0 if
 * they are colinear. The return value corresponds to twice the signed area
 * of the triangle formed by a, b and c, as computed from (Springel, 2010):
 * @f[
 *   \begin{vmatrix}
 *     1 & a_x & a_y \\
 *     1 & b_x & b_y \\
 *     1 & c_x & c_y
 *   \end{vmatrix}
 * @f]
 *
 * Due to roundoff error, this test can be inconsistent with other geometric
 * tests within the same floating point coordinate basis.
 *
 * @param ax, ay First point.
 * @param bx, by Second point.
 * @param cx, cy Third point.
 * @return Signed double area of the triangle formed by a, b and c.
 */
inline static int geometry2d_orient(double ax, double ay, double bx,
                                       double by, double cx, double cy) {

  double det_left = (ax - cx) * (by - cy);
  double det_right = (ay - cy) * (bx - cx);

  double err_bound;
  double result = det_left - det_right;

  if (det_left > 0.0) {
    if (det_right <= 0.0) {
      return sgn(result);
    } else {
      err_bound = det_left + det_right;
    }
  } else if (det_left < 0.0) {
    if (det_right >= 0.0) {
      return sgn(result);
    } else {
      err_bound = -det_left - det_right;
    }
  } else {
    return sgn(result);
  }

  // actual value is smaller, but this will do
  //  err_bound *= 1.e-10;
  err_bound *= DBL_EPSILON * 4;

  if ((result >= err_bound) || (-result >= err_bound)) {
    return sgn(result);
  }

  return 0;
}

/**
 * @brief Arbitrary exact alternative for geometry_orient2d().
 *
 * This function calculates exactly the same thing as the non-exact version, but
 * does so in an integer coordinate basis, using arbitrarily large integers.
 * Since we are not interested in the value of the final result, but only in its
 * sign, this function only returns -1, 0, or 1, depending on the sign of the
 * signed integer triangle area.
 *
 * @param g Geometry object (containing temporary variables that will be used).
 * @param ax, ay, bx, by, cx, cy Integer coordinates of the three points.
 */
inline static int geometry2d_orient_exact(
        struct geometry2d *restrict g, unsigned long int ax, unsigned long int ay,
        unsigned long int bx, unsigned long int by, unsigned long int cx,
        unsigned long int cy) {

  /* store the input coordinates into the temporary large integer variables */
  mpz_set_ui(g->aix, ax);
  mpz_set_ui(g->aiy, ay);

  mpz_set_ui(g->bix, bx);
  mpz_set_ui(g->biy, by);

  mpz_set_ui(g->cix, cx);
  mpz_set_ui(g->ciy, cy);

  /* compute large integer relative coordinates */
  mpz_sub(g->s1x, g->aix, g->cix);
  mpz_sub(g->s1y, g->aiy, g->ciy);

  mpz_sub(g->s2x, g->bix, g->cix);
  mpz_sub(g->s2y, g->biy, g->ciy);

  /* now compute the result using the same 2 steps as the non-exact test */
  mpz_set_ui(g->result, 0);
  mpz_addmul(g->result, g->s1x, g->s2y);
  mpz_submul(g->result, g->s1y, g->s2x);

  /* evaluate the sign of result and return */
  return mpz_sgn(g->result);
}

inline static int geometry2d_orient_adaptive(struct geometry2d* restrict g,
        const unsigned long* al,
        const unsigned long* bl,
        const unsigned long* cl,
        const double* ad, const double* bd,
        const double* cd) {

  int result = geometry2d_orient(ad[0], ad[1], bd[0], bd[1], cd[0], cd[1]);

  if (result == 0) {
    result =
            geometry2d_orient_exact(g, al[0], al[1], bl[0], bl[1], cl[0], cl[1]);
  }

  return result;
}

/**
 * @brief Non-exact 2D in-circle test.
 *
 * This test determines whether the point d is inside the circle through the
 * points a, b and c. Assuming a, b and c are positively oriented (positive
 * return value of geometry2d_orient()), a negative return value means d is
 * outside the circle and a positive value it is inside. A return value of 0
 * means d lies on the circle.
 *
 * This test computes the following (Springel, 2010):
 * @f[
 *   \begin{vmatrix}
 *     1 & a_x & a_y & a_x^2 + a_y^2 \\
 *     1 & b_x & b_y & b_x^2 + b_y^2 \\
 *     1 & c_x & c_y & c_x^2 + c_y^2 \\
 *     1 & d_x & d_y & d_x^2 + d_y^2
 *   \end{vmatrix}
 * @f]
 *
 * Due to roundoff error, this test can be inconsistent with other geometric
 * tests within the same floating point coordinate basis.
 *
 * @param ax, ay First point.
 * @param bx, by Second point.
 * @param cx, cy Third point.
 * @param dx, dy Fourth point.
 */
inline static int geometry2d_in_sphere(double ax, double ay, double bx,
                                       double by, double cx, double cy,
                                       double dx, double dy) {

  /* Compute relative coordinates */
  const double adx = ax - dx;
  const double ady = ay - dy;

  const double bdx = bx - dx;
  const double bdy = by - dy;

  const double cdx = cx - dx;
  const double cdy = cy - dy;

  /* Compute intermediate variables */
  const double bdxcdy = bdx * cdy;
  const double cdxbdy = cdx * bdy;
  const double adnrm2 = adx * adx + ady * ady;

  const double cdxady = cdx * ady;
  const double adxcdy = adx * cdy;
  const double bdnrm2 = bdx * bdx + bdy * bdy;

  const double adxbdy = adx * bdy;
  const double bdxady = bdx * ady;
  const double cdnrm2 = cdx * cdx + cdy * cdy;

  /* Compute errorbound */
  double errbound = (fabs(bdxcdy) + fabs(cdxbdy)) * adnrm2 +
                    (fabs(cdxady) + fabs(adxcdy)) * bdnrm2 +
                    (fabs(adxbdy) + fabs(bdxady)) * cdnrm2;
  // actual value is smaller, but this will do
  //  errbound *= 1.e-10;
  errbound *= DBL_EPSILON * 11;

  /* Compute result */
  const double result = adnrm2 * (bdxcdy - cdxbdy) +
                        bdnrm2 * (cdxady - adxcdy) + cdnrm2 * (adxbdy - bdxady);
  if (result < -errbound || result > errbound) {
    return sgn(result);
  }
  return 0;
}

/**
 * @brief Arbitrary exact alternative for geometry2d_in_sphere().
 *
 * This function calculates exactly the same thing as the non-exact version, but
 * does so in an integer coordinate basis, using arbitrarily large integers.
 * Since we are not interested in the value of the final result, but only in its
 * sign, this function only returns -1, 0, or 1, depending on the sign of the
 * signed integer triangle area.
 *
 * @param g Geometry object (containing temporary variables that will be used).
 * @param ax, ay, bx, by, cx, cy, dx, dy Integer coordinates of the four points.
 */
inline static int geometry2d_in_sphere_exact(
        struct geometry2d *restrict g, unsigned long int ax, unsigned long int ay,
        unsigned long int bx, unsigned long int by, unsigned long int cx,
        unsigned long int cy, unsigned long int dx, unsigned long int dy) {

  /* copy the coordinate values into the large integer temporary variables */
  mpz_set_ui(g->aix, ax);
  mpz_set_ui(g->aiy, ay);

  mpz_set_ui(g->bix, bx);
  mpz_set_ui(g->biy, by);

  mpz_set_ui(g->cix, cx);
  mpz_set_ui(g->ciy, cy);

  mpz_set_ui(g->dix, dx);
  mpz_set_ui(g->diy, dy);

  /* compute the relative coordinates using large integers */
  mpz_sub(g->s1x, g->aix, g->dix);
  mpz_sub(g->s1y, g->aiy, g->diy);

  mpz_sub(g->s2x, g->bix, g->dix);
  mpz_sub(g->s2y, g->biy, g->diy);

  mpz_sub(g->s3x, g->cix, g->dix);
  mpz_sub(g->s3y, g->ciy, g->diy);

  /* compute the result using the same 3 steps as in the non-exact version */
  mpz_set_ui(g->result, 0);

  /* accumulate temporary terms in tmp1 and tmp2 and update result */
  mpz_mul(g->tmp1, g->s2x, g->s3y);
  mpz_submul(g->tmp1, g->s3x, g->s2y);
  mpz_mul(g->tmp2, g->s1x, g->s1x);
  mpz_addmul(g->tmp2, g->s1y, g->s1y);
  mpz_addmul(g->result, g->tmp1, g->tmp2);

  mpz_mul(g->tmp1, g->s3x, g->s1y);
  mpz_submul(g->tmp1, g->s1x, g->s3y);
  mpz_mul(g->tmp2, g->s2x, g->s2x);
  mpz_addmul(g->tmp2, g->s2y, g->s2y);
  mpz_addmul(g->result, g->tmp1, g->tmp2);

  mpz_mul(g->tmp1, g->s1x, g->s2y);
  mpz_submul(g->tmp1, g->s2x, g->s1y);
  mpz_mul(g->tmp2, g->s3x, g->s3x);
  mpz_addmul(g->tmp2, g->s3y, g->s3y);
  mpz_addmul(g->result, g->tmp1, g->tmp2);

  /* evaluate the sign of the result and return */
  return mpz_sgn(g->result);
}

inline static int geometry2d_in_sphere_adaptive(
        struct geometry2d *restrict g, const unsigned long *al,
        const unsigned long *bl, const unsigned long *cl, const unsigned long *dl,
        const double *ad, const double *bd, const double *cd, const double *dd) {

  int result = geometry2d_in_sphere(ad[0], ad[1], bd[0], bd[1], cd[0], cd[1],
                                    dd[0], dd[1]);

  if (result == 0) {
    result = geometry2d_in_sphere_exact(g, al[0], al[1], bl[0], bl[1], cl[0],
                                        cl[1], dl[0], dl[1]);
  }

  return result;
}

/**
 * @brief Compute the coordinates of the circumcenter of the triangle
 * (v0, v1, v2).
 *
 * @param v0x, v0y, v1x, v1y, v2x, v2y Coordinates of the corners of the
 * triangle.
 * @param circumcenter (Returned)
 */
static inline void geometry2d_compute_circumcenter(double v0x, double v0y,
                                                   double v1x, double v1y,
                                                   double v2x, double v2y,
                                                   double *circumcenter) {
  double ax = v1x - v0x;
  double ay = v1y - v0y;
  double bx = v2x - v0x;
  double by = v2y - v0y;

  double D = 2. * (ax * by - ay * bx);
  double a2 = ax * ax + ay * ay;
  double b2 = bx * bx + by * by;
  circumcenter[0] = v0x + (by * a2 - ay * b2) / D;
  circumcenter[1] = v0y + (ax * b2 - bx * a2) / D;
}


inline static void geometry2d_compute_centroid_triangle(double ax, double ay,
                                                        double bx, double by,
                                                        double cx, double cy,
                                                        double* result) {
  result[0] = (ax + bx + cx) / 3;
  result[1] = (ay + by + cy) / 3;
}


/**
 * @brief Compute the volume and centroid of the triangle through the given 3
 * points.
 *
 * @param ax, ay, bx, by, cx, cy Point coordinates.
 * @param result Centroid of the triangle.
 * @return Volume of the triangle.
 */
static inline double geometry2d_compute_centroid_volume_triangle(
        double ax, double ay, double bx, double by, double cx, double cy,
        double *result) {

  geometry2d_compute_centroid_triangle(ax, ay, bx, by, cx, cy, result);

  double s10x = bx - ax;
  double s10y = by - ay;

  double s20x = cx - ax;
  double s20y = cy - ay;

  return 0.5 * fabs(s10x * s20y - s20x * s10y);
}

/**
 * @brief Compute the midpoint and surface area of the face with the given
 * vertex_indices.
 *
 * @param ax, ay, bx, by Face vertex_indices.
 * @param result Midpoint of the face.
 * @return Surface area of the face.
 */
static inline double geometry2d_compute_midpoint_area_face(double ax, double ay,
                                                           double bx, double by,
                                                           double *result) {

  result[0] = 0.5 * (ax + bx);
  result[1] = 0.5 * (ay + by);

  double sx = bx - ax;
  double sy = by - ay;

  return sqrt(sx * sx + sy * sy);
}

inline static int geometry2d_test_line_segment_intersection(
        const double *restrict p1, const double *restrict p2,
        const double *restrict p3, const double *restrict p4) {
  double t =
          (p1[0] - p3[0]) * (p3[1] - p4[1]) - (p1[1] - p3[1]) * (p3[0] - p4[0]);
  double u =
          (p2[0] - p1[0]) * (p1[1] - p3[1]) - (p2[1] - p1[1]) * (p1[0] - p3[0]);

  double denominator =
          (p1[0] - p2[0]) * (p3[1] - p4[1]) - (p1[1] - p2[1]) * (p3[0] - p4[0]);

  int test1 = (fabs(t) <= fabs(denominator)) && (sgn(t) == sgn(denominator));
  int test2 = (fabs(u) <= fabs(denominator)) && (sgn(u) == sgn(denominator));

  return test1 && test2;
}

#endif /* SWIFT_GEOMETRY_H */
