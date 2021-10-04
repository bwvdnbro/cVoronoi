//
// Created by yuyttenh on 18/06/2021.
//

#include <stdlib.h>

#include "delaunay.h"
#include "geometry3d.h"

inline static void get_rand_vec(double *v) {
  v[0] = (double) rand() / INT_MAX + 1.;
  v[1] = (double) rand() / INT_MAX + 1.;
  v[2] = (double) rand() / INT_MAX + 1.;
}

inline static void test_circumcenter() {
  double p0[3], p1[3], p2[3], p3[3], c[3];
  unsigned long p0ul[3], p1ul[3], p2ul[3], p3ul[3];

  geometry3d_compute_circumcenter_relative_non_exact(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, c);
  assert(c[0] == 0.5 && c[1] == 0.5 && c[2] == 0.5);

  struct geometry3d g;
  geometry3d_init(&g);

  for (int l = 0; l < 10000; l++) {
    get_rand_vec(p0);
    get_rand_vec(p1);
    get_rand_vec(p2);
    get_rand_vec(p3);

    for (int i = 0; i < 3; i++) {
      p0ul[i] = delaunay_double_to_int(p0[i]);
      p1ul[i] = delaunay_double_to_int(p1[i]);
      p2ul[i] = delaunay_double_to_int(p2[i]);
      p3ul[i] = delaunay_double_to_int(p3[i]);
    }
    double ce[3];
    geometry3d_compute_circumcenter_relative_exact(&g, p0ul[0], p0ul[1], p0ul[2], p1ul[0], p1ul[1], p1ul[2], p2ul[0],
                                                   p2ul[1], p2ul[2], p3ul[0], p3ul[1], p3ul[2], ce);
    geometry3d_compute_circumcenter_relative_adaptive(&g, p0, p1, p2, p3, p0ul, p1ul, p2ul, p3ul, c);
    assert(fabs(ce[0] - c[0]) / ce[0] < 1e-10 && fabs(ce[1] - c[1]) / ce[1] < 1e-10 && fabs(ce[2] - c[2]) / ce[2] < 1e-10);

    c[0] += p0[0];
    c[1] += p0[1];
    c[2] += p0[2];

    double radius = geometry3d_compute_circumradius_adaptive(&g, p0, p1, p2, p3, p0ul, p1ul, p2ul, p3ul, 1.);

    double d0 = sqrt(
            (c[0] - p0[0]) * (c[0] - p0[0]) + (c[1] - p0[1]) * (c[1] - p0[1]) + (c[2] - p0[2]) * (c[2] - p0[2]));
    double d1 = sqrt(
            (c[0] - p1[0]) * (c[0] - p1[0]) + (c[1] - p1[1]) * (c[1] - p1[1]) + (c[2] - p1[2]) * (c[2] - p1[2]));
    double d2 = sqrt(
            (c[0] - p2[0]) * (c[0] - p2[0]) + (c[1] - p2[1]) * (c[1] - p2[1]) + (c[2] - p2[2]) * (c[2] - p2[2]));
    double d3 = sqrt(
            (c[0] - p3[0]) * (c[0] - p3[0]) + (c[1] - p3[1]) * (c[1] - p3[1]) + (c[2] - p3[2]) * (c[2] - p3[2]));

    assert(fabs(d0 - d1) / d0 < 1e-10 && fabs(d0 - d2) / d0 < 1e-10 && fabs(d0 - d3) / d0 < 1e-10 &&
           fabs(d0 - radius) / d0 < 1e-10);
  }

  geometry3d_destroy(&g);
}

inline static void test_area() {
  double area = geometry3d_compute_area_triangle(0, 0, 0, 1, 0, 0, 0, 0, 1);
  if (area != 0.5) {
    abort();
  }
  area = geometry3d_compute_area_triangle(0, 0, 0, 0, 1, 0, 0, 0, 1);
  if (area != 0.5) {
    abort();
  }
  area = geometry3d_compute_area_triangle(0, 0, 0, 1, 0, 0, 0, 1, 0);
  if (area != 0.5) {
    abort();
  }
  area = geometry3d_compute_area_triangle(1, 0, 0, 0, 1, 0, 0, 0, 1);
  if (area != sqrt(3.) / 2.) {
    abort();
  }
}

inline static void test_orientation() {
  struct geometry3d g;
  geometry3d_init(&g);

  unsigned long al[3] = {0, 0, 0};
  unsigned long bl[3] = {0, 0, 1};
  unsigned long cl[3] = {0, 1, 0};
  unsigned long dl[3] = {1, 0, 0};

  double ad[3] = {0, 0, 0};
  double bd[3] = {0, 0, 1};
  double cd[3] = {0, 1, 0};
  double dd[3] = {1, 0, 0};

  int orientation = geometry3d_orient_adaptive(&g, al, bl, cl, dl, ad, bd, cd, dd);
  if (orientation != 1) {
    abort();
  }

  geometry3d_destroy(&g);
}

inline static void test_volume() {
  double V =
          geometry3d_compute_volume_tetrahedron(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1);
  if (V != 1. / 6.) {
    abort();
  }

  V = geometry3d_compute_volume_tetrahedron(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, -0.5,
                                            1);
  if (V != 1. / 6.) {
    abort();
  }

  struct delaunay d;
  struct hydro_space hs;
  double dim[3] = {1, 1, 1};
  hydro_space_init(&hs, dim);
  delaunay_init(&d, &hs, 8, 16);

  delaunay_add_local_vertex(&d, 0, 0, 0, 0);
  delaunay_add_local_vertex(&d, 1, 0.4, 0, 1);
  delaunay_add_local_vertex(&d, 2, 0, 0.54, 0);
  delaunay_add_local_vertex(&d, 3, .2, 0.678, 1);

  V = 0.;
  for (int i = 4; i < d.tetrahedron_index; i++) {
    struct tetrahedron *t = &d.tetrahedra[i];
    if (!t->active) continue;

    int v0 = t->vertices[0];
    int v1 = t->vertices[1];
    int v2 = t->vertices[2];
    int v3 = t->vertices[3];

    V += geometry3d_compute_volume_tetrahedron(
            d.vertices[3 * v0], d.vertices[3 * v0 + 1], d.vertices[3 * v0 + 2],
            d.vertices[3 * v1], d.vertices[3 * v1 + 1], d.vertices[3 * v1 + 2],
            d.vertices[3 * v2], d.vertices[3 * v2 + 1], d.vertices[3 * v2 + 2],
            d.vertices[3 * v3], d.vertices[3 * v3 + 1], d.vertices[3 * v3 + 2]);
  }
  if (fabs(V - 121.5) / V > 0.0000000001) {
    abort();
  }
  delaunay_destroy(&d);
}

inline static void test_centroid_exact() {
  unsigned long p0[3] = {0, 0, 0};
  unsigned long p1[3] = {0, 0, 0};
  unsigned long p2[3] = {0, 20, 0};
  unsigned long p3[3] = {0, 0, 29};
  p1[0] -= 1;
  p2[0] -= 2;
  p3[0] -= 3;

  unsigned long centroid[3];
  geometry3d_compute_centroid_tetrahedron_exact(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], p3[0],
                                                p3[1], p3[2], centroid);
  assert(centroid[0] == (0ul - 1) / 4 * 3 + 2 && centroid[1] == 5 && centroid[2] == 7);
}

inline static void test_ray_triangle_intersection() {
  double p0[3] = {1.1, 1.1, 1.1};
  double p1[3] = {1.8, 1.1, 1.1};
  double p2[3] = {1.1, 1.5, 1.1};
  double p3[3] = {1.1, 1.1, 1.3};
  double p4[3] = {1.8, 1.67, 1.43};

  unsigned long p0ul[3], p1ul[3], p2ul[3], p3ul[3], p4ul[3];
  for (int i = 0; i < 3; i++) {
    p0ul[i] = delaunay_double_to_int(p0[i]);
    p1ul[i] = delaunay_double_to_int(p1[i]);
    p2ul[i] = delaunay_double_to_int(p2[i]);
    p3ul[i] = delaunay_double_to_int(p3[i]);
    p4ul[i] = delaunay_double_to_int(p4[i]);
  }

  struct ray r;
  ray_init(&r, p0, p4, p0ul, p4ul);

  double distance;
  int intersects = geometry3d_ray_triangle_intersect_non_exact(&r, p1, p2, p3, &distance);
  assert(intersects);

  struct geometry3d g;
  geometry3d_init(&g);
  double distance_exact;
  int intersects_exact = geometry3d_ray_triangle_intersect_exact(&g, &r, p1ul, p2ul, p3ul, &distance_exact);
  assert(intersects_exact);

  double distance_adaptive;
  int intersects_adaptive = geometry3d_ray_triangle_intersect(&g, &r, p1, p2, p3, p1ul, p2ul, p3ul, &distance_adaptive);
  assert(intersects_adaptive);

  assert(fabs(distance - distance_exact) / distance < 1e-10);

  /* Now again for random points */
  for (int i = 0; i < 10000; i++) {
    get_rand_vec(p0);
    get_rand_vec(p1);
    get_rand_vec(p2);
    get_rand_vec(p3);
    get_rand_vec(p4);

    for (int j = 0; j < 3; j++) {
      p0ul[j] = delaunay_double_to_int(p0[j]);
      p1ul[j] = delaunay_double_to_int(p1[j]);
      p2ul[j] = delaunay_double_to_int(p2[j]);
      p3ul[j] = delaunay_double_to_int(p3[j]);
      p4ul[j] = delaunay_double_to_int(p4[j]);
    }

    ray_init(&r, p0, p4, p0ul, p4ul);

    intersects = geometry3d_ray_triangle_intersect_non_exact(&r, p1, p2, p3, &distance);
    intersects_exact = geometry3d_ray_triangle_intersect_exact(&g, &r, p1ul, p2ul, p3ul, &distance_exact);
    intersects_adaptive = geometry3d_ray_triangle_intersect(&g, &r, p1, p2, p3, p1ul, p2ul, p3ul, &distance_adaptive);

    assert(intersects_exact == intersects);
    assert((isinf(distance) && isinf(distance_exact) && isinf(distance_adaptive))
            || (fabs(distance - distance_exact) / fabs(distance) < 1e-9 && fabs(distance - distance_adaptive) / fabs(distance) < 1e-9));
  }
  geometry3d_destroy(&g);
}

/**
 * @brief Test for functions of geometry3d.h
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 */
int main(int argc, char **argv) {
  test_circumcenter();
  test_area();
  test_volume();
  test_orientation();
  test_ray_triangle_intersection();
  test_centroid_exact();
}