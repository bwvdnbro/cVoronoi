//
// Created by yuyttenh on 18/06/2021.
//

#include <stdlib.h>

#include "delaunay.h"
#include "geometry3d.h"

inline static void test_circumcenter() {
  double circumcenter[3];
  geometry3d_compute_circumcenter(0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1,
                                  circumcenter);
  if (circumcenter[0] != 0.5 || circumcenter[1] != 0.5 ||
      circumcenter[2] != 0.5) {
    abort();
  }
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
}