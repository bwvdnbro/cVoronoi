//
// Created by yuyttenh on 18/06/2021.
//

#include <stdlib.h>

#include "geometry3d.h"

inline static void test_circumcenter(){
  double circumcenter[3];
  geometry3d_compute_circumcenter(0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, circumcenter);
  if (circumcenter[0] != 0.5 || circumcenter[1] != 0.5 || circumcenter[2] != 0.5) {
    abort();
  }
}

inline static void test_area(){
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
  if (area != sqrt(3.)/2.) {
    abort();
  }
}

inline static void test_volume() {
  double V = geometry3d_compute_volume_tetrahedron(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1);
  if (V != 1. / 6.) {
    abort();
  }
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