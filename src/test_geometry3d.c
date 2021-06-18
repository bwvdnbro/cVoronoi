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

/**
 * @brief Test for functions of geometry3d.h
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 */
int main(int argc, char **argv) {
  test_circumcenter();
}