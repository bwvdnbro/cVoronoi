//
// Created by yuyttenh on 21/06/2021.
//

#include <delaunay.h>

inline static void test_cube() {
  struct delaunay d;
  struct hydro_space hs;
  double dim[3] = {1, 1, 1};
  hydro_space_init(&hs, dim);
  delaunay_init(&d, &hs, 8, 16);

  delaunay_add_local_vertex(&d, 0, 0, 0, 0);
  delaunay_add_local_vertex(&d, 1, 0, 0, 1);
  delaunay_add_local_vertex(&d, 2, 0, 1, 0);
  delaunay_add_local_vertex(&d, 3, 0, 1, 1);
  delaunay_add_local_vertex(&d, 4, 1, 0, 0);
  delaunay_add_local_vertex(&d, 5, 1, 0, 1);
  delaunay_add_local_vertex(&d, 6, 1, 1, 0);
  delaunay_add_local_vertex(&d, 7, 1, 1, 1);

  delaunay_print_tessellation(&d, "cube.txt");

  delaunay_destroy(&d);
}

int main() {
  test_cube();
}

