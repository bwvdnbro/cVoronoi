//
// Created by yuyttenh on 08/06/2021.
//

#ifndef CVORONOI_VORONOI3D_H
#define CVORONOI_VORONOI3D_H

struct voronoi {
  // TODO
  double *cell_centroid;
};

inline static void voronoi_init(struct voronoi *restrict v,
                                struct delaunay *restrict d) {
  // TODO
}

inline static void voronoi_destroy(struct voronoi *restrict v) {
  // TODO
}

inline static void voronoi_print_grid(const struct voronoi *v,
                                      const char *filename) {
  // TODO
}

#endif  // CVORONOI_VORONOI3D_H
