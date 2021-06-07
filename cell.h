//
// Created by yuyttenh on 07/06/2021.
//

#ifndef CVORONOI_CELL_H

#include "delaunay.h"
#include "hilbert.h"
#include "hydro_space.h"
#include "sort.h"
#include "voronoi.h"

/**
 * @brief Generate a random uniform double in the range [0, 1].
 *
 * This function makes use of C's rand() and is therefore not suited for
 * accurate statistical sampling. It will however generate sufficiently
 * random numbers for any code test.
 *
 * @return Random uniform double in the range [0, 1].
 */
static inline double get_random_uniform_double() {
  return ((double)rand()) / ((double)RAND_MAX);
}

struct cell {
  int count;

  double *vertices;

  unsigned long *hilbert_keys;

  /* Arg-sort indices for directions: x, y, xyp and xym and the hilbert key */
  int *r_sort_lists[5];

  struct hydro_space hs;

  struct delaunay *d;

  struct voronoi *v;
  int voronoi_active;
};

void cell_init(struct cell *c, const int *count, const double pert, const double *dim) {
  hydro_space_init(&c->hs, dim);
  c->count = count[0] * count[1] * count[2];

  /* slightly randomized vertices */
  c->vertices = (double *)malloc(2 * c->count * sizeof(double));
  int index = 0;
  for (int ix = 0; ix < count[0]; ++ix) {
    for (int iy = 0; iy < count[1]; ++iy) {
      c->vertices[index] = (ix + 0.5 + pert * get_random_uniform_double()) *
                           (dim[0] / (double)count[0]);
      ++index;
      c->vertices[index] = (iy + 0.5 + pert * get_random_uniform_double()) *
                           (dim[1] / (double)count[1]);
      ++index;
    }
  }

  /* hilbert keys */
  c->hilbert_keys = (unsigned long *)malloc(c->count * sizeof(unsigned long));
  for (int i = 0; i < c->count; i++) {
    unsigned long bits[2];
    bits[0] = c->vertices[2 * i] / c->hs.side[0] * (1ul << 32);
    bits[1] = c->vertices[2 * i + 1] / c->hs.side[1] * (1ul << 32);
    c->hilbert_keys[i] = hilbert_get_key_2d(bits, 64);
  }

  /* sorting arrays */
  for (int i = 0; i < 5; i++) {
    c->r_sort_lists[i] = (int *)malloc(c->count * sizeof(int));
    for (int j = 0; j < c->count; j++) {
      c->r_sort_lists[i][j] = j;
    }
  }
  qsort_r(c->r_sort_lists[0], c->count, sizeof(int), sort_x_comp, c->vertices);
  qsort_r(c->r_sort_lists[1], c->count, sizeof(int), sort_y_comp, c->vertices);
  qsort_r(c->r_sort_lists[2], c->count, sizeof(int), sort_xyp_comp,
          c->vertices);
  qsort_r(c->r_sort_lists[3], c->count, sizeof(int), sort_xym_comp,
          c->vertices);
  qsort_r(c->r_sort_lists[4], c->count, sizeof(int), sort_h_comp, c->vertices);

  delaunay_init(c->d, &c->hs, c->count, 10 * c->count);
  c->voronoi_active = 0;
}

void cell_destroy(struct cell *c) {
  free(c->vertices);
  free(c->hilbert_keys);
  for (int i = 0; i < 5; i++) {
    free(c->r_sort_lists[i]);
  }
  delaunay_destroy(c->d);
  if (c->voronoi_active) {
    voronoi_destroy(c->v);
  }
}

void cell_construct_local_delaunay(struct cell *c) {
  /* Add the local vertices, one by one, in Hilbert order. */
  for (int i = 0; i < c->count; ++i) {
    int j = c->r_sort_lists[4][i];
    delaunay_add_vertex(c->d, c->vertices[2 * j], c->vertices[2 * j + 1]);
  }

  /* we are done adding the original vertices. We need to consolidate the
     indices of the original vertices within the Delaunay tessellation, so
     that the tessellation knows that vertices added from now on are ghosts.
   */
  delaunay_consolidate(c->d);
}

void cell_make_delaunay_periodic(struct cell *c) {
  /* Add ghosts to impose the periodic boundaries. These ghosts will be
     periodic copies of the original vertices that ensure that the incomplete
     cells at the boundaries of the simulation box have the right shape.
     Within SWIFT, a similar technique needs to be used to guarantee that
     cells near the boundary of a SWIFT-cell correctly "feel" the cells in the
     neighbouring SWIFT-cells. Springel (2010) provides a criterion for
     completeness, which is based on the radius of the circumcircles of the
     triangles that connect to original vertices. Starting from an initial
     search radius, we will iteratively add ghost vertices and increase this
     search radius until all relevant circumcircle radii are smaller than the
     search radius. This mechanism can very easily be combined with SWIFT's
     existing neighbour search algorithms. */

  /* Initial search radius. We use twice the average inter-particle
     separation, but other values would also work. */
  double r = 2. * c->hs.side[0] * c->hs.side[1] / c->count;
  /* add ghosts for the positive horizontal boundary */
  int i = 0;
  /* update the search radii for all triangles in the tessellation and count
       the number of triangles with circumcircles larger than the current search
       radius */
  int count = delaunay_update_search_radii(c->d, r);
  printf("count: %i\n", count);
  while (count > 0) {
    /* we do not want to add the same ghost twice (this causes the incremental
       construction algorithm to crash), so we need to keep track of the
       previous search radius, so that we can only add the ghosts that have
       not been added before */
    double old_r = r;
    /* now gradually increase the search radius */
    r *= 1.5;
    /* and repeat the additions as above */
    i = 0;
    int vi = c->r_sort_lists[0][i];
    while ((c->vertices[2 * vi] >= old_r) && (c->vertices[2 * vi] < r)) {
      delaunay_add_vertex(
          c->d, c->vertices[2 * vi] + c->hs.anchor[0] + c->hs.side[0],
          c->vertices[2 * vi + 1]);
      ++i;
      if (i == c->count) break;
      vi = c->r_sort_lists[0][i];
    }
    i = c->count - 1;
    vi = c->r_sort_lists[0][i];
    while ((c->hs.anchor[0] + c->hs.side[0] - c->vertices[2 * vi] >= old_r) &&
           (c->hs.anchor[0] + c->hs.side[0] - c->vertices[2 * vi] < r)) {
      delaunay_log("x: %g, old_r: %g, r: %g",
                   c->hs.anchor[0] + c->hs.side[0] - c->vertices[2 * vi],
                   old_r, r);
      delaunay_add_vertex(
          c->d, c->vertices[2 * vi] - c->hs.anchor[0] + c->hs.side[0],
          c->vertices[2 * vi + 1]);
      --i;
      if (i == -1) break;
      vi = c->r_sort_lists[0][i];
    }
    i = 0;
    vi = c->r_sort_lists[1][i];
    while ((c->vertices[2 * vi + 1] >= old_r) &&
           (c->vertices[2 * vi + 1] < r)) {
      delaunay_add_vertex(
          c->d, c->vertices[2 * vi],
          c->vertices[2 * vi + 1] + c->hs.anchor[1] + c->hs.side[1]);
      ++i;
      if (i == c->count) break;
      vi = c->r_sort_lists[1][i];
    }
    i = c->count - 1;
    vi = c->r_sort_lists[1][i];
    while ((c->hs.anchor[1] + c->hs.side[1] - c->vertices[2 * vi + 1] >=
            old_r) &&
           (c->hs.anchor[1] + c->hs.side[1] - c->vertices[2 * vi + 1] < r)) {
      delaunay_add_vertex(
          c->d, c->vertices[2 * vi],
          c->vertices[2 * vi + 1] - c->hs.anchor[1] + c->hs.side[1]);
      --i;
      if (i == -1) break;
      vi = c->r_sort_lists[1][i];
    }
    i = 0;
    vi = c->r_sort_lists[2][i];
    while ((c->vertices[2 * vi] + c->vertices[2 * vi + 1] >= old_r) &&
           (c->vertices[2 * vi] + c->vertices[2 * vi + 1] < r)) {
      delaunay_add_vertex(
          c->d, c->vertices[2 * vi] + c->hs.anchor[0] + c->hs.side[0],
          c->vertices[2 * vi + 1] + c->hs.anchor[1] + c->hs.side[1]);
      ++i;
      if (i == c->count) break;
      vi = c->r_sort_lists[2][i];
    }
    i = c->count - 1;
    vi = c->r_sort_lists[2][i];
    while ((c->hs.anchor[0] + c->hs.side[0] - c->vertices[2 * vi] +
                c->hs.anchor[1] + c->hs.side[1] - c->vertices[2 * vi + 1] >=
            old_r) &&
           (c->hs.anchor[0] + c->hs.side[0] - c->vertices[2 * vi] +
                c->hs.anchor[1] + c->hs.side[1] - c->vertices[2 * vi + 1] <
            r)) {
      delaunay_add_vertex(
          c->d, c->vertices[2 * vi] - c->hs.anchor[0] + c->hs.side[0],
          c->vertices[2 * vi + 1] - c->hs.anchor[1] + c->hs.side[1]);
      --i;
      if (i == -1) break;
      vi = c->r_sort_lists[2][i];
    }
    i = 0;
    vi = c->r_sort_lists[3][i];
    while ((c->vertices[2 * vi] - c->hs.anchor[1] + c->hs.side[1] +
                c->vertices[2 * vi + 1] >=
            old_r) &&
           (c->vertices[2 * vi] - c->hs.anchor[1] + c->hs.side[1] +
                c->vertices[2 * vi + 1] <
            r)) {
      delaunay_add_vertex(
          c->d, c->vertices[2 * vi] + c->hs.anchor[0] + c->hs.side[0],
          c->vertices[2 * vi + 1] - c->hs.anchor[1] + c->hs.side[1]);
      ++i;
      if (i == c->count) break;
      vi = c->r_sort_lists[3][i];
    }
    i = c->count - 1;
    vi = c->r_sort_lists[3][i];
    while ((c->hs.anchor[0] + c->hs.side[0] - c->vertices[2 * vi] -
                c->vertices[2 * vi + 1] >=
            old_r) &&
           (c->hs.anchor[0] + c->hs.side[0] - c->vertices[2 * vi] -
                c->vertices[2 * vi + 1] <
            r)) {
      delaunay_add_vertex(
          c->d, c->vertices[2 * vi] - c->hs.anchor[0] + c->hs.side[0],
          c->vertices[2 * vi + 1] + c->hs.anchor[1] + c->hs.side[1]);
      --i;
      if (i == -1) break;
      vi = c->r_sort_lists[3][i];
    }
    /* update the search radii to the new value and count the number of larger
       circumcircles. */
    count = delaunay_update_search_radii(c->d, r);
    printf("count: %i\n", count);
  }
}

void cell_construct_voronoi(struct cell *c) {
  c->voronoi_active = 1;
  voronoi_init(c->v, c->d);
}

#define CVORONOI_CELL_H

#endif  // CVORONOI_CELL_H
