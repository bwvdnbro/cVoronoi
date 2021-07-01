//
// Created by yuyttenh on 07/06/2021.
//

/**
 * @file cell.h
 *
 * @brief Simple cell struct which stores a Delaunay tesselation and the
 * corresponding Voronoi mesh.
 *
 * Emulates some of the basic functionality of a SWIFT cell. Can relax the
 * Voronoi mesh.
 */

#ifndef CVORONOI_CELL_H
#define CVORONOI_CELL_H

#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#include "delaunay.h"
#include "dimensionality.h"
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

/*! @brief Struct containing the necessary information to build a delaunay and
 * vornonoi tesselation for a given set of vertices.
 */
struct cell {
  /*! @brief Number of vertices in cell */
  int count;

  /*! @brief Array of vertices */
  double *vertices;

  /*! @brief Array of hilbert keys of vertices */
  unsigned long *hilbert_keys;

  /*! @brief Arg-sort indices for directions: x, y, xyp and xym and the hilbert
   * key */
  int *r_sort_lists[5];

  /*! @brief Simulation volume */
  struct hydro_space hs;

  /*! @brief Delaunay triangulation, empty upon initialization */
  struct delaunay d;

  /*! @brief Voronoi tesselation, empty upon initialization */
  struct voronoi v;

  /*! @brief Flag indication whether or not the current cell has already
   * constructed its voronoi tesselation */
  int voronoi_active;
};

/*! @brief Calculate the hilbert keys of the vertices
 *
 * @param c Cell containing the vertices
 */
static inline void cell_update_hilbert_keys(struct cell *c) {
  for (int i = 0; i < c->count; i++) {
#if defined(DIMENSIONALITY_2D)
    unsigned long bits[2];
    int nbits = 32;
    bits[0] =
        (c->vertices[3 * i] - c->hs.anchor[0]) / c->hs.side[0] * (1ul << nbits);
    bits[1] = (c->vertices[3 * i + 1] - c->hs.anchor[1]) / c->hs.side[1] *
              (1ul << nbits);
#else
    unsigned long bits[3];
    int nbits = 21;
    bits[0] =
        (c->vertices[3 * i] - c->hs.anchor[0]) / c->hs.side[0] * (1ul << nbits);
    bits[1] = (c->vertices[3 * i + 1] - c->hs.anchor[1]) / c->hs.side[1] *
              (1ul << nbits);
    bits[2] = (c->vertices[3 * i + 2] - c->hs.anchor[2]) / c->hs.side[2] *
              (1ul << nbits);
#endif
    c->hilbert_keys[i] = hilbert_get_key(bits, nbits);
  }
}

/*! @brief Update the arg-sort arrays in the various directions
 *
 * @param c Cell containing the vertices to be sorted.
 */
static inline void cell_update_sorts(struct cell *c) {
  qsort_r(c->r_sort_lists[0], c->count, sizeof(int), sort_x_comp, c->vertices);
  qsort_r(c->r_sort_lists[1], c->count, sizeof(int), sort_y_comp, c->vertices);
  qsort_r(c->r_sort_lists[2], c->count, sizeof(int), sort_xyp_comp,
          c->vertices);
  qsort_r(c->r_sort_lists[3], c->count, sizeof(int), sort_xym_comp,
          c->vertices);
  qsort_r(c->r_sort_lists[4], c->count, sizeof(int), sort_h_comp, c->vertices);
}

/*! @brief Initialize a new cell with slightly randomized vertices
 *
 * @param c Pointer to cell to be initialized
 * @param count Number of vertices to add to cell
 * @param pert Relative scale of the perturbations
 * @param dim Dimensions of the simulation volume
 */
static inline void cell_init(struct cell *c, const int *count,
                             const double pert, const double *dim) {
  hydro_space_init(&c->hs, dim);
  c->count = count[0] * count[1] * count[2];

  /* slightly randomized vertices */
  c->vertices = (double *)malloc(3 * c->count * sizeof(double));
  int index = 0;
  for (int ix = 0; ix < count[0]; ++ix) {
    for (int iy = 0; iy < count[1]; ++iy) {
      for (int iz = 0; iz < count[2]; ++iz) {
        c->vertices[index] =
            (ix + 0.5 + pert * (get_random_uniform_double() - 0.5)) *
            (dim[0] / (double)count[0]);
        ++index;
        c->vertices[index] =
            (iy + 0.5 + pert * (get_random_uniform_double() - 0.5)) *
            (dim[1] / (double)count[1]);
        ++index;
        c->vertices[index] =
            (iz + 0.5 + pert * (get_random_uniform_double() - 0.5)) *
            (dim[2] / (double)count[2]);
        ++index;
      }
    }
  }

  /* hilbert keys */
  c->hilbert_keys = (unsigned long *)malloc(c->count * sizeof(unsigned long));
  cell_update_hilbert_keys(c);

  /* sorting arrays */
  for (int i = 0; i < 5; i++) {
    c->r_sort_lists[i] = (int *)malloc(c->count * sizeof(int));
    for (int j = 0; j < c->count; j++) {
      c->r_sort_lists[i][j] = j;
    }
  }
  cell_update_sorts(c);

  delaunay_init(&c->d, &c->hs, c->count, 10 * c->count);
  c->voronoi_active = 0;
}

/*! @brief Clean up cell
 *
 * @param c pointer to cell to be freed
 */
static inline void cell_destroy(struct cell *c) {
  free(c->vertices);
  free(c->hilbert_keys);
  for (int i = 0; i < 5; i++) {
    free(c->r_sort_lists[i]);
  }
  delaunay_destroy(&c->d);
  if (c->voronoi_active) {
    voronoi_destroy(&c->v);
  }
}

/*! @brief Construct the delaunay triangulation of all the local vertices (no
 * periodic boundaries).
 *
 * @param c Pointer to cell containing the vertices to add to the delaunay
 * triangulation.
 */
static inline void cell_construct_local_delaunay(struct cell *c) {
  /* Add the local vertices, one by one, in Hilbert order. */
  for (int i = 0; i < c->count; ++i) {
    int j = c->r_sort_lists[4][i];
    delaunay_add_local_vertex(&c->d, j, c->vertices[3 * j],
                              c->vertices[3 * j + 1], c->vertices[3 * j + 2]);
  }

  /* we are done adding the original vertices. We need to consolidate the
     indices of the original vertices within the Delaunay tessellation, so
     that the tessellation knows that vertices added from now on are ghosts.
   */
  delaunay_consolidate(&c->d);
}

/*! @brief Impose periodic boundaries by adding the necessary ghost vertices
 *
 * Add ghosts to impose the periodic boundaries. These ghosts will be
 * periodic copies of the original vertices that ensure that the incomplete
 * cells at the boundaries of the simulation box have the right shape.
 * Within SWIFT, a similar technique needs to be used to guarantee that
 * cells near the boundary of a SWIFT-cell correctly "feel" the cells in the
 * neighbouring SWIFT-cells. Springel (2010) provides a criterion for
 * completeness, which is based on the radius of the circumcircles of the
 * triangles that connect to original vertices. Starting from an initial
 * search radius, we will iteratively add ghost vertices and increase this
 * search radius until all relevant circumcircle radii are smaller than the
 * search radius. This mechanism can very easily be combined with SWIFT's
 * existing neighbour search algorithms.

 * @param c Cell containing the consolidated delaunay triangulation.
 */
static inline void cell_make_delaunay_periodic(struct cell *c) {
#ifdef DIMENSIONALITY_2D
  /* Initial search radius. We use twice the average inter-particle
     separation, but other values would also work. */
  double r = 2. * c->hs.side[0] * c->hs.side[1] / c->count;
  double old_r = -DBL_MAX;
  double sqrt2 = sqrt(2.);

  /* update the search radii for all triangles in the tessellation and count
     the number of triangles with circumcircles larger than the current search
     radius */
  int count = delaunay_update_search_radii(&c->d, old_r);
//  printf("count: %i\n", count);
  while (count > 0) {
    /* add ghosts for the positive horizontal boundary */
    int i = 0;
    int vi = c->r_sort_lists[0][i];
    while (c->vertices[3 * vi] < r) {
      if (i == c->count) break;
      ++i;
      if (c->vertices[3 * vi] < old_r) {
        vi = c->r_sort_lists[0][i];
        continue;
      }
      delaunay_add_new_vertex(
          &c->d, c->vertices[3 * vi] + c->hs.anchor[0] + c->hs.side[0],
          c->vertices[3 * vi + 1]);
      vi = c->r_sort_lists[0][i];
    }
    /* add ghosts for the negative horizontal boundary */
    i = c->count - 1;
    vi = c->r_sort_lists[0][i];
    while (c->hs.anchor[0] + c->hs.side[0] - c->vertices[3 * vi] < r) {
      if (i == -1) break;
      --i;
      if (c->hs.anchor[0] + c->hs.side[0] - c->vertices[3 * vi] < old_r) {
        vi = c->r_sort_lists[0][i];
        continue;
      }
      delaunay_log("x: %g, old_r: %g, r: %g",
                   c->hs.anchor[0] + c->hs.side[0] - c->vertices[3 * vi], old_r,
                   r);
      delaunay_add_new_vertex(
          &c->d, c->vertices[3 * vi] - (c->hs.anchor[0] + c->hs.side[0]),
          c->vertices[3 * vi + 1]);
      vi = c->r_sort_lists[0][i];
    }
    /* add ghosts for the positive vertical boundary */
    i = 0;
    vi = c->r_sort_lists[1][i];
    while (c->vertices[3 * vi + 1] < r) {
      if (i == c->count) break;
      ++i;
      if (c->vertices[3 * vi + 1] < old_r) {
        vi = c->r_sort_lists[1][i];
        continue;
      }
      delaunay_add_new_vertex(
          &c->d, c->vertices[3 * vi],
          c->vertices[3 * vi + 1] + c->hs.anchor[1] + c->hs.side[1]);
      vi = c->r_sort_lists[1][i];
    }
    /* add ghosts for the negative vertical boundary */
    i = c->count - 1;
    vi = c->r_sort_lists[1][i];
    while (c->hs.anchor[1] + c->hs.side[1] - c->vertices[3 * vi + 1] < r) {
      if (i == -1) break;
      --i;
      if (c->hs.anchor[1] + c->hs.side[1] - c->vertices[3 * vi + 1] < old_r) {
        vi = c->r_sort_lists[1][i];
        continue;
      }
      delaunay_add_new_vertex(
          &c->d, c->vertices[3 * vi],
          c->vertices[3 * vi + 1] - (c->hs.anchor[1] + c->hs.side[1]));
      vi = c->r_sort_lists[1][i];
    }
    /* add ghosts for the positive x=y diagonal (top right) corner */
    i = 0;
    vi = c->r_sort_lists[2][i];
    while (c->vertices[3 * vi] + c->vertices[3 * vi + 1] < r * sqrt2) {
      if (i == c->count) break;
      ++i;
      if (c->vertices[3 * vi] + c->vertices[3 * vi + 1] < old_r * sqrt2) {
        vi = c->r_sort_lists[2][i];
        continue;
      }
      delaunay_add_new_vertex(
          &c->d, c->vertices[3 * vi] + c->hs.anchor[0] + c->hs.side[0],
          c->vertices[3 * vi + 1] + c->hs.anchor[1] + c->hs.side[1]);
      vi = c->r_sort_lists[2][i];
    }
    /* add ghosts for the negative x=y diagonal (bottom left) corner */
    i = c->count - 1;
    vi = c->r_sort_lists[2][i];
    while (c->hs.anchor[0] + c->hs.side[0] - c->vertices[3 * vi] +
               c->hs.anchor[1] + c->hs.side[1] - c->vertices[3 * vi + 1] <
           r * sqrt2) {
      if (i == -1) break;
      --i;
      if (c->hs.anchor[0] + c->hs.side[0] - c->vertices[3 * vi] +
              c->hs.anchor[1] + c->hs.side[1] - c->vertices[3 * vi + 1] <
          old_r * sqrt2) {
        vi = c->r_sort_lists[2][i];
        continue;
      }
      delaunay_add_new_vertex(
          &c->d, c->vertices[3 * vi] - (c->hs.anchor[0] + c->hs.side[0]),
          c->vertices[3 * vi + 1] - (c->hs.anchor[1] + c->hs.side[1]));
      vi = c->r_sort_lists[2][i];
    }
    /* add ghosts for the positive x=-y diagonal (bottom right) corner */
    i = 0;
    vi = c->r_sort_lists[3][i];
    while (c->vertices[3 * vi] - (c->hs.anchor[1] + c->hs.side[1]) +
               c->vertices[3 * vi + 1] <
           r * sqrt2) {
      if (i == c->count) break;
      ++i;
      if (c->vertices[3 * vi] - (c->hs.anchor[1] + c->hs.side[1]) +
              c->vertices[3 * vi + 1] <
          old_r * sqrt2) {
        vi = c->r_sort_lists[3][i];
        continue;
      }
      delaunay_add_new_vertex(
          &c->d, c->vertices[3 * vi] + c->hs.anchor[0] + c->hs.side[0],
          c->vertices[3 * vi + 1] - (c->hs.anchor[1] + c->hs.side[1]));
      vi = c->r_sort_lists[3][i];
    }
    /* add ghosts for the negative x=-y diagonal (top left) corner */
    i = c->count - 1;
    vi = c->r_sort_lists[3][i];
    while (c->hs.anchor[0] + c->hs.side[0] - c->vertices[3 * vi] -
               c->vertices[3 * vi + 1] <
           r * sqrt2) {
      if (i == -1) break;
      --i;
      if (c->hs.anchor[0] + c->hs.side[0] - c->vertices[3 * vi] -
              c->vertices[3 * vi + 1] <
          old_r * sqrt2) {
        vi = c->r_sort_lists[3][i];
        continue;
      }
      delaunay_add_new_vertex(
          &c->d, c->vertices[3 * vi] - (c->hs.anchor[0] + c->hs.side[0]),
          c->vertices[3 * vi + 1] + c->hs.anchor[1] + c->hs.side[1]);
      vi = c->r_sort_lists[3][i];
    }
    /* update the search radii to the new value and count the number of larger
       circumcircles. */
    count = delaunay_update_search_radii(&c->d, r);
//    printf("count: %i\n", count);

    /* we do not want to add the same ghost twice (this causes the incremental
       construction algorithm to crash), so we need to keep track of the
       previous search radius, so that we can only add the ghosts that have
       not been added before */
    old_r = r;
    /* now gradually increase the search radius and repeat */
    r *= 1.25;
  }
#else
  for (int i = 0; i < 27; i++) {
    if (i == 13) continue; /* skip (0, 0, 0) */
    double shiftx = ((i % 3) - 1.) * c->hs.side[0];
    double shifty = ((i / 3) % 3 - 1.) * c->hs.side[1];
    double shiftz = ((i / 9) % 3 - 1.) * c->hs.side[2];

    for (int k = 0; k < c->count; ++k) {
      int l = c->r_sort_lists[4][k];
      delaunay_add_new_vertex(&c->d, c->vertices[3 * l] + shiftx,
                              c->vertices[3 * l + 1] + shifty,
                              c->vertices[3 * l + 2] + shiftz);
    }
  }
#ifdef DELAUNAY_CHECKS
  for (int i = c->d.vertex_start; i < c->d.vertex_end; i++){
    double search_radius = delaunay_get_search_radius(&c->d, i);
    double radius = delaunay_get_radius(&c->d, c->d.vertex_tetrahedron_links[i]);
    delaunay_assert(search_radius >= 2. * radius);
    delaunay_assert(search_radius < c->hs.side[0]);
  }
#endif
#endif
}

/*! @brief Construct the voronoi grid from this cells delaunay triangulation
 *
 * @param c The cell containing the delaunay triangulation
 */
static inline void cell_construct_voronoi(struct cell *c) {
  c->voronoi_active = 1;
  voronoi_init(&c->v, &c->d);
}

/*! @brief Relax this cells vertices by moving them to the centroids of their
 * corresponding voronoi faces (Lloyds relaxation).
 *
 * After moving the vertices, the hilbert keys and sort lists are updated and
 * the tessellations are rebuild.
 *
 * @param c The cell containing the voronoi tessellation
 */
static inline void cell_lloyd_relax_vertices(struct cell *c) {
  if (!c->voronoi_active) {
    fprintf(stderr, "Voronoi tesselation is uninitialized!\n");
    abort();
  }

  for (int i = 0; i < c->count; ++i) {
    int j = c->r_sort_lists[4][i];
#if defined(DIMENSIONALITY_2D)
    c->vertices[3 * j] = c->v.cells[i].centroid[0];
    c->vertices[3 * j + 1] = c->v.cells[i].centroid[1];
#else
    c->vertices[3 * j] = c->v.cells[i].centroid[0];
    c->vertices[3 * j + 1] = c->v.cells[i].centroid[1];
    c->vertices[3 * j + 2] = c->v.cells[i].centroid[2];
#endif
  }
  /* Destroy existing tesselations */
  delaunay_destroy(&c->d);
  voronoi_destroy(&c->v);
  /* Update sorts */
  cell_update_hilbert_keys(c);
  cell_update_sorts(c);
  /* Rebuild tesselations */
  delaunay_init(&c->d, &c->hs, c->count, 10 * c->count);
  cell_construct_local_delaunay(c);
  cell_make_delaunay_periodic(c);
  cell_construct_voronoi(c);
}

/*! @brief Print this cells voronoi and delaunay tessellations
 *
 * @param c The cell containing the tessellations
 * @param vor_file_name Filename to write the voronoi tessellation to
 * @param del_file_name Filename to write the delaunay tessellation to
 */
static inline void cell_print_tesselations(const struct cell *c,
                                           const char *vor_file_name,
                                           const char *del_file_name) {
  if (!c->voronoi_active) {
    fprintf(stderr, "Voronoi tesselation is uninitialized!\n");
    abort();
  }

  voronoi_print_grid(&c->v, vor_file_name);
  delaunay_print_tessellation(&c->d, del_file_name);
}

#endif  // CVORONOI_CELL_H
