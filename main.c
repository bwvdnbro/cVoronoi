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
 * @file main.c
 *
 * @brief Main program entry point.
 *
 * This prototype is meant as a test and demonstration of a 2D Voronoi grid
 * construction algorithm. The algorithm constructs the Voronoi grid through
 * its dual Delaunay tessellation, using the incremental construction
 * technique described by Springel (2010).
 *
 * This prototype acts as an exploratory study for the future implementation of
 * a 2D Voronoi grid in the simulation code SWIFT.
 *
 * The prototype code can be compiled using
 *   gcc -std=gnu99 -o test main.c -lgmp -lm
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

/* Since I didn't feel like implementing my own sorting routines, I use the
   non-standard sorting function qsort_r(). The define below supresses a
   compiler warning caused by the use of this non-standard function. */
#define _GNU_SOURCE

#include "delaunay.h"
#include "hilbert.h"
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

/**
 * @brief Comparison function for two double precision floating point values.
 *
 * @param a First value.
 * @param b Second value.
 * @return -1 if a < b, 0 if a == b, +1 if a > b.
 */
int compare_double(const double a, const double b) {
  if (a < b) {
    return -1;
  } else {
    if (a > b) {
      return 1;
    } else {
      return 0;
    }
  }
}

/**
 * @brief Comparison function for two unsigned long values.
 *
 * @param a First value.
 * @param b Second value.
 * @return -1 if a < b, 0 if a == b, +1 if a > b.
 */
int compare_unsigned_long(const unsigned long a, const unsigned long b) {
  if (a < b) {
    return -1;
  } else {
    if (a > b) {
      return 1;
    } else {
      return 0;
    }
  }
}

/**
 * @brief Sorting function used to sort vertices along the horizontal direction.
 *
 * @param a First index.
 * @param b Second index.
 * @param x Vertex array to sort.
 * @return Return value of compare_double() for the x-coordinates of vertices a
 * and b.
 */
int sort_x_comp(const void *a, const void *b, void *x) {
  int ai = *(int *)a;
  int bi = *(int *)b;
  double *vertices = (double *)x;
  double ax = vertices[2 * ai];
  double bx = vertices[2 * bi];
  return compare_double(ax, bx);
}

/**
 * @brief Sorting function used to sort vertices along the vertical direction.
 *
 * @param a First index.
 * @param b Second index.
 * @param x Vertex array to sort.
 * @return Return value of compare_double() for the y-coordinates of vertices a
 * and b.
 */
int sort_y_comp(const void *a, const void *b, void *x) {
  int ai = *(int *)a;
  int bi = *(int *)b;
  double *vertices = (double *)x;
  double ay = vertices[2 * ai + 1];
  double by = vertices[2 * bi + 1];
  return compare_double(ay, by);
}

/**
 * @brief Sorting function used to sort vertices along the first diagonal.
 *
 * The first diagonal is the line x = y.
 *
 * @param a First index.
 * @param b Second index.
 * @param x Vertex array to sort.
 * @return Return value of compare_double() for the diagonal sum value (x+y) of
 * the x- and y-coordinates of vertices a and b.
 */
int sort_xyp_comp(const void *a, const void *b, void *x) {
  int ai = *(int *)a;
  int bi = *(int *)b;
  double *vertices = (double *)x;
  double ax = vertices[2 * ai];
  double ay = vertices[2 * ai + 1];
  double bx = vertices[2 * bi];
  double by = vertices[2 * bi + 1];
  return compare_double(ax + ay, bx + by);
}

/**
 * @brief Sorting function used to sort vertices along the second diagonal.
 *
 * The second diagonal is the line x = -y.
 *
 * @param a First index.
 * @param b Second index.
 * @param x Vertex array to sort.
 * @return Return value of compare_double() for the diagonal difference value
 * (x-y) of the x- and y-coordinates of vertices a and b.
 */
int sort_xym_comp(const void *a, const void *b, void *x) {
  int ai = *(int *)a;
  int bi = *(int *)b;
  double *vertices = (double *)x;
  double ax = vertices[2 * ai];
  double ay = vertices[2 * ai + 1];
  double bx = vertices[2 * bi];
  double by = vertices[2 * bi + 1];
  return compare_double(ax - ay, bx - by);
}

/**
 * @brief Sorting function used to sort vertices on their Hilbert key.
 *
 * @param a First index.
 * @param b Second index.
 * @param x Hilbert key array to sort.
 * @return Return value of compare_unsigned_long() for the hilbert key of
 * vertices a and b.
 */
int sort_h_comp(const void *a, const void *b, void *x) {
  int ai = *(int *)a;
  int bi = *(int *)b;
  unsigned long *keys = (unsigned long *)x;
  unsigned long ah = keys[ai];
  unsigned long bh = keys[bi];
  return compare_unsigned_long(ah, bh);
}

/**
 * @brief Auxiliary function used to print an arg-sorted list of vertices to a
 * file.
 *
 * @param order Array containing the indices that sort the vertex array.
 * @param vertices Vertex array.
 * @param N Size of the order and vertex arrays.
 * @param file_name Name of the output file.
 */
static inline void print_list(int *order, double *vertices, int N,
                              const char *file_name) {
  FILE *file = fopen(file_name, "w");

  for (int i = 0; i < N; ++i) {
    int index = order[i];
    fprintf(file, "%g\t%g\n", vertices[2 * index], vertices[2 * index + 1]);
  }

  fclose(file);
}

/**
 * @brief Check if the given vertex is part of a special animated path.
 *
 * @param i Vertex index.
 * @param loop Loop index.
 * @return 1 if the vertex is part of a special path, 0 otherwise.
 */
static inline int path_vertex(int i, int loop) {
  return loop > 9 && (i == 12 || i == 41 || i == 67 || i == 73);
}

/**
 * @brief Move the given vertex along the circle it is on.
 *
 * @param p Vertex coordinates.
 * @param dphi Angular direction.
 */
static inline void move_circle(double *p, double dphi) {
  double x = p[0] - 5.e9;
  double y = p[1] - 5.e9;
  double r2 = x * x + y * y;
  double r = sqrt(r2);
  double phi = atan2(y, x);
  phi += dphi;
  p[0] = r * cos(phi) + 5.e9;
  p[1] = r * sin(phi) + 5.e9;
}

/**
 * @brief Update the positions of vertices that follow a special path.
 *
 * @param loop Loop index.
 * @param vertices Vertex coordinates.
 */
static inline void update_paths(int loop, double *vertices) {
  if (loop > 9 && loop < 25) {
    vertices[2 * 12] += 1.e8;
    vertices[2 * 12 + 1] += 1.e8;

    vertices[2 * 41] -= 1.e8;
    vertices[2 * 41 + 1] += 1.e8;

    vertices[2 * 67] += 1.e8;
    vertices[2 * 67 + 1] -= 1.e8;

    vertices[2 * 73] -= 1.e8;
    vertices[2 * 73 + 1] -= 1.e8;
  } else if (loop >= 25) {
    move_circle(vertices + 2 * 12, 0.01 * M_PI);
    move_circle(vertices + 2 * 41, 0.02 * M_PI);
    move_circle(vertices + 2 * 67, -0.01 * M_PI);
    move_circle(vertices + 2 * 73, -0.02 * M_PI);
  }
}

/**
 * @brief Main program entry point.
 */
int main() {

  /* seed the random generator with the most random seed ever */
  srand(42);

  /* dimensions of the simulation box. We purposefully use large numbers to
     make sure the internal conversion routines are properly tested. */
  /*  double dim[3] = {1.e10, 1.e10, 1.e10};*/
  double dim[3] = {1., 1., 1.};

  /* define the number of vertices and set up the vertex array. */
  const int nvert = 64;
  double *vertices = (double *)malloc(2 * nvert * sizeof(double));

#ifndef REGULAR_GRID
  FILE *ifile = fopen("generators-241.txt", "r");
  for (int i = 0; i < nvert; ++i) {
    fscanf(ifile, "%lg\t%lg\n", &vertices[2 * i], &vertices[2 * i + 1]);
  }
  fclose(ifile);
  /* initialize random uniform vertex positions */
/*  for (int i = 0; i < nvert; ++i) {
    vertices[2 * i] = 1.e10 * get_random_uniform_double();
    vertices[2 * i + 1] = 1.e10 * get_random_uniform_double();
  }*/
#else
  /* initialize a regular grid of vertex positions with small random
     perturbations. */
  int index = 0;
  for (int ix = 0; ix < 20; ++ix) {
    for (int iy = 0; iy < 20; ++iy) {
      vertices[index] =
          (ix + 0.5 + 0.1 * get_random_uniform_double()) * 0.05e10;
      ++index;
      vertices[index] =
          (iy + 0.5 + 0.1 * get_random_uniform_double()) * 0.05e10;
      ++index;
    }
  }
#endif

  /* initialize the hydro space. This is a simple copy of the box dimensions,
     but then wrapped in the same way as it will be in the SWIFT code. */
  struct hydro_space hs;
  hydro_space_init(&hs, dim);

  /* set up the Hilbert key array */
  unsigned long *hkeys = (unsigned long *)malloc(nvert * sizeof(unsigned long));

  /* set up index arrays for sorting */
  int *sortx = (int *)malloc(nvert * sizeof(int));
  int *sorty = (int *)malloc(nvert * sizeof(int));
  int *sortxyp = (int *)malloc(nvert * sizeof(int));
  int *sortxym = (int *)malloc(nvert * sizeof(int));
  int *sorth = (int *)malloc(nvert * sizeof(int));
  for (int i = 0; i < nvert; ++i) {
    sortx[i] = i;
    sorty[i] = i;
    sortxyp[i] = i;
    sortxym[i] = i;
    sorth[i] = i;
  }

  for (int loop = 0; loop < 100; ++loop) {

    /* calculate/update Hilbert keys for all vertices */
    for (int i = 0; i < nvert; ++i) {
      unsigned long bits[2];
      bits[0] = vertices[2 * i] * 1.e-10 * (1ul << 32);
      bits[1] = vertices[2 * i + 1] * 1.e-10 * (1ul << 32);
      hkeys[i] = hilbert_get_key_2d(bits, 64);
    }

    /* sort the vertices. We don't actually sort the vertices themselves, but
       arg-sort them in 4 different directions: along the horizontal and
       vertical direction and along the two diagonals. This way, we can very
       easily determine which vertices are closest to a side or corner of the
       simulation box, which we will use to add ghost particles required to
       impose the periodic boundaries on the grid. */
    qsort_r(sortx, nvert, sizeof(int), sort_x_comp, vertices);
    qsort_r(sorty, nvert, sizeof(int), sort_y_comp, vertices);
    qsort_r(sortxyp, nvert, sizeof(int), sort_xyp_comp, vertices);
    qsort_r(sortxym, nvert, sizeof(int), sort_xym_comp, vertices);
    /* we also sort on Hilbert key, for faster insertion during the initial
       loop */
    qsort_r(sorth, nvert, sizeof(int), sort_h_comp, hkeys);

    /* we print out the sorted vertices to check everything workes as expected
     */
    if (loop == 0) {
      print_list(sortx, vertices, nvert, "sortx.txt");
      print_list(sorty, vertices, nvert, "sorty.txt");
      print_list(sortxyp, vertices, nvert, "sortxyp.txt");
      print_list(sortxym, vertices, nvert, "sortxym.txt");
    }

    /* initialize the delaunay tessellation structure, with initially space for
       100 vertices and 100 triangles.
       This will already create the extra vertices and triangles required to
       make the incremental construction algorithm work. */
    struct delaunay d;
    delaunay_init(&d, &hs, 100, 100);

    /* now add the vertices, one by one, in Hilbert order. */
    for (int i = 0; i < nvert; ++i) {
      int j = sorth[i];
      delaunay_add_vertex(&d, vertices[2 * j], vertices[2 * j + 1]);
    }

    /* we are done adding the original vertices. We need to consolidate the
       indices of the original vertices within the Delaunay tessellation, so
       that the tessellation knows that vertices added from now on are ghosts.
     */
    delaunay_consolidate(&d);

    /* we are now done with the local tessellation. However, we still need to
       add ghosts to impose the periodic boundaries. These ghosts will be
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
    double r = 2. * dim[0] / sqrt(nvert);
    /* add ghosts for the positive horizontal boundary */
    int i = 0;
    int vi = sortx[i];
    while (vertices[2 * vi] < r) {
      delaunay_add_vertex(&d, vertices[2 * vi] + dim[0], vertices[2 * vi + 1]);
      ++i;
      if (i == nvert) break;
      vi = sortx[i];
    }
    /* add ghosts for the negative horizontal boundary */
    i = nvert - 1;
    vi = sortx[i];
    while (dim[0] - vertices[2 * vi] < r) {
      delaunay_log("x: %g, r: %g", dim[0] - vertices[2 * vi], r);
      delaunay_add_vertex(&d, vertices[2 * vi] - dim[0], vertices[2 * vi + 1]);
      --i;
      if (i == -1) break;
      vi = sortx[i];
    }
    /* add ghosts for the positive vertical boundary */
    i = 0;
    vi = sorty[i];
    while (vertices[2 * vi + 1] < r) {
      delaunay_add_vertex(&d, vertices[2 * vi], vertices[2 * vi + 1] + dim[1]);
      ++i;
      if (i == nvert) break;
      vi = sorty[i];
    }
    /* add ghosts for the negative vertical boundary */
    i = nvert - 1;
    vi = sorty[i];
    while (dim[1] - vertices[2 * vi + 1] < r) {
      delaunay_add_vertex(&d, vertices[2 * vi], vertices[2 * vi + 1] - dim[1]);
      --i;
      if (i == -1) break;
      vi = sorty[i];
    }
    /* add ghosts for the positive x=y diagonal (top right) corner */
    i = 0;
    vi = sortxyp[i];
    while (vertices[2 * vi] + vertices[2 * vi + 1] < r) {
      delaunay_add_vertex(&d, vertices[2 * vi] + dim[0],
                          vertices[2 * vi + 1] + dim[1]);
      ++i;
      if (i == nvert) break;
      vi = sortxyp[i];
    }
    /* add ghosts for the negative x=y diagonal (bottom left) corner */
    i = nvert - 1;
    vi = sortxyp[i];
    while (dim[0] - vertices[2 * vi] + dim[1] - vertices[2 * vi + 1] < r) {
      delaunay_add_vertex(&d, vertices[2 * vi] - dim[0],
                          vertices[2 * vi + 1] - dim[1]);
      --i;
      if (i == -1) break;
      vi = sortxyp[i];
    }
    /* add ghosts for the positive x=-y diagonal (bottom right) corner */
    i = 0;
    vi = sortxym[i];
    while (vertices[2 * vi] - dim[1] + vertices[2 * vi + 1] < r) {
      delaunay_add_vertex(&d, vertices[2 * vi] + dim[0],
                          vertices[2 * vi + 1] - dim[1]);
      ++i;
      if (i == nvert) break;
      vi = sortxym[i];
    }
    /* add ghosts for the negative x=-y diagonal (top left) corner */
    i = nvert - 1;
    vi = sortxym[i];
    while (dim[0] - vertices[2 * vi] - vertices[2 * vi + 1] < r) {
      delaunay_add_vertex(&d, vertices[2 * vi] - dim[0],
                          vertices[2 * vi + 1] + dim[1]);
      --i;
      if (i == -1) break;
      vi = sortxym[i];
    }
    /* update the search radii for all triangles in the tessellation and count
       the number of triangles with circumcircles larger than the current search
       radius */
    int count = delaunay_update_search_radii(&d, r);
    printf("count: %i\n", count);
    /* now repeat the above until all circumcircles are smaller than the search
       radius */
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
      vi = sortx[i];
      while ((vertices[2 * vi] >= old_r) && (vertices[2 * vi] < r)) {
        delaunay_add_vertex(&d, vertices[2 * vi] + dim[0],
                            vertices[2 * vi + 1]);
        ++i;
        if (i == nvert) break;
        vi = sortx[i];
      }
      i = nvert - 1;
      vi = sortx[i];
      while ((dim[0] - vertices[2 * vi] >= old_r) &&
             (dim[0] - vertices[2 * vi] < r)) {
        delaunay_log("x: %g, old_r: %g, r: %g", dim[0] - vertices[2 * vi],
                     old_r, r);
        delaunay_add_vertex(&d, vertices[2 * vi] - dim[0],
                            vertices[2 * vi + 1]);
        --i;
        if (i == -1) break;
        vi = sortx[i];
      }
      i = 0;
      vi = sorty[i];
      while ((vertices[2 * vi + 1] >= old_r) && (vertices[2 * vi + 1] < r)) {
        delaunay_add_vertex(&d, vertices[2 * vi],
                            vertices[2 * vi + 1] + dim[1]);
        ++i;
        if (i == nvert) break;
        vi = sorty[i];
      }
      i = nvert - 1;
      vi = sorty[i];
      while ((dim[1] - vertices[2 * vi + 1] >= old_r) &&
             (dim[1] - vertices[2 * vi + 1] < r)) {
        delaunay_add_vertex(&d, vertices[2 * vi],
                            vertices[2 * vi + 1] - dim[1]);
        --i;
        if (i == -1) break;
        vi = sorty[i];
      }
      i = 0;
      vi = sortxyp[i];
      while ((vertices[2 * vi] + vertices[2 * vi + 1] >= old_r) &&
             (vertices[2 * vi] + vertices[2 * vi + 1] < r)) {
        delaunay_add_vertex(&d, vertices[2 * vi] + dim[0],
                            vertices[2 * vi + 1] + dim[1]);
        ++i;
        if (i == nvert) break;
        vi = sortxyp[i];
      }
      i = nvert - 1;
      vi = sortxyp[i];
      while ((dim[0] - vertices[2 * vi] + dim[1] - vertices[2 * vi + 1] >=
              old_r) &&
             (dim[0] - vertices[2 * vi] + dim[1] - vertices[2 * vi + 1] < r)) {
        delaunay_add_vertex(&d, vertices[2 * vi] - dim[0],
                            vertices[2 * vi + 1] - dim[1]);
        --i;
        if (i == -1) break;
        vi = sortxyp[i];
      }
      i = 0;
      vi = sortxym[i];
      while ((vertices[2 * vi] - dim[1] + vertices[2 * vi + 1] >= old_r) &&
             (vertices[2 * vi] - dim[1] + vertices[2 * vi + 1] < r)) {
        delaunay_add_vertex(&d, vertices[2 * vi] + dim[0],
                            vertices[2 * vi + 1] - dim[1]);
        ++i;
        if (i == nvert) break;
        vi = sortxym[i];
      }
      i = nvert - 1;
      vi = sortxym[i];
      while ((dim[0] - vertices[2 * vi] - vertices[2 * vi + 1] >= old_r) &&
             (dim[0] - vertices[2 * vi] - vertices[2 * vi + 1] < r)) {
        delaunay_add_vertex(&d, vertices[2 * vi] - dim[0],
                            vertices[2 * vi + 1] + dim[1]);
        --i;
        if (i == -1) break;
        vi = sortxym[i];
      }
      /* update the search radii to the new value and count the number of larger
         circumcircles. */
      count = delaunay_update_search_radii(&d, r);
      printf("count: %i\n", count);
    }

    /* the Delaunay tessellation is now complete. Print it out for visual
       inspection. */
    char filename[50];
    sprintf(filename, "test%03i.txt", loop);
    delaunay_print_tessellation(&d, filename);

    /* Convert the Delaunay tessellation into a Voronoi grid. */
    struct voronoi v;
    voronoi_init(&v, &d);

    /* Get rid of the Delaunay tessellation. */
    delaunay_destroy(&d);

    /* sanity check(s) on the Voronoi grid */
    voronoi_check_grid(&v);

    /* Now print the Voronoi grid for visual inspection. */
    sprintf(filename, "vtest%03i.txt", loop);
    voronoi_print_grid(&v, filename);

    return 0;

    /* apply Lloyd's algorithm to regularise the grid */
    for (int i = 0; i < nvert; ++i) {
      int j = sorth[i];
      if (!path_vertex(j, loop)) {
        vertices[2 * j] = v.cell_centroid[2 * i];
        vertices[2 * j + 1] = v.cell_centroid[2 * i + 1];
        /* we also add a small circular movement to the grid */
        move_circle(vertices + 2 * j, 0.004 * M_PI);
      }
    }
    /* update the positions of vertices on special paths */
    update_paths(loop, vertices);
    /* impose periodic boundaries (we apply this after moving vertices on
       special paths to make sure these also satisfy the periodic boundaries) */
    for (int i = 0; i < nvert; ++i) {
      if (vertices[2 * i] < 0.) {
        vertices[2 * i] += dim[0];
      }
      if (vertices[2 * i] >= dim[0]) {
        vertices[2 * i] -= dim[0];
      }
      if (vertices[2 * i + 1] < 0.) {
        vertices[2 * i + 1] += dim[1];
      }
      if (vertices[2 * i + 1] >= dim[1]) {
        vertices[2 * i + 1] -= dim[1];
      }
    }

    /* Get rid of the Voronoi grid. */
    voronoi_destroy(&v);
  }

  /* Get rid of the sort arrays and the vertex array. */
  free(sortx);
  free(sorty);
  free(sortxyp);
  free(sortxym);
  free(sorth);
  free(hkeys);
  free(vertices);

  /* We are done. */
  return 0;
}
