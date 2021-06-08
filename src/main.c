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

#include "cell.h"

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
  int count[3] = {10, 10, 1};
  double dim[3] = {2., 2., 2.};
  struct cell c;
  cell_init(&c, count, 1., dim);
  cell_construct_local_delaunay(&c);
  cell_make_delaunay_periodic(&c);
  cell_construct_voronoi(&c);

  /* Now print the Voronoi grid for visual inspection. */
  char vor_filename[50];
  char del_filename[50];
  sprintf(vor_filename, "vtest.txt");
  sprintf(del_filename, "test.txt");
  cell_print_tesselations(&c, vor_filename, del_filename);

  /* Lloyd's relaxation */
  for (int loop = 1; loop <= 5; ++loop) {
    printf("Relaxation loop %i\n", loop);
    cell_lloyd_relax_vertices(&c);
    sprintf(vor_filename, "vtest%03i.txt", loop);
    sprintf(del_filename, "test%03i.txt", loop);
    cell_print_tesselations(&c, vor_filename, del_filename);
  }

  /* cleanup */
  cell_destroy(&c);
  return 0;
}
