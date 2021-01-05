/*******************************************************************************
 * This file is part of cVoronoi.
 * Copyright (c) 2021 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file test_hilbert.c
 *
 * @brief Test for the Hilbert key routines.
 *
 * This test can be compiled using
 *   gcc -std=gnu99 -o test_hilbert test_hilbert.c
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#define _GNU_SOURCE

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "hilbert.h"

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
 * @brief Test for the Hilbert key functions.
 *
 * This program takes one optional extra command line argument: the number of
 * points to randomly generate (default: 100). It then moves on to randomly
 * generate these points in 2D, computes a Hilbert key for each point, arg-sorts
 * the points on this key and outputs the point coordinates and their keys to
 * the standard output.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 */
int main(int argc, char **argv) {

  int nvert = 100;
  if (argc > 1) {
    nvert = atoi(argv[1]);
  }

  double *x = (double *)malloc(nvert * 2 * sizeof(double));
  for (int i = 0; i < 2 * nvert; ++i) {
    x[i] = get_random_uniform_double();
  }

  unsigned long *keys = (unsigned long *)malloc(nvert * sizeof(unsigned long));
  for (int i = 0; i < nvert; ++i) {
    unsigned long bits[2];
    bits[0] = x[2 * i] * (1ul << 32);
    bits[1] = x[2 * i + 1] * (1ul << 32);
    keys[i] = hilbert_get_key_2d(bits, 64);
  }

  int *sorth = (int *)malloc(nvert * sizeof(int));
  for (int i = 0; i < nvert; ++i) {
    sorth[i] = i;
  }

  qsort_r(sorth, nvert, sizeof(int), sort_h_comp, keys);

  for (int i = 0; i < nvert; ++i) {
    int j = sorth[i];
    printf("%.3f\t%.3f\t%lu\n", x[2 * j], x[2 * j + 1], keys[j]);
  }

  return 0;
}
