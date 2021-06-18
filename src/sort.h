//
// Created by yuyttenh on 07/06/2021.
//

#ifndef CVORONOI_SORT_H

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
 * @brief Sorting function used to sort vertex_indices along the horizontal direction.
 *
 * @param a First index.
 * @param b Second index.
 * @param x Vertex array to sort.
 * @return Return value of compare_double() for the x-coordinates of vertex_indices a
 * and b.
 */
int sort_x_comp(const void *a, const void *b, void *x) {
  int ai = *(int *)a;
  int bi = *(int *)b;
  double *vertices = (double *)x;
  double ax = vertices[3 * ai];
  double bx = vertices[3 * bi];
  return compare_double(ax, bx);
}

/**
 * @brief Sorting function used to sort vertex_indices along the vertical direction.
 *
 * @param a First index.
 * @param b Second index.
 * @param x Vertex array to sort.
 * @return Return value of compare_double() for the y-coordinates of vertex_indices a
 * and b.
 */
int sort_y_comp(const void *a, const void *b, void *x) {
  int ai = *(int *)a;
  int bi = *(int *)b;
  double *vertices = (double *)x;
  double ay = vertices[3 * ai + 1];
  double by = vertices[3 * bi + 1];
  return compare_double(ay, by);
}

/**
 * @brief Sorting function used to sort vertex_indices along the first diagonal.
 *
 * The first diagonal is the line x = y.
 *
 * @param a First index.
 * @param b Second index.
 * @param x Vertex array to sort.
 * @return Return value of compare_double() for the diagonal sum value (x+y) of
 * the x- and y-coordinates of vertex_indices a and b.
 */
int sort_xyp_comp(const void *a, const void *b, void *x) {
  int ai = *(int *)a;
  int bi = *(int *)b;
  double *vertices = (double *)x;
  double ax = vertices[3 * ai];
  double ay = vertices[3 * ai + 1];
  double bx = vertices[3 * bi];
  double by = vertices[3 * bi + 1];
  return compare_double(ax + ay, bx + by);
}

/**
 * @brief Sorting function used to sort vertex_indices along the second diagonal.
 *
 * The second diagonal is the line x = -y.
 *
 * @param a First index.
 * @param b Second index.
 * @param x Vertex array to sort.
 * @return Return value of compare_double() for the diagonal difference value
 * (x-y) of the x- and y-coordinates of vertex_indices a and b.
 */
int sort_xym_comp(const void *a, const void *b, void *x) {
  int ai = *(int *)a;
  int bi = *(int *)b;
  double *vertices = (double *)x;
  double ax = vertices[3 * ai];
  double ay = vertices[3 * ai + 1];
  double bx = vertices[3 * bi];
  double by = vertices[3 * bi + 1];
  return compare_double(ax - ay, bx - by);
}

/**
 * @brief Sorting function used to sort vertex_indices on their Hilbert key.
 *
 * @param a First index.
 * @param b Second index.
 * @param x Hilbert key array to sort.
 * @return Return value of compare_unsigned_long() for the hilbert key of
 * vertex_indices a and b.
 */
int sort_h_comp(const void *a, const void *b, void *x) {
  int ai = *(int *)a;
  int bi = *(int *)b;
  unsigned long *keys = (unsigned long *)x;
  unsigned long ah = keys[ai];
  unsigned long bh = keys[bi];
  return compare_unsigned_long(ah, bh);
}

#define CVORONOI_SORT_H

#endif  // CVORONOI_SORT_H
