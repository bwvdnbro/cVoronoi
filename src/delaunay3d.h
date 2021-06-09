//
// Created by yuyttenh on 08/06/2021.
//

#ifndef CVORONOI_DELAUNAY3D_H
#define CVORONOI_DELAUNAY3D_H

#include "geometry.h"
#include "hydro_space.h"
#include "tetrahedron.h"

struct delaunay {

  /*! @brief Anchor of the simulation volume. */
  double anchor[3];

  /*! @brief Inverse side length of the simulation volume. */
  double inverse_side;

  /*! @brief Vertex positions. This array is a copy of the array defined in
   *  main() and we probably want to get rid of it in a SWIFT implementation. */
  double* vertices;

#ifdef DELAUNAY_NONEXACT
  /*! @brief Vertex positions, rescaled to the range 1-2. Kept here in case we
   *  want to adopt hybrid geometrical checks (floating point checks when safe,
   *  integer checks when there is a risk of numerical error leading to the
   *  wrong result) to speed things up. */
  double* rescaled_vertices;
#endif

  /*! @brief Integer vertices. These are the vertex coordinates that are
   *  actually used during the incremental construction. */
  unsigned long int* integer_vertices;

  /*! @brief Vertex search radii. For every vertex, this array contains twice
   *  the radius of the largest circumsphere of the tetrahedra that vertex is
   *  part of. */
  double* search_radii;

  /*! @brief Next available index within the vertex array. Corresponds to the
   *  actual size of the vertex array. */
  int vertex_index;

  /*! @brief Current size of the vertex array in memory. If vertex_size matches
   *  vertex_index, the memory buffer is full and needs to be expanded. */
  int vertex_size;

  /*! @brief Begin index of the normal vertices. This skips the 3 auxiliary
   *  vertices required for the incremental construction algorithm. */
  int vertex_start;

  /*! @brief End index of the normal vertices. This variable is set by calling
   *  delaunay_consolidate() and contains the offset of the ghost vertices
   *  within the vertex array. */
  int vertex_end;

  /*! @brief Triangles that make up the tessellation. */
  struct tetrahedron* tetrahedra;

  /*! @brief Next available index within the triangle array. Corresponds to the
   *  actual size of the triangle array. */
  int tetrahedron_index;

  /*! @brief Current size of the triangle array in memory. If triangle_size
   *  matches triangle_index, the memory buffer is full and needs to be
   *  expanded. */
  int tetrahedron_size;

  /*! @brief Lifo queue of tetrahedra that need checking during the incremental
   *  construction algorithm. After a new vertex has been added, all new
   *  tetrahedra are added to this queue and are tested to see if the
   *  Delaunay criterion (empty circumcircles) still holds. New triangles
   *  created when flipping invalid triangles are also added to the queue. */
  int* tetrahedron_queue;

  /*! @brief Next available index in the tetrahedron_queue. Determines both the
   * actual size of the tetrahedron_queue as the first element that will be
   * popped from the tetrahedron_queue. */
  int tetrahedron_q_index;

  /*! @brief Current size of the tetrahedron_queue in memory. If
   * tetrahedron_q_size matches tetrahedron_q_index, the memory buffer
   * is full and needs to be expanded. */
  int tetrahedron_q_size;

  /*! @brief Lifo queue of free spots in the tetrahedra array. Sometimes 3
   * tetrahedra can be split into 2 new ones. This leaves a free spot in the
   * array.
   */
  int* free_indices_queue;

  /*! @brief Next available index in de free_indices_queue. Determines both the
   * actual size of the free_indices_queue as the first element that will be
   * popped from the free_indices_queue. */
  int free_indices_q_index;

  /*! @brief Current size of the free_indices_queue in memory. If
   * free_indices_q_size matches free_indices_q_index, the memory buffer
   * is full and needs to be expanded. */
  int free_indices_q_size;

  /*! @brief Array of tetrahedra containing the current vertex */
  int* tetrahedra_containing_vertex;

  /*! @brief Next index in the array of tetrahedra containing the current
   * vertex. This is also the number of tetrahedra containing the current
   * vertex. */
  int tetrahedra_containing_vertex_index;

  /*! @brief Size of the array of tetrahedra containing the current vertex in
   * memory. */
  int tetrahedra_containing_vertex_size;

  /*! @brief Index of the last triangle that was accessed. Used as initial
   *  guess for the triangle that contains the next vertex that will be added.
   *  If vertices are added in some sensible order (e.g. in Peano-Hilbert curve
   *  order) then this will greatly speed up the algorithm. */
  int last_tetrahedron;

  /*! @brief Geometry variables. Auxiliary variables used by the exact integer
   *  geometry tests that need to be stored in between tests, since allocating
   *  and deallocating them for every test is too expensive. */
  struct geometry geometry;
};

/* Forward declarations */
inline static void delaunay_check_tessellation(struct delaunay* d);
inline static int delaunay_new_vertex(struct delaunay* restrict d, double x,
                                      double y, double z);
inline static void delaunay_add_vertex(struct delaunay* restrict d, int v);
inline static int delaunay_new_tetrahedron(struct delaunay* restrict d);
inline static int delaunay_free_indices_queue_pop(struct delaunay* restrict d);
inline static int delaunay_tetrahedron_queue_pop(struct delaunay* restrict d);
inline static void delaunay_append_tetrahedron_containing_vertex(
    struct delaunay* d, int t);
inline static int delaunay_find_tetrahedra_containing_vertex(struct delaunay* d,
                                                             const int v);
inline static void delaunay_one_to_four_flip(struct delaunay* d, int v, int t);
inline static void delaunay_two_to_six_flip(struct delaunay* d, int v, int* t);
inline static void delaunay_n_to_2n_flip(struct delaunay* d, int v, int* t,
                                         int n);
inline static void delaunay_check_tetrahedra(struct delaunay* d);
inline static void delaunay_check_tetrahedron(struct delaunay* d, int t);

/**
 * @brief Initialize the Delaunay tessellation.
 *
 * This function allocates memory for all arrays that make up the tessellation
 * and initializes the variables used for bookkeeping.
 *
 * It then sets up a large tetrahedron that contains the entire simulation box
 * and additional buffer space to deal with boundary ghost vertices, and 4
 * additional dummy tetrahedron that provide valid neighbours for the 4 sides of
 * this tetrahedron (these dummy tetrahedra themselves have an invalid tip
 * vertex and are therefore simply placeholders).
 *
 * @param d Delaunay tesselation.
 * @param hs Spatial extents of the simulation box.
 * @param vertex_size Initial size of the vertex array.
 * @param tetrahedron_size Initial size of the tetrahedra array.
 */
inline static void delaunay_init(struct delaunay* restrict d,
                                 const struct hydro_space* restrict hs,
                                 int vertex_size, int tetrahedron_size) {
  /* allocate memory for the vertex arrays */
  d->vertices = (double*)malloc(vertex_size * 3 * sizeof(double));
#ifdef DELAUNAY_NONEXACT
  d->rescaled_vertices = (double*)malloc(vertex_size * 3 * sizeof(double));
#endif
  d->integer_vertices =
      (unsigned long int*)malloc(vertex_size * 3 * sizeof(unsigned long int));
  d->search_radii = (double*)malloc(vertex_size * sizeof(double));
  d->vertex_size = vertex_size;
  /* set vertex start and end (indicating where the local vertices start and
   * end)*/
  d->vertex_start = 0;
  d->vertex_end = vertex_size;
  /* we add the dummy vertices behind the local vertices and before the ghost
   * vertices (see below) */
  d->vertex_index = vertex_size;

  /* allocate memory for the tetrahedra array */
  d->tetrahedra = (struct tetrahedron*)malloc(tetrahedron_size *
                                              sizeof(struct tetrahedron));
  d->tetrahedron_index = 0;
  d->tetrahedron_size = tetrahedron_size;

  /* allocate memory for the tetrahedron_queue (note that the tetrahedron_queue
     size of 10 was chosen arbitrarily, and a proper value should be chosen
     based on performance measurements) */
  d->tetrahedron_queue = (int*)malloc(10 * sizeof(int));
  d->tetrahedron_q_index = 0;
  d->tetrahedron_q_size = 10;

  /* allocate memory for the free_indices_queue */
  d->free_indices_queue = (int*)malloc(10 * sizeof(int));
  d->tetrahedron_q_index = 0;
  d->tetrahedron_q_size = 10;

  /* allocate memory for the array of tetrahedra containing the current vertex
   * Initial size is 10, may grow a lot for degenerate cases... */
  d->tetrahedra_containing_vertex = (int*)malloc(10 * sizeof(int));
  d->tetrahedra_containing_vertex_index = 0;
  d->tetrahedra_containing_vertex_size = 10;

  /* determine the size of a box large enough to accommodate the entire
   * simulation volume and all possible ghost vertices required to deal with
   * boundaries. Note that we convert the generally rectangular box to a
   * square. */
  double box_anchor[3] = {hs->anchor[0] - hs->side[0],
                          hs->anchor[1] - hs->side[1],
                          hs->anchor[2] - hs->side[2]};
  /* Notice we have to take box_side rather large, because we want to fit the
   * cell and all neighbouring cells inside the first tetrahedron. This comes at
   * a loss of precision in the integer arithmetic, though... A better solution
   * would possibly be to start from 5 tetrahedra forming a cube (box_side would
   * have to be 3 in that case). */
  double box_side = fmax(hs->side[0], hs->side[1]);
  box_side = 9 * fmax(box_side, hs->side[2]);
  /* store the anchor and inverse side_length for the conversion from box
     coordinates to rescaled (integer) coordinates */
  d->anchor[0] = box_anchor[0];
  d->anchor[1] = box_anchor[1];
  d->anchor[2] = box_anchor[2];
  /* the 1.e-13 makes sure converted values are in the range [1, 2[ instead of
   * [1,2] (unlike Springel, 2010) */
  d->inverse_side = (1. - 1.e-13) / box_side;

  /* initialise the structure used to perform exact geometrical tests */
  geometry_init(&d->geometry);

  /* set up vertices for large initial tetrahedron */
  int v0 = delaunay_new_vertex(d, d->anchor[0], d->anchor[1], d->anchor[2]);
  delaunay_log("Creating vertex %i: %g %g %g", v0, d->anchor[0], d->anchor[1],
               d->anchor[2]);
  int v1 = delaunay_new_vertex(d, d->anchor[0] + box_side, d->anchor[1],
                               d->anchor[2]);
  delaunay_log("Creating vertex %i: %g %g %g", v1, d->anchor[0] + box_side,
               d->anchor[1], d->anchor[2]);
  int v2 = delaunay_new_vertex(d, d->anchor[0], d->anchor[1] + box_side,
                               d->anchor[2]);
  delaunay_log("Creating vertex %i: %g %g %g", v2, d->anchor[0],
               d->anchor[1] + box_side, d->anchor[2]);
  int v3 = delaunay_new_vertex(d, d->anchor[0], d->anchor[1],
                               d->anchor[2] + box_side);
  delaunay_log("Creating vertex %i: %g %g %g", v3, d->anchor[0], d->anchor[1],
               d->anchor[2] + box_side);
  /* Create initial large tetrahedron and 4 dummy neighbours */
  int dummy0 = delaunay_new_tetrahedron(d); /* opposite of v0 */
  int dummy1 = delaunay_new_tetrahedron(d); /* opposite of v1 */
  int dummy2 = delaunay_new_tetrahedron(d); /* opposite of v2 */
  int dummy3 = delaunay_new_tetrahedron(d); /* opposite of v3 */
  int first_tetrahedron = delaunay_new_tetrahedron(d);
  delaunay_log("Creating dummy tetrahedron %i: %i %i %i %i", dummy0, v1, v2, v3,
               -1);
  tetrahedron_init(&d->tetrahedra[dummy0], v1, v2, v3, -1);
  tetrahedron_swap_neighbour(&d->tetrahedra[dummy0], 3, first_tetrahedron, 0);
  delaunay_log("Creating dummy tetrahedron %i: %i %i %i %i", dummy1, v2, v0, v3,
               -1);
  tetrahedron_init(&d->tetrahedra[dummy1], v2, v0, v3, -1);
  tetrahedron_swap_neighbour(&d->tetrahedra[dummy1], 3, first_tetrahedron, 1);
  delaunay_log("Creating dummy tetrahedron %i: %i %i %i %i", dummy2, v3, v0, v1,
               -1);
  tetrahedron_init(&d->tetrahedra[dummy2], v3, v0, v1, -1);
  tetrahedron_swap_neighbour(&d->tetrahedra[dummy2], 3, first_tetrahedron, 2);
  delaunay_log("Creating dummy tetrahedron %i: %i %i %i %i", dummy3, v0, v2, v1,
               -1);
  tetrahedron_init(&d->tetrahedra[dummy3], v0, v2, v1, -1);
  tetrahedron_swap_neighbour(&d->tetrahedra[dummy3], 3, first_tetrahedron, 3);
  delaunay_log("Creating first tetrahedron %i: %i %i %i %i", first_tetrahedron,
               v0, v1, v2, v3);
  tetrahedron_init(&d->tetrahedra[first_tetrahedron], v0, v1, v2, v3);
  tetrahedron_swap_neighbours(&d->tetrahedra[first_tetrahedron], dummy0, dummy1,
                              dummy2, dummy3, 3, 3, 3, 3);

  /* TODO: Setup vertex-tetrahedra links... */

  /* Set pointer to last created tetrahedron, this will be our first guess when
   * adding new points to the tesselation */
  d->last_tetrahedron = first_tetrahedron;

  /* Perform sanity checks */
  delaunay_check_tessellation(d);
  delaunay_log("Passed post init check");
}

inline static void delaunay_destroy(struct delaunay* restrict d) {
  free(d->vertices);
#ifdef DELAUNAY_NONEXACT
  free(d->rescaled_vertices);
#endif
  free(d->integer_vertices);
  free(d->search_radii);
  free(d->tetrahedra);
  free(d->tetrahedron_queue);
  geometry_destroy(&d->geometry);
}

inline static int delaunay_new_tetrahedron(struct delaunay* restrict d) {
  /* check whether there is a free spot somewhere in the array */
  int index = delaunay_free_indices_queue_pop(d);
  if (index >= 0) {
    return index;
  }
  /* Else: check that we still have space for tetrahedrons available */
  if (d->tetrahedron_index == d->tetrahedron_size) {
    d->tetrahedron_size <<= 1;
    d->tetrahedra = (struct tetrahedron*)realloc(
        d->tetrahedra, d->tetrahedron_size * sizeof(struct tetrahedron));
  }
  /* return and then increase */
  return d->tetrahedron_index++;
}

inline static void delaunay_init_vertex(struct delaunay* restrict d,
                                        const int v, double x, double y,
                                        double z) {
  /* store a copy of the vertex coordinates (we should get rid of this for
     SWIFT) */
  d->vertices[3 * v] = x;
  d->vertices[3 * v + 1] = y;
  d->vertices[3 * v + 2] = z;

  /* compute the rescaled coordinates. We do this because floating point values
     in the range [1,2[ all have the same exponent (0), which guarantees that
     their mantissas form a linear sequence */
  double rescaled_x = 1. + (x - d->anchor[0]) * d->inverse_side;
  double rescaled_y = 1. + (y - d->anchor[1]) * d->inverse_side;
  double rescaled_z = 1. + (z - d->anchor[2]) * d->inverse_side;

  delaunay_assert(rescaled_x >= 1.);
  delaunay_assert(rescaled_x < 2.);
  delaunay_assert(rescaled_y >= 1.);
  delaunay_assert(rescaled_y < 2.);
  delaunay_assert(rescaled_z >= 1.);
  delaunay_assert(rescaled_z < 2.);

#ifdef DELAUNAY_NONEXACT
  /* store a copy of the rescaled coordinates to apply non-exact tests */
  d->rescaled_vertices[3 * d->vertex_index] = rescaled_x;
  d->rescaled_vertices[3 * d->vertex_index + 1] = rescaled_y;
  d->rescaled_vertices[3 * d->vertex_index + 2] = rescaled_z;
#endif

  /* convert the rescaled coordinates to integer coordinates and store these */
  d->integer_vertices[3 * v] = delaunay_double_to_int(rescaled_x);
  d->integer_vertices[3 * v + 1] = delaunay_double_to_int(rescaled_y);
  d->integer_vertices[3 * v + 2] = delaunay_double_to_int(rescaled_z);

  /* TODO: initialise the variables that keep track of the link between vertices
   * and tetrahedra. We use negative values so that we can later detect missing
   * links. */
  /*d->vertex_triangles[v] = -1;
  d->vertex_triangle_index[v] = -1;*/

  /* initialise the search radii to the largest possible value */
  d->search_radii[v] = DBL_MAX;

  delaunay_log("Initialized new vertex with index %i", v);
}

/**
 * @brief Add a new vertex with the given coordinates.
 *
 * This function first makes sure there is sufficient memory to store the
 * vertex and all its properties. It then initializes the vertex
 *
 * @param d Delaunay tessellation.
 * @param x Horizontal coordinate of the vertex.
 * @param y Vertical coordinate of the vertex.
 * @param z Z position of the vertex.
 * @return Index of the new vertex within the vertex array.
 */
inline static int delaunay_new_vertex(struct delaunay* restrict d, double x,
                                      double y, double z) {

  /* check the size of the vertex arrays against the allocated memory size */
  if (d->vertex_index == d->vertex_size) {
    /* dynamically grow the size of the arrays with a factor 2 */
    d->vertex_size <<= 1;
    d->vertices =
        (double*)realloc(d->vertices, d->vertex_size * 3 * sizeof(double));
#ifdef DELAUNAY_NONEXACT
    d->rescaled_vertices = (double*)realloc(
        d->rescaled_vertices, d->vertex_size * 3 * sizeof(double));
#endif
    d->integer_vertices = (unsigned long int*)realloc(
        d->integer_vertices, d->vertex_size * 3 * sizeof(unsigned long int));
    d->search_radii =
        (double*)realloc(d->search_radii, d->vertex_size * sizeof(double));
  }

  delaunay_init_vertex(d, d->vertex_index, x, y, z);

  /* return the vertex index and then increase it by 1.
     After this operation, vertex_index will correspond to the size of the
     vertex arrays and is also the index of the next vertex that will be
     created. */
  return d->vertex_index++;
}

/**
 * @brief Add a local (non ghost) vertex at the given index.
 * @param d Delaunay tessellation
 * @param v Index to add vertex at
 * @param x, y, z Position of vertex
 */
inline static void delaunay_add_local_vertex(struct delaunay* restrict d, int v,
                                             double x, double y, double z) {
  delaunay_assert(v < d->vertex_end && d->vertex_start <= v);
  delaunay_init_vertex(d, v, x, y, z);
  delaunay_log("Adding local vertex with position %g %g %g", x, y, z);
  delaunay_add_vertex(d, v);
}

/**
 * @brief Add a new (ghost) vertex.
 * @param d Delaunay tessellation
 * @param x, y, z Position of vertex
 */
inline static void delaunay_add_new_vertex(struct delaunay* restrict d,
                                           double x, double y, double z) {
  int v = delaunay_new_vertex(d, x, y, z);
  delaunay_log("Created new vertex with position %g %g %g", x, y, z);
  delaunay_add_vertex(d, v);
}

/**
 * @brief Finalize adding a new vertex to the tessellation.
 *
 * This function locates the tetrahedron in the current tessellation that
 * contains the new vertex. Depending on the case (see below) new tetrahedra are
 * added to the tessellation and some are removed.
 *
 * Cases:
 * 1. Point is fully inside a tetrahedron. In this case, this tetrahedron will
 *    be replaced by 4 new tetrahedra.
 * 2. Point is on face between two tetrahedra. In this case, these two will be
 *    replaced by 6 new ones.
 * 3. Point is on edge of N tetrahedra. These N tetrahedra will be replaced by
 *    2N new ones.
 *
 * @param d Delaunay tessellation.
 * @param v Index of new vertex
 */
inline static void delaunay_add_vertex(struct delaunay* restrict d, int v) {
    int number_of_tetrahedra = delaunay_find_tetrahedra_containing_vertex(d, v);

  if (number_of_tetrahedra == 1) {
    // normal case: split 'tetrahedra[0]' into 4 new tetrahedra
    delaunay_one_to_four_flip(d, v, d->tetrahedra_containing_vertex[0]);
  } else if (number_of_tetrahedra == 2) {
    // point on face: replace the 2 tetrahedra with 6 new ones
    delaunay_two_to_six_flip(d, v, d->tetrahedra_containing_vertex);
  } else if (number_of_tetrahedra > 2) {
    // point on edge: replace the N tetrahedra with 2N new ones
    delaunay_n_to_2n_flip(d, v, d->tetrahedra_containing_vertex,
                          number_of_tetrahedra);
  } else {
    fprintf(stderr, "Unknown case of number of tetrahedra: %i!\n",
            number_of_tetrahedra);
    abort();
  }

  /* Now check all tetrahedra in de queue */
  delaunay_check_tetrahedra(d);

  /* perform sanity checks if enabled */
  delaunay_check_tessellation(d);
  delaunay_log("Passed checks after inserting vertex %i", v);
}

inline static int delaunay_find_tetrahedra_containing_vertex(
    struct delaunay* restrict d, const int v) {
  /* Before we do anything: reset the index in the array of tetrahedra
   * containing the current vertex */
  d->tetrahedra_containing_vertex_index = 0;
  // TODO
  return 1;
  return d->tetrahedra_containing_vertex_index;
}

inline static void delaunay_one_to_four_flip(struct delaunay* d, int v, int t) {
  // TODO
}

inline static void delaunay_two_to_six_flip(struct delaunay* d, int v, int* t) {
  // TODO
}

inline static void delaunay_n_to_2n_flip(struct delaunay* d, int v, int* t,
                                         int n) {
  // TODO
}

/**
 * @brief Check the Delaunay criterion for tetrahedra in the queue until the
 * queue is empty.
 * @param d Delaunay triangulation
 */
inline static void delaunay_check_tetrahedra(struct delaunay* d) {
  int t = delaunay_tetrahedron_queue_pop(d);
  while (t >= 0) {
    delaunay_check_tetrahedron(d, t);
    t = delaunay_tetrahedron_queue_pop(d);
  }
}

/**
 * @brief Check the Delaunay criterion for the given tetrahedron.
 *
 * Per convention, we assume this check was triggered by inserting the final
 * vertex of this tetrahedron, so only one check is required.
 * If this check fails, this function also performs the necessary flips. All
 * new tetrahedra created by this function are also pushed to the queue for
 * checking.
 * @param d Delaunay triangulation
 * @param t The tetrahedron to check.
 */
inline static void delaunay_check_tetrahedron(struct delaunay* d, int t) {
  // TODO
}

/**
 * @brief Add the given index to the queue of free indices in the tetrahedra
 * array.
 *
 * @param d Delaunay tessellation.
 * @param idx New free index.
 */
inline static void delaunay_free_index_enqueue(struct delaunay* restrict d,
                                               int idx) {
  /* make sure there is sufficient space in the queue */
  if (d->free_indices_q_index == d->free_indices_q_size) {
    /* there isn't: increase the size of the queue with a factor 2. */
    d->free_indices_q_size <<= 1;
    d->free_indices_queue = (int*)realloc(d->free_indices_queue,
                                          d->free_indices_q_size * sizeof(int));
  }
  delaunay_log("Enqueuing free index %i", idx);
  /* add the triangle to the queue and advance the queue index */
  d->free_indices_queue[d->free_indices_q_index] = idx;
  ++d->free_indices_q_index;
}

/**
 * @brief Pop a free index in the tetrahedra array from the end of the queue.
 *
 * If no more indices are queued, this function returns a negative value.
 *
 * Note that the returned index is effectively removed from the queue
 * and will be overwritten by subsequent calls to delaunay_free_index_enqueue().
 *
 * @param d Delaunay tessellation.
 * @return Index of the next triangle to test, or -1 if the queue is empty.
 */
inline static int delaunay_free_indices_queue_pop(struct delaunay* restrict d) {
  if (d->free_indices_q_index > 0) {
    --d->free_indices_q_index;
    return d->free_indices_queue[d->free_indices_q_index];
  } else {
    return -1;
  }
}

/**
 * @brief Add the given tetrahedron to the queue of tetrahedra that need
 * checking.
 *
 * @param d Delaunay tessellation.
 * @param t Tetrahedron index.
 */
inline static void delaunay_tetrahedron_enqueue(struct delaunay* restrict d,
                                                int t) {
  /* make sure there is sufficient space in the queue */
  if (d->tetrahedron_q_index == d->tetrahedron_q_size) {
    /* there isn't: increase the size of the queue with a factor 2. */
    d->tetrahedron_q_size <<= 1;
    d->tetrahedron_queue =
        (int*)realloc(d->tetrahedron_queue, d->tetrahedron_size * sizeof(int));
  }
  delaunay_log("Enqueuing tetrahedron %i", t);
  /* add the triangle to the queue and advance the queue index */
  d->tetrahedron_queue[d->tetrahedron_q_index] = t;
  ++d->tetrahedron_q_index;
}

/**
 * @brief Pop the next tetrahedron to check from the end of the queue.
 *
 * If no more tetrahedrons are queued, this function returns a negative value.
 * Note that the returned triangle index is effectively removed from the queue
 * and will be overwritten by subsequent calls to
 * delaunay_tetrahedron_enqueue().
 *
 * @param d Delaunay tessellation.
 * @return Index of the next triangle to test, or -1 if the queue is empty.
 */
inline static int delaunay_tetrahedron_queue_pop(struct delaunay* restrict d) {
  if (d->tetrahedron_q_index > 0) {
    --d->tetrahedron_q_index;
    return d->tetrahedron_queue[d->tetrahedron_q_index];
  } else {
    return -1;
  }
}

/**
 * @brief Append a tetrahedron containing the current vertex to the array.
 * @param d Delaunay tesselation
 * @param t The tetrahedron containing the current vertex.
 */
inline static void delaunay_append_tetrahedron_containing_vertex(
    struct delaunay* d, const int t) {
  /* Enough space? */
  if (d->tetrahedra_containing_vertex_index ==
      d->tetrahedra_containing_vertex_size) {
    d->tetrahedra_containing_vertex_size <<= 1;
    d->tetrahedra_containing_vertex =
        (int*)realloc(d->tetrahedra_containing_vertex,
                      d->tetrahedra_containing_vertex_size * sizeof(int));
  }
  /* Append and increase index for next tetrahedron */
  d->tetrahedra_containing_vertex[d->tetrahedra_containing_vertex_index] = t;
  d->tetrahedra_containing_vertex_index++;
}

inline static void delaunay_consolidate(struct delaunay* restrict d) {
  /* perform a consistency test if enabled */
  delaunay_check_tessellation(d);
}

inline static void delaunay_print_tessellation(
    const struct delaunay* restrict d, const char* file_name) {
  // TODO
}

inline static void delaunay_check_tessellation(struct delaunay* restrict d) {
#ifndef DELAUNAY_CHECKS
  /* No expensive checks will be performed */
  return;
#endif

  /* loop over all non-dummy tetrahedra */
  for (int t0 = 4; t0 < d->tetrahedron_index; t0++) {
    int vt0_0 = d->tetrahedra[t0].vertices[0];
    int vt0_1 = d->tetrahedra[t0].vertices[1];
    int vt0_2 = d->tetrahedra[t0].vertices[2];
    int vt0_3 = d->tetrahedra[t0].vertices[3];
    /* loop over neighbours */
    for (int i = 0; i < 4; i++) {
      int t_ngb = d->tetrahedra[t0].neighbours[i];
      /* check neighbour relations */
      int idx_in_ngb0 = d->tetrahedra[t0].index_in_neighbour[i];
      if (d->tetrahedra[t_ngb].neighbours[idx_in_ngb0] != t0) {
        fprintf(stderr, "Wrong neighbour!\n");
        fprintf(stderr, "Tetrahedron %i: %i %i %i %i\n", t0, vt0_0, vt0_1,
                vt0_2, vt0_3);
        fprintf(
            stderr, "\tNeighbours: %i %i %i %i\n",
            d->tetrahedra[t0].neighbours[0], d->tetrahedra[t0].neighbours[1],
            d->tetrahedra[t0].neighbours[2], d->tetrahedra[t0].neighbours[3]);
        fprintf(stderr, "\tIndex in neighbour: %i %i %i %i\n",
                d->tetrahedra[t0].index_in_neighbour[0],
                d->tetrahedra[t0].index_in_neighbour[1],
                d->tetrahedra[t0].index_in_neighbour[2],
                d->tetrahedra[t0].index_in_neighbour[3]);
        fprintf(
            stderr, "Neighbour triangle %i: %i %i %i %i\n", t_ngb,
            d->tetrahedra[t_ngb].vertices[0], d->tetrahedra[t_ngb].vertices[1],
            d->tetrahedra[t_ngb].vertices[2], d->tetrahedra[t_ngb].vertices[3]);
        fprintf(stderr, "\tNeighbours: %i %i %i %i\n",
                d->tetrahedra[t_ngb].neighbours[0],
                d->tetrahedra[t_ngb].neighbours[1],
                d->tetrahedra[t_ngb].neighbours[2],
                d->tetrahedra[t_ngb].neighbours[3]);
        fprintf(stderr, "\tIndex in neighbour: %i %i %i\n",
                d->tetrahedra[t_ngb].index_in_neighbour[0],
                d->tetrahedra[t_ngb].index_in_neighbour[1],
                d->tetrahedra[t_ngb].index_in_neighbour[2],
                d->tetrahedra[t_ngb].index_in_neighbour[3]);
        abort();
      }
      /* TODO: check in_sphere criterion for vertex of t_ngb not shared with t0
       */
    }
  }
}

#endif  // CVORONOI_DELAUNAY3D_H
