//
// Created by yuyttenh on 29/06/2021.
//

#ifndef CVORONOI_TETRAHEDRON_VERTEX_QUEUE_H
#define CVORONOI_TETRAHEDRON_VERTEX_QUEUE_H

/**
 * @brief FIFO queue used during the voronoi construction.
 */
struct tetrahedron_vertex_queue {
  /*! Array of indices (in Delaunay tesselation) of next neighbouring Delaunay
   * vertices to check */
  int *vertex_indices;
  /*! Array of indices (in Delaunay tesselation) of tetrahedra linked next
   * neighbouring Delaunay vertices to check */
  int *tetrahedron_indices;
  /*! Array of indices (in respective tetrahedra) of next neighbouring Delaunay
   * vertices to check */
  int *vertex_tetrahedron_indices;
  /*! Current allocated size for the queue */
  int size;
  /*! Current index of the next element to pop */
  int index_start;
  /*! Current index to push the next new element at */
  int index_end;
};

/**
 * @brief Initialize an empty queue with initial size of 10
 * @param t tetrahedron_vertex_queue
 */
inline static void tetrahedron_vertex_queue_init(
    struct tetrahedron_vertex_queue *t) {
  t->vertex_indices = (int *)malloc(10 * sizeof(int));
  t->tetrahedron_indices = (int *)malloc(10 * sizeof(int));
  t->vertex_tetrahedron_indices = (int *)malloc(10 * sizeof(int));
  t->size = 10;
  t->index_end = 0;
  t->index_start = 0;
}

/**
 * @brief Free up memory allocated by the queue (3 arrays).
 * @param t tetrahedron_vertex_queue
 */
inline static void tetrahedron_vertex_queue_destroy(
    struct tetrahedron_vertex_queue *t) {
  free(t->vertex_indices);
  free(t->tetrahedron_indices);
  free(t->vertex_tetrahedron_indices);
}

/**
 * @brief Reset the queue without reallocating memory.
 *
 * This resets the start and end indices, effectively deleting all items from
 * the queue as they will be overwritten when new elements are added.
 *
 * @param t tetrahedron_vertex_queue
 */
inline static void tetrahedron_vertex_queue_reset(
    struct tetrahedron_vertex_queue *t) {
  t->index_start = 0;
  t->index_end = 0;
}

/**
 * @brief Push a new tuple of (vertex index, associated tetrahedron index and
 * vertex index in associated tetrahedron) to the end of the queue.
 *
 * @param t tetrahedron_vertex_queue
 * @param t_idx Index of a tetrahedron containing the Delaunay vertex.
 * @param v_idx Index of the Delaunay vertex (in the Delaunay tesselation).
 * @param v_idx_in_t Index of Delaunay vertex in its associated tetrahedron
 */
inline static void tetrahedron_vertex_queue_push(
    struct tetrahedron_vertex_queue *t, int t_idx, int v_idx, int v_idx_in_t) {
  if (t->index_end == t->size) {
    t->size <<= 1;
    t->vertex_indices =
        (int *)realloc(t->vertex_indices, t->size * sizeof(int));
    t->tetrahedron_indices =
        (int *)realloc(t->tetrahedron_indices, t->size * sizeof(int));
    t->vertex_tetrahedron_indices =
        (int *)realloc(t->vertex_tetrahedron_indices, t->size * sizeof(int));
  }
  t->vertex_indices[t->index_end] = v_idx;
  t->tetrahedron_indices[t->index_end] = t_idx;
  t->vertex_tetrahedron_indices[t->index_end] = v_idx_in_t;
  t->index_end++;
}

/**
 * Pop a (vertex index, associated tetrahedron index and vertex index in
 * associated tetrahedron) tuple from the front of the queue.
 *
 * @param t tetrahedron_vertex_queue
 * @param t_idx (return) Index of a tetrahedron containing the Delaunay vertex.
 * @param v_idx (return) Index of the Delaunay vertex (in the Delaunay
 *                       tesselation).
 * @param v_idx_in_t (return) Index of Delaunay vertex in its associated
 *                            tetrahedron
 */
inline static void tetrahedron_vertex_queue_pop(
    struct tetrahedron_vertex_queue *t, int *t_idx, int *v_idx,
    int *v_idx_in_t) {
  if (t->index_start < t->index_end) {
    *t_idx = t->tetrahedron_indices[t->index_start];
    *v_idx = t->vertex_indices[t->index_start];
    *v_idx_in_t = t->vertex_tetrahedron_indices[t->index_start];
    t->index_start++;
  } else {
    *t_idx = -1;
    *v_idx = -1;
    *v_idx_in_t = -1;
  }
}

/**
 * @brief Check whether queue is empty or not.
 *
 * @param t tetrahedron_vertex_queue
 * @return 1 for empty queue, 0 for nonempty queue.
 */
inline static int tetrahedron_vertex_queue_is_empty(
    struct tetrahedron_vertex_queue *t) {
  return t->index_start == t->index_end;
}

#endif  // CVORONOI_TETRAHEDRON_VERTEX_QUEUE_H
