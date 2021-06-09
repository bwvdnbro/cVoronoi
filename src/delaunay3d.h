//
// Created by yuyttenh on 08/06/2021.
//

#ifndef CVORONOI_DELAUNAY3D_H
#define CVORONOI_DELAUNAY3D_H

#include "geometry.h"
#include "hydro_space.h"

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
  struct tetrahedon* tetrahedas;

  /*! @brief Next available index within the triangle array. Corresponds to the
   *  actual size of the triangle array. */
  int tetrahedon_index;

  /*! @brief Current size of the triangle array in memory. If triangle_size
   *  matches triangle_index, the memory buffer is full and needs to be
   *  expanded. */
  int tetrahedon_size;

  /*! @brief Queue of triangles that need checking during the incremental
   *  construction algorithm. After a new vertex has been added, all new
   *  triangles are added to this queue and are tested to see if the
   *  Delaunay criterion (empty circumcircles) still holds. New triangles
   *  created when flipping invalid triangles are also added to the queue. */
  int* queue;

  /*! @brief Next available index in the queue. Determines both the actual size
   *  of the queue as the first element that will be popped from the queue. */
  int queue_index;

  /*! @brief Current size of the queue in memory. If queue_size matches
   *  queue_index, the memory buffer is full and needs to be expanded. */
  int queue_size;

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

inline static void delaunay_init(struct delaunay* restrict d,
                                 const struct hydro_space* restrict hs,
                                 int vertex_size, int triangle_size) {
  // TODO
}

inline static void delaunay_destroy(struct delaunay* restrict d) {
  // TODO
}

inline static void delaunay_add_local_vertex(struct delaunay* restrict d, int v,
                                             double x, double y, double z) {
  // TODO
}

inline static void delaunay_consolidate(struct delaunay* restrict d) {
  // TODO
}

inline static void delaunay_print_tessellation(
    const struct delaunay* restrict d, const char* file_name) {
  // TODO
}

#endif  // CVORONOI_DELAUNAY3D_H
