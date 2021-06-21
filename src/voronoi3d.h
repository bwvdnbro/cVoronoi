//
// Created by yuyttenh on 08/06/2021.
//

#ifndef CVORONOI_VORONOI3D_H
#define CVORONOI_VORONOI3D_H

/**
 * @brief Voronoi interface.
 *
 * An interface is a connection between two neighbouring Voronoi cells. It is
 * completely defined by the indices of the generators that generate the two
 * neighbouring cells, a surface area and a midpoint position.
 */
struct voronoi_pair {
  /*! Pointer to particle corresponding to the generator on the left of the
   * interface (always a particle within the local cell). */
  int left;

  /*! Pointer to particle corresponding to the generator on the right of the
   * interface (can be a local particle, but also a particle in a
   * neighbouring cell). */
  int right;

  struct cell *right_cell;

  /*! Surface area of the interface. */
  double surface_area;

  /*! Midpoint of the interface. */
  double midpoint[3];

#ifdef VORONOI_STORE_CONNECTIONS
  /*! Vertices of the interface. */
  double *vertices;

  /*! Index to put the next vertex at. This is also the size of the
   * vertex_indices array. */
  int vertex_index;

  /*! current allocated size for the vertex_indices */
  int vertex_size;
#endif
};

inline static void voronoi_pair_init(struct voronoi_pair *pair,
                                     struct cell *restrict c,
                                     int left_part_pointer,
                                     int right_part_pointer, double *vertices,
                                     int n_vertices) {
  pair->right_cell = c;
  pair->left = left_part_pointer;
  pair->right = right_part_pointer;

  pair->surface_area =
      geometry3d_compute_centroid_area(vertices, n_vertices, pair->midpoint);

#ifdef VORONOI_STORE_CONNECTIONS
  pair->vertices = (double *)malloc(3 * n_vertices * sizeof(double));
  for (int i = 0; i < n_vertices; i++) {
    pair->vertices[3 * i] = vertices[3 * i];
    pair->vertices[3 * i + 1] = vertices[3 * i + 1];
    pair->vertices[3 * i + 2] = vertices[3 * i + 2];
  }
#endif
}

inline static void voronoi_pair_destroy(struct voronoi_pair *pair) {
#ifdef VORONOI_STORE_CONNECTIONS
  free(pair->vertices);
#endif
}

struct tetrahedron_vertex_queue {
  int *vertex_indices;
  int *tetrahedron_indices;
  int *vertex_tetrahedron_indices;
  int size;
  int index_start;
  int index_end;
};

inline static void tetrahedron_vertex_queue_init(
    struct tetrahedron_vertex_queue *t) {
  t->vertex_indices = (int *)malloc(10 * sizeof(int));
  t->tetrahedron_indices = (int *)malloc(10 * sizeof(int));
  t->vertex_tetrahedron_indices = (int *)malloc(10 * sizeof(int));
  t->size = 10;
  t->index_end = 0;
  t->index_start = 0;
}

inline static void tetrahedron_vertex_queue_destroy(
    struct tetrahedron_vertex_queue *t) {
  free(t->vertex_indices);
  free(t->tetrahedron_indices);
  free(t->vertex_tetrahedron_indices);
}

inline static void tetrahedron_vertex_queue_reset(
    struct tetrahedron_vertex_queue *t) {
  t->index_start = 0;
  t->index_end = 0;
}

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

inline static int tetrahedron_vertex_queue_is_empty(
    struct tetrahedron_vertex_queue *t) {
  return t->index_start == t->index_end;
}

/**
 * @brief Voronoi cell.
 *
 * A cell stores geometrical information about a Voronoi cell: its volume and
 * the location of its centroid.
 */
struct voronoi_cell {
  /*! Cell volume. */
  double volume;

  /*! Cell centroid. */
  double centroid[3];

#ifdef VORONOI_STORE_GENERATORS
  /*! Position of the cell generator. */
  double generator[3];
#endif

#ifdef VORONOI_STORE_CELL_STATS
  /*! Number of faces of this cell. */
  int nface;
#endif
};

struct voronoi {
  /*! @brief Voronoi cells. */
  struct voronoi_cell *cells;

  /*! @brief Number of cells. */
  int number_of_cells;

  /*! @brief Voronoi cell pairs. We store these per (SWIFT) cell, i.e. pairs[0]
   *  contains all pairs that are completely contained within this cell, while
   *  pairs[1] corresponds to pairs crossing the boundary between this cell and
   *  the cell with coordinates that are lower in all coordinate directions (the
   *  cell to the left, front, bottom, sid=0), and so on. */
  struct voronoi_pair *pairs[2];

  /*! @brief Current number of pairs per cell index. */
  int pair_index[2];

  /*! @brief Allocated number of pairs per cell index. */
  int pair_size[2];
};

/* Forward declarations */
inline static int double_cmp(double double1, double double2,
                             unsigned long precision);
inline static int voronoi_new_face(struct voronoi *v, int sid,
                                   struct cell *restrict c,
                                   int left_part_pointer,
                                   int right_part_pointer, double *vertices,
                                   int n_vertices);
inline static void voronoi_check_grid(struct voronoi *restrict v);

inline static void voronoi_init(struct voronoi *restrict v,
                                struct delaunay *restrict d) {
  delaunay_assert(d->vertex_end > 0);

  /* the number of cells equals the number of non-ghost and non-dummy
     vertex_indices in the Delaunay tessellation */
  v->number_of_cells = d->vertex_end - d->vertex_start;
  /* allocate memory for the voronoi cells */
  v->cells = (struct voronoi_cell *)malloc(v->number_of_cells *
                                           sizeof(struct voronoi_cell));
  /* Allocate memory to store voronoi vertices (will be freed at end) */
  double *voronoi_vertices =
      (double *)malloc(3 * (d->tetrahedron_index - 4) * sizeof(double));

  /* loop over the tetrahedra in the Delaunay tessellation and compute the
     midpoints of their circumspheres. These happen to be the vertices of
     the Voronoi grid (because they are the points of equal distance to 3
     generators, while the Voronoi edges are the lines of equal distance to 2
     generators) */
  for (int i = 0; i < d->tetrahedron_index - 4; i++) {
    struct tetrahedron *t = &d->tetrahedra[i + 4];
    /* Skip inactive (deleted, but not yet replaced) tetrahedra */
    if (!t->active) continue;

    int v0 = t->vertices[0];
    int v1 = t->vertices[1];
    int v2 = t->vertices[2];
    int v3 = t->vertices[3];
    voronoi_assert(v0 >= 0 && v1 >= 0 && v2 >= 0 && v3 >= 0);

    /* if the triangle is not linked to a non-ghost, non-dummy vertex, it is not
     * a grid vertex and we can skip it. */
    if (v0 >= v->number_of_cells && v1 >= v->number_of_cells &&
        v2 >= v->number_of_cells && v3 >= v->number_of_cells) {
      continue;
    }

    /* Extract coordinates from the Delaunay vertices (generators)
     * FUTURE NOTE: In swift we should read this from the particles themselves!
     * */
    double v0x, v0y, v0z, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z;
    if (v0 < d->vertex_end || v0 >= d->ghost_offset) {
      v0x = d->vertices[3 * v0];
      v0y = d->vertices[3 * v0 + 1];
      v0z = d->vertices[3 * v0 + 2];
    } else {
      /* This could mean that a neighbouring cell of this grids cell is empty!
       * Or that we did not add all the necessary ghost vertex_indices to the
       * delaunay tesselation. */
      voronoi_error(
          "Vertex is part of triangle with Dummy vertex! This could mean that "
          "one of the neighbouring cells is empty.");
    }
    if (v1 < d->vertex_end || v1 >= d->ghost_offset) {
      v1x = d->vertices[3 * v1];
      v1y = d->vertices[3 * v1 + 1];
      v1z = d->vertices[3 * v1 + 2];
    } else {
      voronoi_error(
          "Vertex is part of triangle with Dummy vertex! This could mean that "
          "one of the neighbouring cells is empty.");
    }
    if (v2 < d->vertex_end || v2 >= d->ghost_offset) {
      v2x = d->vertices[3 * v2];
      v2y = d->vertices[3 * v2 + 1];
      v2z = d->vertices[3 * v2 + 2];
    } else {
      voronoi_error(
          "Vertex is part of triangle with Dummy vertex! This could mean that "
          "one of the neighbouring cells is empty.");
    }
    if (v3 < d->vertex_end || v3 >= d->ghost_offset) {
      v3x = d->vertices[3 * v3];
      v3y = d->vertices[3 * v3 + 1];
      v3z = d->vertices[3 * v3 + 2];
    } else {
      voronoi_error(
          "Vertex is part of triangle with Dummy vertex! This could mean that "
          "one of the neighbouring cells is empty.");
    }

    geometry3d_compute_circumcenter(v0x, v0y, v0z, v1x, v1y, v1z, v2x, v2y, v2z,
                                    v3x, v3y, v3z, &voronoi_vertices[3 * i]);
#ifdef VORONOI_CHECKS
    const double cx = voronoi_vertices[3 * i];
    const double cy = voronoi_vertices[3 * i + 1];
    const double cz = voronoi_vertices[3 * i + 2];

    const double r0 = (cx - v0x) * (cx - v0x) + (cy - v0y) * (cy - v0y) +
                      (cz - v0z) * (cz - v0z);
    const double r1 = (cx - v1x) * (cx - v1x) + (cy - v1y) * (cy - v1y) +
                      (cz - v1z) * (cz - v1z);
    const double r2 = (cx - v2x) * (cx - v2x) + (cy - v2y) * (cy - v2y) +
                      (cz - v2z) * (cz - v2z);
    const double r3 = (cx - v3x) * (cx - v3x) + (cy - v3y) * (cy - v3y) +
                      (cz - v3z) * (cz - v3z);
    voronoi_assert(double_cmp(r0, r1, 1e10) && double_cmp(r0, r2, 1e10) &&
                   double_cmp(r0, r3, 1e10));
#endif
  } /* loop over the Delaunay tetrahedra and compute the circumcenters */

  /* Allocate memory for the voronoi pairs (faces). */
  for (int i = 0; i < 2; ++i) {
    v->pairs[i] =
        (struct voronoi_pair *)malloc(10 * sizeof(struct voronoi_pair));
    v->pair_index[i] = 0;
    v->pair_size[i] = 10;
  }

  /* Allocate memory for the neighbour flags and initialize them to 0 */
  int *neighbour_flags = (int *)malloc(d->vertex_index * sizeof(int));
  for (int i = 0; i < d->vertex_index; i++) {
    neighbour_flags[i] = 0;
  }

  /* Allocate a tetrahedron_vertex_queue */
  struct tetrahedron_vertex_queue queue;
  tetrahedron_vertex_queue_init(&queue);

  /* The size of the array used to temporarily store the vertices of the voronoi
   * faces in */
  int face_vertices_size = 10;
  /* Temporary array to store face vertices in */
  double *face_vertices =
      (double *)malloc(3 * face_vertices_size * sizeof(double));

  /* loop over all cell generators, and hence over all non-ghost, non-dummy
     Delaunay vertex_indices */
  for (int gen_idx_in_d = 0; gen_idx_in_d < v->number_of_cells;
       gen_idx_in_d++) {
    /* First reset the tetrahedron_vertex_queue */
    tetrahedron_vertex_queue_reset(&queue);
    /* Set the flag of the central generator so that we never pick it as
     * possible neighbour */
    neighbour_flags[gen_idx_in_d] = 1;

    /* Create a new voronoi cell for this generator */
    struct voronoi_cell *this_cell = &v->cells[gen_idx_in_d - d->vertex_start];
    this_cell->volume = 0.;
    this_cell->centroid[0] = 0.;
    this_cell->centroid[1] = 0.;
    this_cell->centroid[2] = 0.;
    int nface = 0;

    /* get the generator position, we use it during centroid/volume
       calculations */
    voronoi_assert(gen_idx_in_d < d->vertex_end);
    double ax = d->vertices[3 * gen_idx_in_d];
    double ay = d->vertices[3 * gen_idx_in_d + 1];
    double az = d->vertices[3 * gen_idx_in_d + 2];

#ifdef VORONOI_STORE_GENERATORS
    this_cell->generator[0] = ax;
    this_cell->generator[1] = ay;
    this_cell->generator[2] = az;
#endif

    /* Get a tetrahedron containing the central generator */
    int t_idx = d->vertex_tetrahedron_links[gen_idx_in_d];
    int gen_idx_in_t = d->vertex_tetrahedron_index[gen_idx_in_d];

    /* Pick another vertex (generator) from this tetrahedron and add it to the
     * queue */
    int other_v_idx_in_t = (gen_idx_in_t + 1) % 4;
    struct tetrahedron *t = &d->tetrahedra[t_idx];
    int other_v_idx_in_d = t->vertices[other_v_idx_in_t];
    tetrahedron_vertex_queue_push(&queue, t_idx, other_v_idx_in_d,
                                  other_v_idx_in_t);
    /* update flag of the other vertex */
    neighbour_flags[other_v_idx_in_d] = 1;

    while (!tetrahedron_vertex_queue_is_empty(&queue)) {
      /* with each delaunay edge corresponds a voronoi face */
      nface++;

      /* Pop the next axis vertex and corresponding tetrahedron from the queue
       */
      int first_t_idx, axis_idx_in_d, axis_idx_in_t;
      tetrahedron_vertex_queue_pop(&queue, &first_t_idx, &axis_idx_in_d,
                                   &axis_idx_in_t);
      voronoi_assert(axis_idx_in_d >= 0 && (axis_idx_in_d < d->vertex_end ||
                                            axis_idx_in_d >= d->ghost_offset));
      struct tetrahedron *first_t = &d->tetrahedra[first_t_idx];

      /* Get a non axis vertex from first_t */
      int non_axis_idx_in_first_t = (axis_idx_in_t + 1) % 4;
      if (first_t->vertices[non_axis_idx_in_first_t] == gen_idx_in_d) {
        non_axis_idx_in_first_t = (non_axis_idx_in_first_t + 1) % 4;
      }
      int non_axis_idx_in_d = first_t->vertices[non_axis_idx_in_first_t];

      if (!neighbour_flags[non_axis_idx_in_d]) {
        /* Add this vertex and tetrahedron to the queue and update its flag */
        tetrahedron_vertex_queue_push(&queue, first_t_idx, non_axis_idx_in_d,
                                      non_axis_idx_in_first_t);
        neighbour_flags[non_axis_idx_in_d] |= 1;
      }

      /* Get a neighbouring tetrahedron of first_t sharing the axis */
      int cur_t_idx = first_t->neighbours[non_axis_idx_in_first_t];
      struct tetrahedron *cur_t = &d->tetrahedra[cur_t_idx];
      int prev_t_idx_in_cur_t =
          first_t->index_in_neighbour[non_axis_idx_in_first_t];

      /* Get a neighbouring tetrahedron of cur_t that is not first_t, sharing
       * the same axis */
      int next_t_idx_in_cur_t = (prev_t_idx_in_cur_t + 1) % 4;
      while (cur_t->vertices[next_t_idx_in_cur_t] == gen_idx_in_d ||
             cur_t->vertices[next_t_idx_in_cur_t] == axis_idx_in_d) {
        next_t_idx_in_cur_t = (next_t_idx_in_cur_t + 1) % 4;
      }
      int next_t_idx = cur_t->neighbours[next_t_idx_in_cur_t];

      /* Get the next non axis vertex and add it to the queue if necessary */
      int next_non_axis_idx_in_d = cur_t->vertices[next_t_idx_in_cur_t];
      if (!neighbour_flags[next_non_axis_idx_in_d]) {
        tetrahedron_vertex_queue_push(&queue, cur_t_idx, next_non_axis_idx_in_d,
                                      next_t_idx_in_cur_t);
        neighbour_flags[next_non_axis_idx_in_d] |= 1;
      }

      /* Get the coordinates of the voronoi vertex of the new face */
      int vor_vertex0_idx = first_t_idx - 4;
      face_vertices[0] = voronoi_vertices[3 * vor_vertex0_idx];
      face_vertices[1] = voronoi_vertices[3 * vor_vertex0_idx + 1];
      face_vertices[2] = voronoi_vertices[3 * vor_vertex0_idx + 2];
      int face_vertices_index = 3;

      /* Loop around the axis */
      while (next_t_idx != first_t_idx) {
        /* Get the coordinates of the voronoi vertex corresponding to cur_t and
         * next_t */
        if (face_vertices_index + 6 > face_vertices_size) {
          face_vertices_size <<= 1;
          face_vertices = (double *)realloc(
              face_vertices, 3 * face_vertices_size * sizeof(double));
        }
        const int vor_vertex1_idx = cur_t_idx - 4;
        face_vertices[face_vertices_index] =
            voronoi_vertices[3 * vor_vertex1_idx];
        face_vertices[face_vertices_index + 1] =
            voronoi_vertices[3 * vor_vertex1_idx + 1];
        face_vertices[face_vertices_index + 2] =
            voronoi_vertices[3 * vor_vertex1_idx + 2];
        const int vor_vertex2_idx = next_t_idx - 4;
        face_vertices[face_vertices_index + 3] =
            voronoi_vertices[3 * vor_vertex2_idx];
        face_vertices[face_vertices_index + 4] =
            voronoi_vertices[3 * vor_vertex2_idx + 1];
        face_vertices[face_vertices_index + 5] =
            voronoi_vertices[3 * vor_vertex2_idx + 2];
        face_vertices_index += 6;

        /* Update cell volume and triangle_centroid */
        double tetrahedron_centroid[3];
        const double V = geometry3d_compute_centroid_volume_tetrahedron(
            ax, ay, az, face_vertices[0], face_vertices[1], face_vertices[2],
            face_vertices[face_vertices_index - 6],
            face_vertices[face_vertices_index - 5],
            face_vertices[face_vertices_index - 4],
            face_vertices[face_vertices_index - 3],
            face_vertices[face_vertices_index - 2],
            face_vertices[face_vertices_index - 1], tetrahedron_centroid);
        this_cell->volume += V;
        this_cell->centroid[0] += V * tetrahedron_centroid[0];
        this_cell->centroid[1] += V * tetrahedron_centroid[1];
        this_cell->centroid[2] += V * tetrahedron_centroid[2];

        /* Update variables */
        prev_t_idx_in_cur_t = cur_t->index_in_neighbour[next_t_idx_in_cur_t];
        cur_t_idx = next_t_idx;
        cur_t = &d->tetrahedra[cur_t_idx];
        next_t_idx_in_cur_t = (prev_t_idx_in_cur_t + 1) % 4;
        while (cur_t->vertices[next_t_idx_in_cur_t] == gen_idx_in_d ||
               cur_t->vertices[next_t_idx_in_cur_t] == axis_idx_in_d) {
          next_t_idx_in_cur_t = (next_t_idx_in_cur_t + 1) % 4;
        }
        next_t_idx = cur_t->neighbours[next_t_idx_in_cur_t];
        /* Get the next non axis vertex and add it to the queue if necessary */
        next_non_axis_idx_in_d = cur_t->vertices[next_t_idx_in_cur_t];
        if (!neighbour_flags[next_non_axis_idx_in_d]) {
          tetrahedron_vertex_queue_push(
              &queue, cur_t_idx, next_non_axis_idx_in_d, next_t_idx_in_cur_t);
          neighbour_flags[next_non_axis_idx_in_d] |= 1;
        }
      }
      if (axis_idx_in_d < d->vertex_end) {
        /* Store faces only once */
        if (gen_idx_in_d < axis_idx_in_d) {
          voronoi_new_face(v, 0, NULL, gen_idx_in_d, axis_idx_in_d,
                           face_vertices, face_vertices_index);
        }
      } else { /* axis_idx_in_d >= d->ghost_offset */
        voronoi_new_face(v, 1, NULL, gen_idx_in_d, axis_idx_in_d, face_vertices,
                         face_vertices_index);
      }
    }
    this_cell->centroid[0] /= this_cell->volume;
    this_cell->centroid[1] /= this_cell->volume;
    this_cell->centroid[2] /= this_cell->volume;
#ifdef VORONOI_STORE_CELL_STATS
    this_cell->nface = nface;
#endif
    /* reset flags for all neighbours of this cell */
    neighbour_flags[gen_idx_in_d] = 0;
    for (int i = 0; i < queue.index_end; i++) {
      voronoi_assert(queue.vertex_indices[i] < d->vertex_index)
          neighbour_flags[queue.vertex_indices[i]] = 0;
    }
#ifdef VORONOI_CHECKS
    for (int i = 0; i < d->vertex_index; i++) {
      voronoi_assert(neighbour_flags[i] == 0);
    }
#endif
  }
  free(voronoi_vertices);
  free(neighbour_flags);
  free(face_vertices);
  tetrahedron_vertex_queue_destroy(&queue);
  voronoi_check_grid(v);
}

inline static int voronoi_new_face(struct voronoi *v, int sid,
                                   struct cell *restrict c,
                                   int left_part_pointer,
                                   int right_part_pointer, double *vertices,
                                   int n_vertices) {
  if (v->pair_index[sid] == v->pair_size[sid]) {
    v->pair_size[sid] <<= 1;
    v->pairs[sid] = (struct voronoi_pair *)realloc(
        v->pairs[sid], v->pair_size[sid] * sizeof(struct voronoi_pair));
  }
  /* Initialize pair */
  struct voronoi_pair *this_pair = &v->pairs[sid][v->pair_index[sid]];
  voronoi_pair_init(this_pair, c, left_part_pointer, right_part_pointer,
                    vertices, n_vertices);
  /* return and then increase */
  return v->pair_index[sid]++;
}

inline static void voronoi_destroy(struct voronoi *restrict v) {
  free(v->cells);
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < v->pair_index[i]; j ++) {
      voronoi_pair_destroy(&v->pairs[i][j]);
    }
    free(v->pairs[i]);
  }
}

inline static void voronoi_check_grid(struct voronoi *restrict v) {
  double total_volume = 0.;
  for (int i = 0; i < v->number_of_cells; i++) {
    total_volume += v->cells[i].volume;
  }
  fprintf(stderr, "Total volume: %g", total_volume);
}

inline static void voronoi_print_grid(const struct voronoi *v,
                                      const char *filename) {
  // TODO
}

inline static int double_cmp(double double1, double double2,
                             unsigned long precision) {
  long long1, long2;
  if (double1 > 0)
    long1 = (long)(double1 * precision + .5);
  else
    long1 = (long)(double1 * precision - .5);
  if (double2 > 0)
    long2 = (long)(double2 * precision + .5);
  else
    long2 = (long)(double2 * precision - .5);
  return (long1 == long2);
}

#endif  // CVORONOI_VORONOI3D_H
