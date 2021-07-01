# cVoronoi

### TODO: Delaunay3D

- [X] check in-sphere criterion in check_tesselation
- [X] print and plot tesselation
- [x] Flips for degenerate cases
- [X] Periodic boundary conditions (implemented in dummy way by simpy copying the entire domain to all directions)
- [X] Vertex tetrahedra links (needed for Voronoi construction)
- [X] Delaunay search radius
- [X] Use queue structs
- [ ] better PBC
- [ ] Non-exact tests

### TODO: Voronoi3D

- [X] Refactor voronoi2d code to match current version in swift
- [X] Construction
- [X] Documentation
- [ ] Check with CmacIonize
- [ ] Better visualisation

#### Construction:

Correspondences:

- With every vertex of the delaunay tesselation corresponds a voronoi cell.
- With every edge of the delaunay tesselation corresponds a voronoi face.
- With every tetrahedron of the delaunay tesselation corresponds a voronoi vertex.

Algorithm:

1. Calculate the centroid of each tetrahedron of the Delaunay tesselation and store those centroids in a temporary
   array.
1. For each non-dummy, non-ghost generator:
    1. Create an array of flags to check whether a vertex was already considered as neighbour of the current generator
    1. Get a tetrahedron connected to this generator
    1. For all the other vertex_indices in this tetrahedron,
        - loop around the tetrahedra sharing the edge connecting both generators and create the corresponding face from
          the centroids.
        - Add pairs of other vertex_indices from the encountered tetrahedra and pointers to their corresponding
          tetrahedra to the queue of vertex_indices/edges to check, if they were not added before (check and update
          flags).
    1. Get a new vertex/tetrahedron pair from the queue and repeat step iii.
    1. Repeat until queue is empty
    1. Reset all flags (we can also build a second queue/list of the indices of all the vertex_indices that were checked
       during the construction of a face, to avoid having to reset ___all___ the flags of the ___entire___ tesselation
       to false. __Both are managed with one FIFO queue__).
          
