cmake_minimum_required(VERSION 3.17)
project(cVoronoi C)

set(CMAKE_C_STANDARD 99)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=leak -fsanitize=address -fsanitize=undefined -g")

find_package(GMP REQUIRED)

include_directories(src)

# Main program #
add_executable(cVoronoi src/main.c)
target_link_libraries(cVoronoi gmp m)

# Tests #
add_executable(testHilbert test/test_hilbert.c)

add_executable(testGeometry3D test/test_geometry3d.c)
target_link_libraries(testGeometry3D gmp m)

add_executable(testDelaunay test/test_delaunay.c)
target_link_libraries(testDelaunay gmp m)

add_executable(testQueues test/test_queues.c)
