cmake_minimum_required(VERSION 3.17)
project(cVoronoi C)

set(CMAKE_C_STANDARD 99)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=leak -fsanitize=address -fsanitize=undefined -g")

find_package(GMP REQUIRED)

include_directories(src)

add_executable(cVoronoi src/main.c)
target_link_libraries(cVoronoi gmp m)

add_executable(testHilbert src/test_hilbert.c)

add_executable(testGeometry3D src/test_geometry3d.c)
target_link_libraries(testGeometry3D m)
