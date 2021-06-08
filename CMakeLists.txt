cmake_minimum_required(VERSION 3.17)
project(cVoronoi C)

set(CMAKE_C_STANDARD 99)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")
find_package(GMP REQUIRED)

include_directories(.)

add_executable(cVoronoi main.c)

target_link_libraries(cVoronoi gmp m)
