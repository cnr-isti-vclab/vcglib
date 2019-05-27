project(metro CXX)
cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)

include_directories("../../")
add_executable(metro metro.cpp ../../wrap/ply/plylib.cpp)
