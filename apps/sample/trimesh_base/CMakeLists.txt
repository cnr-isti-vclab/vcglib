cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11
project (trimesh_base)
add_executable(trimesh_base trimesh_base.cpp)
include_directories(../../..)
include_directories(../../../eigenlib)
