project (embree_sample)
cmake_minimum_required (VERSION 3.10)
set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_STANDARD 14)
set(SOURCES embree_sample.cpp)
set(./build ${PROJECT_BINARY_DIR}/bin)

add_executable(embree_sample ${SOURCES})

FIND_PACKAGE(embree 3.0 REQUIRED)

INCLUDE_DIRECTORIES(${EMBREE_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(embree_sample ${EMBREE_LIBRARY})

INCLUDE_DIRECTORIES(../../..)
INCLUDE_DIRECTORIES(../../../eigenlib)

