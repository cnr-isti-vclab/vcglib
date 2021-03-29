cmake_minimum_required(VERSION 3.13)
project(trimesh_clustering)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_clustering.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_clustering
	${SOURCES})

target_link_libraries(
	trimesh_clustering
	PUBLIC vcglib
	)
