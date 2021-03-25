cmake_minimum_required(VERSION 3.13)
project(trimesh_voronoiclustering)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_voronoiclustering.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_voronoiclustering
	${SOURCES})

target_link_libraries(
	trimesh_voronoiclustering
	PUBLIC
		vcglib
	)
