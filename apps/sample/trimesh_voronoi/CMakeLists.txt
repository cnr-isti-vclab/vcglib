cmake_minimum_required(VERSION 3.13)
project(trimesh_voronoi)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_voronoi.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_voronoi
	${SOURCES})

target_link_libraries(
	trimesh_voronoi
	PUBLIC
		vcglib
	)
