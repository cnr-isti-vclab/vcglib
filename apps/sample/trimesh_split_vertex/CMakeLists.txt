cmake_minimum_required(VERSION 3.13)
project(trimesh_split_vertex)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_split_vertex.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_split_vertex
	${SOURCES})

target_link_libraries(
	trimesh_split_vertex
	PUBLIC
		vcglib
	)
