cmake_minimum_required(VERSION 3.13)
project(trimesh_kdtree)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_kdtree.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_kdtree
	${SOURCES})

target_link_libraries(
	trimesh_kdtree
	PUBLIC
		vcglib
	)
