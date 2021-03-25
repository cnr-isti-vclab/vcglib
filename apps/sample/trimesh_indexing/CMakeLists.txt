cmake_minimum_required(VERSION 3.13)
project(trimesh_indexing)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_indexing.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_indexing
	${SOURCES})

target_link_libraries(
	trimesh_indexing
	PUBLIC
		vcglib
	)
