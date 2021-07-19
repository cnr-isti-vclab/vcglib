cmake_minimum_required(VERSION 3.13)
project(trimesh_refine)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_refine.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_refine
	${SOURCES})

target_link_libraries(
	trimesh_refine
	PUBLIC
		vcglib
	)
