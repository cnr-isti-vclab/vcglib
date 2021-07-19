cmake_minimum_required(VERSION 3.13)
project(trimesh_hole)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_hole.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_hole
	${SOURCES})

target_link_libraries(
	trimesh_hole
	PUBLIC
		vcglib
	)
