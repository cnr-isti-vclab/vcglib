cmake_minimum_required(VERSION 3.13)
project(trimesh_cylinder_clipping)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_cylinder_clipping.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_cylinder_clipping
	${SOURCES})

target_link_libraries(
	trimesh_cylinder_clipping
	PUBLIC
		vcglib
	)
