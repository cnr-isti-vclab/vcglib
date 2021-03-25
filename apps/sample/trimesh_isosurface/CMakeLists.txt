cmake_minimum_required(VERSION 3.13)
project(trimesh_isosurface)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_isosurface.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_isosurface
	${SOURCES})

target_link_libraries(
	trimesh_isosurface
	PUBLIC
		vcglib
	)
