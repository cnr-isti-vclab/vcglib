cmake_minimum_required(VERSION 3.13)
project(trimesh_color)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_color.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_color
	${SOURCES})

target_link_libraries(
	trimesh_color
	PUBLIC
		vcglib
	)
