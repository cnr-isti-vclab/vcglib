cmake_minimum_required(VERSION 3.13)
project(trimesh_copy)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_copy.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_copy
	${SOURCES})

target_link_libraries(
	trimesh_copy
	PUBLIC
		vcglib
	)
