cmake_minimum_required(VERSION 3.13)
project(trimesh_align_pair)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_align_pair.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_align_pair
	${SOURCES})

target_link_libraries(
	trimesh_align_pair
	PUBLIC
		vcglib
	)
