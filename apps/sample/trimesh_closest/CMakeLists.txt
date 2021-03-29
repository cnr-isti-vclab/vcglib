cmake_minimum_required(VERSION 3.13)
project(trimesh_closest)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_closest.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_closest
	${SOURCES})

target_link_libraries(
	trimesh_closest
	PUBLIC
		vcglib
	)
