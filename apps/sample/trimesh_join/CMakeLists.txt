cmake_minimum_required(VERSION 3.13)
project(trimesh_join)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_join.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_join
	${SOURCES})

target_link_libraries(
	trimesh_join
	PUBLIC
		vcglib
	)
