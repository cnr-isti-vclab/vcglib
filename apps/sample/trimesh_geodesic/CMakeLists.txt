cmake_minimum_required(VERSION 3.13)
project(trimesh_geodesic)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_geodesic.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_geodesic
	${SOURCES})

target_link_libraries(
	trimesh_geodesic
	PUBLIC
		vcglib
	)
