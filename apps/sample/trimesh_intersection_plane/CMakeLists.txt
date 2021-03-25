cmake_minimum_required(VERSION 3.13)
project(trimesh_intersection_plane)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_intersection_plane.cpp)
endif()

add_executable(trimesh_intersection_plane
	${SOURCES})

target_link_libraries(
	trimesh_intersection_plane
	PUBLIC
		vcglib
	)
