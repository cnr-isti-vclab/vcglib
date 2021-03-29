cmake_minimum_required(VERSION 3.13)
project(trimesh_intersection_mesh)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_intersection_mesh.cpp)
endif()

add_executable(trimesh_intersection_mesh
	${SOURCES})

target_link_libraries(
	trimesh_intersection_mesh
	PUBLIC
		vcglib
	)
