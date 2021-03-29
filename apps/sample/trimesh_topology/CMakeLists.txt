cmake_minimum_required(VERSION 3.13)
project(trimesh_topology)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_topology.cpp)
endif()

add_executable(trimesh_topology
	${SOURCES})

target_link_libraries(
	trimesh_topology
	PUBLIC
		vcglib
	)
