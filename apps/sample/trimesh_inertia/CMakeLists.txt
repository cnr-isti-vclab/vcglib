cmake_minimum_required(VERSION 3.13)
project(trimesh_inertia)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_inertia.cpp)
endif()

add_executable(trimesh_inertia
	${SOURCES})

target_link_libraries(
	trimesh_inertia
	PUBLIC
		vcglib
	)
