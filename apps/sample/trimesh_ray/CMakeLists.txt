cmake_minimum_required(VERSION 3.13)
project(trimesh_ray)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_ray.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_ray
	${SOURCES})

target_link_libraries(
	trimesh_ray
	PUBLIC
		vcglib
	)
