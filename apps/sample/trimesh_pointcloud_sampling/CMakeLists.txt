cmake_minimum_required(VERSION 3.13)
project(trimesh_pointcloud_sampling)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_pointcloud_sampling.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_pointcloud_sampling
	${SOURCES})

target_link_libraries(
	trimesh_pointcloud_sampling
	PUBLIC
		vcglib
	)
