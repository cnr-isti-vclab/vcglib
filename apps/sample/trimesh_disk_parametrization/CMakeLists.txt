cmake_minimum_required(VERSION 3.13)
project(trimesh_disk_parametrization)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_disk_parametrization.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_disk_parametrization
	${SOURCES})

target_link_libraries(
	trimesh_disk_parametrization
	PUBLIC
		vcglib
	)
