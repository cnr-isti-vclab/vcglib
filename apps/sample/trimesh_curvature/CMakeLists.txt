cmake_minimum_required(VERSION 3.13)
project(trimesh_curvature)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_curvature.cpp)
endif()

add_executable(trimesh_curvature
	${SOURCES})

target_link_libraries(
	trimesh_curvature
	PUBLIC
		vcglib
	)
