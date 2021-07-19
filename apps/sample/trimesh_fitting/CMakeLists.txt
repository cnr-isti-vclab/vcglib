cmake_minimum_required(VERSION 3.13)
project(trimesh_fitting)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_fitting.cpp)
endif()

add_executable(trimesh_fitting
	${SOURCES})

target_link_libraries(
	trimesh_fitting
	PUBLIC
		vcglib
	)
