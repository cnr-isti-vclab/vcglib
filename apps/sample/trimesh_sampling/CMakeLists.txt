cmake_minimum_required(VERSION 3.13)
project(trimesh_sampling)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_sampling.cpp)
endif()

add_executable(trimesh_sampling
	${SOURCES})

target_link_libraries(
	trimesh_sampling
	PUBLIC
		vcglib
	)
