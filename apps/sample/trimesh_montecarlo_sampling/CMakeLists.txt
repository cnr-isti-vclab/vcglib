cmake_minimum_required(VERSION 3.13)
project(trimesh_montecarlo_sampling)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_montecarlo_sampling.cpp)
endif()

add_executable(trimesh_montecarlo_sampling
	${SOURCES})

target_link_libraries(
	trimesh_montecarlo_sampling
	PUBLIC
		vcglib
	)
