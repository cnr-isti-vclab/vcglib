cmake_minimum_required(VERSION 3.13)
project(trimesh_voronoisampling)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_voronoisampling.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_voronoisampling
	${SOURCES})

target_link_libraries(
	trimesh_voronoisampling
	PUBLIC
		vcglib
	)
