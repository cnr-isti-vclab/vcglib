cmake_minimum_required(VERSION 3.13)
project(trimesh_voronoiatlas)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_voronoiatlas.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_voronoiatlas
	${SOURCES})

target_link_libraries(
	trimesh_voronoiatlas
	PUBLIC
		vcglib
	)
