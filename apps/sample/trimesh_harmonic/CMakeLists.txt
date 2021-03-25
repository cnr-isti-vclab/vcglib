cmake_minimum_required(VERSION 3.13)
project(trimesh_harmonic)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_harmonic.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_harmonic
	${SOURCES})

target_link_libraries(
	trimesh_harmonic
	PUBLIC
		vcglib
	)
