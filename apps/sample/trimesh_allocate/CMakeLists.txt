cmake_minimum_required(VERSION 3.13)
project(trimesh_allocate)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_allocate.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_allocate
	${SOURCES})

target_link_libraries(
	trimesh_allocate
	PUBLIC
		vcglib
	)
