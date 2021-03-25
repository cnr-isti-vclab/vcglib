cmake_minimum_required(VERSION 3.13)
project(trimesh_normal)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_normal.cpp)
endif()

add_executable(trimesh_normal
	${SOURCES})

target_link_libraries(
	trimesh_normal
	PUBLIC
		vcglib
	)
