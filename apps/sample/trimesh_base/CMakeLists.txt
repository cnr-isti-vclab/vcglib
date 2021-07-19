cmake_minimum_required(VERSION 3.13)
project(trimesh_base)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_base.cpp)
endif()

add_executable(trimesh_base
	${SOURCES})

target_link_libraries(
	trimesh_base
	PUBLIC
		vcglib
	)
