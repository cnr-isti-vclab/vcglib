cmake_minimum_required(VERSION 3.13)
project(trimesh_select)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_select.cpp)
endif()

add_executable(trimesh_select
	${SOURCES})

target_link_libraries(
	trimesh_select
	PUBLIC
		vcglib
	)
