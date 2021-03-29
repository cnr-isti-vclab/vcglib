cmake_minimum_required(VERSION 3.13)
project(trimesh_create)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_create.cpp)
endif()

add_executable(trimesh_create
	${SOURCES})

target_link_libraries(
	trimesh_create
	PUBLIC
		vcglib
	)
