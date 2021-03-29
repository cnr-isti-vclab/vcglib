cmake_minimum_required(VERSION 3.13)
project(trimesh_pointmatching)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_pointmatching.cpp)
endif()

add_executable(trimesh_pointmatching
	${SOURCES})

target_link_libraries(
	trimesh_pointmatching
	PUBLIC
		vcglib
	)
