cmake_minimum_required(VERSION 3.13)
project(trimesh_texture_clean)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_texture_clean.cpp)
endif()

add_executable(trimesh_texture_clean
	${SOURCES})

target_link_libraries(
	trimesh_texture_clean
	PUBLIC
		vcglib
	)
