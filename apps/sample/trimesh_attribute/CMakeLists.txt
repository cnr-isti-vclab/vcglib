cmake_minimum_required(VERSION 3.13)
project(trimesh_attribute)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_attribute.cpp)
endif()

add_executable(trimesh_attribute
	${SOURCES})

target_link_libraries(
	trimesh_attribute
	PUBLIC
		vcglib
	)
