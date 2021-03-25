cmake_minimum_required(VERSION 3.13)
project(trimesh_attribute_saving)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_attribute_saving.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_attribute_saving
	${SOURCES})

target_link_libraries(
	trimesh_attribute_saving
	PUBLIC
		vcglib
	)
