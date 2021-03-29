cmake_minimum_required(VERSION 3.13)
project(trimesh_optional)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_optional.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_optional
	${SOURCES})

target_link_libraries(
	trimesh_optional
	PUBLIC
		vcglib
	)
