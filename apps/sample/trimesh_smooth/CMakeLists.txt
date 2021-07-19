cmake_minimum_required(VERSION 3.13)
project(trimesh_smooth)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_smooth.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_smooth
	${SOURCES})

target_link_libraries(
	trimesh_smooth
	PUBLIC
		vcglib
	)
