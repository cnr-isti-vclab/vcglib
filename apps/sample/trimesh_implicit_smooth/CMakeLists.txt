cmake_minimum_required(VERSION 3.13)
project(trimesh_implicit_smooth)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_implicit_smooth.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_implicit_smooth
	${SOURCES})

target_link_libraries(
	trimesh_implicit_smooth
	PUBLIC
		vcglib
	)
