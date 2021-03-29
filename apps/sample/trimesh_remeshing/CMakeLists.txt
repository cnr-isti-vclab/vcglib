cmake_minimum_required(VERSION 3.13)
project(trimesh_remeshing)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_remeshing.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_remeshing
	${SOURCES})

target_link_libraries(
	trimesh_remeshing
	PUBLIC
		vcglib
	)
