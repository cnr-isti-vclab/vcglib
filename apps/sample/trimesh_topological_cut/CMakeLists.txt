cmake_minimum_required(VERSION 3.13)
project(trimesh_topological_cut)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_topological_cut.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_topological_cut
	${SOURCES})

target_link_libraries(
	trimesh_topological_cut
	PUBLIC
		vcglib
	)
