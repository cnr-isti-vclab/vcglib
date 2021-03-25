cmake_minimum_required(VERSION 3.13)
project(polygonmesh_base)

if (VCG_HEADER_ONLY)
	set(SOURCES
		polygonmesh.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(polygonmesh_base
	${SOURCES})

target_link_libraries(
	polygonmesh_base
	PUBLIC
		vcglib
	)
