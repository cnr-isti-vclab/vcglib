cmake_minimum_required(VERSION 3.13)
project(polygonmesh_optimize)

if (VCG_HEADER_ONLY)
	set(SOURCES
		polygonmesh_optimize.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(polygonmesh_optimize
	${SOURCES})

target_link_libraries(
	polygonmesh_optimize
	PUBLIC
		vcglib
	)
