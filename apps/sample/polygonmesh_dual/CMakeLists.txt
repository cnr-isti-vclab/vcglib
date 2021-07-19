cmake_minimum_required(VERSION 3.13)
project(polygonmesh_dual)

if (VCG_HEADER_ONLY)
	set(SOURCES
		polygonmesh_dual.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(polygonmesh_dual
	${SOURCES})

target_link_libraries(
	polygonmesh_dual
	PUBLIC
		vcglib
	)
