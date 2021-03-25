cmake_minimum_required(VERSION 3.13)
project(polygonmesh_smooth)

if (VCG_HEADER_ONLY)
	set(SOURCES
		polygonmesh_smooth.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(polygonmesh_smooth
	${SOURCES})

target_link_libraries(
	polygonmesh_smooth
	PUBLIC
		vcglib
	)
