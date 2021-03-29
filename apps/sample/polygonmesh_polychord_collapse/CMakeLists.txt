cmake_minimum_required(VERSION 3.13)
project(polygonmesh_polychord_collapse)

if (VCG_HEADER_ONLY)
	set(SOURCES
		polygonmesh_polychord_collapse.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(polygonmesh_polychord_collapse
	${SOURCES})

target_link_libraries(
	polygonmesh_polychord_collapse
	PUBLIC
		vcglib
	)
