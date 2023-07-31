cmake_minimum_required(VERSION 3.13)
project(trimesh_geodesic_heat)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_geodesic_heat.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_geodesic_heat
	${SOURCES})

target_link_libraries(
	trimesh_geodesic_heat
	PUBLIC
		vcglib
	)