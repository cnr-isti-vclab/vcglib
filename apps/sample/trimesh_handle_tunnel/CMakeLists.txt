cmake_minimum_required(VERSION 3.13)
project(trimesh_handle_tunnel)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_handle_tunnel.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_handle_tunnel
	${SOURCES})

target_link_libraries(
	trimesh_handle_tunnel
	PUBLIC
		vcglib
	)
