cmake_minimum_required(VERSION 3.13)
project(trimesh_ball_pivoting)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_ball_pivoting.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_ball_pivoting
	${SOURCES})

target_link_libraries(
	trimesh_ball_pivoting
	PUBLIC
		vcglib
	)
