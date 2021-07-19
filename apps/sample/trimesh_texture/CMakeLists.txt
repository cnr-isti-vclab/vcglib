cmake_minimum_required(VERSION 3.13)
project(trimesh_texture)

find_package(
	Qt5
	COMPONENTS OpenGL Svg)

if (TARGET Qt5::OpenGL AND TARGET Qt5::Svg)

	if (VCG_HEADER_ONLY)
		set(SOURCES
			trimesh_texture.cpp
			${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp
			${VCG_INCLUDE_DIRS}/wrap/qt/Outline2ToQImage.cpp
			${VCG_INCLUDE_DIRS}/wrap/qt/outline2_rasterizer.cpp)
	endif()

	add_executable(trimesh_texture
		${SOURCES})

	target_link_libraries(
		trimesh_texture
		PUBLIC
			vcglib
			Qt5::OpenGL
			Qt5::Svg
		)
else()
	message(
		STATUS "VCG examples - Skipping trimesh_texture example")
endif()
