cmake_minimum_required(VERSION 3.13)
project(space_rasterized_packer)

find_package(
	Qt5
	COMPONENTS OpenGL Svg)

if (TARGET Qt5::OpenGL AND TARGET Qt5::Svg)
	if (VCG_HEADER_ONLY)
		set(SOURCES
			space_rasterized_packer.cpp
			${VCG_INCLUDE_DIRS}/wrap/qt/Outline2ToQImage.cpp
			${VCG_INCLUDE_DIRS}/wrap/qt/outline2_rasterizer.cpp)
	endif()

	add_executable(space_rasterized_packer
		${SOURCES})

	target_link_libraries(
		space_rasterized_packer
		PUBLIC
			vcglib
			Qt5::OpenGL
			Qt5::Svg
		)
else()
	message(
		STATUS "VCG examples - Skipping space_rasterized_packer example")
endif()
