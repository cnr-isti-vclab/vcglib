cmake_minimum_required(VERSION 3.13)
project(colorspace)

set(SOURCES
	colorspace.cpp)

add_executable(colorspace
	${SOURCES})

target_link_libraries(
	colorspace
	PUBLIC
		vcglib
	)
