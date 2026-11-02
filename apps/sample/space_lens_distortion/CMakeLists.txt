cmake_minimum_required(VERSION 3.13)
project(space_lens_distortion)

set(SOURCES
	space_lens_distortion.cpp)

add_executable(space_lens_distortion
	${SOURCES})

target_link_libraries(
	space_lens_distortion
	PUBLIC
		vcglib
	)
