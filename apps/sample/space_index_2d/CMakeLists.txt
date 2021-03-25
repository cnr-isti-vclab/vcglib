cmake_minimum_required(VERSION 3.13)
project(space_index_2d)

if (VCG_HEADER_ONLY)
	set(SOURCES
		space_index_2d.cpp)
endif()

add_executable(space_index_2d
	${SOURCES})

target_link_libraries(
	space_index_2d
	PUBLIC
		vcglib
	)
