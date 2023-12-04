cmake_minimum_required(VERSION 3.13)
project(nanoply_demo)

if (VCG_HEADER_ONLY)
	set(SOURCES
		main.cpp)
endif()

add_executable(nanoply_demo
	${SOURCES})

target_link_libraries(
	nanoply_demo
	PUBLIC
		vcglib
	)
