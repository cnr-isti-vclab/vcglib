cmake_minimum_required(VERSION 3.13)
project(trimesh_base)
add_subdirectory(vcglib)
set(CMAKE_CXX_STANDARD 14)
set(SOURCES
	simple_main.cpp)

add_executable(simple_main
	${SOURCES})

target_link_libraries(
	simple_main
	PUBLIC
		vcglib
	)
