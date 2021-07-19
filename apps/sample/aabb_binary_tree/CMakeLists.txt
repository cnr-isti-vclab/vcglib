cmake_minimum_required(VERSION 3.13)
project(aabb_binary_tree)

set(SOURCES
	aabb_binary_tree.cpp)

add_executable(aabb_binary_tree
	${SOURCES})

target_link_libraries(
	aabb_binary_tree
	PUBLIC
		vcglib
	)
