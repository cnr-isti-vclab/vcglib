cmake_minimum_required(VERSION 3.13)
project(space_shot_resection)

set(SOURCES
	space_shot_resection.cpp)

add_executable(space_shot_resection
	${SOURCES})

target_link_libraries(
	space_shot_resection
	PUBLIC
		vcglib
	)
