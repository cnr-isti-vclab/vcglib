cmake_minimum_required(VERSION 3.13)
project (metro)

if (VCG_HEADER_ONLY)
	set(SOURCES
		metro.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(metro
	${SOURCES})

target_link_libraries(
	metro
	PUBLIC
		vcglib
	)
