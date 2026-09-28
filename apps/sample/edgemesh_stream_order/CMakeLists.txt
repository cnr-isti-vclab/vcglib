cmake_minimum_required(VERSION 3.13)
project(edgemesh_stream_order)

if (VCG_HEADER_ONLY)
	set(SOURCES
		edgemesh_stream_order.cpp)
endif()

add_executable(edgemesh_stream_order
	${SOURCES})

target_link_libraries(
	edgemesh_stream_order
	PUBLIC
		vcglib
	)
