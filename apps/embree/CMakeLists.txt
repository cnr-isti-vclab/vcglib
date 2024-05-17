cmake_minimum_required (VERSION 3.13)
project (embree_sample)

set(VCG_INCLUDE_DIRS ../../../vcglib)
add_subdirectory(${VCG_INCLUDE_DIRS} ${CMAKE_BINARY_DIR}/vcglib)

set(SOURCES 
			./embree_sample.cpp
			${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp
			)

FIND_PACKAGE(embree 4.3)

if(embree_FOUND)
	add_executable(embree_sample ${SOURCES})

	target_link_libraries(
		embree_sample
		PUBLIC
			embree
			vcglib)
endif()


