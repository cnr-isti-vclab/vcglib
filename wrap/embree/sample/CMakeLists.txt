project (embree_sample)
cmake_minimum_required (VERSION 2.8.4)
add_subdirectory (embreeTest)

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 14)
set(SOURCES embree_sample.cpp)
add_executable(embree_sample ${SOURCES})

target_link_libraries(
	embree_sample
	PUBLIC
		vcglib
	)
#find_package(embree 3 REQUIRED)
