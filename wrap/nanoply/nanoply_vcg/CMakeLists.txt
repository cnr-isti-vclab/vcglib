cmake_minimum_required(VERSION 3.13)

if ("cxx_std_17" IN_LIST CMAKE_CXX_COMPILE_FEATURES)

project(nanoply_vcg)

if (MSVC)
	add_compile_options(/bigobj)
endif()

if (VCG_HEADER_ONLY)
	set(SOURCES
		main.cpp)
endif()

add_executable(nanoply_vcg
	${SOURCES})
target_compile_features(nanoply_vcg PRIVATE cxx_std_17)

target_link_libraries(
	nanoply_vcg
	PUBLIC
		vcglib
	)

endif()
