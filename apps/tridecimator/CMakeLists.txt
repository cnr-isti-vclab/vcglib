cmake_minimum_required(VERSION 3.13)
project (tridecimator)

if (VCG_HEADER_ONLY)
	set(SOURCES
		tridecimator.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(tridecimator
	${SOURCES})

target_link_libraries(
	tridecimator
	PUBLIC
		vcglib
	)
