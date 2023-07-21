cmake_minimum_required(VERSION 3.13)
project(trimesh_geodesic_heat)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_geodesic_heat.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_geodesic_heat
	${SOURCES})

target_link_libraries(
	trimesh_geodesic_heat
	PUBLIC
		vcglib
	)


# cmake_minimum_required(VERSION 3.13)
# project(heatmethod)
# set(CMAKE_CXX_STANDARD 14)

# # include the vcg library
# set(VCGLIB_DIR "../../../../vcglib" CACHE INTERNAL "Path to vcglib")
# add_subdirectory(${VCGLIB_DIR} subproject/vcglib)
# include_directories(vcglib)
# link_libraries(vcglib)

# # include the Eigen library
# find_package(Eigen3 REQUIRED)
# include_directories(${EIGEN3_INCLUDE_DIR})
# link_libraries(${EIGEN3_LIBRARIES})

# # executables
# add_executable(trimesh_geodesic_heat trimesh_geodesic_heat.cpp ${VCGLIB_DIR}/wrap/ply/plylib.cpp)