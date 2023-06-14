cmake_minimum_required(VERSION 3.13)
project(heatmethod)
set(CMAKE_CXX_STANDARD 14)

# include the vcg library
set(VCGLIB_DIR "../../../../vcglib" CACHE INTERNAL "Path to vcglib")
add_subdirectory(${VCGLIB_DIR} subproject/vcglib)
include_directories(vcglib)
link_libraries(vcglib)

# include the Eigen library
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
link_libraries(${EIGEN3_LIBRARIES})

# executables
add_executable(test_heatmethod test_heatmethod.cpp ${VCGLIB_DIR}/wrap/ply/plylib.cpp)
add_executable(compute_geodesic compute_geodesic.cpp ${VCGLIB_DIR}/wrap/ply/plylib.cpp)



