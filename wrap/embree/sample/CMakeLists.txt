project (embree_sample)
cmake_minimum_required (VERSION 3.10)
#set(CMAKE_TOOLCHAIN_FILE "C:/Users/super/Desktop/vcpkg-master/vcpkg-master/scripts/buildsystems/vcpkg.cmake")
set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_STANDARD 14)
set(SOURCES embree_sample.cpp)
set(./build ${PROJECT_BINARY_DIR}/bin)

add_executable(embree_sample ${SOURCES})

### Eigen
set(Eigen3_DIR "E:/vcpkg/vcpkg-master/vcpkg-master/packages/eigen3_x64-windows/share/eigen3")
find_package(Eigen3 CONFIG REQUIRED)
    target_link_libraries(embree_sample PRIVATE Eigen3::Eigen)

set(embree_DIR "E:/vcpkg/vcpkg-master/vcpkg-master/packages/embree3_x64-windows-release/share/embree")
find_package(embree 3 REQUIRED)
    target_link_libraries(embree_sample PRIVATE embree)
list(APPEND VCG_HEADERS wrap/EmbreeAdaptor.h) 

find_path(VCGLIB_INCLUDE_DIRS "img/img.h")
    target_include_directories(embree_sample PRIVATE ${VCGLIB_INCLUDE_DIRS})


file(COPY E:/vcpkg/vcpkg-master/vcpkg-master/packages/tbb_x64-windows-release/bin/tbb.dll
     DESTINATION E:/UniversitàMagistrale/secondoSemestre/3DgeometricModelingProcessing/vcglib/wrap/embree/sample/build/Debug)

file(COPY E:/vcpkg/vcpkg-master/vcpkg-master/packages/embree3_x64-windows-release/bin/embree3.dll
     DESTINATION E:/UniversitàMagistrale/secondoSemestre/3DgeometricModelingProcessing/vcglib/wrap/embree/sample/build/Debug)

