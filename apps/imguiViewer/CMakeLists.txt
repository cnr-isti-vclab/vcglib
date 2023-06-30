cmake_minimum_required(VERSION 3.1)
project(imguiViewer)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)

if(WIN32)
   message("Windows specific stuff")
elseif(APPLE)
   message("Apple specific stuff")
   add_compile_definitions(GL_SILENCE_DEPRECATION)
elseif(UNIX)
    message("Linux specific stuff")
endif()


find_package(GLFW3 REQUIRED)
find_package(GLEW 2.0 REQUIRED)
find_package(EIGEN3 REQUIRED NO_MODULE)

set(IMGUI_DIR "../../../github/imgui" CACHE INTERNAL "Path to imgui")
set(VCGLIB_DIR "../.." CACHE INTERNAL "Path to vcglib")
set(EIGEN_DIR  = "${VCGLIB_DIR}/eigenlib" CACHE INTERNAL "")


add_executable(${PROJECT_NAME} main.cpp               
               ${VCGLIB_DIR}/wrap/gui/trackball.cpp
               ${VCGLIB_DIR}/wrap/gui/trackmode.cpp
               ${VCGLIB_DIR}/wrap/ply/plylib.cpp
               ${IMGUI_DIR}/imgui.cpp
               ${IMGUI_DIR}/imgui_demo.cpp
               ${IMGUI_DIR}/imgui_draw.cpp
               ${IMGUI_DIR}/imgui_tables.cpp
               ${IMGUI_DIR}/imgui_widgets.cpp
               ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
               ${IMGUI_DIR}/backends/imgui_impl_opengl2.cpp
               FileBrowser/ImGuiFileBrowser.cpp)


target_include_directories(${PROJECT_NAME} PRIVATE ${IMGUI_DIR})
target_include_directories(${PROJECT_NAME} PRIVATE ${IMGUI_DIR}/backends)
target_include_directories(${PROJECT_NAME} PRIVATE ${GLEW_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PRIVATE ${VCGLIB_DIR})
target_include_directories(${PROJECT_NAME} PRIVATE ${EIGEN_DIR})
target_link_directories(${PROJECT_NAME} PRIVATE ${GLEW_LIBRARY_DIRS})
target_link_directories(${PROJECT_NAME} PRIVATE ${GLFW3_LIBRARY_DIRS})
target_link_libraries(${PROJECT_NAME}
    GLEW::GLEW
    glfw
)
