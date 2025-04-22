
# Open terminal
# cd path_to_the_folder/PIC3D/build/
# cmake -DCMAKE_BUILD_TYPE=Debug .. or cmake -DCMAKE_BUILD_TYPE=Release ..
# make
# ./gbc

cmake_minimum_required(VERSION 3.3)

project(gbc)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

message(STATUS "C++ standard set to: ${CMAKE_CXX_STANDARD}")
# 禁用警告信息
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -Wextra -std=c++11 ")

# 使用eigen库
set(EIGEN3_INCLUDE_DIR "/usr/include/eigen3")
include_directories(${EIGEN3_INCLUDE_DIR})
find_package(Eigen3 REQUIRED)


# 添加 libigl 的路径
set(LIBIGL_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/ThirdPartyLib/libigl/include")
include_directories(${LIBIGL_INCLUDE_DIR})
add_subdirectory("${CMAKE_SOURCE_DIR}/ThirdPartyLib/libigl" libigl)

add_subdirectory(BarycentricCoordinates)
add_subdirectory(Basic)
add_subdirectory("${CMAKE_SOURCE_DIR}/ThirdPartyLib/LBC" LBC)

# 查找 OpenMP
find_package(OpenMP REQUIRED)

set(SRC main.cpp
)

add_executable(gbc ${SRC})

target_link_libraries(gbc PUBLIC BarycentricCoordinates)
target_link_libraries(gbc PUBLIC Basic)
target_link_libraries(gbc PUBLIC LBC)

# Link TetGen with libigl
target_link_libraries(gbc PRIVATE igl::core tetgen)

# 将 OpenMP 链接到目标
if(OpenMP_CXX_FOUND)
    target_link_libraries(gbc PUBLIC OpenMP::OpenMP_CXX)
endif()
