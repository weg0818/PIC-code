# Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2017.

# Open terminal
# cd path_to_the_folder/2d/bin/
# cmake -DCMAKE_BUILD_TYPE=Debug .. or cmake -DCMAKE_BUILD_TYPE=Release ..
# make
# ./gbc

cmake_minimum_required(VERSION 3.3)

project(gbc)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")


# 禁用警告信息
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -Wextra -std=c++11 ")

# 使用eigen库
set(EIGEN3_INCLUDE_DIR "/usr/include/eigen3")
include_directories(${EIGEN3_INCLUDE_DIR})
find_package(Eigen3 REQUIRED)

# 使用boost库
set(BOOST_INCLUDE_DIR "/usr/include/boost")
include_directories(${BOOST_INCLUDE_DIR})
find_package(Boost REQUIRED)

find_package(OpenMP REQUIRED)  # 查找 OpenMP

set(SRC 
        extra/BarycentricCoordinatesR2.hpp
        extra/SegmentCoordinatesR2.hpp
        extra/TriangleCoordinatesR2.hpp
        extra/VertexExpressionsR2.hpp
        extra/VertexR2.hpp
        extra/Halfedge.hpp
        extra/Face.hpp
        extra/MeshR2.hpp
        extra/TriangulatorR2.hpp
        extra/IsolinerR2.hpp
        extra/TestCoordinatesR2.hpp
        extra/math.cpp
        extra/math.h
        extra/LoadandSave.h
        extra/LoadandSave.cpp
        extra/libs/gaussinteg/gauss_legendre.hpp
        extra/libs/gaussinteg/gauss_legendre.cpp
        extra/libs/triangle/triangle.hpp
        extra/libs/triangle/triangle.cpp
        extra/libs/lbcsolver/DataSetup.hpp
        extra/libs/lbcsolver/DataTypes.hpp
        extra/libs/lbcsolver/LbcSolver.hpp
        coords/MeanValueR2.hpp
        coords/MaximumEntropyR2.hpp
        coords/HarmonicR2.hpp
        coords/LocalR2.hpp
        coords/PointwiseIterativeR2.h
        coords/PointwiseIterativeR2.cpp
        coords/IteractiveR2.cpp
        coords/IteractiveR2.h
        coords/MaximumLikelihoodR2.h
        coords/MaximumLikelihoodR2.cpp
        coords/PointwiseIterativeR2_Scaled.h
        coords/PointwiseIterativeR2_Scaled.cpp
        coords/PositiveMeanValueR2.cpp
        coords/PositiveMeanValueR2.h
        coords/WachspressR2.hpp

        main.cpp
)

add_executable(gbc ${SRC})

if(OpenMP_CXX_FOUND)  # 如果找到了 OpenMP
    target_link_libraries(gbc OpenMP::OpenMP_CXX)  # 链接 OpenMP 库
endif()
# 链接 CGAL 库
# target_link_libraries(gbc CGAL::CGAL)
