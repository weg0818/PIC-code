project(BarycentricCoordinates)

file(GLOB _SRCS *.cpp *.h *.hpp)

add_library(${PROJECT_NAME} ${_SRCS})

target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)

target_include_directories(${PROJECT_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})