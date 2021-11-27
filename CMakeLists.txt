cmake_minimum_required(VERSION 3.0.2)
project(Eigen_Quadratic)

set(CMAKE_BUILD_TYPE Release)

find_package(OsqpEigen REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
)

file(GLOB SRC 
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/locomotion/*.cpp
)

add_executable(Eigen_Quadratic ${SRC})
target_link_libraries(Eigen_Quadratic 
    OsqpEigen::OsqpEigen
    ${EIGEN_LIBRARIES}
)