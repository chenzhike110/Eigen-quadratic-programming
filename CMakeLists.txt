cmake_minimum_required(VERSION 3.0.2)
project(Eigen_Quadratic)

set(CMAKE_BUILD_TYPE Release)

find_package(OsqpEigen REQUIRED)
find_package(Eigen3 REQUIRED)

add_executable(Eigen_Quadratic main.cpp)
target_link_libraries(Eigen_Quadratic 
    OsqpEigen::OsqpEigen
    ${EIGEN_LIBRARIES}
)