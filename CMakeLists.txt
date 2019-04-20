cmake_minimum_required(VERSION 3.10)
project(pose_slam_2d_benchmark)

set(CMAKE_CXX_STANDARD 17)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(pose_slam_2d_benchmark main.cpp pose_graph.h se2.h se2.cpp)
target_link_libraries (pose_slam_2d_benchmark Eigen3::Eigen)
