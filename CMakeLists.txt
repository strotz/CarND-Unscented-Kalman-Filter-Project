project(UnscentedKF)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++0x)

enable_testing()

add_subdirectory(src)
add_subdirectory(test)
