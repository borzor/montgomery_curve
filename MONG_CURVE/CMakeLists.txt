cmake_minimum_required(VERSION 3.5)

project(mong-curve LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(NTL REQUIRED)

add_executable(mong-curve main.cpp parameters.h mong_params.h computations.h tests.h)
target_link_libraries(mong-curve ${NTL_LIBRARY})

