cmake_minimum_required(VERSION 3.20)
project(advancedProgramming)

set(CMAKE_CXX_STANDARD 23)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp" )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast" )

add_executable(quickpi quickpi.cpp)
add_executable(quickpi_512 quickpi_512.cpp)