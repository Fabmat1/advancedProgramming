cmake_minimum_required(VERSION 3.20)
project(advancedProgramming)

set(CMAKE_CXX_STANDARD 23)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp" )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast" )

add_executable(bonus_DM Bonus_DM.cpp)