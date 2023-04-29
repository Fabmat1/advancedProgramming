cmake_minimum_required(VERSION 3.20)
project(advancedProgramming)

set(CMAKE_CXX_STANDARD 23)

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/output )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/output )

add_subdirectory("${PROJECT_SOURCE_DIR}/Ex1" "${PROJECT_SOURCE_DIR}/Ex1/Ex1_output")