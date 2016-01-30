cmake_minimum_required(VERSION 2.8.9)

project(SemiclassicalMonteCarlo)

file(GLOB SOURCES "src/*.cpp")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(scmc ${SOURCES})