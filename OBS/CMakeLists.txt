cmake_minimum_required(VERSION 3.18)

project(shape_functions)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(shape_functions
    src/solver.cpp
    src/bindings.cpp
)

target_include_directories(shape_functions PRIVATE
    include
    src
)

install(TARGETS shape_functions DESTINATION .)