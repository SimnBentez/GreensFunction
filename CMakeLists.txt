cmake_minimum_required(VERSION 3.15)
project(GreensFunction)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(pybind11 REQUIRED)

# Incluye la carpeta core para headers
include_directories(src/core)

pybind11_add_module(greensfunction
    src/bindings/module.cpp
    src/core/green_kernel.cpp
    src/core/beam_solver.cpp
)