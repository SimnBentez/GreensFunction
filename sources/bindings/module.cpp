#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include "beam_solver.hpp"

namespace py = pybind11;

PYBIND11_MODULE(greensfunction, m) {
    m.doc() = "Green Function Beam Solver";

    py::class_<gf::BeamSolver>(m, "BeamSolver")
        .def(py::init<double, double, double>())
        .def("solve", &gf::BeamSolver::solve,
             py::arg("x"),
             py::arg("q"),
             py::arg("n_points") = 1000);
}