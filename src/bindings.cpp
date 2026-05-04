#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <cmath>

#include "solver.h"

namespace py = pybind11;

// evaluar polinomio
double eval_poly(const std::vector<double>& c, double x) {
    double s = 0;
    for (int k = 0; k < c.size(); ++k)
        s += c[k] * std::pow(x, k);
    return s;
}

py::dict compute_shape_functions(int N, double L,
                                 std::vector<double> EI,
                                 std::vector<double> AsG) {

    Material mat;
    mat.EI = EI;
    mat.AsG = AsG;

    Solver solver(N, L, mat);

    auto psi2 = solver.solveShape(1,0,0,0);
    auto psi3 = solver.solveShape(0,1,0,0);
    auto psi5 = solver.solveShape(0,0,1,0);
    auto psi6 = solver.solveShape(0,0,0,1);

    py::dict result;

    // devolver funciones evaluables
    result["psi2_v"] = py::cpp_function([psi2](double x){
        return eval_poly(psi2.v, x);
    });

    result["psi3_v"] = py::cpp_function([psi3](double x){
        return eval_poly(psi3.v, x);
    });

    result["psi5_v"] = py::cpp_function([psi5](double x){
        return eval_poly(psi5.v, x);
    });

    result["psi6_v"] = py::cpp_function([psi6](double x){
        return eval_poly(psi6.v, x);
    });

    result["psi2_theta"] = py::cpp_function([psi2](double x){
        return eval_poly(psi2.theta, x);
    });

    result["psi3_theta"] = py::cpp_function([psi3](double x){
        return eval_poly(psi3.theta, x);
    });

    result["psi5_theta"] = py::cpp_function([psi5](double x){
        return eval_poly(psi5.theta, x);
    });

    result["psi6_theta"] = py::cpp_function([psi6](double x){
        return eval_poly(psi6.theta, x);
    });

    return result;
}

// módulo Python
PYBIND11_MODULE(shape_functions, m) {
    m.def("compute_shape_functions", &compute_shape_functions,
          "Compute FEM shape functions");
}