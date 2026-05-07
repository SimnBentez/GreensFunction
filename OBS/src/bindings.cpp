#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>

#include "solver.h"

namespace py = pybind11;

// Evaluación de polinomios mediante Horner
double eval_poly(const std::vector<double>& c, double x) {

    double result = 0.0;

    for (int i = static_cast<int>(c.size()) - 1; i >= 0; --i) {
        result = result * x + c[i];
    }

    return result;
}

// Construcción de funciones evaluables desde Python
py::cpp_function make_eval_function(
    const std::vector<double>& coeffs
) {
    return py::cpp_function(
        [coeffs](double x) {
            return eval_poly(coeffs, x);
        }
    );
}

py::dict compute_shape_functions(
    int N,
    double L,
    std::vector<double> EI,
    std::vector<double> AsG
) {

    Material mat;
    mat.EI = EI;
    mat.AsG = AsG;

    Solver solver(N, L, mat);

    auto psi2 = solver.solveShape(1,0,0,0);
    auto psi3 = solver.solveShape(0,1,0,0);
    auto psi5 = solver.solveShape(0,0,1,0);
    auto psi6 = solver.solveShape(0,0,0,1);

    py::dict result;

    // Desplazamientos
    result["psi2_v"] = make_eval_function(psi2.v);
    result["psi3_v"] = make_eval_function(psi3.v);
    result["psi5_v"] = make_eval_function(psi5.v);
    result["psi6_v"] = make_eval_function(psi6.v);

    // Rotaciones
    result["psi2_theta"] = make_eval_function(psi2.theta);
    result["psi3_theta"] = make_eval_function(psi3.theta);
    result["psi5_theta"] = make_eval_function(psi5.theta);
    result["psi6_theta"] = make_eval_function(psi6.theta);

    return result;
}

PYBIND11_MODULE(shape_functions, m) {

    m.doc() = "Shape functions based on recurrence expansions";

    m.def(
        "compute_shape_functions",
        &compute_shape_functions,
        "Compute shape functions"
    );
}