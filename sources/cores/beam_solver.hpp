#pragma once
#include "green_kernel.hpp"
#include <functional>

namespace gf {

class BeamSolver {
public:
    BeamSolver(double E, double I, double L);

    double solve(double x,
                 std::function<double(double)> q,
                 int n_points = 1000) const;

private:
    GreenKernel kernel_;
    double L_;
};

}