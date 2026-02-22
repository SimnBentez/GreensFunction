#include "beam_solver.hpp"

namespace gf {

BeamSolver::BeamSolver(double E, double I, double L)
    : kernel_(E, I, L), L_(L) {}

double BeamSolver::solve(double x,
                         std::function<double(double)> q,
                         int n_points) const
{
    double h = L_ / n_points;
    double integral = 0.0;

    for (int i = 0; i <= n_points; ++i) {
        double xi = i * h;
        double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;
        integral += weight * kernel_.evaluate(x, xi) * q(xi);
    }

    return integral * h;
}

}