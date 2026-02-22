#include "green_kernel.hpp"
#include <cmath>

namespace gf {

GreenKernel::GreenKernel(double E, double I, double L)
    : E_(E), I_(I), L_(L) {}

double GreenKernel::evaluate(double x, double xi) const {
    // Simplified symmetric kernel for testing
    double dx = std::abs(x - xi);
    return (dx * dx * dx) / (6.0 * E_ * I_);
}

}