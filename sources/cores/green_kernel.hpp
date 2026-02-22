#pragma once
#include <functional>

namespace gf {

class GreenKernel {
public:
    GreenKernel(double E, double I, double L);

    double evaluate(double x, double xi) const;

private:
    double E_;
    double I_;
    double L_;
};

}