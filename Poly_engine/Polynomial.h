#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>

class Polynomial {

public:

    explicit Polynomial(
        const std::vector<double>& coefficients
    );

    double evaluate(double x) const;

    std::vector<double> evaluate_range(
        double x0,
        double xN,
        int num_points = 100
    ) const;

private:

    std::vector<double> coeffs;
};

#endif