#include "polynomial.h"

Polynomial::Polynomial(
    const std::vector<double>& coefficients
)
    : coeffs(coefficients) {}

double Polynomial::evaluate(double x) const {

    double result = 0.0;

    for (
        int i = static_cast<int>(coeffs.size()) - 1;
        i >= 0;
        --i
    ) {
        result = result * x + coeffs[i];
    }

    return result;
}

std::vector<double> Polynomial::evaluate_range(
    double x0,
    double xN,
    int num_points
) const {

    std::vector<double> results;

    results.reserve(num_points + 1);

    double step =
        (xN - x0) / static_cast<double>(num_points);

    for (int i = 0; i <= num_points; ++i) {

        double x = x0 + i * step;

        results.push_back(evaluate(x));
    }

    return results;
}