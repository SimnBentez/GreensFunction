#include <iostream>
#include <vector>
#include "polynomial.h"

int main() {
    // Ejemplo: P(x) = 1 + 2x + 3x² + 4x³
    std::vector<double> coeffs = {1.0, 2.0, 3.0, 4.0};

    Polynomial poly(coeffs);

    double x0 = 0.0;
    double xN = 5.0;
    int points = 100;

    std::vector<double> values = poly.evaluate_range(x0, xN, points);

    for (int i = 0; i <= points; ++i) { // imprime algunos
        std::cout << values[i] << std::endl;
    }

    return 0;
}