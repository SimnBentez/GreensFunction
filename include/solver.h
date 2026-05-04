#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include "coeff.h"
#include "material.h"

struct ShapeFunctions {
    std::vector<double> v;
    std::vector<double> theta;
};

class Solver {
public:
    Solver(int N, double L, const Material& mat);

    ShapeFunctions solveShape(double ui, double ti, double vj, double tj);

private:
    int N;
    double L;
    Material mat;

    std::vector<Coeff> theta;
    std::vector<Coeff> v;

    void computeTheta();
    void computeV();
    Coeff evaluate(const std::vector<Coeff>& coeffs);
};

#endif