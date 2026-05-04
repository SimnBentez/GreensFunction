#include <iostream>
#include "solver.h"

double eval(const std::vector<double>& c, double x) {
    double s=0;
    for (int k=0;k<c.size();k++) s+=c[k]*pow(x,k);
    return s;
}

int main() {

    int N=10;
    double L=1.0;

    Material mat;
    mat.EI = {1.0,0.0};
    mat.AsG = {1.0,0.0};

    Solver solver(N,L,mat);

    auto psi2 = solver.solveShape(1,0,0,0);
    auto psi3 = solver.solveShape(0,1,0,0);
    auto psi5 = solver.solveShape(0,0,1,0);
    auto psi6 = solver.solveShape(0,0,0,1);

    double x=0.5;

    std::cout<<"psi2(x)="<<eval(psi2.v,x)<<std::endl;

    return 0;
}