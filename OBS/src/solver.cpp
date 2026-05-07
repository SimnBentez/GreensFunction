#include "solver.h"

#include <cmath>
#include <stdexcept>

Solver::Solver(
    int N_,
    double L_,
    const Material& mat_
)
    : N(N_), L(L_), mat(mat_) {

    theta.resize(N + 1);
    v.resize(N + 1);
}

void Solver::computeTheta() {

    Coeff FY(0,0,1,0);
    Coeff M(0,0,0,1);

    theta[0] = {0,1,0,0};

    theta[1] =
        (-1.0 / mat.getEI(0)) * M;

    theta[2] =
        (1.0 / (2.0 * mat.getEI(0))) *
        (((-1.0) * mat.getEI(1)) * theta[1] + FY);

    for (int k = 2; k < N; ++k) {

        Coeff sum;

        for (int k1 = 0; k1 < k; ++k1) {

            double factor =
                static_cast<double>(k1 + 1)
                * mat.getEI(k - k1);

            sum =
                sum + (factor * theta[k1 + 1]);
        }

        theta[k + 1] =
            (-1.0 / ((k + 1) * mat.getEI(0))) * sum;
    }
}

void Solver::computeV() {

    Coeff FY(0,0,1,0);

    v[0] = {1,0,0,0};

    v[1] =
        theta[0]
        + ((-1.0 / mat.getAsG(0)) * FY);

    for (int k = 1; k < N; ++k) {

        v[k + 1] =
            (1.0 / (k + 1)) * theta[k];

        for (int k1 = 0; k1 < k; ++k1) {

            Coeff term =
                ((double)(k1 + 1) * v[k1 + 1])
                + ((-1.0) * theta[k1]);

            v[k + 1] =
                v[k + 1]
                - (1.0 / ((k + 1) * mat.getAsG(0)))
                * (mat.getAsG(k - k1) * term);
        }
    }
}

Coeff Solver::evaluate(
    const std::vector<Coeff>& coeffs
) {

    Coeff result;

    double powerL = 1.0;

    for (int k = 0; k <= N; ++k) {

        result =
            result + (powerL * coeffs[k]);

        powerL *= L;
    }

    return result;
}

ShapeFunctions Solver::solveShape(
    double ui,
    double ti,
    double vj,
    double tj
) {

    computeTheta();
    computeV();

    Coeff vL = evaluate(v);
    Coeff thetaL = evaluate(theta);

    double a11 = vL.FY;
    double a12 = vL.M;

    double a21 = thetaL.FY;
    double a22 = thetaL.M;

    double b1 =
        vj - (vL.ui * ui + vL.ti * ti);

    double b2 =
        tj - (thetaL.ui * ui + thetaL.ti * ti);

    double det =
        a11 * a22 - a12 * a21;

    constexpr double EPS = 1e-14;

    if (std::abs(det) < EPS) {
        throw std::runtime_error(
            "Singular system detected in solveShape()"
        );
    }

    double FY =
        (b1 * a22 - b2 * a12) / det;

    double M =
        (a11 * b2 - a21 * b1) / det;

    std::vector<double> v_num(N + 1);
    std::vector<double> theta_num(N + 1);

    for (int k = 0; k <= N; ++k) {

        v_num[k] =
            v[k].ui * ui
            + v[k].ti * ti
            + v[k].FY * FY
            + v[k].M * M;

        theta_num[k] =
            theta[k].ui * ui
            + theta[k].ti * ti
            + theta[k].FY * FY
            + theta[k].M * M;
    }

    return {v_num, theta_num};
}