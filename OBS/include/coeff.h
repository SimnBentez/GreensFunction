#ifndef COEFF_H
#define COEFF_H

struct Coeff {
    double ui, ti, FY, M;

    Coeff(double a=0, double b=0, double c=0, double d=0)
        : ui(a), ti(b), FY(c), M(d) {}
};

inline Coeff operator+(const Coeff& a, const Coeff& b) {
    return {a.ui+b.ui, a.ti+b.ti, a.FY+b.FY, a.M+b.M};
}

inline Coeff operator-(const Coeff& a, const Coeff& b) {
    return {a.ui-b.ui, a.ti-b.ti, a.FY-b.FY, a.M-b.M};
}

inline Coeff operator*(double s, const Coeff& a) {
    return {s*a.ui, s*a.ti, s*a.FY, s*a.M};
}

#endif