#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>

struct Material {
    std::vector<double> EI;
    std::vector<double> AsG;

    double getEI(int k) const {
        return k < EI.size() ? EI[k] : 0.0;
    }

    double getAsG(int k) const {
        return k < AsG.size() ? AsG[k] : 0.0;
    }
};

#endif