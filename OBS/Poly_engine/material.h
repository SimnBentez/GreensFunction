#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>

struct Material {

    std::vector<double> EI;
    std::vector<double> AsG;

    double getEI(int k) const {

        if (k < 0 || k >= static_cast<int>(EI.size())) {
            return 0.0;
        }

        return EI[k];
    }

    double getAsG(int k) const {

        if (k < 0 || k >= static_cast<int>(AsG.size())) {
            return 0.0;
        }

        return AsG[k];
    }
};

#endif