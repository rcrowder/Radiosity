#pragma once

#include "scene.h"
#include <vector>

namespace Radiosity {
    /*
     * Progressive Refinement Radiosity (Cohen et al. 1988).
     * F is the flat N×N form-factor matrix from Hemicube::computeFormFactors().
     */
    void solve(Scene& scene, const std::vector<float>& F,
               int maxIters = 2000, float threshold = 1e-4f);
}
