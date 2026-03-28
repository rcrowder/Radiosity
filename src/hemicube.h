#pragma once

#include "scene.h"
#include <vector>

namespace Hemicube {
    /*
     * Compute the N×N form-factor matrix for all patch pairs using the
     * hemicube method (Cohen & Greenberg 1985).
     * Returns a flat row-major vector: result[i*N + j] = F_ij.
     */
    [[nodiscard]]
    std::vector<float> computeFormFactors(const Scene& scene, int resolution = 64);
}
