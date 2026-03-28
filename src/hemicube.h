#pragma once

#include "scene.h"
#include <vector>

namespace Hemicube {
    /*
     * Compute the N×N form-factor matrix for all patch pairs using the
     * GPU FBO hemicube method (Cohen & Greenberg 1985).
     *
     * For each source patch i, five 90° perspective renders into an offscreen
     * GL_RGB8 FBO (top + four side faces) produce an ID-colour buffer.  Pixels
     * are decoded and accumulated with precomputed delta-weight tables,
     * replacing the O(n²·res²) CPU ray-casting inner loop with GPU rasterization.
     *
     * Requires an active OpenGL context — call after glutCreateWindow().
     * Returns a flat row-major vector: result[i*N + j] = F_ij.
     */
    [[nodiscard]]
    std::vector<float> computeFormFactors(const Scene& scene, int resolution = 64);
}
