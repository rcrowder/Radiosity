#include <cstdio>
#include "radiosity.h"

/*
 * Progressive Refinement Radiosity  (Cohen et al. 1988)
 * ------------------------------------------------------
 * Each iteration shoots unshot energy from the patch with the maximum
 * area-weighted unshot luminance:
 *
 *   ΔB_j  = ρ_j · F(i*→j) · (A_i/A_j) · unshot_i*
 *   B_j  += ΔB_j,  unshot_j += ΔB_j,  unshot_i* = 0
 */

namespace Radiosity {

static float luminance(const Vec3& v)
{
    // CIE 1931 luminance weights
    return 0.2126f*v.x + 0.7152f*v.y + 0.0722f*v.z;
}

void solve(Scene& scene, const std::vector<float>& F,
           int maxIters, float threshold)
{
    const int n = static_cast<int>(scene.patches.size());

    for (int iter = 0; iter < maxIters; ++iter) {

        // 1. Find patch with maximum unshot energy (luminance × area)
        int   best   = -1;
        float best_e = -1.f;
        for (int i = 0; i < n; ++i) {
            float e = luminance(scene.patches[i].unshot) * scene.patches[i].area;
            if (e > best_e) { best_e = e; best = i; }
        }

        if (best < 0 || best_e < threshold) {
            std::printf("Radiosity converged after %d iterations (unshot energy=%.6f)\n",
                        iter, best_e);
            return;
        }

        // 2. Shoot unshot energy from `best` to all other patches
        Patch& pi = scene.patches[best];
        for (int j = 0; j < n; ++j) {
            if (j == best) continue;
            Patch& pj  = scene.patches[j];
            float  fij = F[static_cast<std::size_t>(best) * n + j];
            if (fij < 1e-9f) continue;

            float area_ratio = (pj.area > 1e-9f) ? pi.area / pj.area : 0.f;
            Vec3  delta      = pj.reflectance * pi.unshot * (fij * area_ratio);
            pj.radiosity += delta;
            pj.unshot    += delta;
        }

        // 3. Clear shooter's unshot radiosity
        pi.unshot = {};

        //if ((iter + 1) % 10 == 0)
        //    std::printf("  iter %4d  best_energy=%.5f\n", iter + 1, best_e);
    }
}

} // namespace Radiosity
