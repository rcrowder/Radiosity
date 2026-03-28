#include <cmath>
#include <vector>
#include "hemicube.h"

/*
 * Hemicube Form-Factor Computation  (Cohen & Greenberg 1985)
 * ----------------------------------------------------------
 * For each patch i we place a hemicube at its centre aligned with its normal.
 * Five faces (top + four sides) are ray-cast against all other patches.
 * Each pixel contributes a weighted delta form-factor according to the
 * analytical hemicube kernel; row i is then normalised so weights sum to 1.
 */

namespace Hemicube {

// ---------- ray–quad intersection (Möller–Trumbore, two-triangle split) ------
// Returns t > 0 on hit, -1 otherwise.
static float rayQuad(Vec3 orig, Vec3 dir, const Vec3 verts[4])
{
    // Triangle 0: v0,v1,v2
    {
        Vec3  e1 = verts[1] - verts[0], e2 = verts[2] - verts[0];
        Vec3  h  = dir.cross(e2);
        float a  = e1.dot(h);
        if (std::abs(a) > 1e-7f) {
            float f  = 1.f / a;
            Vec3  s  = orig - verts[0];
            float u  = f * s.dot(h);
            if (u >= 0.f && u <= 1.f) {
                Vec3  q  = s.cross(e1);
                float vv = f * dir.dot(q);
                if (vv >= 0.f && u + vv <= 1.f) {
                    float t = f * e2.dot(q);
                    if (t > 1e-4f) return t;
                }
            }
        }
    }
    // Triangle 1: v0,v2,v3
    {
        Vec3  e1 = verts[2] - verts[0], e2 = verts[3] - verts[0];
        Vec3  h  = dir.cross(e2);
        float a  = e1.dot(h);
        if (std::abs(a) > 1e-7f) {
            float f  = 1.f / a;
            Vec3  s  = orig - verts[0];
            float u  = f * s.dot(h);
            if (u >= 0.f && u <= 1.f) {
                Vec3  q  = s.cross(e1);
                float vv = f * dir.dot(q);
                if (vv >= 0.f && u + vv <= 1.f) {
                    float t = f * e2.dot(q);
                    if (t > 1e-4f) return t;
                }
            }
        }
    }
    return -1.f;
}

// Build an orthonormal frame (right, up) perpendicular to n.
static void makeFrame(Vec3 n, Vec3& right, Vec3& up)
{
    Vec3 tmp = (std::abs(n.x) < 0.9f) ? Vec3{1,0,0} : Vec3{0,1,0};
    right = n.cross(tmp).normalized();
    up    = right.cross(n);
}

// Delta form-factor weight for the top face pixel at (u,v) in [-1,1].
// dF ∝ 1/(u²+v²+1)²
static float weightTop(float u, float v)
{
    float r2 = u*u + v*v + 1.0f;
    return 1.0f / (r2 * r2);
}

// Delta form-factor weight for a side face pixel; w ≥ 0 faces the horizon.
static float weightSide(float v, float w)
{
    if (w < 0.f) return 0.f;
    float r2 = 1.0f + v*v + w*w;
    return w / (r2 * r2);
}

std::vector<float> computeFormFactors(const Scene& scene, int res)
{
    const int n = static_cast<int>(scene.patches.size());
    std::vector<float> F(static_cast<std::size_t>(n * n), 0.f);

    const float dp = 2.0f / static_cast<float>(res);

    for (int i = 0; i < n; ++i) {
        const Patch& pi = scene.patches[i];
        Vec3 orig = pi.center;
        Vec3 N    = pi.normal;
        Vec3 R, U;
        makeFrame(N, R, U);

        float total_weight = 0.f;

        // Top face (+N direction): full res×res grid
        for (int py = 0; py < res; ++py) {
            for (int px = 0; px < res; ++px) {
                float u    = -1.f + (px + 0.5f) * dp;
                float v    = -1.f + (py + 0.5f) * dp;
                float w_px = weightTop(u, v);
                Vec3  dir  = (R*u + U*v + N).normalized();

                float t_min = 1e30f;
                int   hit   = -1;
                for (int j = 0; j < n; ++j) {
                    if (j == i) continue;
                    float t = rayQuad(orig, dir, scene.patches[j].verts);
                    if (t > 0.f && t < t_min) { t_min = t; hit = j; }
                }
                if (hit >= 0) F[i*n + hit] += w_px;
                total_weight += w_px;
            }
        }

        // Four side faces (+R, -R, +U, -U): upper half only (res × res/2)
        const Vec3 sideN[4]  = { R,    -R,   U,    -U  };
        const Vec3 sideT1[4] = { N,     N,   N,     N  };
        const Vec3 sideT2[4] = { U,     U,   R,     R  };

        for (int face = 0; face < 4; ++face) {
            for (int py = 0; py < res / 2; ++py) {
                for (int px = 0; px < res; ++px) {
                    float t1   = -1.f + (px + 0.5f) * dp;
                    float t2   =        (py + 0.5f) * dp;
                    float w_px = weightSide(t1, t2);
                    Vec3  dir  = (sideT1[face]*t1 + sideT2[face]*t2 + sideN[face]).normalized();

                    float t_min = 1e30f;
                    int   hit   = -1;
                    for (int j = 0; j < n; ++j) {
                        if (j == i) continue;
                        float t = rayQuad(orig, dir, scene.patches[j].verts);
                        if (t > 0.f && t < t_min) { t_min = t; hit = j; }
                    }
                    if (hit >= 0) F[i*n + hit] += w_px;
                    total_weight += w_px;
                }
            }
        }

        // Normalise row i
        if (total_weight > 1e-9f) {
            float inv = 1.f / total_weight;
            for (int j = 0; j < n; ++j)
                F[i*n + j] *= inv;
        }
    }

    return F;
}

} // namespace Hemicube
