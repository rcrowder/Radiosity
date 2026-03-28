#pragma once
#include "scene.h"

namespace Teapot {

// Tessellate the Utah/Newell teapot (32 bicubic Bézier patches) and add it
// as radiosity scene geometry.
//
// Parameters:
//   scene   — target scene; faces[] and patches[] are extended in-place.
//   origin  — 3-D position of the teapot base centre (bottom of the body).
//   scale   — uniform scale applied to the canonical teapot coordinates.
//             The canonical teapot fits in a ≈ 2×2×3 unit bounding box;
//             scale=0.065 gives a ~130 mm tall teapot at Cornell scale.
//   divs    — tessellation subdivisions per Bézier patch side (default 6).
//             Each of the 32 patches produces divs² quads, so the total
//             patch count is 32 × divs².
//   reflectance — surface reflectance of the teapot material.
//   is_box  — passed to Face::is_box (set true so the teapot is included in
//             the box-geometry rendering pass and polygon-offset correctly).
void addToScene(Scene& scene,
                Vec3  origin,
                float scale,
                Vec3  reflectance,
                int   divs   = 6,
                bool  is_box = true);

} // namespace Teapot
