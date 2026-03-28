#include "teapot.h"

/*
 * Utah Teapot — captured triangle mesh
 * ----------------------------------------
 * Triangle mesh for 3 136 triangles captured from glutSolidTeapot(1.0) via
 * GL_FEEDBACK with an orthographic projection, preserving object-space
 * coordinates.
 *
 * GLUT teapot (size = 1.0) coordinate system:
 *   +Y  up
 *   Body bottom  y ~ -0.7875   (GLUT_Y_MIN)
 *   Body top     y ~ +0.7875
 *   Total height ~ 1.575 = Newell canonical height (3.15) / 2
 *   Spout on +X side; handle on -X side; body rotationally symmetric in Z.
 *
 * addToScene() applies a normalisation factor of 2.0 so that the externally
 * supplied `scale` parameter has the same meaning as in the original Bezier-
 * patch implementation (scale = 0.054 -> teapot ~ 94 mm tall in scene).
 */

// Include the auto-generated vertex table (TEAPOT_TRIS[][3], TEAPOT_NTRIS).
#include "teapot_tris.inc"

namespace Teapot {

void addToScene(Scene&  scene,
                Vec3    origin,
                float   scale,
                Vec3    reflectance,
                int     /*divs*/,   // unused; tessellation fixed by capture
                bool    /*is_box*/)
{
    // Body bottom is at y = -0.7875 in GLUT object space.  Shift so that the
    // bottom lands exactly at origin.y.  Multiply by NORM so that
    // scale == 0.054 gives the same physical height as the old Bezier code.
    constexpr float Y_MIN = -0.7875f;   // GLUT body bottom y
    constexpr float NORM  =  2.0f;      // GLUT height 1.575 -> Newell height 3.15

    const Vec3 black{};

    // Assign a grouping face_id from the current face count.  We do NOT push
    // a Face into scene.faces — doing so would cause meshUniform / DiscMesh to
    // re-mesh the registration quad into a flat horizontal disc and clear the
    // teapot triangle patches.  The Shepard interpolator in render.cpp builds
    // its by_face map directly from patches (not from faces[]), so a Face entry
    // is not required.
    const int face_id = static_cast<int>(scene.faces.size());

    // Transform: GLUT object coords -> scene coords
    //   x_scene =  v[0]          * NORM * scale + origin.x
    //   y_scene = (v[1] - Y_MIN) * NORM * scale + origin.y
    //   z_scene =  v[2]          * NORM * scale + origin.z
    const float ks = NORM * scale;
    auto xfm = [&](const float* v) -> Vec3 {
        return {
             v[0]          * ks + origin.x,
            (v[1] - Y_MIN) * ks + origin.y,
             v[2]          * ks + origin.z
        };
    };

    // Emit each captured triangle as a degenerate quad (v3 = v2).
    // With v3 = v2:
    //   normal = (v1-v0) x (v2-v0)   (correct triangle normal)
    //   area   = 0.5 * |(v2-v0) x (v1-v0)|  (correct triangle area)
    // Both computed correctly by Scene::addPatch.
    for (int i = 0; i < TEAPOT_NTRIS; ++i) {
        const float* a = TEAPOT_TRIS[i * 3 + 0];
        const float* b = TEAPOT_TRIS[i * 3 + 1];
        const float* c = TEAPOT_TRIS[i * 3 + 2];

        const Vec3 v0 = xfm(a);
        const Vec3 v1 = xfm(b);
        const Vec3 v2 = xfm(c);

        // Skip degenerate triangles (zero area after transform)
        if ((v1 - v0).cross(v2 - v0).length() < 1e-10f) continue;

        scene.addPatch(v0, v1, v2, v2, black, reflectance, face_id);
    }
}

} // namespace Teapot
