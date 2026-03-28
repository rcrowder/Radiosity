#include <cmath>
#include <unordered_map>
#include <array>
#include <cstring>
#include <cstdio>
#include <vector>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include "render.h"

/*
 * OpenGL renderer with Reinhard tone-mapping + gamma correction.
 * True Gouraud shading: per-vertex colours computed by averaging the
 * radiosity of every patch that shares each vertex (Cohen et al. 1988).
 */

namespace Render {

static float tonemap(float x)
{
    x = x / (1.0f + x);
    return std::pow(x < 0.f ? 0.f : x, 1.0f / 2.2f);
}

// ---------------------------------------------------------------------------
// Vertex-radiosity interpolation — Shepard's method (1968)
//
// Reference: D. Shepard, "A two-dimensional interpolation function for
//   irregularly spaced data", Proc. 23rd ACM National Conference, 1968,
//   pp. 517–524.  DOI:10.1145/800186.810616
//
// Goal: assign a smooth radiosity colour to every mesh vertex for Gouraud
// interpolation, using the solved per-patch radiosity values.
//
// Why not simple per-patch averaging?
//   1. Triangle patches (v3 == v2) would double-count vertex 2 if iterated
//      over all 4 corners, biasing that vertex's colour toward the triangle's
//      own radiosity and sharpening shadow boundaries.
//   2. Fan-triangle outer-tip vertices (from disc-mesh polygon decomposition)
//      may belong to only one patch, so the colour would equal that patch's
//      radiosity with no blending — causing abrupt jumps across critical lines.
//
// Shepard's scattered-data interpolation:
//   Given a set of sample points {x_j} with values {B_j}, the interpolated
//   value at query point x is:
//
//     B(x) = Σ_j  w_j(x) · B_j  /  Σ_j  w_j(x)
//
//   with weights:
//
//     w_j(x) = A_j / ( ||x - x_j||²  +  ε² )
//
//   where A_j is the area of patch j and ε is a smoothing radius.
//   Multiplying by area A_j makes larger (better-sampled) patches dominate
//   and prevents tiny disc-mesh slivers from causing colour spikes.
//   The ε² term (modified Shepard, Hardy 1971) prevents division-by-zero
//   when the query vertex coincides with a patch centre, and controls the
//   spatial extent of each patch's influence.
//
// Scope: interpolation is done separately per face (patches grouped by
// face_id).  The VKey encodes position, surface normal, AND face_id, so
// corner vertices shared between two different faces are always treated as
// independent — even when the faces are co-planar and co-normal (e.g. the
// flush area light and its surrounding ceiling strips both have normal -Y
// and share 4 corner positions; without face_id in the key, whichever face
// is iterated first would claim those VKeys, and the light would render with
// ceiling radiosities instead of its own emission).
// ---------------------------------------------------------------------------
struct VKey {
    int x, y, z;       // quantised position (1/10000 units)
    int nx, ny, nz;    // quantised normal   (1/1000 — normals are unit length)
    int face_id;       // logical face index — disambiguates co-planar co-normal faces
    bool operator==(const VKey& o) const {
        return x==o.x && y==o.y && z==o.z
            && nx==o.nx && ny==o.ny && nz==o.nz
            && face_id==o.face_id;
    }
};
struct VKeyHash {
    std::size_t operator()(const VKey& k) const {
        std::size_t h = 2166136261u;
        auto mix = [&](int v){ h ^= static_cast<std::size_t>(v); h *= 16777619u; };
        mix(k.x); mix(k.y); mix(k.z);
        mix(k.nx); mix(k.ny); mix(k.nz);
        mix(k.face_id);
        return h;
    }
};

static VKey toKey(const Vec3& v, const Vec3& n, int face_id)
{
    return { static_cast<int>(std::round(v.x * 10000)),
             static_cast<int>(std::round(v.y * 10000)),
             static_cast<int>(std::round(v.z * 10000)),
             static_cast<int>(std::round(n.x * 1000)),
             static_cast<int>(std::round(n.y * 1000)),
             static_cast<int>(std::round(n.z * 1000)),
             face_id };
}

// Assign a Shepard-interpolated radiosity colour to every patch vertex.
// See the block comment above for the full algorithm derivation and references.
static std::unordered_map<VKey, Vec3, VKeyHash>
buildVertexColours(const Scene& scene)
{
    // Group patches by face for intra-face interpolation.
    std::unordered_map<int, std::vector<const Patch*>> by_face;
    by_face.reserve(32);
    for (const auto& p : scene.patches)
        by_face[p.face_id].push_back(&p);

    std::unordered_map<VKey, Vec3, VKeyHash> result;
    result.reserve(scene.patches.size() * 4);

    // ε² = (0.04)²: smoothing radius ≈ 4 % of a unit-wide face (~22 mm at
    // real Cornell scale).  Follows Hardy (1971) modified-Shepard convention
    // to avoid singularities when a vertex coincides with a patch centre.
    static constexpr float SMOOTH_EPS2 = 0.04f * 0.04f;

    for (const auto& p : scene.patches) {
        // Avoid double-counting v3==v2 for triangle patches.
        bool is_tri = (p.verts[2].x == p.verts[3].x &&
                       p.verts[2].y == p.verts[3].y &&
                       p.verts[2].z == p.verts[3].z);
        int nv = is_tri ? 3 : 4;

        const auto& fps = by_face.at(p.face_id);

        for (int i = 0; i < nv; ++i) {
            VKey key = toKey(p.verts[i], p.normal, p.face_id);
            if (result.count(key)) continue; // already computed

            // Shepard (1968): B(x) = Σ w_j·B_j / Σ w_j,  w_j = A_j/(dist²+ε²)
            Vec3  wsum{};
            float wtotal = 0.f;
            for (const Patch* q : fps) {
                float dx = q->center.x - p.verts[i].x;
                float dy = q->center.y - p.verts[i].y;
                float dz = q->center.z - p.verts[i].z;
                float w  = q->area / (dx*dx + dy*dy + dz*dz + SMOOTH_EPS2);
                wsum.x  += q->radiosity.x * w;
                wsum.y  += q->radiosity.y * w;
                wsum.z  += q->radiosity.z * w;
                wtotal  += w;
            }
            if (wtotal > 0.f)
                result[key] = wsum * (1.f / wtotal);
        }
    }
    return result;
}

// ---------------------------------------------------------------------------

void init()
{
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glClearColor(0.f, 0.f, 0.f, 1.f);
    glShadeModel(GL_SMOOTH);   // Gouraud interpolation across triangles
}

void reshape(int w, int h)
{
    if (h == 0) h = 1;
    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // FOV from Cornell dataset: 2*atan(film_half/focal_len) = 2*atan(12.5/35) ≈ 39.3°
    gluPerspective(39.3, static_cast<double>(w) / h, 0.01, 10.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    // Camera position/direction from Cornell dataset, normalised by 555.
    // Eye at (278, 273, -800)/555 looking in +z toward the back wall.
    gluLookAt(278.0/555, 273.0/555, -800.0/555,
              278.0/555, 273.0/555,  559.2/555,
              0.0, 1.0, 0.0);
}

static bool g_wireframe = false;
void setWireframe(bool enable) { g_wireframe = enable; }

void draw(const Scene& scene)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Build per-vertex averaged radiosity (Gouraud colours) each frame.
    // For a static solved scene this could be cached, but the cost is trivial
    // compared to the solve itself.
    auto vc = buildVertexColours(scene);

    auto vertColour = [&](const Vec3& v, const Vec3& n, int fid) {
        auto it = vc.find(toKey(v, n, fid));
        if (it != vc.end()) {
            glColor3f(tonemap(it->second.x),
                      tonemap(it->second.y),
                      tonemap(it->second.z));
        }
    };

    // Fill mode: two triangles per quad patch.
    auto drawPatchFill = [&](const Patch& p) {
        vertColour(p.verts[0], p.normal, p.face_id); glVertex3fv(&p.verts[0].x);
        vertColour(p.verts[1], p.normal, p.face_id); glVertex3fv(&p.verts[1].x);
        vertColour(p.verts[2], p.normal, p.face_id); glVertex3fv(&p.verts[2].x);
        vertColour(p.verts[0], p.normal, p.face_id); glVertex3fv(&p.verts[0].x);
        vertColour(p.verts[2], p.normal, p.face_id); glVertex3fv(&p.verts[2].x);
        vertColour(p.verts[3], p.normal, p.face_id); glVertex3fv(&p.verts[3].x);
    };

    // Wireframe mode: draw only the true patch boundary as explicit edges
    // (GL_LINES pairs).  A triangle patch has v3==v2; that gives 3 edges.
    // A quad patch gives 4 edges.  Using GL_LINES avoids the internal
    // v0→v2 diagonal that GL_TRIANGLES+GL_LINE mode would otherwise draw
    // across every patch, including unwanted "spokes" in fan-decomposed polygons.
    auto drawPatchWire = [&](const Patch& p) {
        bool is_tri = (p.verts[2].x == p.verts[3].x &&
                       p.verts[2].y == p.verts[3].y &&
                       p.verts[2].z == p.verts[3].z);
        // Edge v0→v1
        vertColour(p.verts[0], p.normal, p.face_id); glVertex3fv(&p.verts[0].x);
        vertColour(p.verts[1], p.normal, p.face_id); glVertex3fv(&p.verts[1].x);
        // Edge v1→v2
        vertColour(p.verts[1], p.normal, p.face_id); glVertex3fv(&p.verts[1].x);
        vertColour(p.verts[2], p.normal, p.face_id); glVertex3fv(&p.verts[2].x);
        if (is_tri) {
            // Edge v2→v0 (close triangle)
            vertColour(p.verts[2], p.normal, p.face_id); glVertex3fv(&p.verts[2].x);
            vertColour(p.verts[0], p.normal, p.face_id); glVertex3fv(&p.verts[0].x);
        } else {
            // Edge v2→v3
            vertColour(p.verts[2], p.normal, p.face_id); glVertex3fv(&p.verts[2].x);
            vertColour(p.verts[3], p.normal, p.face_id); glVertex3fv(&p.verts[3].x);
            // Edge v3→v0
            vertColour(p.verts[3], p.normal, p.face_id); glVertex3fv(&p.verts[3].x);
            vertColour(p.verts[0], p.normal, p.face_id); glVertex3fv(&p.verts[0].x);
        }
    };

    const int n  = static_cast<int>(scene.patches.size());
    const int bs = scene.box_start;

    if (g_wireframe) {
        // Draw patch outlines only — no polygon-fill, no depth offset.
        glBegin(GL_LINES);
        for (int i = 0; i < n; ++i)
            drawPatchWire(scene.patches[i]);
        glEnd();
    } else {
        // Pass 1: environment geometry (floor, ceiling, walls, light).
        glBegin(GL_TRIANGLES);
        for (int i = 0; i < bs; ++i)
            drawPatchFill(scene.patches[i]);
        glEnd();

        // Pass 2: box geometry, pulled slightly toward the camera via polygon
        // offset so box sides win the depth test at the y=0 seam with the floor,
        // eliminating the z-fighting artefact that makes boxes appear to float.
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(-1.0f, -1.0f);
        glBegin(GL_TRIANGLES);
        for (int i = bs; i < n; ++i)
            drawPatchFill(scene.patches[i]);
        glEnd();
        glDisable(GL_POLYGON_OFFSET_FILL);
    }
}

// ---------------------------------------------------------------------------
// Screenshot: write current framebuffer to a binary PPM file.
// OpenGL's pixel origin is bottom-left, so we flip rows vertically.
// ---------------------------------------------------------------------------
void saveScreenshot(const char* filename, int w, int h)
{
    std::vector<unsigned char> pixels(static_cast<std::size_t>(w * h * 3));
    glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());

    std::FILE* f = std::fopen(filename, "wb");
    if (!f) {
        std::fprintf(stderr, "saveScreenshot: cannot open '%s' for writing\n", filename);
        return;
    }
    std::fprintf(f, "P6\n%d %d\n255\n", w, h);
    // Write rows top-to-bottom (flip the bottom-left GL origin)
    for (int y = h - 1; y >= 0; --y)
        std::fwrite(&pixels[static_cast<std::size_t>(y * w * 3)], 3,
                    static_cast<std::size_t>(w), f);
    std::fclose(f);
    std::printf("Screenshot saved to '%s' (%dx%d)\n", filename, w, h);
}

} // namespace Render

