#include <cmath>
#include <unordered_map>
#include <array>
#include <cstring>
#include <cstdio>
#include <vector>
#include <GL/gl.h>
#include <GL/glu.h>
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
// Vertex-radiosity averaging
// The key combines quantised position AND quantised normal so that corner
// vertices shared between two differently-oriented faces (e.g. wall/ceiling)
// are treated as independent vertices for each face.  This prevents the dark
// ceiling corners from bleeding into the top of the walls (Gouraud artefact).
// ---------------------------------------------------------------------------
struct VKey {
    int x, y, z;       // quantised position (1/10000 units)
    int nx, ny, nz;    // quantised normal   (1/1000 — normals are unit length)
    bool operator==(const VKey& o) const {
        return x==o.x && y==o.y && z==o.z
            && nx==o.nx && ny==o.ny && nz==o.nz;
    }
};
struct VKeyHash {
    std::size_t operator()(const VKey& k) const {
        std::size_t h = 2166136261u;
        auto mix = [&](int v){ h ^= static_cast<std::size_t>(v); h *= 16777619u; };
        mix(k.x); mix(k.y); mix(k.z);
        mix(k.nx); mix(k.ny); mix(k.nz);
        return h;
    }
};

struct RadSum { Vec3 sum{}; int count{}; };

static VKey toKey(const Vec3& v, const Vec3& n)
{
    return { static_cast<int>(std::round(v.x * 10000)),
             static_cast<int>(std::round(v.y * 10000)),
             static_cast<int>(std::round(v.z * 10000)),
             static_cast<int>(std::round(n.x * 1000)),
             static_cast<int>(std::round(n.y * 1000)),
             static_cast<int>(std::round(n.z * 1000)) };
}

// Build the vertex→averaged-radiosity table from current patch data.
static std::unordered_map<VKey, Vec3, VKeyHash>
buildVertexColours(const Scene& scene)
{
    std::unordered_map<VKey, RadSum, VKeyHash> acc;
    acc.reserve(scene.patches.size() * 4);

    for (const auto& p : scene.patches) {
        for (int i = 0; i < 4; ++i) {
            auto& entry = acc[toKey(p.verts[i], p.normal)];
            entry.sum.x += p.radiosity.x;
            entry.sum.y += p.radiosity.y;
            entry.sum.z += p.radiosity.z;
            entry.count++;
        }
    }

    std::unordered_map<VKey, Vec3, VKeyHash> result;
    result.reserve(acc.size());
    for (auto& [key, rs] : acc) {
        float inv = 1.f / rs.count;
        result[key] = rs.sum * inv;
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
    glPolygonMode(GL_FRONT_AND_BACK, g_wireframe ? GL_LINE : GL_FILL);

    // Build per-vertex averaged radiosity (Gouraud colours) each frame.
    // For a static solved scene this could be cached, but the cost is trivial
    // compared to the solve itself.
    auto vc = buildVertexColours(scene);

    auto vertColour = [&](const Vec3& v, const Vec3& n) {
        auto it = vc.find(toKey(v, n));
        if (it != vc.end()) {
            glColor3f(tonemap(it->second.x),
                      tonemap(it->second.y),
                      tonemap(it->second.z));
        }
    };

    // Helper: emit both triangles of a quad patch.
    auto drawPatch = [&](const Patch& p) {
        vertColour(p.verts[0], p.normal); glVertex3fv(&p.verts[0].x);
        vertColour(p.verts[1], p.normal); glVertex3fv(&p.verts[1].x);
        vertColour(p.verts[2], p.normal); glVertex3fv(&p.verts[2].x);
        vertColour(p.verts[0], p.normal); glVertex3fv(&p.verts[0].x);
        vertColour(p.verts[2], p.normal); glVertex3fv(&p.verts[2].x);
        vertColour(p.verts[3], p.normal); glVertex3fv(&p.verts[3].x);
    };

    const int n  = static_cast<int>(scene.patches.size());
    const int bs = scene.box_start;

    // Pass 1: environment geometry (floor, ceiling, walls, light).
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < bs; ++i)
        drawPatch(scene.patches[i]);
    glEnd();

    // Pass 2: box geometry, pulled slightly toward the camera via polygon
    // offset so box sides win the depth test at the y=0 seam with the floor,
    // eliminating the z-fighting artefact that makes boxes appear to float.
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(-1.0f, -1.0f);
    glBegin(GL_TRIANGLES);
    for (int i = bs; i < n; ++i)
        drawPatch(scene.patches[i]);
    glEnd();
    glDisable(GL_POLYGON_OFFSET_FILL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // restore default
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

