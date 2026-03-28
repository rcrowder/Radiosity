// GPU FBO Hemicube — replaces the CPU Möller-Trumbore ray-caster with OpenGL
// rasterization.  For each source patch i, five 90° perspective renders are
// made into an offscreen RGB8 FBO (top + four side faces of the hemicube);
// each other patch is drawn with its index encoded as a 24-bit RGB colour.
// Reading back the FBO pixels and summing precomputed delta-weights replicates
// Cohen & Greenberg 1985 while moving the O(n²·res²) work to the GPU.
//
// Requires an active OpenGL context (call after glutCreateWindow).
#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#include "hemicube.h"
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <vector>

namespace Hemicube {

// Encode 0-based patch index j as 24-bit RGB (1-based so that all-zeros
// black = background sentinel, not patch 0).
static void encodeIdx(int j, GLubyte& r, GLubyte& g, GLubyte& b)
{
    unsigned idx = static_cast<unsigned>(j + 1);
    r = static_cast<GLubyte>( idx        & 0xFFu);
    g = static_cast<GLubyte>((idx >>  8) & 0xFFu);
    b = static_cast<GLubyte>((idx >> 16) & 0xFFu);
}

// Decode a readback pixel to a 0-based patch index; -1 = background.
static int decodePixel(GLubyte r, GLubyte g, GLubyte b)
{
    int idx = static_cast<int>(r)
            | (static_cast<int>(g) << 8)
            | (static_cast<int>(b) << 16);
    return idx - 1;
}

// GPU vertex layout: position (12 B) + patch-ID colour (3 B) + 1 B pad = 16 B.
struct GpuVertex {
    float   x, y, z;
    GLubyte r, g, b, _pad;
};
static_assert(sizeof(GpuVertex)      == 16, "GpuVertex size mismatch");
static_assert(offsetof(GpuVertex, r) == 12, "GpuVertex colour offset mismatch");

// Build an orthonormal frame (right, up) perpendicular to n.
static void makeFrame(Vec3 n, Vec3& right, Vec3& up)
{
    Vec3 tmp = (std::abs(n.x) < 0.9f) ? Vec3{1,0,0} : Vec3{0,1,0};
    right = n.cross(tmp).normalized();
    up    = right.cross(n);
}


std::vector<float> computeFormFactors(const Scene& scene, int res)
{
    const int n = static_cast<int>(scene.patches.size());
    std::vector<float> F(static_cast<std::size_t>(n) * n, 0.f);

    // ---- Build VBO: 6 vertices (2 × GL_TRIANGLES) per patch ----
    // All 6 vertices share the same ID colour; GL_SMOOTH with uniform
    // per-patch colours produces no interpolation so pixel == vertex colour.
    std::vector<GpuVertex> gpuVerts;
    gpuVerts.reserve(static_cast<std::size_t>(n) * 6);
    for (int j = 0; j < n; ++j) {
        GLubyte r, g, b;
        encodeIdx(j, r, g, b);
        const Vec3* v = scene.patches[j].verts;
        // Triangle 0: v0–v1–v2
        gpuVerts.push_back({v[0].x, v[0].y, v[0].z, r, g, b, 0});
        gpuVerts.push_back({v[1].x, v[1].y, v[1].z, r, g, b, 0});
        gpuVerts.push_back({v[2].x, v[2].y, v[2].z, r, g, b, 0});
        // Triangle 1: v0–v2–v3
        gpuVerts.push_back({v[0].x, v[0].y, v[0].z, r, g, b, 0});
        gpuVerts.push_back({v[2].x, v[2].y, v[2].z, r, g, b, 0});
        gpuVerts.push_back({v[3].x, v[3].y, v[3].z, r, g, b, 0});
    }
    GLuint vbo = 0;
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(gpuVerts.size() * sizeof(GpuVertex)),
                 gpuVerts.data(), GL_STATIC_DRAW);

    // ---- Build offscreen FBO: res×res GL_RGB8 colour + depth24 ----
    GLuint colorTex = 0, depthRB = 0, fbo = 0;
    glGenTextures(1, &colorTex);
    glBindTexture(GL_TEXTURE_2D, colorTex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, res, res, 0,
                 GL_RGB, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);

    glGenRenderbuffers(1, &depthRB);
    glBindRenderbuffer(GL_RENDERBUFFER, depthRB);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, res, res);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);

    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                           GL_TEXTURE_2D, colorTex, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                              GL_RENDERBUFFER, depthRB);
    {
        GLenum st = glCheckFramebufferStatus(GL_FRAMEBUFFER);
        if (st != GL_FRAMEBUFFER_COMPLETE) {
            std::fprintf(stderr, "Hemicube: FBO incomplete (status 0x%x)\n", st);
            std::exit(1);
        }
    }

    // ---- GL state for ID rendering ----
    glViewport(0, 0, res, res);
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);
    glDisable(GL_DITHER);
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    // Uniform per-patch vertex colours → no interpolation under GL_SMOOTH.
    glShadeModel(GL_SMOOTH);
    // Disable sRGB conversion so ID bytes survive readback unmodified.
#ifdef GL_FRAMEBUFFER_SRGB
    glDisable(GL_FRAMEBUFFER_SRGB);
#endif

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glVertexPointer(3, GL_FLOAT,         sizeof(GpuVertex),
                    reinterpret_cast<const void*>(offsetof(GpuVertex, x)));
    glColorPointer (3, GL_UNSIGNED_BYTE, sizeof(GpuVertex),
                    reinterpret_cast<const void*>(offsetof(GpuVertex, r)));

    // Projection: 90° square FOV; near=1e-3 clips the self-patch (which lies
    // at z≈0 in camera space for all 5 views).
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(90.0, 1.0, 1e-3, 100.0);
    glMatrixMode(GL_MODELVIEW);

    // ---- Precompute delta-weight tables ----
    // top_w[py*res+px]  — Cohen & Greenberg top-face kernel: 1/(u²+v²+1)²
    // side_w[py*res+px] — side-face kernel for the upper half of the image:
    //   t2 = screen-y in (0,1] maps to the elevation (N-component), i.e.
    //   pixel row (half+py) in the readback buffer (see below).
    const float dp   = 2.0f / static_cast<float>(res);
    const int   half = res / 2;

    std::vector<float> top_w(static_cast<std::size_t>(res * res));
    for (int py = 0; py < res; ++py)
        for (int px = 0; px < res; ++px) {
            float u  = -1.f + (px + 0.5f) * dp;
            float v  = -1.f + (py + 0.5f) * dp;
            float r2 = u*u + v*v + 1.f;
            top_w[py * res + px] = 1.f / (r2 * r2);
        }

    // For each side face the camera up-vector is N, so the upper half of the image
    // (rows [half, res-1] in the readback) has positive N-component (hemisphere).
    // side_w[py*res+px] carries the weight for readback row (half+py):
    //   t2 = (py+0.5)*dp  corresponds to the screen-y of row (half+py).
    std::vector<float> side_w(static_cast<std::size_t>(half * res));
    for (int py = 0; py < half; ++py)
        for (int px = 0; px < res; ++px) {
            float t1 = -1.f + (px + 0.5f) * dp;
            float t2 =        (py + 0.5f) * dp;  // > 0
            float r2 = 1.f + t1*t1 + t2*t2;
            side_w[py * res + px] = t2 / (r2 * r2);
        }

    std::vector<GLubyte> pixels(static_cast<std::size_t>(res * res * 3));

    // ---- Progress bar ----
    const int BAR_WIDTH = 50;
    int last_filled = -1;

    // ---- Main loop: one hemicube per source patch ----
    for (int i = 0; i < n; ++i) {
        int filled = (i + 1) * BAR_WIDTH / n;
        if (filled != last_filled) {
            last_filled = filled;
            int pct = (i + 1) * 100 / n;
            std::printf("\r  [");
            for (int b = 0; b < BAR_WIDTH; ++b)
                std::putchar(b < filled ? '#' : '.');
            std::printf("] %3d%%  (%d/%d patches)", pct, i + 1, n);
            std::fflush(stdout);
        }

        const Patch& pi = scene.patches[i];
        Vec3 N = pi.normal;
        Vec3 c = pi.center;
        Vec3 R, U;
        makeFrame(N, R, U);

        float total_w = 0.f;

        // 5 hemicube face cameras, all with up=N (verified: camera y-axis = N
        // for all five views, so the upper half of the rendered image always
        // corresponds to directions with positive N-component).
        // For top face: all res rows are live (start_row=0, rows=res).
        // For side faces: only upper half rows (start_row=half, rows=half).
        struct View { Vec3 look, up; int rows, start_row; const float* wt; };
        const View views[5] = {
            {  N,  U,    res,  0,    top_w.data()  },  // top (+N)
            {  R,  N,    half, half, side_w.data() },  // +R
            { -R,  N,    half, half, side_w.data() },  // -R
            {  U,  N,    half, half, side_w.data() },  // +U (tangent-2)
            { -U,  N,    half, half, side_w.data() },  // -U
        };

        for (const auto& view : views) {
            Vec3 at = c + view.look;
            glLoadIdentity();
            gluLookAt(c.x,  c.y,  c.z,
                      at.x, at.y, at.z,
                      view.up.x, view.up.y, view.up.z);

            glClearColor(0.f, 0.f, 0.f, 1.f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            // Draw all patches except i (skip source patch via two ranges).
            if (i > 0)
                glDrawArrays(GL_TRIANGLES, 0, i * 6);
            if (i + 1 < n)
                glDrawArrays(GL_TRIANGLES, (i + 1) * 6, (n - i - 1) * 6);

            // Synchronous readback.
            glReadPixels(0, 0, res, res, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());

            // Accumulate from the live rows (start_row..start_row+rows-1).
            const float* wt        = view.wt;
            const int    rows      = view.rows;
            const int    start_row = view.start_row;
            for (int py = 0; py < rows; ++py) {
                int py_full = start_row + py;
                for (int px = 0; px < res; ++px) {
                    const std::size_t pidx =
                        static_cast<std::size_t>((py_full * res + px) * 3);
                    int j = decodePixel(pixels[pidx],
                                        pixels[pidx + 1],
                                        pixels[pidx + 2]);
                    if (j < 0 || j >= n) continue;
                    float w = wt[py * res + px];
                    F[static_cast<std::size_t>(i) * n + j] += w;
                    total_w += w;
                }
            }
        }

        // Normalise row i so form factors sum to ≤ 1.
        if (total_w > 1e-9f) {
            float inv = 1.f / total_w;
            for (int j = 0; j < n; ++j)
                F[static_cast<std::size_t>(i) * n + j] *= inv;
        }
    }

    // Finish progress bar.
    std::printf("\r  [");
    for (int b = 0; b < BAR_WIDTH; ++b) std::putchar('#');
    std::printf("] 100%%  (%d/%d patches)\n", n, n);

    // ---- Reciprocity enforcement pass ----
    // The hemicube rasterises 64×64 pixels per face.  Patches smaller than one
    // pixel in solid angle register zero samples and end up with col_sum ≈ 0
    // even though their row_sum ≈ 1 (their own hemicube works fine).  This
    // violates reciprocity (A_i·F_ij = A_j·F_ji) and leaves those patches
    // permanently dark.
    //
    // Fix: for every pair (i,j) where F[i][j] > 0 but F[j][i] == 0, derive
    // F[j][i] via reciprocity:  F[j][i] = F[i][j] * A_i / A_j.
    // After filling in missing reciprocal terms, renormalise each row to
    // ≤ 1 so the energy-conservation budget is not violated.
    {
        for (int i = 0; i < n; ++i) {
            const float Ai = scene.patches[i].area;
            for (int j = i + 1; j < n; ++j) {
                float& fij = F[static_cast<std::size_t>(i) * n + j];
                float& fji = F[static_cast<std::size_t>(j) * n + i];
                if (fij > 0.f && fji == 0.f) {
                    const float Aj = scene.patches[j].area;
                    if (Aj > 1e-10f)
                        fji = fij * Ai / Aj;
                } else if (fji > 0.f && fij == 0.f) {
                    const float Aj = scene.patches[j].area;
                    if (Ai > 1e-10f)
                        fij = fji * Aj / Ai;
                }
            }
        }
        // Renormalise each row so F[i][*] sums to ≤ 1.
        for (int i = 0; i < n; ++i) {
            float row_sum = 0.f;
            for (int j = 0; j < n; ++j)
                row_sum += F[static_cast<std::size_t>(i) * n + j];
            if (row_sum > 1.f) {
                float inv = 1.f / row_sum;
                for (int j = 0; j < n; ++j)
                    F[static_cast<std::size_t>(i) * n + j] *= inv;
            }
        }
    }

    // ---- Restore default GL state ----
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    // Viewport will be reset by Render::reshape via the GLUT resize callback.

    // ---- Free GPU resources ----
    glDeleteBuffers(1, &vbo);
    glDeleteTextures(1, &colorTex);
    glDeleteRenderbuffers(1, &depthRB);
    glDeleteFramebuffers(1, &fbo);

    return F;
}

} // namespace Hemicube

