#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <GL/glut.h>
#include "scene.h"
#include "hemicube.h"
#include "radiosity.h"
#include "render.h"
#include "discmesh.h"

static Scene g_scene;
static bool  g_save_next_frame = true;   // auto-save on the first rendered frame

static void display()
{
    Render::draw(g_scene);
    if (g_save_next_frame) {
        Render::saveScreenshot("output.ppm",
                               glutGet(GLUT_WINDOW_WIDTH),
                               glutGet(GLUT_WINDOW_HEIGHT));
        g_save_next_frame = false;
    }
    glutSwapBuffers();
}

static void reshape(int w, int h)
{
    Render::reshape(w, h);
}

static void keyboard(unsigned char key, int /*x*/, int /*y*/)
{
    if (key == 27 || key == 'q') std::exit(0);   // ESC or q to quit
    if (key == 's') {                            // 's' to save a screenshot
        Render::saveScreenshot("output.ppm",
                               glutGet(GLUT_WINDOW_WIDTH),
                               glutGet(GLUT_WINDOW_HEIGHT));
    }
    if (key == 'w') {                            // 'w' to toggle wireframe
        static bool wf = false;
        wf = !wf;
        Render::setWireframe(wf);
        glutPostRedisplay();
    }
}

int main(int argc, char** argv)
{
    // Parse flags before passing argc/argv to GLUT
    bool use_disc_mesh  = false;
    bool skip_radiosity = false;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--disc-mesh") == 0 ||
            std::strcmp(argv[i], "-d") == 0)
            use_disc_mesh = true;
        if (std::strcmp(argv[i], "--no-radiosity") == 0 ||
            std::strcmp(argv[i], "-n") == 0)
            skip_radiosity = true;
    }

    // Build scene
    std::printf("Building Cornell Box scene...\n");
    g_scene.buildCornellBox();
    std::printf("  %zu logical faces\n", g_scene.faces.size());

    if (use_disc_mesh) {
        std::printf("Meshing: discontinuity meshing (D0 critical lines, Lischinski 1992)\n");
        DiscMesh::apply(g_scene);
    } else {
        std::printf("Meshing: uniform grid (W=10, B=4)\n");
        g_scene.meshUniform(10, 4);
    }
    std::printf("  %zu patches\n", g_scene.patches.size());

    // Initialise OpenGL now — computeFormFactors uses FBO and needs a GL context.
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(640, 480);
    glutCreateWindow(use_disc_mesh
        ? "Radiosity — disc-mesh (Lischinski 1992)"
        : "Radiosity — uniform mesh (Cohen 1988)");
    Render::init();

    if (skip_radiosity) {
        // No solve: use reflectance as flat colour so non-emissive surfaces are
        // visible (otherwise radiosity = emission = 0 → all black).
        for (auto& p : g_scene.patches)
            if (p.radiosity.x == 0.f && p.radiosity.y == 0.f && p.radiosity.z == 0.f)
                p.radiosity = p.reflectance;
        std::printf("Radiosity skipped (--no-radiosity); rendering geometry only.\n");
    } else {
        // Compute form factors via GPU FBO hemicube (requires active GL context).
        constexpr int hemi_res = 64;
        std::printf("Computing form factors (hemicube res=%d, GPU FBO)...\n", hemi_res);
        auto F = Hemicube::computeFormFactors(g_scene, hemi_res);

        // Solve radiosity with progressive refinement
        std::printf("Solving radiosity (progressive refinement)...\n");
        Radiosity::solve(g_scene, F);

        // Print per-face averaged radiosity, grouped by face_id.
        std::printf("Per-face average radiosity:\n");
        const int nf = static_cast<int>(g_scene.faces.size());
        std::vector<Vec3> face_sum(nf, Vec3{});
        std::vector<int>  face_cnt(nf, 0);
        for (const auto& p : g_scene.patches) {
            if (p.face_id >= 0 && p.face_id < nf) {
                face_sum[p.face_id] = face_sum[p.face_id] + p.radiosity;
                face_cnt[p.face_id]++;
            }
        }
        for (int fi = 0; fi < nf; ++fi) {
            if (face_cnt[fi] == 0) continue;
            Vec3 avg = face_sum[fi] * (1.f / face_cnt[fi]);
            std::printf("  %-14s  B=(%.3f, %.3f, %.3f)  (%d patches)\n",
                        g_scene.faces[fi].name.c_str(), avg.x, avg.y, avg.z, face_cnt[fi]);
        }
    }

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);

    std::printf("\nPress ESC or 'q' to quit.  's' = screenshot  'w' = wireframe\n");
    std::printf("Meshing mode:    %s\n",
                use_disc_mesh ? "discontinuity (--disc-mesh / -d)" : "uniform (default)");
    std::printf("Radiosity solve: %s\n",
                skip_radiosity ? "skipped (--no-radiosity / -n)" : "enabled");
    glutMainLoop();
    return 0;
}
