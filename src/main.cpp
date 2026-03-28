#include <cstdio>
#include <cstdlib>
#include <GL/glut.h>
#include "scene.h"
#include "hemicube.h"
#include "radiosity.h"
#include "render.h"

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
    // Build scene
    std::printf("Building Cornell Box scene...\n");
    g_scene.buildCornellBox();
    std::printf("  %zu patches\n", g_scene.patches.size());

    // Compute form factors via hemicube
    constexpr int hemi_res = 64;
    std::printf("Computing form factors (hemicube res=%d)...\n", hemi_res);
    auto F = Hemicube::computeFormFactors(g_scene, hemi_res);

    // Print form-factor matrix (useful for verification)
    const int n = static_cast<int>(g_scene.patches.size());
    /*
    std::printf("Form factor matrix F[i][j]:\n");
    for (int i = 0; i < n; ++i) {
        std::printf("  patch %d: ", i);
        float row_sum = 0.f;
        for (int j = 0; j < n; ++j) {
            std::printf("%.3f ", F[i*n+j]);
            row_sum += F[i*n+j];
        }
        std::printf(" (sum=%.3f)\n", row_sum);
    }
    */

    // Solve radiosity with progressive refinement
    std::printf("Solving radiosity (progressive refinement)...\n");
    Radiosity::solve(g_scene, F);

    // Print per-face averaged radiosity.
    // Each logical face maps to a contiguous run of sub-patches produced by
    // addSubdividedQuad (W×W for walls, B×B for box faces, 1 for the light).
    constexpr int W = 10, B = 4;
    struct FaceDesc { const char* name; int count; };
    const FaceDesc faces[] = {
        {"Floor",        W*W}, {"Ceiling",      W*W}, {"Back wall",    W*W},
        {"Right(green)", W*W}, {"Left(red)",    W*W}, {"Light",          1},
        {"Short top",    B*B}, {"Short front",  B*B}, {"Short right",  B*B},
        {"Short back",   B*B}, {"Short left",   B*B},
        {"Tall top",     B*B}, {"Tall front",   B*B}, {"Tall right",   B*B},
        {"Tall back",    B*B}, {"Tall left",    B*B},
    };
    std::printf("Per-face average radiosity:\n");
    int idx = 0;
    for (const auto& f : faces) {
        Vec3 avg{};
        for (int k = 0; k < f.count && idx + k < n; ++k)
            avg = avg + g_scene.patches[idx + k].radiosity;
        avg = avg * (1.f / f.count);
        std::printf("  %-14s  B=(%.3f, %.3f, %.3f)\n", f.name, avg.x, avg.y, avg.z);
        idx += f.count;
    }

    // Open GL window
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(640, 480);
    glutCreateWindow("Radiosity — Cohen et al. 1985/1988");

    Render::init();
    Render::reshape(640, 480);

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);

    std::printf("\nPress ESC or 'q' to quit.  Press 's' to save a screenshot (output.ppm).\n");
    std::printf("A screenshot will also be saved automatically on the first frame.\n");
    glutMainLoop();
    return 0;
}
