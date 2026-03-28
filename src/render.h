#pragma once

#include "scene.h"

namespace Render {
    void init();
    void reshape(int w, int h);
    void draw(const Scene& scene);
    // Toggle wireframe overlay: when enabled, patch edges are drawn over the
    // shaded scene, useful for inspecting mesh density and discontinuity lines.
    void setWireframe(bool enable);
    // Save the current framebuffer as a NetPBM (PPM) image file.
    void saveScreenshot(const char* filename, int w, int h);
}
