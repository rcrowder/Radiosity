#pragma once

#include "scene.h"

namespace Render {
    void init();
    void reshape(int w, int h);
    void draw(const Scene& scene);
    void setWireframe(bool enable);
    // Save the current framebuffer as a NetPBM (PPM) image file.
    void saveScreenshot(const char* filename, int w, int h);
}
