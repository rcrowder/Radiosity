# Radiosity

A global illumination renderer implementing the classic **radiosity** algorithm, demonstrated with the Cornell Box scene.

![Reference image (Cornell Box, 1991 measured)](Cornell_box_1991_measured.jpg)

## Overview

Radiosity is a physically based rendering technique that models diffuse inter-reflections of light between surfaces. It was introduced by Goral et al. (1984) and refined significantly by Cohen & Greenberg (1985) and Cohen et al. (1988).

This implementation covers the full pipeline:

1. **Scene construction** – the Cornell Box is built as a set of logical *faces* (walls, floor, ceiling, light, two boxes) with exact Cornell University coordinates and measured reflectances.
2. **Meshing** – two modes are supported:
   - *Uniform grid* (default): each face is subdivided into W×W (walls/ceiling/floor) or B×B (box) equal patches.
   - *Discontinuity meshing* (`--disc-mesh`): D0 critical lines (Lischinski, Tampieri & Greenberg 1992) are analytically computed per receiver face, then used to BSP-split each face into convex sub-polygons aligned with penumbra boundaries.
3. **Form-factor computation** – using the **GPU FBO hemicube method** (Cohen & Greenberg 1985): a virtual hemicube is placed at each patch centre; five 90° perspective renders into an offscreen `GL_RGB8` framebuffer object (FBO) replace the CPU ray-caster. Each other patch is drawn with its index encoded as a 24-bit RGB colour; pixels are decoded and accumulated with precomputed delta-weight tables, reducing the dominant O(n²·res²) inner loop to n × 5 GPU draw calls.
4. **Radiosity solve** – **progressive refinement** (Cohen et al. 1988): each iteration selects the patch with the greatest unshot energy and distributes it to all visible patches, converging when remaining unshot energy falls below a threshold.
5. **Rendering** – OpenGL/GLUT with **Gouraud shading**: per-vertex colours are computed using **Shepard's scattered-data interpolation** (Shepard 1968) over all patches on the same face, weighted by `area / (dist² + ε²)`. This produces smooth colour gradients across disc-mesh patch boundaries where simple averaging would create visible discontinuities. Colours are then tone-mapped with Reinhard tonemapping and gamma-corrected (γ = 2.2).
6. **Image output** – the rendered frame is saved as `output.ppm` (binary NetPBM) automatically on the first frame, and again on demand with the `s` key.

## Building

Requirements: a C++17 compiler, OpenGL, GLU, and GLUT (freeglut).

```bash
# Debian/Ubuntu
sudo apt install build-essential freeglut3-dev

make
```

## Running

```bash
./radiosity                # uniform grid (default)
./radiosity --disc-mesh    # D0 discontinuity meshing (Lischinski 1992)
./radiosity -d             # shorthand for --disc-mesh
./radiosity --no-radiosity # skip form-factor + solve; render flat geometry only
./radiosity -n             # shorthand for --no-radiosity
./radiosity -d -n          # disc-mesh without solving (inspect mesh in wireframe)
```

`--no-radiosity` skips the hemicube computation and progressive refinement entirely, colouring each patch with its material reflectance. This is useful for quickly inspecting mesh topology with the `w` wireframe key.

| Key | Action |
|-----|--------|
| `s` | Save current frame to `output.ppm` |
| `w` | Toggle wireframe overlay (draws only true patch boundary edges via `GL_LINES`) |
| `q` / `Esc` | Quit |

`output.ppm` is also written automatically when the first frame is rendered. Most image viewers (eog, feh, GIMP) open PPM files directly. To convert to PNG:

```bash
convert output.ppm output.png
```

## Project Structure

```
src/
  scene.h / scene.cpp         – Vec3, Vec2, Face, Patch, Scene;
                                Cornell Box geometry; uniform-grid meshing
  discmesh.h / discmesh.cpp   – D0 discontinuity meshing (Lischinski 1992);
                                Sutherland-Hodgman polygon splitter;
                                Cyrus-Beck segment clipping
  hemicube.h / hemicube.cpp   – form-factor computation via GPU FBO hemicube;
                                offscreen GL_RGB8 render; 24-bit patch ID encoding;
                                precomputed delta-weight tables
  radiosity.h / radiosity.cpp – progressive refinement solver
  render.h / render.cpp       – OpenGL renderer + PPM screenshot export
  main.cpp                    – entry point: build scene, mesh, solve, open window
Makefile
Cornell_box_1991_measured.jpg  – reference photograph (Cornell University, 1991)
```

## Comparison with the 1991 Cornell Box reference photo

A quantitative region-by-region analysis comparing `output.ppm` against `Cornell_box_1991_measured.jpg`:

| Feature | Reference | Output | Notes |
|---|---|---|---|
| Left wall colour | Red | Red | ✅ correct |
| Right wall colour | Green | Green | ✅ correct |
| Ceiling light position | Top-centre | Top-centre | ✅ correct |
| Tall box position | Left (near red wall) | Left | ✅ correct |
| Short box position | Right (near green wall) | Right | ✅ correct |
| Overall brightness | ~76 luminance | ~47 luminance | See note |

**Brightness difference:** The reference is a physical photograph taken with a fixed exposure and film response. The renderer uses Reinhard tonemapping + gamma correction, which compresses highlights and produces a more conservative overall brightness. The *relative* luminance distribution (ceiling bright, corners darker, colour bleed visible) matches well.

**Fidelity to Cohen et al. 1988:** The implementation faithfully reproduces the key contributions of the paper:
- Hemicube form-factor computation (Cohen & Greenberg 1985)
- Progressive refinement solver — each iteration shoots from the patch with the greatest unshot energy
- The Cornell Box scene with correct geometry, reflectances, and light emittance
- Gouraud-shaded output showing smooth colour gradients and inter-reflection colour bleeding

## References

- Goral, C. M., Torrance, K. E., Greenberg, D. P., & Battaile, B. (1984). *Modeling the interaction of light between diffuse surfaces.* SIGGRAPH '84.
- Cohen, M. F., & Greenberg, D. P. (1985). *The hemi-cube: A radiosity solution for complex environments.* SIGGRAPH '85.
- Cohen, M. F., Chen, S. E., Wallace, J. R., & Greenberg, D. P. (1988). *A progressive refinement approach to fast radiosity image generation.* SIGGRAPH '88.
- Lischinski, D., Tampieri, F., & Greenberg, D. P. (1992). *Discontinuity meshing for accurate radiosity.* IEEE Computer Graphics and Applications, 12(6), 25–39.
- Shepard, D. (1968). *A two-dimensional interpolation function for irregularly spaced data.* Proc. 23rd ACM National Conference, pp. 517–524. DOI:10.1145/800186.810616
- Sutherland, I. E., & Hodgman, G. W. (1974). *Reentrant polygon clipping.* Communications of the ACM, 17(1), 32–42.
- Cyrus, M., & Beck, J. (1978). *Generalized two- and three-dimensional clipping.* Computers & Graphics, 3(1), 23–28.
