#pragma once

#include <cmath>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Vec3
// ---------------------------------------------------------------------------
struct Vec3 {
    float x{}, y{}, z{};

    constexpr Vec3() = default;
    constexpr Vec3(float x, float y, float z) : x(x), y(y), z(z) {}

    Vec3  operator+(const Vec3& b) const { return {x+b.x, y+b.y, z+b.z}; }
    Vec3  operator-(const Vec3& b) const { return {x-b.x, y-b.y, z-b.z}; }
    Vec3  operator-()              const { return {-x, -y, -z}; }
    Vec3  operator*(float s)       const { return {x*s,   y*s,   z*s};   }
    Vec3  operator*(const Vec3& b) const { return {x*b.x, y*b.y, z*b.z}; } // component-wise
    Vec3& operator+=(const Vec3& b)      { x+=b.x; y+=b.y; z+=b.z; return *this; }

    float dot  (const Vec3& b) const { return x*b.x + y*b.y + z*b.z; }
    Vec3  cross(const Vec3& b) const {
        return {y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x};
    }
    float length()     const { return std::sqrt(x*x + y*y + z*z); }
    Vec3  normalized() const {
        float l = length();
        return l > 1e-9f ? *this * (1.0f/l) : Vec3{};
    }
};

inline Vec3 operator*(float s, const Vec3& v) { return v * s; }

// ---------------------------------------------------------------------------
// Vec2  (face-local UV coordinates used by discontinuity meshing)
//
// u runs along verts[0]→verts[1] (local_u direction).
// v runs perpendicular in the face plane (local_v = normal × local_u).
// The frame (local_u, local_v, normal) is right-handed, so CCW vertex
// order in 3-D maps to CCW order in UV space.
// ---------------------------------------------------------------------------
struct Vec2 {
    float u{}, v{};
    constexpr Vec2() = default;
    constexpr Vec2(float u, float v) : u(u), v(v) {}
    Vec2  operator+(const Vec2& b) const { return {u+b.u, v+b.v}; }
    Vec2  operator-(const Vec2& b) const { return {u-b.u, v-b.v}; }
    Vec2  operator*(float s)       const { return {u*s, v*s}; }
    float dot(const Vec2& b)       const { return u*b.u + v*b.v; }
    float cross(const Vec2& b)     const { return u*b.v - v*b.u; } // 2-D "cross" (scalar)
    float length()                 const { return std::sqrt(u*u + v*v); }
};

// ---------------------------------------------------------------------------
// Face  — one logical surface (quad) before meshing
//
// A Face is a pre-mesh abstraction populated by buildCornellBox().  It
// stores the geometry and material of a single planar quad together with
// an orthonormal tangent frame used by the discontinuity mesher to project
// 3-D critical lines into face-local UV coordinates.
//
// Two-phase build model:
//   Phase 1 — buildCornellBox() fills Scene::faces[] with Face objects.
//   Phase 2 — either meshUniform() or DiscMesh::apply() consumes faces[]
//             and fills Scene::patches[] with Patch objects.
//
// Vertex ordering: verts[0..3] are CCW when viewed from the outward normal
// (the side a ray from inside the box would hit).
// ---------------------------------------------------------------------------
struct Face {
    Vec3        verts[4];       // corners, CCW when viewed from the front normal
    Vec3        normal;
    Vec3        emission;       // zero for non-emissive surfaces
    Vec3        reflectance;
    std::string name;
    bool        is_box{false};  // true for the two Cornell box objects

    // Orthonormal tangent frame in the face plane (computed on construction
    // by addFace()).  local_u × local_v = normal (right-handed).
    Vec3 local_u;   // unit vector along verts[0]→verts[1]
    Vec3 local_v;   // = normal × local_u (unit vector, in-plane, ⊥ local_u)

    // Project a 3-D point (assumed coplanar) to face-local UV.
    Vec2 project(const Vec3& p) const {
        Vec3 d = p - verts[0];
        return {d.dot(local_u), d.dot(local_v)};
    }
    // Unproject face-local UV back to 3-D.
    Vec3 unproject(const Vec2& uv) const {
        return verts[0] + local_u * uv.u + local_v * uv.v;
    }
};

// ---------------------------------------------------------------------------
// Patch  — one meshed element (product of meshing a Face)
// ---------------------------------------------------------------------------
struct Patch {
    Vec3  verts[4];         // corners, CCW order
    Vec3  normal, center;
    float area{};
    Vec3  emission;         // emitted radiance (light sources)
    Vec3  reflectance;      // diffuse reflectance ρ
    Vec3  radiosity;        // accumulated B
    Vec3  unshot;           // unshot B (progressive refinement)
    int   face_id{-1};      // index into Scene::faces (-1 = not set)
};

// ---------------------------------------------------------------------------
// Scene
//
// Holds the two representations of the geometry:
//   faces[]   — logical surfaces, populated once by buildCornellBox()
//   patches[] — meshed elements, populated by meshUniform() or DiscMesh::apply()
//
// box_start is the index of the first patch belonging to the Cornell box
// objects (the two blocks inside the room).  The renderer uses this to
// separate room surfaces from box surfaces for rendering purposes.
// ---------------------------------------------------------------------------
class Scene {
public:
    std::vector<Face>  faces;       // logical surfaces (populated by buildCornellBox)
    std::vector<Patch> patches;     // meshed elements (populated by meshUniform / DiscMesh)
    int box_start{0};               // index of first box patch (set after meshing)

    // Populate faces[] with Cornell Box geometry (does NOT create patches).
    void buildCornellBox();

    // Uniform-grid fallback: subdivide each face into W×W (walls) or B×B (boxes) patches.
    void meshUniform(int W, int B);

    // Low-level patch emitter used by meshing passes (public so DiscMesh can call it).
    void addPatch(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 v3,
                  Vec3 emission, Vec3 reflectance, int face_id = -1);

private:
    void addFace(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 v3,
                 Vec3 emission, Vec3 reflectance, const char* name, bool is_box = false);
    void addSubdividedQuad(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 v3,
                           Vec3 emission, Vec3 reflectance, int divs, int face_id);
};
