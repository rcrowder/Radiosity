#pragma once

#include <cmath>
#include <vector>

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

struct Patch {
    Vec3  verts[4];         // corners, CCW order
    Vec3  normal, center;
    float area{};
    Vec3  emission;         // emitted radiance (light sources)
    Vec3  reflectance;      // diffuse reflectance ρ
    Vec3  radiosity;        // accumulated B
    Vec3  unshot;           // unshot B (progressive refinement)
};

class Scene {
public:
    std::vector<Patch> patches;
    int box_start{0};   // index of first box patch (set by buildCornellBox)

    void buildCornellBox();

private:
    void addPatch(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 v3,
                  Vec3 emission, Vec3 reflectance);
    // Subdivide a quad into divs×divs sub-patches before adding.
    void addSubdividedQuad(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 v3,
                           Vec3 emission, Vec3 reflectance, int divs);
    // Add a rotated box (5 visible faces), each face subdivided.
    void addBox(Vec3 b0, Vec3 b1, Vec3 b2, Vec3 b3, float height,
                Vec3 emission, Vec3 reflectance, int divs);
};
