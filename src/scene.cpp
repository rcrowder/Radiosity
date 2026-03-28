#include "scene.h"

void Scene::addPatch(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 v3,
                     Vec3 emission, Vec3 reflectance)
{
    Patch p;
    p.verts[0] = v0;  p.verts[1] = v1;
    p.verts[2] = v2;  p.verts[3] = v3;
    p.center   = (v0 + v1 + v2 + v3) * 0.25f;
    p.normal   = (v1 - v0).cross(v3 - v0).normalized();
    p.area     = 0.5f * ((v2-v0).cross(v1-v0).length()
                       + (v2-v0).cross(v3-v0).length());
    p.emission    = emission;
    p.reflectance = reflectance;
    p.radiosity   = emission;   // initial B = E
    p.unshot      = emission;
    patches.push_back(p);
}

/*
 * Cornell Box (simplified, one patch per wall face).
 * Coordinates: x,y,z in [0,1]; y is up, z=0 is the open front.
 */
/*
 * Bilinearly subdivide a quad into divs×divs equal sub-patches.
 * v0=bottom-left, v1=bottom-right, v2=top-right, v3=top-left (CCW).
 */
void Scene::addSubdividedQuad(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 v3,
                               Vec3 emission, Vec3 reflectance, int divs)
{
    // Bilinear interpolation across the quad so non-rectangular quads
    // (e.g. rotated box faces) subdivide correctly.
    float inv = 1.f / divs;
    for (int row = 0; row < divs; ++row) {
        for (int col = 0; col < divs; ++col) {
            // Four corners via bilinear lerp: s in [0,1] along v0->v1/v3->v2,
            //                                 t in [0,1] along v0->v3/v1->v2
            auto blerp = [&](float s, float t) {
                Vec3 bot = v0 + (v1 - v0) * s;
                Vec3 top = v3 + (v2 - v3) * s;
                return bot + (top - bot) * t;
            };
            float s0 = col * inv,       s1 = (col + 1) * inv;
            float t0 = row * inv,       t1 = (row + 1) * inv;
            addPatch(blerp(s0,t0), blerp(s1,t0),
                     blerp(s1,t1), blerp(s0,t1),
                     emission, reflectance);
        }
    }
}

void Scene::addBox(Vec3 b0, Vec3 b1, Vec3 b2, Vec3 b3, float height,
                   Vec3 emission, Vec3 reflectance, int divs)
{
    Vec3 t0{b0.x, b0.y + height, b0.z};
    Vec3 t1{b1.x, b1.y + height, b1.z};
    Vec3 t2{b2.x, b2.y + height, b2.z};
    Vec3 t3{b3.x, b3.y + height, b3.z};

    // Top face
    addSubdividedQuad(t0, t1, t2, t3, emission, reflectance, divs);
    // Four side faces (outward winding)
    addSubdividedQuad(b1, b0, t0, t1, emission, reflectance, divs);
    addSubdividedQuad(b2, b1, t1, t2, emission, reflectance, divs);
    addSubdividedQuad(b3, b2, t2, t3, emission, reflectance, divs);
    addSubdividedQuad(b0, b3, t3, t0, emission, reflectance, divs);
}

void Scene::buildCornellBox()
{
    patches.clear();

    constexpr Vec3 black  {0.0f,  0.0f,   0.0f};
    constexpr Vec3 white  {0.78f, 0.78f,  0.78f};
    constexpr Vec3 red    {0.63f, 0.065f, 0.05f};
    constexpr Vec3 green  {0.14f, 0.45f,  0.091f};
    constexpr Vec3 light_e{12.0f, 12.0f,  12.0f};

    // Official Cornell Box dataset — graphics.cornell.edu/online/box/data.html
    // All coordinates are the measured values divided by 555 (room ≈ 555 units wide).
    // Vertex winding taken directly from the dataset; normals point into the room.
    constexpr int   W = 10;          // wall/floor/ceiling subdivisions
    constexpr int   B = 4;           // box face subdivisions
    constexpr float S = 1.f / 555.f; // normalisation factor

    // Floor — white
    addSubdividedQuad({552.8f*S,       0,          0},
                      {  0,            0,          0},
                      {  0,            0,  559.2f*S},
                      {549.6f*S,       0,  559.2f*S}, black, white, W);

    // Ceiling — white
    addSubdividedQuad({556.0f*S, 548.8f*S,          0},
                      {556.0f*S, 548.8f*S,  559.2f*S},
                      {  0,      548.8f*S,  559.2f*S},
                      {  0,      548.8f*S,          0}, black, white, W);

    // Back wall — white
    addSubdividedQuad({549.6f*S,         0, 559.2f*S},
                      {  0,              0, 559.2f*S},
                      {  0,      548.8f*S,  559.2f*S},
                      {556.0f*S, 548.8f*S,  559.2f*S}, black, white, W);

    // Right wall — green  (x = 0 plane)
    addSubdividedQuad({0,       0,         559.2f*S},
                      {0,       0,                0},
                      {0, 548.8f*S,               0},
                      {0, 548.8f*S,       559.2f*S}, black, green, W);

    // Left wall — red  (x ≈ 555/555 plane)
    addSubdividedQuad({552.8f*S,         0,          0},
                      {549.6f*S,         0,  559.2f*S},
                      {556.0f*S, 548.8f*S,   559.2f*S},
                      {556.0f*S, 548.8f*S,          0}, black, red,   W);

    // Light — emissive; offset 1 mm below ceiling to avoid z-fighting
    addPatch({343.0f*S, 548.8f*S - 0.001f, 227.0f*S},
             {343.0f*S, 548.8f*S - 0.001f, 332.0f*S},
             {213.0f*S, 548.8f*S - 0.001f, 332.0f*S},
             {213.0f*S, 548.8f*S - 0.001f, 227.0f*S}, light_e, black);

    // -------------------------------------------------------------------
    box_start = static_cast<int>(patches.size());

    // Short block  (h = 165/555 ≈ 0.297, near green/right wall at x ≈ 0)
    addSubdividedQuad({130.0f*S, 165.0f*S,  65.0f*S},   // top
                      { 82.0f*S, 165.0f*S, 225.0f*S},
                      {240.0f*S, 165.0f*S, 272.0f*S},
                      {290.0f*S, 165.0f*S, 114.0f*S}, black, white, B);
    addSubdividedQuad({290.0f*S,        0, 114.0f*S},
                      {290.0f*S, 165.0f*S, 114.0f*S},
                      {240.0f*S, 165.0f*S, 272.0f*S},
                      {240.0f*S,        0, 272.0f*S}, black, white, B);
    addSubdividedQuad({130.0f*S,        0,  65.0f*S},
                      {130.0f*S, 165.0f*S,  65.0f*S},
                      {290.0f*S, 165.0f*S, 114.0f*S},
                      {290.0f*S,        0, 114.0f*S}, black, white, B);
    addSubdividedQuad({ 82.0f*S,        0, 225.0f*S},
                      { 82.0f*S, 165.0f*S, 225.0f*S},
                      {130.0f*S, 165.0f*S,  65.0f*S},
                      {130.0f*S,        0,  65.0f*S}, black, white, B);
    addSubdividedQuad({240.0f*S,        0, 272.0f*S},
                      {240.0f*S, 165.0f*S, 272.0f*S},
                      { 82.0f*S, 165.0f*S, 225.0f*S},
                      { 82.0f*S,        0, 225.0f*S}, black, white, B);

    // Tall block  (h = 330/555 ≈ 0.595, near red/left wall at x ≈ 1)
    addSubdividedQuad({423.0f*S, 330.0f*S, 247.0f*S},   // top
                      {265.0f*S, 330.0f*S, 296.0f*S},
                      {314.0f*S, 330.0f*S, 456.0f*S},
                      {472.0f*S, 330.0f*S, 406.0f*S}, black, white, B);
    addSubdividedQuad({423.0f*S,        0, 247.0f*S},
                      {423.0f*S, 330.0f*S, 247.0f*S},
                      {472.0f*S, 330.0f*S, 406.0f*S},
                      {472.0f*S,        0, 406.0f*S}, black, white, B);
    addSubdividedQuad({472.0f*S,        0, 406.0f*S},
                      {472.0f*S, 330.0f*S, 406.0f*S},
                      {314.0f*S, 330.0f*S, 456.0f*S},
                      {314.0f*S,        0, 456.0f*S}, black, white, B);
    addSubdividedQuad({314.0f*S,        0, 456.0f*S},
                      {314.0f*S, 330.0f*S, 456.0f*S},
                      {265.0f*S, 330.0f*S, 296.0f*S},
                      {265.0f*S,        0, 296.0f*S}, black, white, B);
    addSubdividedQuad({265.0f*S,        0, 296.0f*S},
                      {265.0f*S, 330.0f*S, 296.0f*S},
                      {423.0f*S, 330.0f*S, 247.0f*S},
                      {423.0f*S,        0, 247.0f*S}, black, white, B);
}
