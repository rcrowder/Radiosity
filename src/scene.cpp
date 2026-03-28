#include "scene.h"

// ---------------------------------------------------------------------------
// addPatch — emit a single quad patch (public so DiscMesh can call it)
// ---------------------------------------------------------------------------
void Scene::addPatch(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 v3,
                     Vec3 emission, Vec3 reflectance, int face_id)
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
    p.face_id     = face_id;
    patches.push_back(p);
}

// ---------------------------------------------------------------------------
// addFace — register a logical surface into faces[]
// ---------------------------------------------------------------------------
void Scene::addFace(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 v3,
                    Vec3 emission, Vec3 reflectance,
                    const char* name, bool is_box)
{
    Face f;
    f.verts[0] = v0; f.verts[1] = v1;
    f.verts[2] = v2; f.verts[3] = v3;
    f.normal     = (v1 - v0).cross(v3 - v0).normalized();
    f.emission   = emission;
    f.reflectance = reflectance;
    f.name       = name;
    f.is_box     = is_box;
    // Build orthonormal tangent frame
    f.local_u = (v1 - v0).normalized();
    f.local_v = f.normal.cross(f.local_u).normalized();
    faces.push_back(f);
}

// ---------------------------------------------------------------------------
// addSubdividedQuad (private) — bilinear subdivision into divs×divs patches
// ---------------------------------------------------------------------------
void Scene::addSubdividedQuad(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 v3,
                               Vec3 emission, Vec3 reflectance,
                               int divs, int face_id)
{
    float inv = 1.f / divs;
    for (int row = 0; row < divs; ++row) {
        for (int col = 0; col < divs; ++col) {
            auto blerp = [&](float s, float t) {
                Vec3 bot = v0 + (v1 - v0) * s;
                Vec3 top = v3 + (v2 - v3) * s;
                return bot + (top - bot) * t;
            };
            float s0 = col * inv,   s1 = (col + 1) * inv;
            float t0 = row * inv,   t1 = (row + 1) * inv;
            addPatch(blerp(s0,t0), blerp(s1,t0),
                     blerp(s1,t1), blerp(s0,t1),
                     emission, reflectance, face_id);
        }
    }
}

// ---------------------------------------------------------------------------
// buildCornellBox — populate faces[] with the official Cornell dataset geometry.
// Does NOT create any patches; call meshUniform() or DiscMesh::apply() after.
// ---------------------------------------------------------------------------
void Scene::buildCornellBox()
{
    faces.clear();
    patches.clear();
    box_start = 0;

    constexpr Vec3 black  {0.0f,  0.0f,   0.0f};
    constexpr Vec3 white  {0.78f, 0.78f,  0.78f};
    constexpr Vec3 red    {0.63f, 0.065f, 0.05f};
    constexpr Vec3 green  {0.14f, 0.45f,  0.091f};
    constexpr Vec3 light_e{12.0f, 12.0f,  12.0f};

    constexpr float S = 1.f / 555.f;

    // ---- Environment ----
    addFace({552.8f*S,       0,          0},
            {  0,            0,          0},
            {  0,            0,  559.2f*S},
            {549.6f*S,       0,  559.2f*S}, black, white, "Floor");

    // Ceiling split into 4 strips framing the light opening (light is coplanar at y=548.8*S)
    // Front strip: full width, z=0..227S
    addFace({  0,      548.8f*S,          0},
            {556.0f*S, 548.8f*S,          0},
            {556.0f*S, 548.8f*S,  227.0f*S},
            {  0,      548.8f*S,  227.0f*S}, black, white, "Ceiling");
    // Back strip: full width, z=332S..559.2S
    addFace({  0,      548.8f*S,  332.0f*S},
            {556.0f*S, 548.8f*S,  332.0f*S},
            {556.0f*S, 548.8f*S,  559.2f*S},
            {  0,      548.8f*S,  559.2f*S}, black, white, "Ceiling");
    // Left strip: x=0..213S, z=227S..332S
    addFace({  0,      548.8f*S,  227.0f*S},
            {213.0f*S, 548.8f*S,  227.0f*S},
            {213.0f*S, 548.8f*S,  332.0f*S},
            {  0,      548.8f*S,  332.0f*S}, black, white, "Ceiling");
    // Right strip: x=343S..556S, z=227S..332S
    addFace({343.0f*S, 548.8f*S,  227.0f*S},
            {556.0f*S, 548.8f*S,  227.0f*S},
            {556.0f*S, 548.8f*S,  332.0f*S},
            {343.0f*S, 548.8f*S,  332.0f*S}, black, white, "Ceiling");

    addFace({549.6f*S,         0, 559.2f*S},
            {  0,              0, 559.2f*S},
            {  0,      548.8f*S,  559.2f*S},
            {556.0f*S, 548.8f*S,  559.2f*S}, black, white, "Back wall");

    addFace({0,       0,         559.2f*S},
            {0,       0,                0},
            {0, 548.8f*S,               0},
            {0, 548.8f*S,       559.2f*S}, black, green, "Right(green)");

    addFace({552.8f*S,         0,          0},
            {549.6f*S,         0,  559.2f*S},
            {556.0f*S, 548.8f*S,   559.2f*S},
            {556.0f*S, 548.8f*S,          0}, black, red, "Left(red)");

    // Light — flush with ceiling (coplanar), no offset
    addFace({343.0f*S, 548.8f*S, 227.0f*S},
            {343.0f*S, 548.8f*S, 332.0f*S},
            {213.0f*S, 548.8f*S, 332.0f*S},
            {213.0f*S, 548.8f*S, 227.0f*S}, light_e, black, "Light");

    // ---- Short block ----
    addFace({130.0f*S, 165.0f*S,  65.0f*S},
            { 82.0f*S, 165.0f*S, 225.0f*S},
            {240.0f*S, 165.0f*S, 272.0f*S},
            {290.0f*S, 165.0f*S, 114.0f*S}, black, white, "Short top", true);
    addFace({290.0f*S,        0, 114.0f*S},
            {290.0f*S, 165.0f*S, 114.0f*S},
            {240.0f*S, 165.0f*S, 272.0f*S},
            {240.0f*S,        0, 272.0f*S}, black, white, "Short front", true);
    addFace({130.0f*S,        0,  65.0f*S},
            {130.0f*S, 165.0f*S,  65.0f*S},
            {290.0f*S, 165.0f*S, 114.0f*S},
            {290.0f*S,        0, 114.0f*S}, black, white, "Short right", true);
    addFace({ 82.0f*S,        0, 225.0f*S},
            { 82.0f*S, 165.0f*S, 225.0f*S},
            {130.0f*S, 165.0f*S,  65.0f*S},
            {130.0f*S,        0,  65.0f*S}, black, white, "Short back", true);
    addFace({240.0f*S,        0, 272.0f*S},
            {240.0f*S, 165.0f*S, 272.0f*S},
            { 82.0f*S, 165.0f*S, 225.0f*S},
            { 82.0f*S,        0, 225.0f*S}, black, white, "Short left", true);

    // ---- Tall block ----
    addFace({423.0f*S, 330.0f*S, 247.0f*S},
            {265.0f*S, 330.0f*S, 296.0f*S},
            {314.0f*S, 330.0f*S, 456.0f*S},
            {472.0f*S, 330.0f*S, 406.0f*S}, black, white, "Tall top", true);
    addFace({423.0f*S,        0, 247.0f*S},
            {423.0f*S, 330.0f*S, 247.0f*S},
            {472.0f*S, 330.0f*S, 406.0f*S},
            {472.0f*S,        0, 406.0f*S}, black, white, "Tall front", true);
    addFace({472.0f*S,        0, 406.0f*S},
            {472.0f*S, 330.0f*S, 406.0f*S},
            {314.0f*S, 330.0f*S, 456.0f*S},
            {314.0f*S,        0, 456.0f*S}, black, white, "Tall right", true);
    addFace({314.0f*S,        0, 456.0f*S},
            {314.0f*S, 330.0f*S, 456.0f*S},
            {265.0f*S, 330.0f*S, 296.0f*S},
            {265.0f*S,        0, 296.0f*S}, black, white, "Tall back", true);
    addFace({265.0f*S,        0, 296.0f*S},
            {265.0f*S, 330.0f*S, 296.0f*S},
            {423.0f*S, 330.0f*S, 247.0f*S},
            {423.0f*S,        0, 247.0f*S}, black, white, "Tall left", true);
}

// ---------------------------------------------------------------------------
// meshUniform — subdivide every face in the faces[] list into a regular grid.
// ---------------------------------------------------------------------------
void Scene::meshUniform(int W, int B)
{
    patches.clear();
    box_start = -1;

    for (int fi = 0; fi < static_cast<int>(faces.size()); ++fi) {
        const Face& f = faces[fi];

        if (f.is_box && box_start < 0)
            box_start = static_cast<int>(patches.size());

        // Light and emissive faces → single patch, no subdivision
        bool emissive = f.emission.x > 0.f || f.emission.y > 0.f || f.emission.z > 0.f;
        if (emissive) {
            addPatch(f.verts[0], f.verts[1], f.verts[2], f.verts[3],
                     f.emission, f.reflectance, fi);
            continue;
        }

        int divs = f.is_box ? B : W;
        addSubdividedQuad(f.verts[0], f.verts[1], f.verts[2], f.verts[3],
                          f.emission, f.reflectance, divs, fi);
    }

    if (box_start < 0)
        box_start = static_cast<int>(patches.size());
}
