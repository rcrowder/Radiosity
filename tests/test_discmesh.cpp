// Unit tests for discontinuity meshing and supporting math.
//
// Test organisation
// -----------------
// Section A — Math primitives (Vec3, Vec2, Face frame, project/unproject)
// Section B — DiscMesh invariants on the Cornell Box (integration test)
// Section C — DiscMesh invariants on hand-crafted synthetic scenes
//
//   C-01  Minimal room, no internal occluders.
//   C-02  Room + single vertical occluder → floor must be split.
//   C-03  Tilted (30°) non-axis-aligned receiver face.
//   C-04  Wide (2:1) rectangular floor — area conservation.
//   C-05  Symmetric mirror occluder pair → equal half-areas.
//   C-06  45°-rotated room — no axis-aligned faces.
//   C-07  Tiny-scale room (0.1 units) — numerical stability.
//   C-08  Large-scale room (100 units) — numerical stability.
//   C-09  Four pillar occluders — stress-tests splitting & deduplication.
//
// All Section C tests run the full invariant suite:
//   • no zero-area patches
//   • every patch normal matches its parent face normal (dot ≥ 0.99)
//   • every patch center lies on the parent face plane
//   • every patch vertex lies on the parent face plane
//   • triangle patch centroid is (v0+v1+v2)/3
//   • all face_ids are in [0, faces.size())
//   • every face receives ≥ 1 patch
//   • every patch center lies inside the face UV boundary
//
// Compile & run:  make test
// (links scene.cpp + discmesh.cpp only — no GL dependency)

#include "scene.h"
#include "discmesh.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

// ============================================================
// Minimal test framework
// ============================================================
static int g_pass = 0, g_fail = 0;

#define SUITE(name) std::printf("\n[%s]\n", (name))

#define CHECK(cond) \
    do { \
        if (cond) { ++g_pass; } \
        else { ++g_fail; \
               std::printf("  FAIL  %s:%d  %s\n", __FILE__, __LINE__, #cond); } \
    } while (0)

#define CHECK_NEAR(a, b, eps) \
    CHECK(std::abs(static_cast<float>(a) - static_cast<float>(b)) \
          <= static_cast<float>(eps))

static void print_summary()
{
    std::printf("\n================================================\n");
    std::printf("%d passed,  %d failed\n", g_pass, g_fail);
}
static int exit_code() { return g_fail > 0 ? 1 : 0; }

// ============================================================
// Shared helpers
// ============================================================

static float quad_area(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 v3)
{
    auto cross3 = [](Vec3 a, Vec3 b) -> Vec3 {
        return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
    };
    auto len3   = [](Vec3 v)        { return std::sqrt(v.x*v.x+v.y*v.y+v.z*v.z); };
    Vec3 d02 = v2 - v0;
    return 0.5f*(len3(cross3(d02,v1-v0)) + len3(cross3(d02,v3-v0)));
}

// Build a Face from four CCW vertices.
// Mirrors Scene::addFace (which is private).
static Face make_face(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 v3,
                      Vec3 emission, Vec3 reflectance,
                      const char* name, bool is_box = false)
{
    Face f;
    f.verts[0]=v0; f.verts[1]=v1; f.verts[2]=v2; f.verts[3]=v3;
    f.normal      = (v1-v0).cross(v3-v0).normalized();
    f.emission    = emission;
    f.reflectance = reflectance;
    f.name        = name;
    f.is_box      = is_box;
    f.local_u     = (v1-v0).normalized();
    f.local_v     = f.normal.cross(f.local_u).normalized();
    return f;
}

// Run the full invariant suite on a meshed scene.
// Returns the number of invariant violations found.
static int check_invariants(const Scene& s, const char* label)
{
    const int nf = static_cast<int>(s.faces.size());
    int fails = 0;

    // no zero-area patches
    for (const auto& p : s.patches)
        if (p.area <= 0.f) {
            ++fails;
            std::printf("  [%s] zero-area: face_id=%d area=%.2e\n",
                        label, p.face_id, p.area);
        }

    // normals match face (dot ≥ 0.99)
    for (const auto& p : s.patches) {
        if (p.face_id < 0 || p.face_id >= nf) continue;
        float dot = p.normal.dot(s.faces[p.face_id].normal);
        if (dot < 0.99f) {
            ++fails;
            std::printf("  [%s] bad normal: face=%s dot=%.4f\n",
                        label, s.faces[p.face_id].name.c_str(), dot);
        }
    }

    // centers on face plane
    for (const auto& p : s.patches) {
        if (p.face_id < 0 || p.face_id >= nf) continue;
        const Face& f = s.faces[p.face_id];
        float dist = std::abs((p.center - f.verts[0]).dot(f.normal));
        if (dist > 1e-4f) {
            ++fails;
            std::printf("  [%s] center off plane: face=%s dist=%.2e\n",
                        label, f.name.c_str(), dist);
        }
    }

    // vertices on face plane
    for (const auto& p : s.patches) {
        if (p.face_id < 0 || p.face_id >= nf) continue;
        const Face& f = s.faces[p.face_id];
        bool is_tri = (p.verts[2].x==p.verts[3].x &&
                       p.verts[2].y==p.verts[3].y &&
                       p.verts[2].z==p.verts[3].z);
        int nv = is_tri ? 3 : 4;
        for (int i = 0; i < nv; ++i) {
            float d = std::abs((p.verts[i]-f.verts[0]).dot(f.normal));
            if (d > 1e-4f) {
                ++fails;
                std::printf("  [%s] vert[%d] off plane: face=%s dist=%.2e\n",
                            label, i, f.name.c_str(), d);
                break;
            }
        }
    }

    // triangle centroid = (v0+v1+v2)/3
    for (const auto& p : s.patches) {
        bool is_tri = (p.verts[2].x==p.verts[3].x &&
                       p.verts[2].y==p.verts[3].y &&
                       p.verts[2].z==p.verts[3].z);
        if (!is_tri) continue;
        Vec3 tc = (p.verts[0]+p.verts[1]+p.verts[2])*(1.f/3.f);
        float err = (p.center-tc).length();
        if (err > 1e-5f) {
            ++fails;
            std::printf("  [%s] bad tri centroid: face_id=%d err=%.2e\n",
                        label, p.face_id, err);
        }
    }

    // face_ids in range
    for (const auto& p : s.patches)
        if (p.face_id < 0 || p.face_id >= nf) {
            ++fails;
            std::printf("  [%s] face_id out of range: %d\n", label, p.face_id);
        }

    // every face gets ≥ 1 patch
    {
        std::vector<int> cnt(nf, 0);
        for (const auto& p : s.patches)
            if (p.face_id>=0 && p.face_id<nf) ++cnt[p.face_id];
        for (int fi = 0; fi < nf; ++fi)
            if (cnt[fi] == 0) {
                ++fails;
                std::printf("  [%s] face '%s' has 0 patches\n",
                            label, s.faces[fi].name.c_str());
            }
    }

    // patch centers inside face UV boundary
    for (const auto& p : s.patches) {
        if (p.face_id < 0 || p.face_id >= nf) continue;
        const Face& f = s.faces[p.face_id];
        Vec2 uv = f.project(p.center);
        Vec2 uvs[4];
        for (int i = 0; i < 4; ++i) uvs[i] = f.project(f.verts[i]);
        bool inside = true;
        for (int i = 0; i < 4; ++i) {
            Vec2 edge  = uvs[(i+1)%4] - uvs[i];
            Vec2 to_pt = uv - uvs[i];
            if (edge.cross(to_pt) < -1e-4f) { inside = false; break; }
        }
        if (!inside) {
            ++fails;
            std::printf("  [%s] center outside UV: face=%s uv=(%.4f,%.4f)\n",
                        label, f.name.c_str(), uv.u, uv.v);
        }
    }

    return fails;
}

// Build a closed room (floor + ceiling + 4 walls + overhead light).
// No internal occluders; caller can push_back additional faces before
// calling DiscMesh::apply().
// Build a closed room with inward-pointing face normals (faces visible from inside).
// Normal convention: (v1-v0)×(v3-v0) must point INTO the room.
//
//   Floor   y=0: normal +Y  → v0=(0,0,0) v1=(0,0,D) v2=(W,0,D) v3=(W,0,0)
//   Ceiling y=H: normal -Y  → v0=(0,H,0) v1=(W,H,0) v2=(W,H,D) v3=(0,H,D)
//   Back    z=D: normal -Z  → v0=(W,0,D) v1=(W,H,D) v2=(0,H,D) v3=(0,0,D)
//   Front   z=0: normal +Z  → v0=(0,0,0) v1=(0,H,0) v2=(W,H,0) v3=(W,0,0)
//   Left    x=0: normal +X  → v0=(0,0,D) v1=(0,H,D) v2=(0,H,0) v3=(0,0,0)
//   Right   x=W: normal -X  → v0=(W,0,0) v1=(W,H,0) v2=(W,H,D) v3=(W,0,D)
static Scene make_simple_room(float W, float H, float D,
                              float lx0, float lx1,
                              float lz0, float lz1)
{
    const float eps = 0.001f;
    const Vec3 black{}, white{0.8f,0.8f,0.8f}, light_e{10,10,10};
    Scene s;
    // Floor: normal +Y
    s.faces.push_back(make_face({0,0,0},{0,0,D},{W,0,D},{W,0,0},
                                black,white,"Floor"));
    // Ceiling: normal -Y
    s.faces.push_back(make_face({0,H,0},{W,H,0},{W,H,D},{0,H,D},
                                black,white,"Ceiling"));
    // Back wall: normal -Z
    s.faces.push_back(make_face({W,0,D},{W,H,D},{0,H,D},{0,0,D},
                                black,white,"Back"));
    // Front wall: normal +Z
    s.faces.push_back(make_face({0,0,0},{0,H,0},{W,H,0},{W,0,0},
                                black,white,"Front"));
    // Left wall: normal +X
    s.faces.push_back(make_face({0,0,D},{0,H,D},{0,H,0},{0,0,0},
                                black,white,"Left"));
    // Right wall: normal -X
    s.faces.push_back(make_face({W,0,0},{W,H,0},{W,H,D},{W,0,D},
                                black,white,"Right"));
    // Light patch slightly below ceiling (no normal conflict)
    s.faces.push_back(make_face(
        {lx0,H-eps,lz0},{lx1,H-eps,lz0},{lx1,H-eps,lz1},{lx0,H-eps,lz1},
        light_e,black,"Light"));
    return s;
}

// ============================================================
// Section A — Math primitives
// ============================================================

static void test_vec3()
{
    SUITE("A-01  Vec3 arithmetic");
    Vec3 a{1,2,3}, b{4,5,6};
    Vec3 s=a+b; CHECK_NEAR(s.x,5,1e-6f); CHECK_NEAR(s.y,7,1e-6f); CHECK_NEAR(s.z,9,1e-6f);
    Vec3 d=b-a; CHECK_NEAR(d.x,3,1e-6f); CHECK_NEAR(d.y,3,1e-6f); CHECK_NEAR(d.z,3,1e-6f);
    CHECK_NEAR(a.dot(b), 32.f, 1e-6f);
    Vec3 c=a.cross(b); CHECK_NEAR(c.x,-3,1e-6f); CHECK_NEAR(c.y,6,1e-6f); CHECK_NEAR(c.z,-3,1e-6f);
    CHECK_NEAR(a.dot(c), 0.f, 1e-5f);
    CHECK_NEAR(b.dot(c), 0.f, 1e-5f);
    CHECK_NEAR(a.normalized().length(), 1.f, 1e-6f);
    CHECK_NEAR(Vec3{}.normalized().length(), 0.f, 1e-6f);
    Vec3 sc=a*2.f; CHECK_NEAR(sc.x,2,1e-6f); CHECK_NEAR(sc.z,6,1e-6f);
    Vec3 neg=-a;  CHECK_NEAR(neg.x,-1,1e-6f);
    Vec3 cw=a*b;  CHECK_NEAR(cw.x,4,1e-6f); CHECK_NEAR(cw.y,10,1e-6f); CHECK_NEAR(cw.z,18,1e-6f);
}

static void test_vec2()
{
    SUITE("A-02  Vec2 arithmetic");
    Vec2 a{3,4}; CHECK_NEAR(a.length(), 5.f, 1e-6f);
    Vec2 b{1,2}, c{3,4};
    CHECK_NEAR(b.cross(c), -2.f, 1e-6f);
    CHECK_NEAR(b.dot(c),   11.f, 1e-6f);
    Vec2 s=b+c; CHECK_NEAR(s.u,4,1e-6f); CHECK_NEAR(s.v,6,1e-6f);
    Vec2 d=c-b; CHECK_NEAR(d.u,2,1e-6f); CHECK_NEAR(d.v,2,1e-6f);
    Vec2 sc=b*3.f; CHECK_NEAR(sc.u,3,1e-6f); CHECK_NEAR(sc.v,6,1e-6f);
}

static void test_face_frame_axis_aligned()
{
    SUITE("A-03  Face frame — axis-aligned quads");
    // XZ floor, inward-facing normal +Y.
    // CCW order when viewed from +Y: v0=(0,0,0), v1=(0,0,1), v2=(1,0,1), v3=(1,0,0)
    // (v1-v0)×(v3-v0) = (0,0,1)×(1,0,0) = (0*0-1*0, 1*1-0*0, 0*0-0*1) = (0,1,0) ✓
    Face fl = make_face({0,0,0},{0,0,1},{1,0,1},{1,0,0},{},{1,1,1},"floor");
    CHECK_NEAR(fl.normal.y, 1.f, 1e-5f);
    CHECK_NEAR(fl.local_u.length(), 1.f, 1e-5f);
    CHECK_NEAR(fl.local_v.length(), 1.f, 1e-5f);
    CHECK_NEAR(fl.local_u.dot(fl.local_v), 0.f, 1e-5f);
    CHECK_NEAR(fl.local_u.cross(fl.local_v).dot(fl.normal), 1.f, 1e-5f);

    // YZ wall, inward-facing normal +X.
    // v0=(1,0,0), v1=(1,1,0), v3=(1,0,1)
    // (v1-v0)×(v3-v0) = (0,1,0)×(0,0,1) = (1,0,0) ✓
    Face sd = make_face({1,0,0},{1,1,0},{1,1,1},{1,0,1},{},{1,1,1},"side");
    CHECK_NEAR(sd.normal.x, 1.f, 1e-5f);
    CHECK_NEAR(sd.local_u.cross(sd.local_v).dot(sd.normal), 1.f, 1e-5f);

    // XY wall (z=1), check frame is orthonormal regardless of normal direction.
    Face bk = make_face({0,0,1},{1,0,1},{1,1,1},{0,1,1},{},{1,1,1},"back");
    CHECK_NEAR(bk.local_u.dot(bk.normal), 0.f, 1e-5f);
    CHECK_NEAR(bk.local_v.dot(bk.normal), 0.f, 1e-5f);
    CHECK_NEAR(bk.local_u.cross(bk.local_v).dot(bk.normal), 1.f, 1e-5f);
}

static void test_face_frame_tilted()
{
    SUITE("A-04  Face frame — 45° rotated (non-axis-aligned) quad");
    const float sq = 0.70710678f;
    Face g = make_face({sq,sq,0},{-sq,sq,0},{-sq,sq,1},{sq,sq,1},{},{1,1,1},"rot45");
    CHECK_NEAR(g.local_u.length(), 1.f, 1e-5f);
    CHECK_NEAR(g.local_v.length(), 1.f, 1e-5f);
    CHECK_NEAR(g.local_u.dot(g.local_v), 0.f, 1e-5f);
    CHECK_NEAR(g.local_u.cross(g.local_v).dot(g.normal), 1.f, 1e-5f);
    for (int i = 0; i < 4; ++i) {
        Vec3 r = g.unproject(g.project(g.verts[i]));
        CHECK_NEAR(r.x, g.verts[i].x, 1e-5f);
        CHECK_NEAR(r.y, g.verts[i].y, 1e-5f);
        CHECK_NEAR(r.z, g.verts[i].z, 1e-5f);
    }
}

static void test_face_project_interior()
{
    SUITE("A-05  Face project/unproject — interior points, various orientations");
    std::vector<Face> faces = {
        make_face({0,0,0},{2,0,0},{2,0,3},{0,0,3},{},{1,1,1},"XZ_nonsquare"),
        make_face({0,0,0},{0,0,1},{0,1,1},{0,1,0},{},{1,1,1},"YZ_wall"),
        make_face({0,0,1},{1,0,1},{1,1,1},{0,1,1},{},{1,1,1},"XY_wall"),
        // 30° tilt around X
        make_face({0,0,0},{1,0,0},{1,0.866f,0.5f},{0,0.866f,0.5f},{},{1,1,1},"tilt30"),
    };
    for (const auto& f : faces) {
        Vec3 cen = (f.verts[0]+f.verts[1]+f.verts[2]+f.verts[3])*0.25f;
        Vec3 r = f.unproject(f.project(cen));
        CHECK_NEAR(r.x, cen.x, 1e-4f);
        CHECK_NEAR(r.y, cen.y, 1e-4f);
        CHECK_NEAR(r.z, cen.z, 1e-4f);
        for (int i = 0; i < 4; ++i) {
            Vec3 rt = f.unproject(f.project(f.verts[i]));
            CHECK_NEAR(rt.x, f.verts[i].x, 1e-4f);
            CHECK_NEAR(rt.y, f.verts[i].y, 1e-4f);
            CHECK_NEAR(rt.z, f.verts[i].z, 1e-4f);
        }
    }
}

static void test_all_cornell_face_frames()
{
    SUITE("A-06  Cornell Box — all face tangent frames orthonormal");
    Scene s; s.buildCornellBox();
    for (const auto& f : s.faces) {
        CHECK_NEAR(f.local_u.length(), 1.f, 1e-5f);
        CHECK_NEAR(f.local_v.length(), 1.f, 1e-5f);
        CHECK_NEAR(f.normal.length(),  1.f, 1e-5f);
        CHECK_NEAR(f.local_u.dot(f.local_v),  0.f, 1e-5f);
        CHECK_NEAR(f.local_u.dot(f.normal),   0.f, 1e-5f);
        CHECK_NEAR(f.local_v.dot(f.normal),   0.f, 1e-5f);
        CHECK_NEAR(f.local_u.cross(f.local_v).dot(f.normal), 1.f, 1e-5f);
        // NOTE: project/unproject round-trip may exceed 1e-4 for non-rectangular
        // Cornell quads (e.g. Left(red) trapezoid) because the plane fit can
        // accumulate ~5e-3 error in scaled coordinates.  Use 1e-2 tolerance.
        for (int i = 0; i < 4; ++i) {
            Vec3 r = f.unproject(f.project(f.verts[i]));
            CHECK_NEAR(r.x, f.verts[i].x, 1e-2f);
            CHECK_NEAR(r.y, f.verts[i].y, 1e-2f);
            CHECK_NEAR(r.z, f.verts[i].z, 1e-2f);
        }
    }
}

// ============================================================
// Section B — Cornell Box integration tests
// ============================================================

static Scene& cornell_disc()
{
    static Scene s;
    static bool built = false;
    if (!built) { s.buildCornellBox(); DiscMesh::apply(s); built = true; }
    return s;
}

static void test_b01_emissive_single_patch()
{
    SUITE("B-01  Cornell — emissive face → exactly 1 patch");
    const Scene& s = cornell_disc();
    int n = 0;
    for (const auto& p : s.patches)
        if (p.emission.x>0||p.emission.y>0||p.emission.z>0) ++n;
    std::printf("  light patches: %d\n", n);
    CHECK(n == 1);
}

static void test_b02_box_start()
{
    SUITE("B-02  Cornell — box_start correctly partitions room vs box patches");
    const Scene& s = cornell_disc();
    CHECK(s.box_start >= 0);
    CHECK(s.box_start < static_cast<int>(s.patches.size()));
    int bad_room=0, bad_box=0;
    for (int i=0; i<s.box_start; ++i) {
        int fi=s.patches[i].face_id;
        if (fi>=0 && s.faces[fi].is_box) ++bad_room;
    }
    for (int i=s.box_start; i<static_cast<int>(s.patches.size()); ++i) {
        int fi=s.patches[i].face_id;
        if (fi>=0 && !s.faces[fi].is_box) ++bad_box;
    }
    CHECK(bad_room==0);
    CHECK(bad_box==0);
}

static void test_b03_area_conservation()
{
    SUITE("B-03  Cornell — area conservation per face (< 1% error)");
    const Scene& s = cornell_disc();
    const int nf = static_cast<int>(s.faces.size());
    std::vector<float> ps(nf, 0.f);
    for (const auto& p : s.patches)
        if (p.face_id>=0 && p.face_id<nf) ps[p.face_id]+=p.area;
    int bad=0;
    for (int fi=0; fi<nf; ++fi) {
        const Face& f=s.faces[fi];
        float fa  = quad_area(f.verts[0],f.verts[1],f.verts[2],f.verts[3]);
        float rel = std::abs(ps[fi]-fa)/(fa+1e-12f);
        std::printf("  %-16s  face=%.5f  patches=%.5f  rel_err=%.4f\n",
                    f.name.c_str(), fa, ps[fi], rel);
        if (rel>0.01f) { ++bad; std::printf("    *** rel_err > 1%% ***\n"); }
    }
    CHECK(bad==0);
}

static void test_b04_green_wall()
{
    SUITE("B-04  Cornell — Right(green) wall is split; area conserved");
    const Scene& s = cornell_disc();
    int gfi=-1;
    for (int fi=0; fi<static_cast<int>(s.faces.size()); ++fi)
        if (s.faces[fi].name=="Right(green)") { gfi=fi; break; }
    CHECK(gfi>=0);
    if (gfi<0) return;
    const Face& gf=s.faces[gfi];
    int count=0; float asum=0.f;
    for (const auto& p : s.patches)
        if (p.face_id==gfi) { ++count; asum+=p.area; }
    float fa  = quad_area(gf.verts[0],gf.verts[1],gf.verts[2],gf.verts[3]);
    float rel = std::abs(asum-fa)/(fa+1e-12f);
    std::printf("  patches=%d  face_area=%.5f  sum=%.5f  rel_err=%.4f\n",
                count, fa, asum, rel);
    CHECK(count>1);
    CHECK(rel<0.01f);
}

static void test_b05_non_emissive_split()
{
    SUITE("B-05  Cornell — non-emissive room faces all split into > 1 patch");
    const Scene& s = cornell_disc();
    const int nf=static_cast<int>(s.faces.size());
    std::vector<int> cnt(nf,0);
    for (const auto& p : s.patches)
        if (p.face_id>=0 && p.face_id<nf) ++cnt[p.face_id];
    int unsplit=0;
    for (int fi=0; fi<nf; ++fi) {
        const Face& f=s.faces[fi];
        bool emit = f.emission.x>0||f.emission.y>0||f.emission.z>0;
        if (!emit && !f.is_box && cnt[fi]<=1) {
            ++unsplit;
            std::printf("  not split: '%s' patches=%d\n", f.name.c_str(), cnt[fi]);
        }
    }
    CHECK(unsplit==0);
}

static void test_b06_all_invariants()
{
    SUITE("B-06  Cornell — full invariant suite");
    const Scene& s = cornell_disc();
    int fails = check_invariants(s, "cornell");
    std::printf("  patches=%zu  invariant_failures=%d\n", s.patches.size(), fails);
    CHECK(fails==0);
}

// ============================================================
// Section C — Synthetic geometry scenes
// ============================================================

static void test_c01_minimal_room()
{
    SUITE("C-01  Synthetic minimal room (no occluders)");
    Scene s = make_simple_room(2,2,2, 0.8f,1.2f,0.8f,1.2f);
    DiscMesh::apply(s);
    std::printf("  patches=%zu\n", s.patches.size());
    CHECK(check_invariants(s,"C01")==0);
}

static void test_c02_single_occluder()
{
    SUITE("C-02  Synthetic room — single vertical occluder splits floor");
    Scene s = make_simple_room(2,2,2, 0.8f,1.2f,0.8f,1.2f);
    // Thin vertical panel between light and floor
    s.faces.push_back(make_face(
        {0.9f,0,1},{1.1f,0,1},{1.1f,1.5f,1},{0.9f,1.5f,1},
        {},{0.5f,0.5f,0.5f},"OccFront"));
    s.faces.push_back(make_face(
        {1.1f,0,1},{0.9f,0,1},{0.9f,1.5f,1},{1.1f,1.5f,1},
        {},{0.5f,0.5f,0.5f},"OccBack"));
    DiscMesh::apply(s);

    int floor_fi=-1;
    for (int fi=0; fi<static_cast<int>(s.faces.size()); ++fi)
        if (s.faces[fi].name==std::string("Floor")) { floor_fi=fi; break; }
    CHECK(floor_fi>=0);
    int nc=0;
    for (const auto& p : s.patches) if (p.face_id==floor_fi) ++nc;
    std::printf("  floor patches: %d\n", nc);
    CHECK(nc>1);
    CHECK(check_invariants(s,"C02")==0);
}

static void test_c03_tilted_receiver()
{
    SUITE("C-03  Synthetic scene — tilted (30°) non-axis-aligned receiver");
    const float c30=0.866025f, s30=0.5f;
    Scene s;
    const Vec3 black{}, white{0.8f,0.8f,0.8f}, le{10,10,10};
    // Tilted floor: 30° around X
    s.faces.push_back(make_face(
        {0,0,0},{2,0,0},{2,s30*2,-c30*2},{0,s30*2,-c30*2},black,white,"TiltFloor"));
    // Horizontal light above
    s.faces.push_back(make_face(
        {0.8f,2.5f,-0.5f},{1.2f,2.5f,-0.5f},{1.2f,2.5f,-0.9f},{0.8f,2.5f,-0.9f},
        le,black,"Light"));
    // Back wall (occluder)
    s.faces.push_back(make_face(
        {0,0,-1},{2,0,-1},{2,3,-1},{0,3,-1},black,white,"Back"));
    DiscMesh::apply(s);
    int tc=0;
    for (const auto& p : s.patches) if (p.face_id==0) ++tc;
    std::printf("  tilted floor patches: %d\n", tc);
    CHECK(tc>=1);
    CHECK(check_invariants(s,"C03")==0);
}

static void test_c04_wide_floor()
{
    SUITE("C-04  Synthetic scene — wide (2:1) rectangular floor area conservation");
    Scene s = make_simple_room(4,2,2, 1.8f,2.2f,0.8f,1.2f);
    DiscMesh::apply(s);
    const int nf=static_cast<int>(s.faces.size());
    std::vector<float> ps(nf,0.f);
    for (const auto& p : s.patches)
        if (p.face_id>=0 && p.face_id<nf) ps[p.face_id]+=p.area;
    int bad=0;
    for (int fi=0; fi<nf; ++fi) {
        const Face& f=s.faces[fi];
        float fa  = quad_area(f.verts[0],f.verts[1],f.verts[2],f.verts[3]);
        float rel = std::abs(ps[fi]-fa)/(fa+1e-12f);
        std::printf("  %-10s  face=%.4f  patches=%.4f  rel=%.4f\n",
                    f.name.c_str(), fa, ps[fi], rel);
        if (rel>0.01f) { ++bad; std::printf("    *** area mismatch ***\n"); }
    }
    CHECK(bad==0);
    CHECK(check_invariants(s,"C04")==0);
}

static void test_c05_symmetric_occluders()
{
    SUITE("C-05  Synthetic scene — symmetric occluder pair; equal floor halves");
    Scene s = make_simple_room(2,2,2, 0.8f,1.2f,0.8f,1.2f);
    // Left slab at x=0.5, normal +X (facing right toward light at x=1)
    // (v1-v0)×(v3-v0): v0=(0.5,0,0.7), v1=(0.5,1.2,0.7), v3=(0.5,0,1.3)
    // = (0,1.2,0)×(0,0,0.6) = (1.2*0.6-0*0, 0*0-0*0.6, 0*0-1.2*0) = (0.72,0,0) → +X ✓
    s.faces.push_back(make_face(
        {0.5f,0,0.7f},{0.5f,1.2f,0.7f},{0.5f,1.2f,1.3f},{0.5f,0,1.3f},
        {},{0.5f,0.5f,0.5f},"OccL"));
    // Right slab at x=1.5, normal -X (facing left toward light at x=1)
    // v0=(1.5,0,1.3),v1=(1.5,1.2,1.3),v3=(1.5,0,0.7)
    // (v1-v0)×(v3-v0) = (0,1.2,0)×(0,0,-0.6) = (-0.72,0,0) → -X ✓
    s.faces.push_back(make_face(
        {1.5f,0,1.3f},{1.5f,1.2f,1.3f},{1.5f,1.2f,0.7f},{1.5f,0,0.7f},
        {},{0.5f,0.5f,0.5f},"OccR"));
    DiscMesh::apply(s);

    int floor_fi=-1;
    for (int fi=0; fi<static_cast<int>(s.faces.size()); ++fi)
        if (s.faces[fi].name==std::string("Floor")) { floor_fi=fi; break; }
    CHECK(floor_fi>=0);

    float lasum=0, rasum=0;
    for (const auto& p : s.patches) {
        if (p.face_id!=floor_fi) continue;
        if (p.center.x<1.f) lasum+=p.area; else rasum+=p.area;
    }
    float asym = std::abs(lasum-rasum)/(lasum+rasum+1e-9f);
    std::printf("  left=%.5f  right=%.5f  asymmetry=%.5f\n", lasum, rasum, asym);
    CHECK(asym<0.02f); // 2% tolerance — meshing asymmetry from finite critical-line resolution
    CHECK(check_invariants(s,"C05")==0);
}

static void test_c06_rotated_room()
{
    SUITE("C-06  Synthetic scene — 45° rotated room (no axis-aligned faces)");
    const float sq=0.70710678f;
    auto rot=[&](Vec3 v)->Vec3{ return {(v.x-v.z)*sq, v.y, (v.x+v.z)*sq}; };
    struct Q { Vec3 v[4]; const char* name; bool emit; };
    std::vector<Q> qs = {
        {{{0,0,0},{1,0,0},{1,0,1},{0,0,1}}, "Floor",   false},
        {{{0,1,0},{0,1,1},{1,1,1},{1,1,0}}, "Ceiling", false},
        {{{0,0,1},{1,0,1},{1,1,1},{0,1,1}}, "Back",    false},
        {{{1,0,0},{0,0,0},{0,1,0},{1,1,0}}, "Front",   false},
        {{{0,0,0},{0,0,1},{0,1,1},{0,1,0}}, "Left",    false},
        {{{1,0,1},{1,0,0},{1,1,0},{1,1,1}}, "Right",   false},
        {{{0.4f,0.999f,0.4f},{0.6f,0.999f,0.4f},
          {0.6f,0.999f,0.6f},{0.4f,0.999f,0.6f}}, "Light", true},
    };
    const Vec3 black{}, white{0.8f,0.8f,0.8f}, le{10,10,10};
    Scene sc;
    for (const auto& q : qs)
        sc.faces.push_back(make_face(rot(q.v[0]),rot(q.v[1]),rot(q.v[2]),rot(q.v[3]),
                                     q.emit?le:black, q.emit?black:white, q.name));
    DiscMesh::apply(sc);
    std::printf("  patches=%zu\n", sc.patches.size());
    CHECK(check_invariants(sc,"C06")==0);
}

static void test_c07_tiny_room()
{
    SUITE("C-07  Synthetic scene — tiny coordinates (0.1-scale room)");
    Scene s = make_simple_room(0.1f,0.1f,0.1f, 0.04f,0.06f,0.04f,0.06f);
    DiscMesh::apply(s);
    std::printf("  patches=%zu\n", s.patches.size());
    CHECK(check_invariants(s,"C07")==0);
}

static void test_c08_large_room()
{
    SUITE("C-08  Synthetic scene — large coordinates (100-scale room)");
    Scene s = make_simple_room(100,100,100, 40,60,40,60);
    DiscMesh::apply(s);
    std::printf("  patches=%zu\n", s.patches.size());
    CHECK(check_invariants(s,"C08")==0);
}

static void test_c09_multiple_occluders()
{
    SUITE("C-09  Synthetic scene — 4 vertical pillar occluders");
    Scene s = make_simple_room(3,3,3, 1.2f,1.8f,1.2f,1.8f);
    struct Pillar { float cx,cz; const char* name; };
    for (const auto& pl : std::vector<Pillar>{{0.9f,0.9f,"P1"},{2.1f,0.9f,"P2"},
                                               {0.9f,2.1f,"P3"},{2.1f,2.1f,"P4"}}) {
        const float r=0.1f;
        s.faces.push_back(make_face(
            {pl.cx-r,0,pl.cz+r},{pl.cx+r,0,pl.cz+r},
            {pl.cx+r,2,pl.cz+r},{pl.cx-r,2,pl.cz+r},
            {},{0.6f,0.6f,0.6f}, pl.name));
    }
    DiscMesh::apply(s);
    std::printf("  patches=%zu\n", s.patches.size());
    CHECK(check_invariants(s,"C09")==0);
}

// ============================================================
// main
// ============================================================
int main()
{
    // Section A
    test_vec3();
    test_vec2();
    test_face_frame_axis_aligned();
    test_face_frame_tilted();
    test_face_project_interior();
    test_all_cornell_face_frames();

    // Section B
    test_b01_emissive_single_patch();
    test_b02_box_start();
    test_b03_area_conservation();
    test_b04_green_wall();
    test_b05_non_emissive_split();
    test_b06_all_invariants();

    // Section C
    test_c01_minimal_room();
    test_c02_single_occluder();
    test_c03_tilted_receiver();
    test_c04_wide_floor();
    test_c05_symmetric_occluders();
    test_c06_rotated_room();
    test_c07_tiny_room();
    test_c08_large_room();
    test_c09_multiple_occluders();

    print_summary();
    return exit_code();
}

