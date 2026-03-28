#include "discmesh.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
#include <utility>

/*
 * Discontinuity Meshing — D0 critical lines only.
 *
 * Algorithm (Lischinski et al. 1992, Section 3):
 *
 *  For each receiver face R:
 *    For each emissive face L (light):
 *      For each other face O (potential occluder, not R, not L):
 *        For each vertex v_l of L and each edge (v0_o,v1_o) of O:
 *          Wedge plane π = plane(v_l, v0_o, v1_o)
 *          Intersect π with receiver plane → 3-D line
 *          Project line endpoints to receiver UV space
 *          Clip to receiver quad boundary
 *          If non-degenerate: record as critical segment
 *        (Repeat with light edges × occluder vertices)
 *      Split the receiver's current polygon list along all critical segs
 *    Emit one Patch per resulting convex polygon
 *
 * Notes:
 *  - Occluder backface culling: only use O if it faces the light
 *  - Degenerate wedge planes (zero-area triangle) are skipped
 *  - Each resulting polygon may have 3..N vertices; we fan-divide into
 *    quads (4-vert) and triangles stored as degenerate quads (v3==v2)
 */

// ============================================================
// Polygon splitter — bidirectional Sutherland-Hodgman (2D)
//
// Sutherland-Hodgman (1974): clip a polygon against a half-plane by
// walking its edge list and emitting vertices/intersections only for
// the kept side.  Standard S-H discards one side; here both halves
// are built simultaneously in a single edge-walk pass:
//   left  ← vertices with sideDist ≥ 0  (positive / left-of-line)
//   right ← vertices with sideDist < 0  (negative / right-of-line)
// Edge–line crossings are appended to both lists.  The sideDist sign
// follows the 2-D cross product: (b-a) × (p-a) > 0 ⟺ p is left of a→b.
// ============================================================

// Minimum 2D polygon area (shoelace) below which the polygon is considered
// degenerate.  Faces are ~1 UV unit wide; this discards polygons covering
// less than 1/(128²) of the face, which would be sub-pixel at any practical
// hemicube resolution.
static constexpr float POLY_AREA_EPS = 6.1e-5f; // 1/(128²) ≈ 6.1e-5

struct Poly2D {
    std::vector<Vec2> verts;
    explicit Poly2D(std::initializer_list<Vec2> il) : verts(il) {}
    explicit Poly2D(std::vector<Vec2> v) : verts(std::move(v)) {}

    // Shoelace (signed) area — positive for CCW, negative for CW.
    float signedArea() const {
        const int n = static_cast<int>(verts.size());
        float a = 0.f;
        for (int i = 0; i < n; ++i) {
            const Vec2& p = verts[i];
            const Vec2& q = verts[(i + 1) % n];
            a += p.u * q.v - q.u * p.v;
        }
        return 0.5f * a;
    }
    // A polygon is empty if it has fewer than 3 vertices or its area is
    // negligible (degenerate sliver from nearly-coincident split lines).
    bool empty() const {
        return verts.size() < 3 || std::abs(signedArea()) < POLY_AREA_EPS;
    }
};

// Signed distance of p from the line a→b (positive = left side)
static float sideDist(Vec2 a, Vec2 b, Vec2 p)
{
    Vec2 ab = b - a;
    Vec2 ap = p - a;
    return ab.cross(ap); // positive when p is left of a→b
}

// Intersection parameter t ∈ [0,1] of segment p→q with line a→b
static float intersectT(Vec2 p, Vec2 q, Vec2 a, Vec2 b)
{
    Vec2 pq = q - p;
    Vec2 ab = b - a;
    float denom = ab.cross(pq);
    if (std::abs(denom) < 1e-10f) return 0.5f; // parallel
    Vec2 pa = a - p;
    return ab.cross(pa) / denom;
}

// Bidirectional Sutherland-Hodgman split: infinite line a→b divides poly into
// 'left' (sideDist ≥ 0) and 'right' (sideDist < 0) halves.
// Returns false if the line does not actually cross the polygon interior
// (one half would have fewer than 3 vertices or negligible area).
// Note: vertices exactly on the split line (sideDist == 0) are assigned to
// the left polygon; both output polygons share those boundary vertices.
//
// Post-pass vertex deduplication: after collecting each half, consecutive
// vertices within VERT_EPS are merged.  This prevents accumulated near-
// coincident vertices from many nearly-parallel critical lines from producing
// zero-area sliver polygons, which would later generate patches with garbage
// normals and permanently-zero radiosity.
static constexpr float VERT_EPS = 1e-5f; // merge vertices closer than this in UV

static bool splitPoly(const Poly2D& poly, Vec2 a, Vec2 b,
                       Poly2D& left, Poly2D& right)
{
    const int n = static_cast<int>(poly.verts.size());
    std::vector<Vec2> lv, rv;

    for (int i = 0; i < n; ++i) {
        Vec2 p = poly.verts[i];
        Vec2 q = poly.verts[(i + 1) % n];
        float dp = sideDist(a, b, p);
        float dq = sideDist(a, b, q);

        if (dp >= 0.f) lv.push_back(p);
        if (dp <= 0.f) rv.push_back(p);

        // If the edge crosses the line, emit the intersection on both sides
        if ((dp > 0.f && dq < 0.f) || (dp < 0.f && dq > 0.f)) {
            float t = intersectT(p, q, a, b);
            Vec2 ix = p + (q - p) * t;
            lv.push_back(ix);
            rv.push_back(ix);
        }
    }

    // Deduplicate consecutive vertices closer than VERT_EPS (wrap-around
    // handled by checking after each push and at the seam v[last] vs v[0]).
    auto dedup = [](std::vector<Vec2>& verts) {
        std::vector<Vec2> out;
        out.reserve(verts.size());
        for (const auto& v : verts) {
            if (!out.empty() && (v - out.back()).length() < VERT_EPS)
                continue;
            out.push_back(v);
        }
        // Seam: last vs first
        while (out.size() >= 2 &&
               (out.back() - out.front()).length() < VERT_EPS)
            out.pop_back();
        verts = std::move(out);
    };
    dedup(lv);
    dedup(rv);

    Poly2D l(std::move(lv)), r(std::move(rv));
    if (l.empty() || r.empty()) return false;

    left  = std::move(l);
    right = std::move(r);
    return true;
}

// ============================================================
// Clip a segment (a,b) so it lies fully inside the convex polygon boundary.
// Returns false if the segment is entirely outside.
//
// Algorithm: Cyrus-Beck (1978) parametric line clipping against a convex
// polygon.  (Liang-Barsky (1984) is a specialisation of Cyrus-Beck for
// axis-aligned rectangles only; the general convex-polygon form is
// attributed to Cyrus & Beck.)
//
// Parametric form: P(t) = a + t·d,  t ∈ [t0, t1], initially [0, 1].
// For each edge e0→e1 we compute t at which P(t) crosses the edge line
// and tighten [t0, t1] to keep only the interior portion.
//
// Normal convention (adjustment vs. textbook Cyrus-Beck):
// The textbook uses the INWARD normal n_in, where n_in·(P-e0) ≥ 0 is
// the "inside" condition.  With that convention, denom = n_in·d > 0
// means the segment is entering the half-plane (t_entry → update t0).
//
// This implementation uses the OUTWARD normal:
//   en = (e1.v-e0.v, -(e1.u-e0.u))  [CW ⊥ of the CCW edge direction]
// which equals −n_in.  The inside condition becomes en·(P-e0) ≤ 0, and
// the entry/exit roles of denom's sign are reversed:
//   denom = en·d < 0  →  entry (update t0 = max)
//   denom = en·d > 0  →  exit  (update t1 = min)
// For the parallel-outside case: num = en·(e0-a); num < 0 means a lies
// strictly outside the half-plane (en points away from the interior), so
// the entire segment is outside.
// ============================================================
static bool clipSegToConvexPoly(Vec2& a, Vec2& b, const Poly2D& boundary)
{
    float t0 = 0.f, t1 = 1.f;
    Vec2 d = b - a;
    const int n = static_cast<int>(boundary.verts.size());
    for (int i = 0; i < n; ++i) {
        Vec2 e0 = boundary.verts[i];
        Vec2 e1 = boundary.verts[(i + 1) % n];
        // Outward normal of edge e0→e1 for a CCW polygon (see note above)
        Vec2 en = {e1.v - e0.v, -(e1.u - e0.u)};
        float denom = en.dot(d);
        float num   = en.dot(e0 - a);
        if (std::abs(denom) < 1e-10f) {
            if (num < 0.f) return false; // parallel and outside
            continue;
        }
        float t = num / denom;
        // With outward normal: denom<0 → entry, denom>0 → exit
        if (denom < 0.f) { if (t > t0) t0 = t; }
        else             { if (t < t1) t1 = t; }
        if (t0 > t1) return false;
    }
    if (t1 - t0 < 1e-6f) return false;
    Vec2 na = a + d * t0;
    Vec2 nb = a + d * t1;
    a = na; b = nb;
    return true;
}

// ============================================================
// Critical line computation (D0) — Lischinski et al. 1992
//
// D0 discontinuities are caused by vertices and edges of light sources
// and occluders.  Each critical line on a receiver surface is the
// projection of a "wedge plane" — a plane through one feature of the
// light and one feature of an occluder (vertex × edge or edge × vertex).
// ============================================================

// Intersect two planes (each defined by a unit normal and a point) to obtain
// a line in 3-D.  Returns {point, direction}; direction is zero for parallel planes.
//
// Derivation: the line direction is d = n1 × n2.  The point on both planes
// that is also in the plane through the origin perpendicular to d satisfies
// the 3×3 system  [n1; n2; d]·x = [n1·p1; n2·p2; 0].  For unit normals the
// 2×2 block reduces to a 2-equation system with determinant 1 − (n1·n2)²,
// giving x = c1·n1 + c2·n2 — the minimum-distance-to-origin point on the line.
static std::pair<Vec3,Vec3> intersectPlanes(Vec3 n1, Vec3 p1, Vec3 n2, Vec3 p2)
{
    Vec3 dir = n1.cross(n2);
    if (dir.length() < 1e-9f) return {{}, {}};

    float d1  = n1.dot(p1);
    float d2  = n2.dot(p2);
    float n12 = n1.dot(n2);
    float det = 1.f - n12 * n12; // = sin²θ between planes; 0 if parallel
    if (std::abs(det) < 1e-9f) return {{}, {}};
    float c1 = (d1 - d2 * n12) / det;
    float c2 = (d2 - d1 * n12) / det;
    Vec3  pt  = n1 * c1 + n2 * c2;
    return {pt, dir.normalized()};
}

// Project the intersection line of the wedge plane and the receiver plane
// to receiver UV, then clip to the receiver quad boundary.
// Returns true and fills seg_a, seg_b if a non-degenerate segment exists.
static bool wedgeToSeg(Vec3 wedge_n, Vec3 wedge_p,
                       const Face& receiver,
                       const Poly2D& recv_boundary,
                       Vec2& seg_a, Vec2& seg_b)
{
    // Intersect wedge plane with receiver plane
    auto [lp, ld] = intersectPlanes(wedge_n, wedge_p, receiver.normal, receiver.verts[0]);
    if (ld.length() < 1e-9f) return false;

    // Pick two points far along the line and project to receiver UV
    constexpr float EXTENT = 10.f;
    Vec2 a = receiver.project(lp - ld * EXTENT);
    Vec2 b = receiver.project(lp + ld * EXTENT);

    if (!clipSegToConvexPoly(a, b, recv_boundary)) return false;
    if ((b - a).length() < 1e-6f) return false;

    seg_a = a;
    seg_b = b;
    return true;
}

// Compute the initial receiver boundary polygon in face-local UV.
static Poly2D receiverBoundary(const Face& f)
{
    return Poly2D({ f.project(f.verts[0]),
                    f.project(f.verts[1]),
                    f.project(f.verts[2]),
                    f.project(f.verts[3]) });
}

// Emit patches from a convex polygon by fan-decomposition from verts[0]:
//   triangle  n==3: one degenerate quad (v2 repeated as v3)
//   quad      n==4: one regular quad
//   n-gon     n>4:  (n-2) triangles, each a fan spoke from verts[0]
// All resulting patches share the face's emission/reflectance and face_id.
// Helper: fix a just-emitted triangle patch's normal and center.
// addPatch computes normal = (v1-v0)×(v3-v0) and center = (v0+v1+v2+v3)/4.
// For disc-mesh patches, both must come from the parent face:
//   • Normal must be the face's outward normal (vertex cross-product winding
//     may produce the opposite sign for non-convex or CW-residual polygons).
//   • Triangles store v3=v2, so center = (v0+v1+2·v2)/4 — biased toward v2
//     and potentially outside the triangle; the true centroid is (v0+v1+v2)/3.
static void fixTriPatch(Scene& scene, const Vec3& f_normal,
                        const Vec3& p0, const Vec3& p1, const Vec3& p2)
{
    Patch& p = scene.patches.back();
    p.normal = f_normal;
    p.center = (p0 + p1 + p2) * (1.f / 3.f);
    (void)p;
}

static void emitPolygon(Scene& scene, const Poly2D& poly, const Face& f, int face_id)
{
    const auto& verts_ref = poly.verts;
    const int n = static_cast<int>(verts_ref.size());
    if (n < 3) return;

    // Guard: skip sub-pixel slivers (should already be filtered by splitPoly,
    // but a second check here ensures nothing slips through).
    float sa = poly.signedArea();
    if (std::abs(sa) < POLY_AREA_EPS) return;

    // Ensure CCW winding (positive signed area) so the 3D patch normal computed
    // by addPatch matches the face's outward normal.  splitPoly preserves winding
    // in the left half but flips it in the right half after deduplication.
    const std::vector<Vec2>* vp = &verts_ref;
    std::vector<Vec2> flipped;
    if (sa < 0.f) {
        flipped.assign(verts_ref.rbegin(), verts_ref.rend());
        vp = &flipped;
    }
    const std::vector<Vec2>& verts = *vp;

    if (n == 3) {
        Vec3 p0 = f.unproject(verts[0]);
        Vec3 p1 = f.unproject(verts[1]);
        Vec3 p2 = f.unproject(verts[2]);
        scene.addPatch(p0, p1, p2, p2, f.emission, f.reflectance, face_id);
        fixTriPatch(scene, f.normal, p0, p1, p2);
    } else if (n == 4) {
        Vec3 p0 = f.unproject(verts[0]);
        Vec3 p1 = f.unproject(verts[1]);
        Vec3 p2 = f.unproject(verts[2]);
        Vec3 p3 = f.unproject(verts[3]);
        scene.addPatch(p0, p1, p2, p3, f.emission, f.reflectance, face_id);
        // Quad normal: force to face normal for the same reasons as above.
        scene.patches.back().normal = f.normal;
    } else {
        // Fan-decompose: triangle fan from verts[0]
        Vec3 p0 = f.unproject(verts[0]);
        for (int i = 1; i + 1 < n; ++i) {
            Vec3 pi  = f.unproject(verts[i]);
            Vec3 pi1 = f.unproject(verts[i + 1]);
            scene.addPatch(p0, pi, pi1, pi1, f.emission, f.reflectance, face_id);
            fixTriPatch(scene, f.normal, p0, pi, pi1);
        }
    }
}

// ============================================================
// Main entry point
// ============================================================

namespace DiscMesh {

void apply(Scene& scene)
{
    scene.patches.clear();
    scene.box_start = -1;

    const int nf = static_cast<int>(scene.faces.size());

    // Identify light faces
    std::vector<int> light_indices;
    for (int i = 0; i < nf; ++i) {
        const Face& f = scene.faces[i];
        if (f.emission.x > 0.f || f.emission.y > 0.f || f.emission.z > 0.f)
            light_indices.push_back(i);
    }

    for (int ri = 0; ri < nf; ++ri) {
        const Face& R = scene.faces[ri];

        if (R.is_box && scene.box_start < 0)
            scene.box_start = static_cast<int>(scene.patches.size());

        // Emissive faces → single patch, no splitting
        bool emissive = R.emission.x > 0.f || R.emission.y > 0.f || R.emission.z > 0.f;
        if (emissive) {
            scene.addPatch(R.verts[0], R.verts[1], R.verts[2], R.verts[3],
                           R.emission, R.reflectance, ri);
            continue;
        }

        // Start with the receiver quad as a single polygon
        std::vector<Poly2D> polys = { receiverBoundary(R) };
        Poly2D recv_boundary = receiverBoundary(R);

        // Collect ALL unique critical segments for this receiver first,
        // then apply them once (deduplication prevents exponential growth).
        std::vector<std::pair<Vec2,Vec2>> all_segs;

        // For each light × each occluder
        for (int li : light_indices) {
            const Face& L = scene.faces[li];
            // Light centroid
            Vec3 Lc = (L.verts[0] + L.verts[1] + L.verts[2] + L.verts[3]) * 0.25f;

            for (int oi = 0; oi < nf; ++oi) {
                if (oi == ri) continue;
                const Face& O = scene.faces[oi];
                // Emissive faces (light sources) are sources, not occluders.
                if (O.emission.x > 0.f || O.emission.y > 0.f || O.emission.z > 0.f) continue;

                // --- Occluder relevance filters ---
                // 1. O must face toward the light (front-face from light's view)
                Vec3 Oc = (O.verts[0] + O.verts[1] + O.verts[2] + O.verts[3]) * 0.25f;
                if (O.normal.dot(Lc - Oc) <= 0.f) continue;

                // 2. O must face toward the receiver (front-face from receiver's view)
                Vec3 Rc = (R.verts[0] + R.verts[1] + R.verts[2] + R.verts[3]) * 0.25f;
                if (O.normal.dot(Rc - Oc) <= 0.f) continue;

                // 3. R must face toward the occluder (receiver can "see" occluder)
                if (R.normal.dot(Oc - Rc) <= 0.f) continue;

                // Collect critical segments for this (L, O, R) triple
                // --- Light vertex × occluder edges ---
                for (int lv = 0; lv < 4; ++lv) {
                    Vec3 vl = L.verts[lv];
                    for (int oe = 0; oe < 4; ++oe) {
                        Vec3 v0 = O.verts[oe];
                        Vec3 v1 = O.verts[(oe + 1) % 4];
                        // Wedge plane through vl, v0, v1
                        Vec3 wn = (v0 - vl).cross(v1 - vl);
                        if (wn.length() < 1e-9f) continue;
                        Vec2 sa, sb;
                        if (wedgeToSeg(wn, vl, R, recv_boundary, sa, sb))
                            all_segs.push_back({sa, sb});
                    }
                }

                // --- Light edges × occluder vertices ---
                for (int le = 0; le < 4; ++le) {
                    Vec3 l0 = L.verts[le];
                    Vec3 l1 = L.verts[(le + 1) % 4];
                    for (int ov = 0; ov < 4; ++ov) {
                        Vec3 vo = O.verts[ov];
                        Vec3 wn = (l0 - vo).cross(l1 - vo);
                        if (wn.length() < 1e-9f) continue;
                        Vec2 sa, sb;
                        if (wedgeToSeg(wn, vo, R, recv_boundary, sa, sb))
                            all_segs.push_back({sa, sb});
                    }
                }
            }
        }

        // Deduplicate critical lines using Hessian normal form.
        // Each segment defines an infinite line; two segments represent the
        // same line if their (angle, signed_dist) pairs agree within LINE_EPS.
        // Without deduplication, near-parallel wedge planes from different
        // (L, O) triples produce nearly-identical lines that each split every
        // sub-polygon in the growing list, causing exponential patch growth.
        constexpr float LINE_EPS = 2e-3f; // tolerance in UV units (face ≈ 1 unit wide)
        std::vector<std::pair<float,float>> used_lines; // (angle, dist) canonical form
        std::vector<std::pair<Vec2,Vec2>> unique_segs;
        for (const auto& [sa, sb] : all_segs) {
            Vec2 d = sb - sa;
            float len = d.length();
            if (len < 1e-6f) continue;
            // Canonical line: normal form n·p = c, with n = (-dy, dx)/len
            // Hessian normal form: unit normal n = (-dv, du)/|d|, offset c = n·sa.
            // Orient n so the angle lies in [0, π) for a unique representation.
            float nx = -d.v / len, ny = d.u / len;
            if (nx < 0.f || (std::abs(nx) < 1e-9f && ny < 0.f)) { nx = -nx; ny = -ny; }
            float c = nx * sa.u + ny * sa.v;
            float angle = std::atan2(ny, nx); // ∈ [0, π)
            // Check against already-used lines
            bool dup = false;
            for (const auto& [ua, uc] : used_lines) {
                if (std::abs(ua - angle) < LINE_EPS && std::abs(uc - c) < LINE_EPS) {
                    dup = true; break;
                }
            }
            if (!dup) {
                used_lines.push_back({angle, c});
                unique_segs.push_back({sa, sb});
            }
        }

        // Apply unique segments to split the polygon list
        for (const auto& [sa, sb] : unique_segs) {
                std::vector<Poly2D> next;
                next.reserve(polys.size() * 2);
                for (auto& poly : polys) {
                    Poly2D lp({});
                    Poly2D rp({});
                    if (splitPoly(poly, sa, sb, lp, rp)) {

                        next.push_back(std::move(lp));
                        next.push_back(std::move(rp));
                    } else {
                        next.push_back(std::move(poly));
                    }
                }
                polys = std::move(next);
        }

        // Emit one patch per resulting convex polygon
        for (const auto& poly : polys)
            emitPolygon(scene, poly, R, ri);
    }

    if (scene.box_start < 0)
        scene.box_start = static_cast<int>(scene.patches.size());

    std::printf("  DiscMesh: %d faces → %zu patches\n",
                nf, scene.patches.size());
    std::fflush(stdout);
}

} // namespace DiscMesh
