#pragma once
#include "scene.h"

/*
 * Discontinuity Meshing  (Lischinski, Tampieri & Greenberg 1992)
 * --------------------------------------------------------------
 * apply() is an alternative to Scene::meshUniform().  It analytically
 * computes D0 critical lines (shadow / penumbra boundaries) on each
 * receiver face, splits each face into convex sub-polygons along those
 * lines, and emits one Scene::Patch per resulting polygon.  The result
 * aligns patch edges with light–shadow transitions, reducing banding
 * artefacts that appear when patch edges cut across soft shadow gradients.
 *
 * Algorithm outline
 * -----------------
 * For each non-emissive receiver face R:
 *   1. Collect all D0 critical segments on R.
 *      For each (light L, occluder O) triple that passes visibility filters:
 *        Build wedge planes from light-vertex × occluder-edges and
 *        light-edges × occluder-vertex combinations.
 *        Intersect each wedge plane with the receiver plane (gives a 3-D line).
 *        Project the line to R's face-local UV; clip to R's quad boundary
 *        using Cyrus-Beck (1978) parametric clipping.
 *   2. Deduplicate segments by Hessian normal form (angle, signed-distance)
 *      to avoid near-identical lines causing exponential polygon growth.
 *   3. Split R's polygon list along every unique segment using a bidirectional
 *      Sutherland-Hodgman (1974) 2-D polygon splitter.
 *   4. Emit one Scene::Patch per resulting convex polygon via fan-decomposition.
 *
 * Occluder relevance filters (applied before wedge computation)
 * -------------------------------------------------------------
 *   - O must face toward the light centroid (front-face from light's view).
 *   - O must face toward the receiver centroid (front-face from R's view).
 *   - R must face toward the occluder centroid (R can "see" O at all).
 * These three dot-product tests eliminate walls and faces that cannot
 * cast or receive shadows relevant to the (L, R) pair.
 *
 * Limitations / scope
 * -------------------
 *   - Only D0 critical lines are computed (vertex × edge and edge × vertex
 *     combinations).  D1 lines (edge × edge, requiring tangent-plane
 *     continuity) are NOT generated; see Lischinski 1992, Section 4.
 *   - Faces are assumed to be convex quads; the algorithm does not handle
 *     non-planar or concave patches.
 *   - The duplicate-line test is O(N²) per receiver; for large scenes a
 *     spatial hash on (angle, distance) would be faster.
 *
 * Suggested future additions
 * --------------------------
 *   - D1 critical lines from edge–edge interactions for smooth surfaces.
 *   - Visibility precomputation (e.g. bounding-box overlap) to cull
 *     (L, O, R) triples before wedge plane construction.
 *   - Adaptive mesh refinement guided by radiosity gradient magnitude
 *     (Lischinski 1992, Section 5), instead of the current one-shot split.
 *   - Support for triangulated or non-quad input geometry.
 *
 * Pre-condition : scene.faces[] must be populated (call buildCornellBox()).
 * Post-condition: scene.patches[] filled; scene.box_start set.
 *                 Emissive faces receive a single unsubdivided patch.
 */
namespace DiscMesh {
    void apply(Scene& scene);
}
