// mesh.cpp
// Rcpp wrapper around CDT constrained Delaunay triangulation

// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
// Suppress -Wunused-parameter in vendored CDT headers (callback interface)
#if defined(__GNUC__) || defined(__clang__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
#include "cdt/CDT.h"
#if defined(__GNUC__) || defined(__clang__)
#  pragma GCC diagnostic pop
#endif
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>

using namespace Rcpp;

// [[Rcpp::export]]
List cpp_triangulate(
    NumericMatrix points,
    Nullable<IntegerMatrix> edges_nullable = R_NilValue,
    double min_dist_tolerance = 0.0
) {
    int n_pts = points.nrow();
    if (n_pts < 3) stop("Need at least 3 points for triangulation");
    if (points.ncol() != 2) stop("points must be an N x 2 matrix");

    // Build CDT vertex list
    std::vector<CDT::V2d<double>> vertices(n_pts);
    for (int i = 0; i < n_pts; i++) {
        vertices[i] = CDT::V2d<double>(points(i, 0), points(i, 1));
    }

    // Build edge constraint list (if provided)
    std::vector<CDT::Edge> constraint_edges;
    if (edges_nullable.isNotNull()) {
        IntegerMatrix edges = as<IntegerMatrix>(edges_nullable);
        int n_edges = edges.nrow();
        for (int i = 0; i < n_edges; i++) {
            // R is 1-based, CDT is 0-based
            CDT::VertInd v1 = static_cast<CDT::VertInd>(edges(i, 0) - 1);
            CDT::VertInd v2 = static_cast<CDT::VertInd>(edges(i, 1) - 1);
            constraint_edges.push_back(CDT::Edge(v1, v2));
        }
    }

    // Create triangulation
    CDT::Triangulation<double> cdt(
        CDT::VertexInsertionOrder::Auto,
        CDT::IntersectingConstraintEdges::TryResolve,
        min_dist_tolerance
    );

    cdt.insertVertices(vertices);
    if (!constraint_edges.empty()) {
        cdt.insertEdges(constraint_edges);
    }
    cdt.eraseOuterTrianglesAndHoles();

    // Extract results
    const auto& tri_verts = cdt.vertices;
    const auto& triangles = cdt.triangles;

    int n_verts = static_cast<int>(tri_verts.size());
    int n_tri = static_cast<int>(triangles.size());

    // Vertices (may include Steiner points added by CDT)
    NumericMatrix out_vertices(n_verts, 2);
    for (int i = 0; i < n_verts; i++) {
        out_vertices(i, 0) = tri_verts[i].x;
        out_vertices(i, 1) = tri_verts[i].y;
    }

    // Triangles (1-based indices for R)
    IntegerMatrix out_triangles(n_tri, 3);
    for (int i = 0; i < n_tri; i++) {
        out_triangles(i, 0) = static_cast<int>(triangles[i].vertices[0]) + 1;
        out_triangles(i, 1) = static_cast<int>(triangles[i].vertices[1]) + 1;
        out_triangles(i, 2) = static_cast<int>(triangles[i].vertices[2]) + 1;
    }

    // Edge list (unique edges from triangles, 1-based)
    std::set<std::pair<int,int>> edge_set;
    for (int i = 0; i < n_tri; i++) {
        for (int j = 0; j < 3; j++) {
            int v1 = out_triangles(i, j);
            int v2 = out_triangles(i, (j + 1) % 3);
            if (v1 > v2) std::swap(v1, v2);
            edge_set.insert({v1, v2});
        }
    }
    IntegerMatrix out_edges(edge_set.size(), 2);
    int idx = 0;
    for (const auto& e : edge_set) {
        out_edges(idx, 0) = e.first;
        out_edges(idx, 1) = e.second;
        idx++;
    }

    return List::create(
        Named("vertices") = out_vertices,
        Named("triangles") = out_triangles,
        Named("edges") = out_edges,
        Named("n_vertices") = n_verts,
        Named("n_triangles") = n_tri,
        Named("n_edges") = static_cast<int>(edge_set.size()),
        Named("n_input_points") = n_pts
    );
}

// ---------------------------------------------------------------------------
// Ruppert refinement: iteratively insert circumcenters of bad triangles
// until all triangles have min angle >= threshold.
// ---------------------------------------------------------------------------

// Circumcenter of triangle (px0,py0), (px1,py1), (px2,py2)
static inline void circumcenter(
    double x0, double y0, double x1, double y1, double x2, double y2,
    double& cx, double& cy
) {
    double ax = x0, ay = y0;
    double bx = x1, by = y1;
    double ccx = x2, ccy = y2;
    double D = 2.0 * (ax * (by - ccy) + bx * (ccy - ay) + ccx * (ay - by));
    if (std::abs(D) < 1e-30) {
        // Degenerate: use centroid
        cx = (x0 + x1 + x2) / 3.0;
        cy = (y0 + y1 + y2) / 3.0;
        return;
    }
    double ax2ay2 = ax * ax + ay * ay;
    double bx2by2 = bx * bx + by * by;
    double cx2cy2 = ccx * ccx + ccy * ccy;
    cx = (ax2ay2 * (by - ccy) + bx2by2 * (ccy - ay) + cx2cy2 * (ay - by)) / D;
    cy = (ax2ay2 * (ccx - bx) + bx2by2 * (ax - ccx) + cx2cy2 * (bx - ax)) / D;
}

// Minimum angle (radians) of a triangle given its three edge lengths
static inline double min_angle_from_lengths(double l0, double l1, double l2) {
    // Law of cosines: cos(A) = (b² + c² - a²) / (2bc)
    // Angle opposite l0
    double c0 = (l1 * l1 + l2 * l2 - l0 * l0) / (2.0 * l1 * l2);
    double c1 = (l0 * l0 + l2 * l2 - l1 * l1) / (2.0 * l0 * l2);
    double c2 = (l0 * l0 + l1 * l1 - l2 * l2) / (2.0 * l0 * l1);
    c0 = std::max(-1.0, std::min(1.0, c0));
    c1 = std::max(-1.0, std::min(1.0, c1));
    c2 = std::max(-1.0, std::min(1.0, c2));
    return std::min({std::acos(c0), std::acos(c1), std::acos(c2)});
}

// Off-center point for Delaunay refinement (Üngör 2009).
// Given a bad triangle, computes the off-center: a point on the perpendicular
// bisector of the shortest edge at distance beta*l from each endpoint.
// Returns the closer of {off-center, circumcenter} to the shortest edge.
// This guarantees the new triangle on the shortest edge has min angle >= theta_min.
static inline void offcenter(
    double x0, double y0, double x1, double y1, double x2, double y2,
    double theta_min_rad,
    double& ox, double& oy
) {
    // Edge vectors from v0
    double xdo = x1 - x0, ydo = y1 - y0;
    double xao = x2 - x0, yao = y2 - y0;

    // Squared edge lengths
    double dodist = xdo * xdo + ydo * ydo;  // |v0-v1|²
    double aodist = xao * xao + yao * yao;  // |v0-v2|²
    double dadist = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);  // |v1-v2|²

    // Circumcenter relative to v0
    double denom = 2.0 * (xdo * yao - xao * ydo);
    double dx, dy;
    if (std::abs(denom) < 1e-30) {
        // Degenerate: use centroid
        ox = (x0 + x1 + x2) / 3.0;
        oy = (y0 + y1 + y2) / 3.0;
        return;
    }
    dx = (yao * dodist - ydo * aodist) / denom;
    dy = (xdo * aodist - xao * dodist) / denom;

    // Off-center constant: 0.475 / tan(theta_min / 2)
    double goodangle = std::cos(theta_min_rad);
    double offconst = (goodangle >= 1.0) ? 0.0
        : 0.475 * std::sqrt((1.0 + goodangle) / (1.0 - goodangle));

    // Find shortest edge and compute off-center on its perpendicular bisector
    if (dodist <= aodist && dodist <= dadist) {
        // Shortest: v0-v1
        double dxoff = 0.5 * xdo - offconst * ydo;
        double dyoff = 0.5 * ydo + offconst * xdo;
        // Use off-center if closer to v0 than circumcenter
        if (dxoff * dxoff + dyoff * dyoff < dx * dx + dy * dy) {
            dx = dxoff; dy = dyoff;
        }
    } else if (aodist <= dadist) {
        // Shortest: v0-v2
        double dxoff = 0.5 * xao + offconst * yao;
        double dyoff = 0.5 * yao - offconst * xao;
        if (dxoff * dxoff + dyoff * dyoff < dx * dx + dy * dy) {
            dx = dxoff; dy = dyoff;
        }
    } else {
        // Shortest: v1-v2
        double xda = x2 - x1, yda = y2 - y1;
        double dxoff = 0.5 * xda - offconst * yda;
        double dyoff = 0.5 * yda + offconst * xda;
        // Off-center relative to v1
        double cc_from_v1_x = dx - xdo, cc_from_v1_y = dy - ydo;
        if (dxoff * dxoff + dyoff * dyoff <
            cc_from_v1_x * cc_from_v1_x + cc_from_v1_y * cc_from_v1_y) {
            dx = xdo + dxoff; dy = ydo + dyoff;
        }
    }

    ox = x0 + dx;
    oy = y0 + dy;
}

// Check if point (px, py) is inside a polygon defined by edges
// Uses ray casting
static inline bool point_in_polygon(
    double px, double py,
    const std::vector<CDT::V2d<double>>& poly_verts,
    const std::vector<CDT::Edge>& poly_edges
) {
    int crossings = 0;
    for (const auto& e : poly_edges) {
        double x0 = poly_verts[e.v1()].x, y0 = poly_verts[e.v1()].y;
        double x1 = poly_verts[e.v2()].x, y1 = poly_verts[e.v2()].y;
        if ((y0 <= py && y1 > py) || (y1 <= py && y0 > py)) {
            double t = (py - y0) / (y1 - y0);
            if (px < x0 + t * (x1 - x0)) {
                crossings++;
            }
        }
    }
    return (crossings % 2) == 1;
}

// Encroachment test: point (px,py) encroaches segment (ax,ay)-(bx,by)
// iff it lies strictly inside the diametral circle of the segment,
// i.e. the angle ∠apb > 90° ⟺ dot product (a-p)·(b-p) < 0.
static inline bool point_encroaches_segment(
    double px, double py,
    double ax, double ay, double bx, double by
) {
    return (ax - px) * (bx - px) + (ay - py) * (by - py) < 0.0;
}

// [[Rcpp::export]]
List cpp_ruppert_refine(
    NumericMatrix points,
    Nullable<IntegerMatrix> edges_nullable,
    double min_angle_deg,
    double max_area,
    int max_steiner,
    double min_dist_tolerance
) {
    int n_pts = points.nrow();
    if (n_pts < 3) stop("Need at least 3 points for triangulation");

    double min_angle_rad = min_angle_deg * M_PI / 180.0;

    // Build mutable vertex and edge lists
    std::vector<CDT::V2d<double>> vertices(n_pts);
    for (int i = 0; i < n_pts; i++) {
        vertices[i] = CDT::V2d<double>(points(i, 0), points(i, 1));
    }

    std::vector<CDT::Edge> constraint_edges;
    if (edges_nullable.isNotNull()) {
        IntegerMatrix edges = as<IntegerMatrix>(edges_nullable);
        for (int i = 0; i < edges.nrow(); i++) {
            constraint_edges.push_back(CDT::Edge(
                static_cast<CDT::VertInd>(edges(i, 0) - 1),
                static_cast<CDT::VertInd>(edges(i, 1) - 1)
            ));
        }
    }

    // Store original constraint edges for encroachment checks
    std::vector<CDT::Edge> original_constraints = constraint_edges;

    // Shewchuk (1997) Ruppert refinement with off-centers (Üngör 2009):
    // Phase 1: split ALL encroached segments before touching bad triangles.
    //          A segment is encroached if any non-endpoint vertex lies inside
    //          its diametral circle (angle ∠apb > 90°).
    // Phase 2: find worst bad triangle, compute off-center.
    //          If it encroaches any segment → split those segments instead.
    //          Otherwise insert the off-center.

    int steiner_count = 0;
    int consecutive_skips = 0;
    double tol2 = min_dist_tolerance * min_dist_tolerance;
    bool has_boundary = !constraint_edges.empty();

    for (int iter = 0; iter < max_steiner * 3; iter++) {
        if (consecutive_skips > 200 || steiner_count >= max_steiner) break;

        // Build CDT with current vertices and constraint edges
        CDT::Triangulation<double> cdt(
            CDT::VertexInsertionOrder::Auto,
            CDT::IntersectingConstraintEdges::TryResolve,
            min_dist_tolerance
        );
        cdt.insertVertices(vertices);
        if (has_boundary) {
            cdt.insertEdges(constraint_edges);
        }
        cdt.eraseOuterTrianglesAndHoles();

        const auto& tv = cdt.vertices;
        const auto& tris = cdt.triangles;
        if (tris.empty()) break;

        // Sync: adopt any vertices CDT added during constraint resolution
        for (size_t vi = vertices.size(); vi < tv.size(); vi++) {
            vertices.push_back(tv[vi]);
        }

        // ---- Phase 1: Split encroached segments (priority) ----
        // All segment encroachments must be resolved before processing
        // bad triangles. This is the key insight from Shewchuk (1997) that
        // makes Ruppert refinement work on constrained domains.
        if (has_boundary) {
            bool found_encroached = false;
            for (size_t si = 0; si < constraint_edges.size(); si++) {
                CDT::VertInd sv1 = constraint_edges[si].v1();
                CDT::VertInd sv2 = constraint_edges[si].v2();
                double ax = vertices[sv1].x, ay = vertices[sv1].y;
                double bx = vertices[sv2].x, by = vertices[sv2].y;

                for (size_t vi = 0; vi < tv.size(); vi++) {
                    if (vi == sv1 || vi == sv2) continue;
                    if (!point_encroaches_segment(tv[vi].x, tv[vi].y,
                                                  ax, ay, bx, by)) continue;

                    // Encroached — split segment at midpoint
                    double mx = (ax + bx) / 2.0;
                    double my = (ay + by) / 2.0;

                    // Check proximity to existing vertices
                    bool too_close = false;
                    CDT::VertInd close_v = 0;
                    for (size_t vj = 0; vj < vertices.size(); vj++) {
                        double ddx = mx - vertices[vj].x;
                        double ddy = my - vertices[vj].y;
                        if (ddx * ddx + ddy * ddy < tol2) {
                            too_close = true;
                            close_v = static_cast<CDT::VertInd>(vj);
                            break;
                        }
                    }

                    if (too_close && close_v != sv1 && close_v != sv2) {
                        // Route through existing nearby vertex (no new Steiner)
                        constraint_edges[si] = CDT::Edge(sv1, close_v);
                        constraint_edges.push_back(CDT::Edge(close_v, sv2));
                    } else if (!too_close) {
                        CDT::VertInd mid = static_cast<CDT::VertInd>(vertices.size());
                        vertices.push_back(CDT::V2d<double>(mx, my));
                        constraint_edges[si] = CDT::Edge(sv1, mid);
                        constraint_edges.push_back(CDT::Edge(mid, sv2));
                        steiner_count++;
                    } else {
                        continue;  // too close to endpoint — try next vertex
                    }
                    found_encroached = true;
                    consecutive_skips = 0;
                    break;
                }
                if (found_encroached) break;
            }
            if (found_encroached) continue;  // rebuild CDT, re-check
        }

        // ---- Phase 2: Refine worst bad triangle ----
        int worst = -1;
        double worst_ma = M_PI;
        for (int t = 0; t < (int)tris.size(); t++) {
            auto q0 = tv[tris[t].vertices[0]];
            auto q1 = tv[tris[t].vertices[1]];
            auto q2 = tv[tris[t].vertices[2]];

            double l0 = std::sqrt((q1.x-q0.x)*(q1.x-q0.x) + (q1.y-q0.y)*(q1.y-q0.y));
            double l1 = std::sqrt((q2.x-q1.x)*(q2.x-q1.x) + (q2.y-q1.y)*(q2.y-q1.y));
            double l2 = std::sqrt((q0.x-q2.x)*(q0.x-q2.x) + (q0.y-q2.y)*(q0.y-q2.y));
            if (l0 < 1e-15 || l1 < 1e-15 || l2 < 1e-15) continue;

            double ma = min_angle_from_lengths(l0, l1, l2);
            double area = 0.5 * std::abs(
                (q1.x-q0.x)*(q2.y-q0.y) - (q2.x-q0.x)*(q1.y-q0.y));

            bool bad = (min_angle_deg > 0 && ma < min_angle_rad) ||
                       (max_area > 0 && area > max_area);
            if (bad && ma < worst_ma) { worst_ma = ma; worst = t; }
        }
        if (worst < 0) break;  // All triangles satisfy constraints

        auto p0 = tv[tris[worst].vertices[0]];
        auto p1 = tv[tris[worst].vertices[1]];
        auto p2 = tv[tris[worst].vertices[2]];

        // Compute off-center (Üngör 2009: closer to shortest edge)
        double ox, oy;
        offcenter(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, min_angle_rad, ox, oy);

        // Check if off-center encroaches any constraint edge
        if (has_boundary) {
            std::vector<size_t> encroached_segs;
            for (size_t e = 0; e < constraint_edges.size(); e++) {
                double eax = vertices[constraint_edges[e].v1()].x;
                double eay = vertices[constraint_edges[e].v1()].y;
                double ebx = vertices[constraint_edges[e].v2()].x;
                double eby = vertices[constraint_edges[e].v2()].y;
                if (point_encroaches_segment(ox, oy, eax, eay, ebx, eby)) {
                    encroached_segs.push_back(e);
                }
            }

            if (!encroached_segs.empty()) {
                // Split ALL encroached segments instead of inserting off-center.
                // Process in reverse index order to preserve earlier indices
                // when modifying constraint_edges in-place.
                std::sort(encroached_segs.rbegin(), encroached_segs.rend());
                bool any_split = false;
                for (size_t eidx : encroached_segs) {
                    if (steiner_count >= max_steiner) break;
                    CDT::VertInd ev1 = constraint_edges[eidx].v1();
                    CDT::VertInd ev2 = constraint_edges[eidx].v2();
                    double mx = (vertices[ev1].x + vertices[ev2].x) / 2.0;
                    double my = (vertices[ev1].y + vertices[ev2].y) / 2.0;

                    bool too_close = false;
                    CDT::VertInd close_v = 0;
                    for (size_t vj = 0; vj < vertices.size(); vj++) {
                        double ddx = mx - vertices[vj].x;
                        double ddy = my - vertices[vj].y;
                        if (ddx * ddx + ddy * ddy < tol2) {
                            too_close = true;
                            close_v = static_cast<CDT::VertInd>(vj);
                            break;
                        }
                    }

                    if (too_close && close_v != ev1 && close_v != ev2) {
                        constraint_edges[eidx] = CDT::Edge(ev1, close_v);
                        constraint_edges.push_back(CDT::Edge(close_v, ev2));
                        any_split = true;
                    } else if (!too_close) {
                        CDT::VertInd mid = static_cast<CDT::VertInd>(vertices.size());
                        vertices.push_back(CDT::V2d<double>(mx, my));
                        constraint_edges[eidx] = CDT::Edge(ev1, mid);
                        constraint_edges.push_back(CDT::Edge(mid, ev2));
                        steiner_count++;
                        any_split = true;
                    }
                }
                if (any_split) { consecutive_skips = 0; continue; }
                consecutive_skips++;
                continue;
            }

            // Off-center doesn't encroach any segment; check domain containment.
            if (!point_in_polygon(ox, oy, vertices, original_constraints)) {
                double gx = (p0.x + p1.x + p2.x) / 3.0;
                double gy = (p0.y + p1.y + p2.y) / 3.0;
                double lo = 0.0, hi = 1.0;
                for (int bs = 0; bs < 20; bs++) {
                    double mid_t = (lo + hi) / 2.0;
                    double tx = gx + mid_t * (ox - gx);
                    double ty = gy + mid_t * (oy - gy);
                    if (point_in_polygon(tx, ty, vertices, original_constraints))
                        lo = mid_t;
                    else
                        hi = mid_t;
                }
                if (lo > 0.01) {
                    ox = gx + lo * (ox - gx);
                    oy = gy + lo * (oy - gy);
                } else {
                    consecutive_skips++;
                    continue;
                }
            }
        }

        // Check proximity to existing vertices
        bool close = false;
        for (size_t vi = 0; vi < vertices.size(); vi++) {
            double ddx = ox - vertices[vi].x, ddy = oy - vertices[vi].y;
            if (ddx * ddx + ddy * ddy < tol2) { close = true; break; }
        }
        if (close) { consecutive_skips++; continue; }

        // Insert off-center
        vertices.push_back(CDT::V2d<double>(ox, oy));
        steiner_count++;
        consecutive_skips = 0;
    }

    // Final triangulation
    CDT::Triangulation<double> cdt_final(
        CDT::VertexInsertionOrder::Auto,
        CDT::IntersectingConstraintEdges::TryResolve,
        min_dist_tolerance
    );
    cdt_final.insertVertices(vertices);
    if (!constraint_edges.empty()) {
        cdt_final.insertEdges(constraint_edges);
    }
    cdt_final.eraseOuterTrianglesAndHoles();

    const auto& final_verts = cdt_final.vertices;
    const auto& final_tris = cdt_final.triangles;

    int n_verts = static_cast<int>(final_verts.size());
    int n_tri = static_cast<int>(final_tris.size());

    NumericMatrix out_vertices(n_verts, 2);
    for (int i = 0; i < n_verts; i++) {
        out_vertices(i, 0) = final_verts[i].x;
        out_vertices(i, 1) = final_verts[i].y;
    }

    IntegerMatrix out_triangles(n_tri, 3);
    for (int i = 0; i < n_tri; i++) {
        out_triangles(i, 0) = static_cast<int>(final_tris[i].vertices[0]) + 1;
        out_triangles(i, 1) = static_cast<int>(final_tris[i].vertices[1]) + 1;
        out_triangles(i, 2) = static_cast<int>(final_tris[i].vertices[2]) + 1;
    }

    std::set<std::pair<int,int>> edge_set;
    for (int i = 0; i < n_tri; i++) {
        for (int j = 0; j < 3; j++) {
            int v1 = out_triangles(i, j);
            int v2 = out_triangles(i, (j + 1) % 3);
            if (v1 > v2) std::swap(v1, v2);
            edge_set.insert({v1, v2});
        }
    }
    IntegerMatrix out_edges(edge_set.size(), 2);
    int idx = 0;
    for (const auto& e : edge_set) {
        out_edges(idx, 0) = e.first;
        out_edges(idx, 1) = e.second;
        idx++;
    }

    return List::create(
        Named("vertices") = out_vertices,
        Named("triangles") = out_triangles,
        Named("edges") = out_edges,
        Named("n_vertices") = n_verts,
        Named("n_triangles") = n_tri,
        Named("n_edges") = static_cast<int>(edge_set.size()),
        Named("n_input_points") = n_pts,
        Named("n_steiner") = steiner_count
    );
}

// Compute FEM matrices (mass C and stiffness G) from a triangulation.
// Uses linear (P1) basis functions on triangles.
// Returns sparse matrices in triplet format (i, j, x) for conversion to dgCMatrix.

// [[Rcpp::export]]
List cpp_fem_matrices(NumericMatrix vertices, IntegerMatrix triangles) {
    int n_verts = vertices.nrow();
    int n_tri = triangles.nrow();

    // Triplet accumulators for C (mass) and G (stiffness)
    std::vector<int> C_i, C_j, G_i, G_j;
    std::vector<double> C_x, G_x;

    for (int t = 0; t < n_tri; t++) {
        // Triangle vertex indices (0-based)
        int v0 = triangles(t, 0) - 1;
        int v1 = triangles(t, 1) - 1;
        int v2 = triangles(t, 2) - 1;

        // Vertex coordinates
        double x0 = vertices(v0, 0), y0 = vertices(v0, 1);
        double x1 = vertices(v1, 0), y1 = vertices(v1, 1);
        double x2 = vertices(v2, 0), y2 = vertices(v2, 1);

        // Triangle area (signed, take absolute)
        double area = 0.5 * std::abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
        if (area < 1e-15) continue;  // degenerate triangle

        // --- Mass matrix C (lumped diagonal: area/3 per vertex) ---
        // Consistent mass: C_ij = area/12 for i!=j, area/6 for i==j
        // We use consistent mass (better accuracy for SPDE)
        int tri_verts[3] = {v0, v1, v2};
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double val = (i == j) ? area / 6.0 : area / 12.0;
                C_i.push_back(tri_verts[i]);
                C_j.push_back(tri_verts[j]);
                C_x.push_back(val);
            }
        }

        // --- Stiffness matrix G ---
        // G_ij = integral(grad_phi_i . grad_phi_j) over triangle
        // For linear elements: grad_phi = (1/2A) * [y_j - y_k, x_k - x_j]
        // where (i, j, k) are cyclic permutations

        // Edge vectors (opposite to each vertex)
        double dx[3], dy[3];
        dx[0] = x2 - x1; dy[0] = y2 - y1;  // opposite v0
        dx[1] = x0 - x2; dy[1] = y0 - y2;  // opposite v1
        dx[2] = x1 - x0; dy[2] = y1 - y0;  // opposite v2

        // G_ij = (1 / (4 * area)) * (dx_i * dx_j + dy_i * dy_j)
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double val = (dx[i] * dx[j] + dy[i] * dy[j]) / (4.0 * area);
                G_i.push_back(tri_verts[i]);
                G_j.push_back(tri_verts[j]);
                G_x.push_back(val);
            }
        }
    }

    // Return triplets (R will assemble into sparse matrix)
    return List::create(
        Named("C_i") = wrap(C_i),
        Named("C_j") = wrap(C_j),
        Named("C_x") = wrap(C_x),
        Named("G_i") = wrap(G_i),
        Named("G_j") = wrap(G_j),
        Named("G_x") = wrap(G_x),
        Named("n") = n_verts
    );
}

// ---------------------------------------------------------------------------
// Parallel FEM assembly using RcppParallel.
// Each thread accumulates triplets for its chunk of triangles,
// then results are merged. Significant speedup for >50K triangles.
// ---------------------------------------------------------------------------

struct FemWorker : public RcppParallel::Worker {
    const RcppParallel::RMatrix<double> vertices;
    const RcppParallel::RMatrix<int> triangles;

    // Thread-local storage: one vector per thread
    // We use a vector of vectors indexed by thread
    std::vector<std::vector<int>>    tl_C_i, tl_C_j, tl_G_i, tl_G_j;
    std::vector<std::vector<double>> tl_C_x, tl_G_x;
    int n_threads;

    FemWorker(NumericMatrix verts, IntegerMatrix tris)
        : vertices(verts), triangles(tris)
    {
        n_threads = std::max(1, (int)tbb::this_task_arena::max_concurrency());
        tl_C_i.resize(n_threads); tl_C_j.resize(n_threads); tl_C_x.resize(n_threads);
        tl_G_i.resize(n_threads); tl_G_j.resize(n_threads); tl_G_x.resize(n_threads);
    }

    void operator()(std::size_t begin, std::size_t end) {
        int tid = tbb::this_task_arena::current_thread_index();
        if (tid < 0 || tid >= n_threads) tid = 0;

        auto& C_i = tl_C_i[tid]; auto& C_j = tl_C_j[tid]; auto& C_x = tl_C_x[tid];
        auto& G_i = tl_G_i[tid]; auto& G_j = tl_G_j[tid]; auto& G_x = tl_G_x[tid];

        for (std::size_t t = begin; t < end; t++) {
            int v0 = triangles(t, 0) - 1;
            int v1 = triangles(t, 1) - 1;
            int v2 = triangles(t, 2) - 1;

            double x0 = vertices(v0, 0), y0 = vertices(v0, 1);
            double x1 = vertices(v1, 0), y1 = vertices(v1, 1);
            double x2 = vertices(v2, 0), y2 = vertices(v2, 1);

            double area = 0.5 * std::abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
            if (area < 1e-15) continue;

            int tri_verts[3] = {v0, v1, v2};

            // Mass matrix
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    double val = (i == j) ? area / 6.0 : area / 12.0;
                    C_i.push_back(tri_verts[i]);
                    C_j.push_back(tri_verts[j]);
                    C_x.push_back(val);
                }
            }

            // Stiffness matrix
            double dx[3], dy[3];
            dx[0] = x2 - x1; dy[0] = y2 - y1;
            dx[1] = x0 - x2; dy[1] = y0 - y2;
            dx[2] = x1 - x0; dy[2] = y1 - y0;

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    double val = (dx[i] * dx[j] + dy[i] * dy[j]) / (4.0 * area);
                    G_i.push_back(tri_verts[i]);
                    G_j.push_back(tri_verts[j]);
                    G_x.push_back(val);
                }
            }
        }
    }
};

// [[Rcpp::export]]
List cpp_fem_matrices_parallel(NumericMatrix vertices, IntegerMatrix triangles) {
    int n_verts = vertices.nrow();
    int n_tri = triangles.nrow();

    FemWorker worker(vertices, triangles);
    RcppParallel::parallelFor(0, n_tri, worker);

    // Merge thread-local results
    std::size_t total_C = 0, total_G = 0;
    for (int t = 0; t < worker.n_threads; t++) {
        total_C += worker.tl_C_i[t].size();
        total_G += worker.tl_G_i[t].size();
    }

    std::vector<int> C_i, C_j, G_i, G_j;
    std::vector<double> C_x, G_x;
    C_i.reserve(total_C); C_j.reserve(total_C); C_x.reserve(total_C);
    G_i.reserve(total_G); G_j.reserve(total_G); G_x.reserve(total_G);

    for (int t = 0; t < worker.n_threads; t++) {
        C_i.insert(C_i.end(), worker.tl_C_i[t].begin(), worker.tl_C_i[t].end());
        C_j.insert(C_j.end(), worker.tl_C_j[t].begin(), worker.tl_C_j[t].end());
        C_x.insert(C_x.end(), worker.tl_C_x[t].begin(), worker.tl_C_x[t].end());
        G_i.insert(G_i.end(), worker.tl_G_i[t].begin(), worker.tl_G_i[t].end());
        G_j.insert(G_j.end(), worker.tl_G_j[t].begin(), worker.tl_G_j[t].end());
        G_x.insert(G_x.end(), worker.tl_G_x[t].begin(), worker.tl_G_x[t].end());
    }

    return List::create(
        Named("C_i") = wrap(C_i),
        Named("C_j") = wrap(C_j),
        Named("C_x") = wrap(C_x),
        Named("G_i") = wrap(G_i),
        Named("G_j") = wrap(G_j),
        Named("G_x") = wrap(G_x),
        Named("n") = n_verts
    );
}

// ---------------------------------------------------------------------------
// P2 (quadratic) FEM matrices.
// 6-node triangles: 3 vertex nodes + 3 edge midpoint nodes.
// Basis functions are quadratic: phi_i = L_i(2L_i - 1) for vertices,
// phi_ij = 4 L_i L_j for midpoints, where L_i are barycentric coords.
//
// The R layer handles midpoint node creation and index mapping.
// This function takes the full vertex set (original + midpoints) and
// the 6-column triangle connectivity.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List cpp_fem_matrices_p2(
    NumericMatrix vertices,
    IntegerMatrix triangles6  // N_tri x 6: [v0, v1, v2, m01, m12, m20]
) {
    int n_verts = vertices.nrow();
    int n_tri = triangles6.nrow();

    std::vector<int> C_i, C_j, G_i, G_j;
    std::vector<double> C_x, G_x;

    // Reference element mass matrix for P2 on the unit triangle.
    // M_ref[i][j] = integral(phi_i * phi_j) over ref triangle (area=0.5)
    // Scaled by 2 (so multiply by area, not area/2, at the end)
    // Symmetric 6x6 matrix, pre-computed from exact quadrature.
    // Nodes: 0=(0,0), 1=(1,0), 2=(0,1), 3=mid(0,1), 4=mid(1,2), 5=mid(2,0)
    static const double M_ref[6][6] = {
        { 6, -1, -1,  0, -4,  0},
        {-1,  6, -1,  0,  0, -4},
        {-1, -1,  6, -4,  0,  0},
        { 0,  0, -4, 32, 16, 16},
        {-4,  0,  0, 16, 32, 16},
        { 0, -4,  0, 16, 16, 32}
    };
    // Scale: multiply by area / 180

    for (int t = 0; t < n_tri; t++) {
        // Vertex nodes (0-based)
        int nodes[6];
        for (int k = 0; k < 6; k++) nodes[k] = triangles6(t, k) - 1;

        int v0 = nodes[0], v1 = nodes[1], v2 = nodes[2];

        double x0 = vertices(v0, 0), y0 = vertices(v0, 1);
        double x1 = vertices(v1, 0), y1 = vertices(v1, 1);
        double x2 = vertices(v2, 0), y2 = vertices(v2, 1);

        double area = 0.5 * std::abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
        if (area < 1e-15) continue;

        // --- Mass matrix ---
        double mass_scale = area / 180.0;
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                double val = M_ref[i][j] * mass_scale;
                if (std::abs(val) > 1e-20) {
                    C_i.push_back(nodes[i]);
                    C_j.push_back(nodes[j]);
                    C_x.push_back(val);
                }
            }
        }

        // --- Stiffness matrix ---
        // For P2, the stiffness integral is:
        //   G_ij = integral(grad_phi_i . grad_phi_j) dA
        // Using 3-point Gauss quadrature on the reference triangle.
        // Quadrature points (barycentric): (2/3, 1/6, 1/6), (1/6, 2/3, 1/6), (1/6, 1/6, 2/3)
        // Weight: 1/3 each (times area)
        double L_qp[3][3] = {
            {2.0/3.0, 1.0/6.0, 1.0/6.0},
            {1.0/6.0, 2.0/3.0, 1.0/6.0},
            {1.0/6.0, 1.0/6.0, 2.0/3.0}
        };
        double w_qp = 1.0 / 3.0;  // weight per quadrature point (area factored separately)

        // Jacobian of the mapping from reference to physical triangle
        double J11 = x1 - x0, J12 = x2 - x0;
        double J21 = y1 - y0, J22 = y2 - y0;
        double detJ = J11 * J22 - J12 * J21;
        if (std::abs(detJ) < 1e-15) continue;

        // Inverse transpose of Jacobian (for gradient transformation)
        double invJt11 =  J22 / detJ, invJt12 = -J21 / detJ;
        double invJt21 = -J12 / detJ, invJt22 =  J11 / detJ;

        // Accumulate stiffness via quadrature
        double G_local[6][6] = {};

        for (int q = 0; q < 3; q++) {
            double L0 = L_qp[q][0], L1 = L_qp[q][1], L2 = L_qp[q][2];

            // Derivatives of P2 basis functions w.r.t. barycentric coords
            // phi0 = L0(2L0-1), phi1 = L1(2L1-1), phi2 = L2(2L2-1)
            // phi3 = 4*L0*L1,   phi4 = 4*L1*L2,   phi5 = 4*L2*L0
            //
            // d/dL0: dphi0 = 4L0-1, dphi3 = 4L1, dphi5 = 4L2
            // d/dL1: dphi1 = 4L1-1, dphi3 = 4L0, dphi4 = 4L2
            // d/dL2: dphi2 = 4L2-1, dphi4 = 4L1, dphi5 = 4L0
            //
            // Reference coords: xi = L1, eta = L2, L0 = 1 - xi - eta
            // d/dxi = d/dL1 - d/dL0
            // d/deta = d/dL2 - d/dL0

            double dL0 = 4.0 * L0 - 1.0;
            double dL1 = 4.0 * L1 - 1.0;
            double dL2 = 4.0 * L2 - 1.0;

            // d(phi)/d(xi) and d(phi)/d(eta) for each of 6 basis functions
            double dphi_dxi[6], dphi_deta[6];

            // phi0 = L0(2L0-1): d/dxi = -dL0, d/deta = -dL0
            dphi_dxi[0]  = -dL0;    dphi_deta[0] = -dL0;
            // phi1 = L1(2L1-1): d/dxi = dL1,  d/deta = 0
            dphi_dxi[1]  =  dL1;    dphi_deta[1] = 0.0;
            // phi2 = L2(2L2-1): d/dxi = 0,    d/deta = dL2
            dphi_dxi[2]  = 0.0;     dphi_deta[2] = dL2;
            // phi3 = 4*L0*L1:   d/dxi = 4(L0 - L1) = 4L0-4L1, d/deta = -4L1
            dphi_dxi[3]  = 4.0*(L0 - L1); dphi_deta[3] = -4.0*L1;
            // phi4 = 4*L1*L2:   d/dxi = 4L2,  d/deta = 4L1
            dphi_dxi[4]  = 4.0*L2;  dphi_deta[4] = 4.0*L1;
            // phi5 = 4*L2*L0:   d/dxi = -4L2, d/deta = 4(L0 - L2)
            dphi_dxi[5]  = -4.0*L2; dphi_deta[5] = 4.0*(L0 - L2);

            // Transform to physical gradients: grad_phi = invJ^T * [dphi/dxi, dphi/deta]
            double grad_x[6], grad_y[6];
            for (int k = 0; k < 6; k++) {
                grad_x[k] = invJt11 * dphi_dxi[k] + invJt21 * dphi_deta[k];
                grad_y[k] = invJt12 * dphi_dxi[k] + invJt22 * dphi_deta[k];
            }

            // Accumulate: G_ij += w * (grad_x_i*grad_x_j + grad_y_i*grad_y_j) * |detJ|
            double wdetJ = w_qp * std::abs(detJ);
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 6; j++) {
                    G_local[i][j] += wdetJ * (grad_x[i]*grad_x[j] + grad_y[i]*grad_y[j]);
                }
            }
        }

        // Push stiffness triplets
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                if (std::abs(G_local[i][j]) > 1e-20) {
                    G_i.push_back(nodes[i]);
                    G_j.push_back(nodes[j]);
                    G_x.push_back(G_local[i][j]);
                }
            }
        }
    }

    return List::create(
        Named("C_i") = wrap(C_i),
        Named("C_j") = wrap(C_j),
        Named("C_x") = wrap(C_x),
        Named("G_i") = wrap(G_i),
        Named("G_j") = wrap(G_j),
        Named("G_x") = wrap(G_x),
        Named("n") = n_verts
    );
}

// ---------------------------------------------------------------------------
// Non-stationary FEM: spatially varying kappa and tau.
// Computes weighted mass and stiffness matrices where the weights
// are per-vertex kappa(s) and tau(s) values interpolated within triangles.
//
// Returns three matrices in triplet form:
//   Ck: kappa-weighted mass   integral(kappa(s)^2 * phi_i * phi_j)
//   Gk: kappa-weighted stiffness  integral(kappa(s)^2 * grad_phi_i . grad_phi_j)
//   Ct: tau-weighted mass     integral(tau(s)^2 * phi_i * phi_j)
//
// These are used to build the non-stationary SPDE precision:
//   Q = Ct * (Ck^2 + 2*Gk + Gk * C0^{-1} * Gk)
// where the kappa and tau fields encode spatially varying range and variance.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List cpp_fem_matrices_nonstationary(
    NumericMatrix vertices,
    IntegerMatrix triangles,
    NumericVector kappa,
    NumericVector tau
) {
    int n_verts = vertices.nrow();
    int n_tri = triangles.nrow();

    // Triplet accumulators
    std::vector<int> Ck_i, Ck_j, Gk_i, Gk_j, Ct_i, Ct_j;
    std::vector<double> Ck_x, Gk_x, Ct_x;

    for (int t = 0; t < n_tri; t++) {
        int v0 = triangles(t, 0) - 1;
        int v1 = triangles(t, 1) - 1;
        int v2 = triangles(t, 2) - 1;

        double x0 = vertices(v0, 0), y0 = vertices(v0, 1);
        double x1 = vertices(v1, 0), y1 = vertices(v1, 1);
        double x2 = vertices(v2, 0), y2 = vertices(v2, 1);

        double area = 0.5 * std::abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
        if (area < 1e-15) continue;

        int tri_verts[3] = {v0, v1, v2};

        // Per-vertex kappa^2 and tau^2
        double k2[3] = {kappa[v0] * kappa[v0], kappa[v1] * kappa[v1], kappa[v2] * kappa[v2]};
        double t2[3] = {tau[v0] * tau[v0], tau[v1] * tau[v1], tau[v2] * tau[v2]};

        // Weighted mass: integral(w(s) * phi_i * phi_j) over triangle
        // Using exact quadrature for P1 * P1 * P1 on a triangle:
        //   integral(phi_i * phi_j * phi_k) = area * (1 + delta_ij + delta_ik + delta_jk) / 60
        // But for w = sum_k w_k phi_k, the integral becomes:
        //   integral(w * phi_i * phi_j) = area * sum_k w_k * I(i,j,k)
        // where I(i,j,k) = area * (1/60) * [2 if all different, 6 if two same, 12 if all same]
        // Simplified: for the symmetric 3x3 element matrix:
        //   M_ij = area/60 * (w_i*(1+delta_ij) + w_j*(1+delta_ij) + w_k*(1-delta_ij))
        // Actually, the exact formula for integral(w * phi_i * phi_j) with linear w:
        //   = area/12 * (w_i + w_j) for i != j
        //   = area/6 * w_i             for i == j  (approximation)
        // More precisely, for w = w0*phi0 + w1*phi1 + w2*phi2:
        //   integral(w * phi_i * phi_i) = area * (2*w_i + w_j + w_k) / 12
        //   integral(w * phi_i * phi_j) = area * (w_i + w_j + 2*w_k) / 24  [where k is the third vertex]
        // But we use a simpler approximation with triangle-average weight:
        //   w_avg = (w0 + w1 + w2) / 3
        //   M_ij = w_avg * (standard mass element)

        double k2_avg = (k2[0] + k2[1] + k2[2]) / 3.0;
        double t2_avg = (t2[0] + t2[1] + t2[2]) / 3.0;

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double base_mass = (i == j) ? area / 6.0 : area / 12.0;

                // Kappa-weighted mass
                Ck_i.push_back(tri_verts[i]);
                Ck_j.push_back(tri_verts[j]);
                Ck_x.push_back(k2_avg * base_mass);

                // Tau-weighted mass
                Ct_i.push_back(tri_verts[i]);
                Ct_j.push_back(tri_verts[j]);
                Ct_x.push_back(t2_avg * base_mass);
            }
        }

        // Kappa-weighted stiffness: integral(kappa^2 * grad_phi_i . grad_phi_j)
        // Use triangle-average kappa^2
        double dx[3], dy[3];
        dx[0] = x2 - x1; dy[0] = y2 - y1;
        dx[1] = x0 - x2; dy[1] = y0 - y2;
        dx[2] = x1 - x0; dy[2] = y1 - y0;

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double val = k2_avg * (dx[i] * dx[j] + dy[i] * dy[j]) / (4.0 * area);
                Gk_i.push_back(tri_verts[i]);
                Gk_j.push_back(tri_verts[j]);
                Gk_x.push_back(val);
            }
        }
    }

    return List::create(
        Named("Ck_i") = wrap(Ck_i), Named("Ck_j") = wrap(Ck_j), Named("Ck_x") = wrap(Ck_x),
        Named("Gk_i") = wrap(Gk_i), Named("Gk_j") = wrap(Gk_j), Named("Gk_x") = wrap(Gk_x),
        Named("Ct_i") = wrap(Ct_i), Named("Ct_j") = wrap(Ct_j), Named("Ct_x") = wrap(Ct_x),
        Named("n") = n_verts
    );
}

// ---------------------------------------------------------------------------
// Barycentric point-in-triangle test (shared helper)
// ---------------------------------------------------------------------------
struct BaryResult {
    bool inside;
    double l0, l1, l2;
};

static inline BaryResult barycentric(
    double px, double py,
    double x0, double y0,
    double x1, double y1,
    double x2, double y2
) {
    BaryResult r{false, 0.0, 0.0, 0.0};
    double denom = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2);
    if (std::abs(denom) < 1e-15) return r;

    r.l0 = ((y1 - y2) * (px - x2) + (x2 - x1) * (py - y2)) / denom;
    r.l1 = ((y2 - y0) * (px - x2) + (x0 - x2) * (py - y2)) / denom;
    r.l2 = 1.0 - r.l0 - r.l1;

    constexpr double tol = -1e-8;
    if (r.l0 >= tol && r.l1 >= tol && r.l2 >= tol) {
        r.l0 = std::max(0.0, r.l0);
        r.l1 = std::max(0.0, r.l1);
        r.l2 = std::max(0.0, r.l2);
        double sum = r.l0 + r.l1 + r.l2;
        r.l0 /= sum; r.l1 /= sum; r.l2 /= sum;
        r.inside = true;
    }
    return r;
}

// ---------------------------------------------------------------------------
// Triangle adjacency + walk-based point location.
//
// Strategy:
//   1. Build triangle adjacency: for each triangle, store its 3 neighbors
//      (or -1 if boundary edge).
//   2. Build a grid of triangle centroids for O(1) initial guess.
//   3. For each query point: grid lookup → walk from initial triangle
//      toward target by crossing the edge with the most negative
//      barycentric coordinate. O(sqrt(n)) expected walk length.
//   4. Remembered-edge: start the next query from the last found triangle
//      (huge win for spatially clustered observation points).
// ---------------------------------------------------------------------------

struct TriLocator {
    int n_tri, n_verts;
    NumericMatrix verts;
    IntegerMatrix tris;

    // Adjacency: adj[t*3 + k] = neighbor triangle across edge k, or -1
    std::vector<int> adj;

    // Grid for initial guess
    double x_min, y_min, cell_w, cell_h;
    int nx, ny;
    std::vector<int> grid_tri;  // one representative triangle per cell

    // Remembered last triangle for sequential queries
    mutable int last_tri;

    TriLocator(NumericMatrix vertices, IntegerMatrix triangles)
        : n_tri(triangles.nrow()), n_verts(vertices.nrow()),
          verts(vertices), tris(triangles), last_tri(0)
    {
        build_adjacency();
        build_grid();
    }

    void build_adjacency() {
        adj.assign(n_tri * 3, -1);

        // Map each directed edge (v1, v2) to (triangle, local_edge_index)
        // Edge k of triangle t is from tris[t, k] to tris[t, (k+1)%3]
        std::unordered_map<int64_t, std::pair<int,int>> edge_map;

        for (int t = 0; t < n_tri; t++) {
            for (int k = 0; k < 3; k++) {
                int v1 = tris(t, k);
                int v2 = tris(t, (k + 1) % 3);
                // Canonical edge key: smaller vertex first
                int lo = std::min(v1, v2);
                int hi = std::max(v1, v2);
                int64_t key = static_cast<int64_t>(lo) * (n_verts + 1) + hi;

                auto it = edge_map.find(key);
                if (it != edge_map.end()) {
                    int other_t = it->second.first;
                    int other_k = it->second.second;
                    adj[t * 3 + k] = other_t;
                    adj[other_t * 3 + other_k] = t;
                    edge_map.erase(it);
                } else {
                    edge_map[key] = {t, k};
                }
            }
        }
    }

    void build_grid() {
        // Bounding box
        x_min = verts(0, 0); double x_max = x_min;
        y_min = verts(0, 1); double y_max = y_min;
        for (int i = 1; i < n_verts; i++) {
            double x = verts(i, 0), y = verts(i, 1);
            if (x < x_min) x_min = x;
            if (x > x_max) x_max = x;
            if (y < y_min) y_min = y;
            if (y > y_max) y_max = y;
        }
        double pad = 1e-10;
        x_min -= pad; x_max += pad;
        y_min -= pad; y_max += pad;

        nx = std::max(1, std::min(1000, (int)std::sqrt((double)n_tri)));
        ny = nx;
        cell_w = (x_max - x_min) / nx;
        cell_h = (y_max - y_min) / ny;

        // For each cell, store the triangle whose centroid is closest
        grid_tri.assign(nx * ny, 0);

        // Place each triangle's centroid into its grid cell
        for (int t = 0; t < n_tri; t++) {
            int v0 = tris(t, 0) - 1, v1 = tris(t, 1) - 1, v2 = tris(t, 2) - 1;
            double cx = (verts(v0, 0) + verts(v1, 0) + verts(v2, 0)) / 3.0;
            double cy = (verts(v0, 1) + verts(v1, 1) + verts(v2, 1)) / 3.0;
            int ix = std::max(0, std::min(nx - 1, (int)((cx - x_min) / cell_w)));
            int iy = std::max(0, std::min(ny - 1, (int)((cy - y_min) / cell_h)));
            grid_tri[iy * nx + ix] = t;
        }

        // Fill empty cells with nearest non-empty cell's triangle
        // (simple flood-fill would be better, but this works)
        for (int iy = 0; iy < ny; iy++) {
            for (int ix = 0; ix < nx; ix++) {
                if (grid_tri[iy * nx + ix] == 0 && n_tri > 1) {
                    // Find nearest non-default cell
                    double best_d = 1e30;
                    for (int t = 0; t < n_tri; t++) {
                        int v0 = tris(t, 0) - 1, v1 = tris(t, 1) - 1, v2 = tris(t, 2) - 1;
                        double cx = (verts(v0, 0) + verts(v1, 0) + verts(v2, 0)) / 3.0;
                        double cy = (verts(v0, 1) + verts(v1, 1) + verts(v2, 1)) / 3.0;
                        double gcx = x_min + (ix + 0.5) * cell_w;
                        double gcy = y_min + (iy + 0.5) * cell_h;
                        double d = (cx - gcx) * (cx - gcx) + (cy - gcy) * (cy - gcy);
                        if (d < best_d) { best_d = d; grid_tri[iy * nx + ix] = t; }
                    }
                }
            }
        }
    }

    int initial_guess(double px, double py) const {
        int ix = std::max(0, std::min(nx - 1, (int)((px - x_min) / cell_w)));
        int iy = std::max(0, std::min(ny - 1, (int)((py - y_min) / cell_h)));
        return grid_tri[iy * nx + ix];
    }

    // Walk from triangle `start` toward point (px, py).
    // Returns triangle index if found, -1 if walked off the mesh.
    int walk(double px, double py, int start) const {
        int current = start;
        int max_steps = n_tri;  // safety bound

        for (int step = 0; step < max_steps; step++) {
            int v0 = tris(current, 0) - 1;
            int v1 = tris(current, 1) - 1;
            int v2 = tris(current, 2) - 1;

            auto br = barycentric(px, py,
                verts(v0, 0), verts(v0, 1),
                verts(v1, 0), verts(v1, 1),
                verts(v2, 0), verts(v2, 1));

            if (br.inside) {
                return current;
            }

            // Walk across the edge with the most negative barycentric coord
            double bary[3] = {br.l0, br.l1, br.l2};
            int worst = 0;
            if (bary[1] < bary[worst]) worst = 1;
            if (bary[2] < bary[worst]) worst = 2;

            // Edge opposite to vertex `worst` is edge index (worst+1)%3
            // Edge k connects vertex k and vertex (k+1)%3
            // So edge opposite vertex 0 is edge 1 (v1→v2), etc.
            int edge_idx = (worst + 1) % 3;
            int neighbor = adj[current * 3 + edge_idx];
            if (neighbor < 0) {
                // Hit boundary — point is outside mesh
                return -1;
            }
            current = neighbor;
        }
        return -1;  // didn't converge
    }
};

// Compute projection matrix A: maps mesh vertices to observation locations.
// Uses jump-and-walk point location with remembered-edge optimization.

// [[Rcpp::export]]
List cpp_projection_matrix(
    NumericMatrix obs_coords,
    NumericMatrix vertices,
    IntegerMatrix triangles
) {
    int n_obs = obs_coords.nrow();
    int n_verts = vertices.nrow();

    TriLocator locator(vertices, triangles);

    std::vector<int> A_i, A_j;
    std::vector<double> A_x;

    for (int obs = 0; obs < n_obs; obs++) {
        double px = obs_coords(obs, 0);
        double py = obs_coords(obs, 1);

        bool found = false;

        // Try 1: walk from last found triangle (remembered-edge)
        int t = locator.walk(px, py, locator.last_tri);

        // Try 2: walk from grid initial guess
        if (t < 0) {
            int guess = locator.initial_guess(px, py);
            t = locator.walk(px, py, guess);
        }

        if (t >= 0) {
            int v0 = triangles(t, 0) - 1;
            int v1 = triangles(t, 1) - 1;
            int v2 = triangles(t, 2) - 1;

            auto br = barycentric(px, py,
                vertices(v0, 0), vertices(v0, 1),
                vertices(v1, 0), vertices(v1, 1),
                vertices(v2, 0), vertices(v2, 1));

            if (br.l0 > 1e-15) { A_i.push_back(obs); A_j.push_back(v0); A_x.push_back(br.l0); }
            if (br.l1 > 1e-15) { A_i.push_back(obs); A_j.push_back(v1); A_x.push_back(br.l1); }
            if (br.l2 > 1e-15) { A_i.push_back(obs); A_j.push_back(v2); A_x.push_back(br.l2); }

            locator.last_tri = t;
            found = true;
        }

        if (!found) {
            // Point outside mesh: project to nearest vertex
            double min_dist = 1e30;
            int nearest = 0;
            for (int v = 0; v < n_verts; v++) {
                double dx = px - vertices(v, 0);
                double dy = py - vertices(v, 1);
                double d = dx * dx + dy * dy;
                if (d < min_dist) { min_dist = d; nearest = v; }
            }
            A_i.push_back(obs);
            A_j.push_back(nearest);
            A_x.push_back(1.0);
        }
    }

    return List::create(
        Named("i") = wrap(A_i),
        Named("j") = wrap(A_j),
        Named("x") = wrap(A_x),
        Named("nrow") = n_obs,
        Named("ncol") = n_verts
    );
}

// ---------------------------------------------------------------------------
// FEM matrices for 3D surface meshes (spherical meshes).
// Uses cross-product for triangle area and local tangent-plane gradients.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List cpp_fem_matrices_3d(NumericMatrix vertices, IntegerMatrix triangles) {
    int n_verts = vertices.nrow();
    int n_tri = triangles.nrow();

    std::vector<int> C_i, C_j, G_i, G_j;
    std::vector<double> C_x, G_x;

    for (int t = 0; t < n_tri; t++) {
        int v0 = triangles(t, 0) - 1;
        int v1 = triangles(t, 1) - 1;
        int v2 = triangles(t, 2) - 1;

        double x0 = vertices(v0, 0), y0 = vertices(v0, 1), z0 = vertices(v0, 2);
        double x1 = vertices(v1, 0), y1 = vertices(v1, 1), z1 = vertices(v1, 2);
        double x2 = vertices(v2, 0), y2 = vertices(v2, 1), z2 = vertices(v2, 2);

        // Edge vectors
        double e1x = x1 - x0, e1y = y1 - y0, e1z = z1 - z0;
        double e2x = x2 - x0, e2y = y2 - y0, e2z = z2 - z0;

        // Cross product = normal * 2 * area
        double nx = e1y * e2z - e1z * e2y;
        double ny = e1z * e2x - e1x * e2z;
        double nz = e1x * e2y - e1y * e2x;

        double area2 = std::sqrt(nx * nx + ny * ny + nz * nz);
        double area = area2 / 2.0;
        if (area < 1e-15) continue;

        // --- Mass matrix C (consistent mass) ---
        int tri_verts[3] = {v0, v1, v2};
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double val = (i == j) ? area / 6.0 : area / 12.0;
                C_i.push_back(tri_verts[i]);
                C_j.push_back(tri_verts[j]);
                C_x.push_back(val);
            }
        }

        // --- Stiffness matrix G (tangent-plane gradients) ---
        // Use local 2D parameterization via edge vectors e1, e2.
        // The metric tensor is:
        //   g11 = e1 . e1,  g12 = e1 . e2,  g22 = e2 . e2
        // The inverse metric times |det g|^{1/2} gives the stiffness.
        double g11 = e1x * e1x + e1y * e1y + e1z * e1z;
        double g12 = e1x * e2x + e1y * e2y + e1z * e2z;
        double g22 = e2x * e2x + e2y * e2y + e2z * e2z;
        double det_g = g11 * g22 - g12 * g12;
        if (det_g < 1e-30) continue;

        // Inverse metric
        double inv_g11 = g22 / det_g;
        double inv_g12 = -g12 / det_g;
        double inv_g22 = g11 / det_g;

        // Gradients of barycentric coordinates in (u, v) parametric space:
        // phi0(u,v) = 1 - u - v  → grad = (-1, -1)
        // phi1(u,v) = u          → grad = (1, 0)
        // phi2(u,v) = v          → grad = (0, 1)
        double du[3] = {-1.0,  1.0, 0.0};
        double dv[3] = {-1.0,  0.0, 1.0};

        // G_ij = integral(grad_phi_i . grad_phi_j) dS
        //      = (grad_phi_i^T * g^{-1} * grad_phi_j) * sqrt(det_g) / 2
        double sqrt_det_g = std::sqrt(det_g);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double dot = (du[i] * inv_g11 + dv[i] * inv_g12) * du[j] +
                             (du[i] * inv_g12 + dv[i] * inv_g22) * dv[j];
                double val = dot * sqrt_det_g / 2.0;
                G_i.push_back(tri_verts[i]);
                G_j.push_back(tri_verts[j]);
                G_x.push_back(val);
            }
        }
    }

    return List::create(
        Named("C_i") = wrap(C_i),
        Named("C_j") = wrap(C_j),
        Named("C_x") = wrap(C_x),
        Named("G_i") = wrap(G_i),
        Named("G_j") = wrap(G_j),
        Named("G_x") = wrap(G_x),
        Named("n") = n_verts
    );
}

// ---------------------------------------------------------------------------
// Projection matrix for 3D surface meshes (spherical meshes).
// Uses 3D barycentric coordinates within each triangle's plane.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List cpp_projection_matrix_3d(
    NumericMatrix obs_coords,
    NumericMatrix vertices,
    IntegerMatrix triangles
) {
    int n_obs = obs_coords.nrow();
    int n_tri = triangles.nrow();
    int n_verts = vertices.nrow();

    std::vector<int> A_i, A_j;
    std::vector<double> A_x;

    for (int obs = 0; obs < n_obs; obs++) {
        double px = obs_coords(obs, 0);
        double py = obs_coords(obs, 1);
        double pz = obs_coords(obs, 2);

        bool found = false;
        double best_dist = 1e30;
        double best_l0 = 0, best_l1 = 0, best_l2 = 0;
        int best_v0 = 0, best_v1 = 0, best_v2 = 0;

        for (int t = 0; t < n_tri; t++) {
            int v0 = triangles(t, 0) - 1;
            int v1 = triangles(t, 1) - 1;
            int v2 = triangles(t, 2) - 1;

            double x0 = vertices(v0, 0), y0 = vertices(v0, 1), z0 = vertices(v0, 2);
            double x1 = vertices(v1, 0), y1 = vertices(v1, 1), z1 = vertices(v1, 2);
            double x2 = vertices(v2, 0), y2 = vertices(v2, 1), z2 = vertices(v2, 2);

            double e1x = x1 - x0, e1y = y1 - y0, e1z = z1 - z0;
            double e2x = x2 - x0, e2y = y2 - y0, e2z = z2 - z0;
            double epx = px - x0, epy = py - y0, epz = pz - z0;

            double d11 = e1x * e1x + e1y * e1y + e1z * e1z;
            double d12 = e1x * e2x + e1y * e2y + e1z * e2z;
            double d22 = e2x * e2x + e2y * e2y + e2z * e2z;
            double dp1 = epx * e1x + epy * e1y + epz * e1z;
            double dp2 = epx * e2x + epy * e2y + epz * e2z;

            double denom = d11 * d22 - d12 * d12;
            if (std::abs(denom) < 1e-30) continue;

            double l1 = (d22 * dp1 - d12 * dp2) / denom;
            double l2 = (d11 * dp2 - d12 * dp1) / denom;
            double l0 = 1.0 - l1 - l2;

            constexpr double tol = -1e-6;
            if (l0 >= tol && l1 >= tol && l2 >= tol) {
                l0 = std::max(0.0, l0);
                l1 = std::max(0.0, l1);
                l2 = std::max(0.0, l2);
                double sum = l0 + l1 + l2;
                l0 /= sum; l1 /= sum; l2 /= sum;

                double proj_x = x0 + l1 * e1x + l2 * e2x;
                double proj_y = y0 + l1 * e1y + l2 * e2y;
                double proj_z = z0 + l1 * e1z + l2 * e2z;
                double d2 = (px - proj_x) * (px - proj_x) +
                            (py - proj_y) * (py - proj_y) +
                            (pz - proj_z) * (pz - proj_z);

                if (d2 < best_dist) {
                    best_dist = d2;
                    best_l0 = l0; best_l1 = l1; best_l2 = l2;
                    best_v0 = v0; best_v1 = v1; best_v2 = v2;
                    found = true;
                }
            }
        }

        if (found) {
            if (best_l0 > 1e-15) { A_i.push_back(obs); A_j.push_back(best_v0); A_x.push_back(best_l0); }
            if (best_l1 > 1e-15) { A_i.push_back(obs); A_j.push_back(best_v1); A_x.push_back(best_l1); }
            if (best_l2 > 1e-15) { A_i.push_back(obs); A_j.push_back(best_v2); A_x.push_back(best_l2); }
        } else {
            double min_dist = 1e30;
            int nearest = 0;
            for (int v = 0; v < n_verts; v++) {
                double dx = px - vertices(v, 0);
                double dy = py - vertices(v, 1);
                double dz = pz - vertices(v, 2);
                double d = dx * dx + dy * dy + dz * dz;
                if (d < min_dist) { min_dist = d; nearest = v; }
            }
            A_i.push_back(obs);
            A_j.push_back(nearest);
            A_x.push_back(1.0);
        }
    }

    return List::create(
        Named("i") = wrap(A_i),
        Named("j") = wrap(A_j),
        Named("x") = wrap(A_x),
        Named("nrow") = n_obs,
        Named("ncol") = n_verts
    );
}
