#include "voronoi.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <algorithm>
#include <cmath>
#include <execution>
#include <limits>
#include <numeric>
#include <ranges>
#include <span>
#include <stdexcept>
#include <vector>

namespace r = std::ranges;
namespace rv = std::ranges::views;

namespace {

    constexpr double k_clip_epsilon = 1e-9;

    using kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
    using cgal_point = kernel::Point_2;
    using vertex_info = size_t;
    using vertex_base = CGAL::Triangulation_vertex_base_with_info_2<vertex_info, kernel>;
    using face_base = CGAL::Triangulation_face_base_2<kernel>;
    using tds = CGAL::Triangulation_data_structure_2<vertex_base, face_base>;
    using delaunay_triangulation = CGAL::Delaunay_triangulation_2<kernel, tds>;
    using vertex_handle = delaunay_triangulation::Vertex_handle;
    using vertex_circulator = delaunay_triangulation::Vertex_circulator;

    struct clip_line {
        cz::point point_on_line;
        cz::point normal;
    };

    bool nearly_equal(double a, double b, double epsilon)
    {
        return std::abs(a - b) <= epsilon;
    }

    bool nearly_equal(const cz::point& a, const cz::point& b, double epsilon)
    {
        return nearly_equal(a.x, b.x, epsilon) &&
            nearly_equal(a.y, b.y, epsilon);
    }

    double signed_area_times_two(const cz::polygon& poly)
    {
        if (poly.size() < 3) {
            return 0.0;
        }

        double area2 = 0.0;

        for (size_t i = 0; i < poly.size(); ++i) {
            const cz::point& a = poly[i];
            const cz::point& b = poly[(i + 1) % poly.size()];
            area2 += a.x * b.y - b.x * a.y;
        }

        return area2;
    }

    bool is_degenerate_polygon(const cz::polygon& poly, double epsilon)
    {
        if (poly.size() < 3) {
            return true;
        }

        return std::abs(signed_area_times_two(poly)) <= epsilon;
    }

    void remove_consecutive_duplicate_vertices(cz::polygon& poly, double epsilon)
    {
        if (poly.empty()) {
            return;
        }

        cz::polygon cleaned;
        cleaned.reserve(poly.size());

        for (const cz::point& p : poly) {
            if (cleaned.empty() || !nearly_equal(cleaned.back(), p, epsilon)) {
                cleaned.push_back(p);
            }
        }

        if (cleaned.size() >= 2 &&
            nearly_equal(cleaned.front(), cleaned.back(), epsilon)) {
            cleaned.pop_back();
        }

        poly = std::move(cleaned);
    }

    void remove_collinear_vertices(cz::polygon& poly, double epsilon)
    {
        if (poly.size() < 3) {
            return;
        }

        bool changed = true;

        while (changed && poly.size() >= 3) {
            changed = false;
            cz::polygon cleaned;
            cleaned.reserve(poly.size());

            for (size_t i = 0; i < poly.size(); ++i) {
                const cz::point& prev = poly[(i + poly.size() - 1) % poly.size()];
                const cz::point& curr = poly[i];
                const cz::point& next = poly[(i + 1) % poly.size()];

                const cz::point a = curr - prev;
                const cz::point b = next - curr;
                const double cross = a.x * b.y - a.y * b.x;

                if (std::abs(cross) <= epsilon &&
                    !(nearly_equal(prev, curr, epsilon) ||
                        nearly_equal(curr, next, epsilon))) {
                    changed = true;
                    continue;
                }

                cleaned.push_back(curr);
            }

            poly = std::move(cleaned);
        }
    }

    void normalize_polygon(cz::polygon& poly, double epsilon)
    {
        remove_consecutive_duplicate_vertices(poly, epsilon);
        remove_collinear_vertices(poly, epsilon);
        remove_consecutive_duplicate_vertices(poly, epsilon);

        if (poly.size() >= 3 && signed_area_times_two(poly) < 0.0) {
            std::reverse(poly.begin(), poly.end());
        }

        if (is_degenerate_polygon(poly, epsilon)) {
            poly.clear();
        }
    }

    bool point_is_inside_half_plane(
        const cz::point& p,
        const clip_line& line,
        double epsilon)
    {
        const cz::point offset = p - line.point_on_line;
        return cz::dot(offset, line.normal) <= epsilon;
    }

    std::optional<cz::point> intersect_segment_with_line(
        const cz::point& a,
        const cz::point& b,
        const clip_line& line,
        double epsilon)
    {
        const cz::point segment = b - a;
        const double denominator = cz::dot(segment, line.normal);

        if (std::abs(denominator) <= epsilon) {
            return std::nullopt;
        }

        const double t = cz::dot(line.point_on_line - a, line.normal) / denominator;
        return a + t * segment;
    }

    void push_if_distinct(cz::polygon& poly, const cz::point& point, double epsilon)
    {
        if (poly.empty() || !nearly_equal(poly.back(), point, epsilon)) {
            poly.push_back(point);
        }
    }

    cz::polygon clip_polygon_against_half_plane(
        const cz::polygon& input,
        const clip_line& line,
        double epsilon)
    {
        if (input.empty()) {
            return {};
        }

        cz::polygon output;
        output.reserve(input.size() + 1);

        for (size_t i = 0; i < input.size(); ++i) {
            const cz::point& current = input[i];
            const cz::point& next = input[(i + 1) % input.size()];

            const bool current_inside = point_is_inside_half_plane(current, line, epsilon);
            const bool next_inside = point_is_inside_half_plane(next, line, epsilon);

            if (current_inside && next_inside) {
                push_if_distinct(output, next, epsilon);
                continue;
            }

            if (current_inside && !next_inside) {
                const std::optional<cz::point> intersection =
                    intersect_segment_with_line(current, next, line, epsilon);

                if (intersection.has_value()) {
                    push_if_distinct(output, *intersection, epsilon);
                }

                continue;
            }

            if (!current_inside && next_inside) {
                const std::optional<cz::point> intersection =
                    intersect_segment_with_line(current, next, line, epsilon);

                if (intersection.has_value()) {
                    push_if_distinct(output, *intersection, epsilon);
                }

                push_if_distinct(output, next, epsilon);
            }
        }

        if (output.size() >= 2 &&
            nearly_equal(output.front(), output.back(), epsilon)) {
            output.pop_back();
        }

        return output;
    }

    clip_line make_voronoi_bisector(const cz::point& site, const cz::point& neighbor)
    {
        const cz::point midpoint{
            (site.x + neighbor.x) * 0.5,
            (site.y + neighbor.y) * 0.5
        };

        const cz::point normal = neighbor - site;
        return { midpoint, normal };
    }

    cz::polygon rect_to_polygon(const cz::rect& bounds)
    {
        return {
            {bounds.min_point.x, bounds.min_point.y},
            {bounds.max_point.x, bounds.min_point.y},
            {bounds.max_point.x, bounds.max_point.y},
            {bounds.min_point.x, bounds.max_point.y}
        };
    }

    bool is_valid_bounds(const cz::rect& bounds)
    {
        return bounds.min_point.x < bounds.max_point.x &&
            bounds.min_point.y < bounds.max_point.y;
    }

    cgal_point to_cgal_point(const cz::point& p)
    {
        return { p.x, p.y };
    }

    delaunay_triangulation build_delaunay(std::span<const cz::point> sites)
    {
        delaunay_triangulation dt;

        for (size_t i = 0; i < sites.size(); ++i) {
            vertex_handle vh = dt.insert(to_cgal_point(sites[i]));
            vh->info() = i;
        }

        return dt;
    }

    std::vector<std::vector<size_t>> build_neighbor_lists(std::span<const cz::point> sites)
    {
        std::vector<std::vector<size_t>> neighbors(sites.size());

        if (sites.empty()) {
            return neighbors;
        }

        delaunay_triangulation dt = build_delaunay(sites);

        for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
            const size_t source_index = vit->info();
            auto& out = neighbors[source_index];

            vertex_circulator circ = dt.incident_vertices(vit);
            if (circ == 0) {
                continue;
            }

            vertex_circulator start = circ;
            do {
                if (!dt.is_infinite(circ)) {
                    out.push_back(circ->info());
                }
                ++circ;
            } while (circ != start);

            std::sort(out.begin(), out.end());
            out.erase(std::unique(out.begin(), out.end()), out.end());
        }

        return neighbors;
    }

    cz::polygon construct_cell_polygon(
        std::span<const cz::point> sites,
        size_t site_index,
        const std::vector<size_t>& neighbors,
        const cz::rect& bounds,
        double epsilon)
    {
        if (site_index >= sites.size()) {
            return {};
        }

        cz::polygon cell = rect_to_polygon(bounds);
        const cz::point& site = sites[site_index];

        for (const size_t neighbor_index : neighbors) {
            if (neighbor_index >= sites.size()) {
                continue;
            }

            const cz::point& neighbor = sites[neighbor_index];

            if (nearly_equal(site, neighbor, epsilon)) {
                continue;
            }

            const clip_line bisector = make_voronoi_bisector(site, neighbor);
            cell = clip_polygon_against_half_plane(cell, bisector, epsilon);

            if (cell.empty()) {
                return {};
            }
        }

        normalize_polygon(cell, epsilon);
        return cell;
    }

} // namespace

std::vector<std::vector<size_t>> cz::to_voronoi_topology(
    std::span<const cz::point> sites,
    const cz::rect& bounds)
{
    if (sites.empty()) {
        return {};
    }

    if (!is_valid_bounds(bounds)) {
        return {};
    }

    return build_neighbor_lists(sites);
}

std::vector<cz::polygon> cz::to_voronoi_polygons(
    std::span<const point> sites,
    const std::vector<std::vector<size_t>>& graph,
    const rect& bounds)
{
    if (graph.empty()) {
        return {};
    }

    std::vector<cz::polygon> result(sites.size());
    std::vector<size_t> indices(sites.size());
    std::iota(indices.begin(), indices.end(), size_t{ 0 });

    std::for_each(
        std::execution::par,
        indices.begin(),
        indices.end(),
        [&](size_t i) {
            result[i] = construct_cell_polygon(
                sites,
                i,
                graph[i],
                bounds,
                k_clip_epsilon
            );
        });

    return result;
}

std::vector<cz::polygon> cz::to_voronoi_polygons(
    std::span<const cz::point> sites,
    const cz::rect& bounds)
{
    auto graph = to_voronoi_topology(sites, bounds);
    return to_voronoi_polygons(sites, graph, bounds);
}

std::vector<cz::point> cz::perform_lloyd_relaxation(
    std::span<const point> sites,
    double min_delta_thresh,
    int max_iterations,
    const rect& bounds)
{
    if (sites.empty()) {
        return {};
    }

    int iter = 0;
    double max_delta = std::numeric_limits<double>::max();
    auto points = sites | r::to<std::vector>();

    while (iter++ < max_iterations && max_delta > min_delta_thresh) {
        auto polygons = to_voronoi_polygons(points, bounds);

        auto centroids = polygons |
            rv::transform(
                [](const auto& p) -> point {
                    if (p.empty()) {
                        throw std::runtime_error("attempted to centroid an empty Voronoi cell");
                    }
                    return cz::centroid(p);
                }
            ) | r::to<std::vector>();

        max_delta = r::max(
            rv::zip(centroids, points) |
            rv::transform([](const auto& pair) -> double {
                const auto& [lhs, rhs] = pair;
                return cz::distance(lhs, rhs);
                })
        );

        points = std::move(centroids);
    }

    return points;
}