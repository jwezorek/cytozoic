#include "voronoi.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Regular_triangulation_face_base_2.h>
#include <CGAL/Regular_triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <algorithm>
#include <cmath>
#include <execution>
#include <limits>
#include <numeric>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace r = std::ranges;
namespace rv = std::ranges::views;

/*------------------------------------------------------------------------------------------------*/

namespace {

    constexpr double k_clip_epsilon = 1e-9;

    using kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
    using cgal_point = kernel::Point_2;
    using cgal_weighted_point = kernel::Weighted_point_2;
    using vertex_info = size_t;

    using delaunay_vertex_base =
        CGAL::Triangulation_vertex_base_with_info_2<vertex_info, kernel>;
    using delaunay_face_base = CGAL::Triangulation_face_base_2<kernel>;
    using delaunay_tds =
        CGAL::Triangulation_data_structure_2<delaunay_vertex_base, delaunay_face_base>;
    using delaunay_triangulation = CGAL::Delaunay_triangulation_2<kernel, delaunay_tds>;
    using delaunay_vertex_handle = delaunay_triangulation::Vertex_handle;
    using delaunay_vertex_circulator = delaunay_triangulation::Vertex_circulator;

    using regular_vertex_base_raw = CGAL::Regular_triangulation_vertex_base_2<kernel>;
    using regular_vertex_base =
        CGAL::Triangulation_vertex_base_with_info_2<vertex_info, kernel, regular_vertex_base_raw>;
    using regular_face_base = CGAL::Regular_triangulation_face_base_2<kernel>;
    using regular_tds =
        CGAL::Triangulation_data_structure_2<regular_vertex_base, regular_face_base>;
    using regular_triangulation = CGAL::Regular_triangulation_2<kernel, regular_tds>;
    using regular_vertex_handle = regular_triangulation::Vertex_handle;
    using regular_vertex_circulator = regular_triangulation::Vertex_circulator;

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

    double squared_norm(const cz::point& p)
    {
        return cz::dot(p, p);
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

    clip_line make_power_bisector(
        const cz::weighted_point& site,
        const cz::weighted_point& neighbor)
    {
        const cz::point normal = neighbor.pt - site.pt;
        const double denom = 2.0 * squared_norm(normal);

        if (denom <= k_clip_epsilon) {
            return { site.pt, { 0.0, 0.0 } };
        }

        const double rhs =
            squared_norm(neighbor.pt) - neighbor.weight -
            squared_norm(site.pt) + site.weight;

        const cz::point point_on_line = normal * (rhs / denom);
        return { point_on_line, normal };
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

    cgal_weighted_point to_cgal_weighted_point(const cz::weighted_point& p)
    {
        return { to_cgal_point(p.pt), p.weight };
    }

    delaunay_triangulation build_delaunay(std::span<const cz::point> sites)
    {
        delaunay_triangulation dt;

        for (size_t i = 0; i < sites.size(); ++i) {
            delaunay_vertex_handle vh = dt.insert(to_cgal_point(sites[i]));
            vh->info() = i;
        }

        return dt;
    }

    regular_triangulation build_regular_triangulation(
        std::span<const cz::weighted_point> sites)
    {
        std::vector<std::pair<cgal_weighted_point, vertex_info>> entries;
        entries.reserve(sites.size());

        for (size_t i = 0; i < sites.size(); ++i) {
            entries.emplace_back(to_cgal_weighted_point(sites[i]), i);
        }

        regular_triangulation rt;
        rt.insert(entries.begin(), entries.end());
        return rt;
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

            delaunay_vertex_circulator circ = dt.incident_vertices(vit);
            if (circ == 0) {
                continue;
            }

            delaunay_vertex_circulator start = circ;
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

    std::vector<std::vector<size_t>> build_power_neighbor_lists(
        std::span<const cz::weighted_point> sites)
    {
        std::vector<std::vector<size_t>> neighbors(sites.size());

        if (sites.empty()) {
            return neighbors;
        }

        regular_triangulation rt = build_regular_triangulation(sites);

        for (auto vit = rt.finite_vertices_begin(); vit != rt.finite_vertices_end(); ++vit) {
            const size_t source_index = vit->info();
            auto& out = neighbors[source_index];

            regular_vertex_circulator circ = rt.incident_vertices(vit);
            if (circ == 0) {
                continue;
            }

            regular_vertex_circulator start = circ;
            do {
                if (!rt.is_infinite(circ)) {
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

    cz::polygon construct_power_cell_polygon(
        std::span<const cz::weighted_point> sites,
        size_t site_index,
        const std::vector<size_t>& neighbors,
        const cz::rect& bounds,
        double epsilon)
    {
        if (site_index >= sites.size()) {
            return {};
        }

        cz::polygon cell = rect_to_polygon(bounds);
        const cz::weighted_point& site = sites[site_index];

        for (const size_t neighbor_index : neighbors) {
            if (neighbor_index >= sites.size()) {
                continue;
            }

            const cz::weighted_point& neighbor = sites[neighbor_index];

            if (nearly_equal(site.pt, neighbor.pt, epsilon)) {
                if (neighbor.weight >= site.weight) {
                    return {};
                }
                continue;
            }

            const clip_line bisector = make_power_bisector(site, neighbor);

            if (std::abs(bisector.normal.x) <= epsilon &&
                std::abs(bisector.normal.y) <= epsilon) {
                continue;
            }

            cell = clip_polygon_against_half_plane(cell, bisector, epsilon);

            if (cell.empty()) {
                return {};
            }
        }

        normalize_polygon(cell, epsilon);
        return cell;
    }

    template<typename Site, typename GetPoint, typename RebindPoint>
    std::vector<Site> perform_lloyd_relaxation_impl(
            std::span<const Site> sites,
            double min_delta_thresh,
            int max_iterations,
            const cz::rect& bounds,
            GetPoint get_point,
            RebindPoint rebind_point) {

        if (sites.empty()) {
            return {};
        }

        int iter = 0;
        double max_delta = std::numeric_limits<double>::max();
        auto input = sites | r::to<std::vector>();

        while (iter++ < max_iterations && max_delta > min_delta_thresh) {
            auto polygons = cz::to_voronoi_polygons(input, bounds);

            auto centroids = polygons | rv::transform(
                    [](const auto& poly) -> cz::point {
                        if (poly.empty()) {
                            throw std::runtime_error(
                                "attempted to centroid an empty Voronoi cell"
                            );
                        }
                        return cz::centroid(poly);
                        }
                ) | r::to<std::vector>();

            max_delta = r::max(
                rv::zip(centroids, input) | rv::transform(
                    [&](const auto& pair) -> double {
                        const auto& [centroid, site] = pair;
                        return cz::distance(centroid, get_point(site));
                    }
                )
            );

            input = rv::zip(centroids, input)  | rv::transform(
                    [&](const auto& pair) -> Site {
                        const auto& [centroid, site] = pair;
                        return rebind_point(site, centroid);
                    }
                ) | r::to<std::vector>();
        }

        return input;
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
    if (sites.empty() || graph.size() != sites.size() || !is_valid_bounds(bounds)) {
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

std::vector<cz::polygon> cz::to_voronoi_polygons( std::span<const cz::point> v, const cz::rect& b) {
    const auto graph = to_voronoi_topology(v, b);
    return to_voronoi_polygons(v, graph, b);
}

std::vector<cz::point> cz::perform_lloyd_relaxation( std::span<const point> sites,
        double min_delta_thresh, int max_iterations, const rect& bounds) {

    return perform_lloyd_relaxation_impl(
        sites,
        min_delta_thresh,
        max_iterations,
        bounds,
        [](const cz::point& p) -> cz::point {
            return p;
        },
        [](const cz::point&, const cz::point& new_pt) -> cz::point {
            return new_pt;
        }
    );
}

std::vector<std::vector<size_t>> cz::to_voronoi_topology(
        std::span<const weighted_point> sites, const rect& bounds) {

    if (sites.empty()) {
        return {};
    }

    if (!is_valid_bounds(bounds)) {
        return {};
    }

    return build_power_neighbor_lists(sites);
}

std::vector<cz::polygon> cz::to_voronoi_polygons( std::span<const weighted_point> sites,
        const std::vector<std::vector<size_t>>& graph, const rect& bounds) {

    if (sites.empty() || graph.size() != sites.size() || !is_valid_bounds(bounds)) {
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
            result[i] = construct_power_cell_polygon(
                sites,
                i,
                graph[i],
                bounds,
                k_clip_epsilon
            );
        });

    return result;
}

std::vector<cz::polygon> cz::to_voronoi_polygons( std::span<const weighted_point> sites,
        const rect& bounds) {
    const auto graph = to_voronoi_topology(sites, bounds);
    return to_voronoi_polygons(sites, graph, bounds);
}

std::vector<cz::point> cz::perform_lloyd_relaxation(
        std::span<const weighted_point> sites,
        double min_delta_thresh,
        int max_iterations,
        const rect& bounds) {

    auto relaxed = perform_lloyd_relaxation_impl(
        sites,
        min_delta_thresh,
        max_iterations,
        bounds,
        [](const cz::weighted_point& p) -> cz::point {
            return p.pt;
        },
        [](const cz::weighted_point& old, const cz::point& new_pt) -> cz::weighted_point {
            return { new_pt, old.weight };
        }
    );

    return relaxed | rv::transform(
            [](const cz::weighted_point& p) -> cz::point {
                return p.pt;
            }
        ) | r::to<std::vector>();
}