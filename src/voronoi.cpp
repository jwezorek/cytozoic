#include "voronoi.hpp"
#include "util.hpp"

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 5055)
#endif

#include <boost/polygon/voronoi.hpp>

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <execution>
#include <limits>
#include <numeric>
#include <optional>
#include <ranges>
#include <span>
#include <vector>

namespace bp = boost::polygon;
namespace r = std::ranges;
namespace rv = std::ranges::views;

namespace cz::detail {

    using int_point = vec2<std::int64_t>;

} // namespace cz::detail

namespace boost::polygon {

    template <>
    struct geometry_concept<cz::detail::int_point> {
        using type = point_concept;
    };

    template <>
    struct point_traits<cz::detail::int_point> {
        using coordinate_type = std::int64_t;

        static coordinate_type get(
            const cz::detail::int_point& point,
            orientation_2d orient)
        {
            return orient == HORIZONTAL ? point.x : point.y;
        }
    };

} // namespace boost::polygon

namespace {

    constexpr double k_coordinate_scale = 1000000.0;
    constexpr double k_clip_epsilon = 1e-9;

    using cz::detail::int_point;

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

    void remove_collinear_vertices(cz::polygon& poly, double epsilon) {
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

    void normalize_polygon(cz::polygon& poly, double epsilon) {
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

    std::vector<int_point> scale_points_to_integer_grid(
        std::span<const cz::point> sites,
        double coordinate_scale)
    {
        std::vector<int_point> result;
        result.reserve(sites.size());

        for (const cz::point& site : sites) {
            result.push_back({
                static_cast<std::int64_t>(std::llround(site.x * coordinate_scale)),
                static_cast<std::int64_t>(std::llround(site.y * coordinate_scale))
                });
        }

        return result;
    }

    bool point_is_inside_half_plane(
        const cz::point& p,
        const clip_line& line,
        double epsilon)
    {
        const cz::point offset = p - line.point_on_line;
        return dot(offset, line.normal) <= epsilon;
    }

    std::optional<cz::point> intersect_segment_with_line(
        const cz::point& a,
        const cz::point& b,
        const clip_line& line,
        double epsilon)
    {
        const cz::point segment = b - a;
        const double denominator = dot(segment, line.normal);

        if (std::abs(denominator) <= epsilon) {
            return std::nullopt;
        }

        const double t = dot(line.point_on_line - a, line.normal) / denominator;
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

    std::vector<std::vector<size_t>> build_neighbor_lists(
        std::span<const cz::point> sites,
        double coordinate_scale)
    {
        std::vector<std::vector<size_t>> neighbors(sites.size());

        if (sites.empty()) {
            return neighbors;
        }

        for (auto& cell_neighbors : neighbors) {
            cell_neighbors.reserve(8);
        }

        const std::vector<int_point> integer_sites =
            scale_points_to_integer_grid(sites, coordinate_scale);

        bp::voronoi_diagram<double> diagram;
        bp::construct_voronoi(integer_sites.begin(), integer_sites.end(), &diagram);

        for (const auto& cell : diagram.cells()) {
            const size_t source_index = cell.source_index();

            if (source_index >= sites.size()) {
                continue;
            }

            const auto* edge = cell.incident_edge();

            if (edge == nullptr) {
                continue;
            }

            const auto* start = edge;

            do {
                const auto* twin = edge->twin();

                if (twin != nullptr) {
                    const auto* adjacent_cell = twin->cell();

                    if (adjacent_cell != nullptr) {
                        const size_t neighbor_index = adjacent_cell->source_index();

                        if (neighbor_index < sites.size() &&
                            neighbor_index != source_index) {
                            neighbors[source_index].push_back(neighbor_index);
                        }
                    }
                }

                edge = edge->next();
            } while (edge != start);
        }

        for (std::vector<size_t>& cell_neighbors : neighbors) {
            std::sort(cell_neighbors.begin(), cell_neighbors.end());
            cell_neighbors.erase(
                std::unique(cell_neighbors.begin(), cell_neighbors.end()),
                cell_neighbors.end());
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

cz::voronoi_diagram cz::construct_voronoi_diagram( 
        std::span<const cz::point> sites, const cz::rect& bounds) 
{
    if (sites.empty()) {
        return {};
    }

    if (!is_valid_bounds(bounds)) {
        return {};
    }

    std::vector<voronoi_cell> result(sites.size());

    const std::vector<std::vector<size_t>> neighbors =
        build_neighbor_lists(sites, k_coordinate_scale);

    std::vector<size_t> indices(sites.size());
    std::iota(indices.begin(), indices.end(), size_t{ 0 });

    std::for_each(
        std::execution::par,
        indices.begin(),
        indices.end(),
        [&](size_t i) {
            result[i].site = sites[i];
            result[i].neighbors = neighbors[i];

            if (neighbors[i].empty()) {
                return;
            }

            result[i].cell = construct_cell_polygon(
                sites,
                i,
                neighbors[i],
                bounds,
                k_clip_epsilon
            );
        });

    return result;
}

std::vector<cz::polygon> cz::to_voronoi_polygons(
    std::span<const point> sites,
    const rect& bounds)
{
    return construct_voronoi_diagram(sites, bounds) |
        rv::transform([](const auto& cell) -> polygon {
        return cell.cell;
            }) |
        r::to<std::vector>();
}

std::vector<cz::point> cz::perform_lloyd_relaxation(
    std::span<const point> sites,
    const rect& bounds,
    double min_delta,
    int max_iterations)
{
    if (sites.empty()) {
        return {};
    }

    int iter = 0;
    double max_delta = std::numeric_limits<double>::max();
    auto points = sites | r::to<std::vector>();

    while (iter++ < max_iterations && max_delta > min_delta) {
        auto voronoi = construct_voronoi_diagram(points, bounds);

        auto centroids = voronoi |
            rv::transform([](const auto& c) -> point {
            if (c.cell.empty()) {
                return c.site;
            }
            return centroid(c.cell);
                }) |
            r::to<std::vector>();

                max_delta = r::max(
                    rv::zip(centroids, points) |
                    rv::transform([](const auto& pair) -> double {
                        const auto& [lhs, rhs] = pair;
                        return distance(lhs, rhs);
                        })
                );

                points = std::move(centroids);
    }

    return points;
}