#include "voronoi.hpp"

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Regular_triangulation_face_base_2.h>
#include <CGAL/Regular_triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <execution>
#include <limits>
#include <map>
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

namespace
{
    constexpr double k_clip_epsilon = 1e-9;
    constexpr long long k_left_side_label = -1;
    constexpr long long k_right_side_label = -2;
    constexpr long long k_bottom_side_label = -3;
    constexpr long long k_top_side_label = -4;

    using kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
    using cgal_point = kernel::Point_2;
    using cgal_weighted_point = kernel::Weighted_point_2;
    using vertex_info = std::size_t;

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

    struct clip_line
    {
        cz::point point_on_line;
        cz::point normal;
    };

    struct point_sites_view
    {
        explicit point_sites_view(std::span<const cz::point> sites) : sites(sites) {}

        std::size_t size() const
        {
            return sites.size();
        }

        cz::point point_at(std::size_t i) const
        {
            return sites[i];
        }

        clip_line bisector(std::size_t a, std::size_t b) const
        {
            const cz::point& site = sites[a];
            const cz::point& neighbor = sites[b];

            const cz::point midpoint{
                (site.x + neighbor.x) * 0.5,
                (site.y + neighbor.y) * 0.5
            };

            const cz::point normal = neighbor - site;
            return { midpoint, normal };
        }

        std::span<const cz::point> sites;
    };

    struct weighted_sites_view
    {
        explicit weighted_sites_view(std::span<const cz::weighted_point> sites) : sites(sites) {}

        std::size_t size() const
        {
            return sites.size();
        }

        cz::point point_at(std::size_t i) const
        {
            return sites[i].pt;
        }

        clip_line bisector(std::size_t a, std::size_t b) const
        {
            const cz::weighted_point& site = sites[a];
            const cz::weighted_point& neighbor = sites[b];

            const cz::point normal = neighbor.pt - site.pt;
            const double denom = 2.0 * cz::dot(normal, normal);

            if (denom <= k_clip_epsilon) {
                return { site.pt, { 0.0, 0.0 } };
            }

            const double rhs =
                cz::dot(neighbor.pt, neighbor.pt) - neighbor.weight -
                cz::dot(site.pt, site.pt) + site.weight;

            const cz::point point_on_line = normal * (rhs / denom);
            return { point_on_line, normal };
        }

        std::span<const cz::weighted_point> sites;
    };

    struct labeled_polygon
    {
        cz::polygon vertices;
        std::vector<long long> edge_labels;
    };

    using vertex_key = std::array<long long, 3>;

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

        for (std::size_t i = 0; i < poly.size(); ++i) {
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

            for (std::size_t i = 0; i < poly.size(); ++i) {
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

        for (std::size_t i = 0; i < input.size(); ++i) {
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

    cz::polygon rect_to_polygon(const cz::rect& bounds)
    {
        return {
            { bounds.min_point.x, bounds.min_point.y },
            { bounds.max_point.x, bounds.min_point.y },
            { bounds.max_point.x, bounds.max_point.y },
            { bounds.min_point.x, bounds.max_point.y }
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

        for (std::size_t i = 0; i < sites.size(); ++i) {
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

        for (std::size_t i = 0; i < sites.size(); ++i) {
            entries.emplace_back(to_cgal_weighted_point(sites[i]), i);
        }

        regular_triangulation rt;
        rt.insert(entries.begin(), entries.end());
        return rt;
    }

    std::vector<std::vector<std::size_t>> build_neighbor_lists(std::span<const cz::point> sites)
    {
        std::vector<std::vector<std::size_t>> neighbors(sites.size());

        if (sites.empty()) {
            return neighbors;
        }

        delaunay_triangulation dt = build_delaunay(sites);

        for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
            const std::size_t source_index = vit->info();
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

    std::vector<std::vector<std::size_t>> build_power_neighbor_lists(
        std::span<const cz::weighted_point> sites)
    {
        std::vector<std::vector<std::size_t>> neighbors(sites.size());

        if (sites.empty()) {
            return neighbors;
        }

        regular_triangulation rt = build_regular_triangulation(sites);

        for (auto vit = rt.finite_vertices_begin(); vit != rt.finite_vertices_end(); ++vit) {
            const std::size_t source_index = vit->info();
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

    template<typename Sites>
    cz::polygon construct_cell_polygon(
        const Sites& sites,
        std::size_t site_index,
        const std::vector<std::size_t>& neighbors,
        const cz::rect& bounds,
        double epsilon)
    {
        if (site_index >= sites.size()) {
            return {};
        }

        cz::polygon cell = rect_to_polygon(bounds);
        const cz::point site = sites.point_at(site_index);

        for (const std::size_t neighbor_index : neighbors) {
            if (neighbor_index >= sites.size()) {
                continue;
            }

            const cz::point neighbor = sites.point_at(neighbor_index);

            if (nearly_equal(site, neighbor, epsilon)) {
                continue;
            }

            const clip_line bisector = sites.bisector(site_index, neighbor_index);

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

    long long classify_boundary_edge(
        const cz::point& midpoint,
        const cz::rect& bounds,
        double epsilon)
    {
        if (nearly_equal(midpoint.x, bounds.min_point.x, epsilon)) {
            return k_left_side_label;
        }

        if (nearly_equal(midpoint.x, bounds.max_point.x, epsilon)) {
            return k_right_side_label;
        }

        if (nearly_equal(midpoint.y, bounds.min_point.y, epsilon)) {
            return k_bottom_side_label;
        }

        if (nearly_equal(midpoint.y, bounds.max_point.y, epsilon)) {
            return k_top_side_label;
        }

        return 0;
    }

    template<typename Sites>
    long long classify_edge_label(
        const Sites& sites,
        std::size_t site_index,
        const std::vector<std::size_t>& neighbors,
        const cz::point& a,
        const cz::point& b,
        const cz::rect& bounds,
        double epsilon)
    {
        const cz::point midpoint = (a + b) * 0.5;

        const long long boundary_label = classify_boundary_edge(midpoint, bounds, epsilon);
        if (boundary_label != 0) {
            return boundary_label;
        }

        double best_score = std::numeric_limits<double>::max();
        std::optional<std::size_t> best_neighbor;

        for (const std::size_t neighbor_index : neighbors) {
            const clip_line bisector = sites.bisector(site_index, neighbor_index);

            const double normal_norm = std::sqrt(cz::dot(bisector.normal, bisector.normal));
            if (normal_norm <= epsilon) {
                continue;
            }

            const double score =
                std::abs(cz::dot(midpoint - bisector.point_on_line, bisector.normal)) /
                normal_norm;

            if (score < best_score) {
                best_score = score;
                best_neighbor = neighbor_index;
            }
        }

        if (!best_neighbor.has_value()) {
            throw std::runtime_error("failed to classify Voronoi edge.");
        }

        return static_cast<long long>(*best_neighbor) + 1;
    }

    template<typename Sites>
    labeled_polygon label_polygon_edges(
        const Sites& sites,
        std::size_t site_index,
        const std::vector<std::size_t>& neighbors,
        const cz::polygon& polygon,
        const cz::rect& bounds,
        double epsilon)
    {
        labeled_polygon result;
        result.vertices = polygon;

        if (polygon.empty()) {
            return result;
        }

        result.edge_labels.reserve(polygon.size());

        for (std::size_t i = 0; i < polygon.size(); ++i) {
            const cz::point& a = polygon[i];
            const cz::point& b = polygon[(i + 1) % polygon.size()];

            result.edge_labels.push_back(
                classify_edge_label(
                    sites,
                    site_index,
                    neighbors,
                    a,
                    b,
                    bounds,
                    epsilon
                )
            );
        }

        return result;
    }

    vertex_key make_vertex_key(long long a, long long b, long long c)
    {
        vertex_key key{ a, b, c };
        std::sort(key.begin(), key.end());
        return key;
    }

    template<typename Sites>
    cz::voronoi_embedding build_embedding_from_polygons(
        const Sites& sites,
        const std::vector<std::vector<std::size_t>>& graph,
        const std::vector<cz::polygon>& polygons,
        const cz::rect& bounds,
        double epsilon)
    {
        cz::voronoi_embedding embedding;
        embedding.cells.resize(polygons.size());

        std::map<vertex_key, std::size_t> vertex_map;

        for (std::size_t cell_index = 0; cell_index < polygons.size(); ++cell_index) {
            const cz::polygon& polygon = polygons[cell_index];

            if (polygon.empty()) {
                continue;
            }

            const labeled_polygon labeled = label_polygon_edges(
                sites,
                cell_index,
                graph[cell_index],
                polygon,
                bounds,
                epsilon
            );

            auto& cell_vertices = embedding.cells[cell_index];
            cell_vertices.reserve(labeled.vertices.size());

            for (std::size_t i = 0; i < labeled.vertices.size(); ++i) {
                const long long prev_label =
                    labeled.edge_labels[(i + labeled.edge_labels.size() - 1) %
                    labeled.edge_labels.size()];
                const long long next_label = labeled.edge_labels[i];
                const long long cell_label = static_cast<long long>(cell_index) + 1;

                const vertex_key key = make_vertex_key(cell_label, prev_label, next_label);

                auto [it, inserted] = vertex_map.emplace(key, embedding.vertices.size());
                if (inserted) {
                    embedding.vertices.push_back(labeled.vertices[i]);
                }

                cell_vertices.push_back(it->second);
            }
        }

        return embedding;
    }

    template<typename Site, typename GetPoint, typename RebindPoint>
    std::vector<Site> perform_lloyd_relaxation_impl(
        std::span<const Site> sites,
        double min_delta_thresh,
        int max_iterations,
        const cz::rect& bounds,
        GetPoint get_point,
        RebindPoint rebind_point)
    {
        if (sites.empty()) {
            return {};
        }

        int iter = 0;
        double max_delta = std::numeric_limits<double>::max();
        auto input = sites | r::to<std::vector>();

        while (iter++ < max_iterations && max_delta > min_delta_thresh) {
            const auto graph = cz::to_voronoi_topology(input, bounds);
            const auto polygons = cz::to_voronoi_polygons(input, graph, bounds);

            auto centroids = polygons
                | rv::transform([](const auto& poly) -> cz::point {
                if (poly.empty()) {
                    throw std::runtime_error(
                        "attempted to centroid an empty Voronoi cell"
                    );
                }

                return cz::centroid(poly);
                    })
                | r::to<std::vector>();

            max_delta = r::max(
                rv::zip(centroids, input)
                | rv::transform([&](const auto& pair) -> double {
                    const auto& [centroid, site] = pair;
                    return cz::distance(centroid, get_point(site));
                    })
            );

            input = rv::zip(centroids, input)
                | rv::transform([&](const auto& pair) -> Site {
                const auto& [centroid, site] = pair;
                return rebind_point(site, centroid);
                    })
                | r::to<std::vector>();
        }

        return input;
    }

} // namespace

/*------------------------------------------------------------------------------------------------*/

std::vector<std::vector<std::size_t>> cz::to_voronoi_topology(
    std::span<const point> sites,
    const rect& bounds)
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
    const std::vector<std::vector<std::size_t>>& graph,
    const rect& bounds)
{
    if (sites.empty() || graph.size() != sites.size() || !is_valid_bounds(bounds)) {
        return {};
    }

    std::vector<cz::polygon> result(sites.size());
    std::vector<std::size_t> indices(sites.size());
    std::iota(indices.begin(), indices.end(), std::size_t{ 0 });

    const point_sites_view site_view{ sites };

    std::for_each(
        std::execution::par,
        indices.begin(),
        indices.end(),
        [&](std::size_t i) {
            result[i] = construct_cell_polygon(
                site_view,
                i,
                graph[i],
                bounds,
                k_clip_epsilon
            );
        }
    );

    return result;
}

cz::voronoi_embedding cz::to_voronoi_embedding(
    std::span<const point> sites,
    const rect& bounds)
{
    if (sites.empty() || !is_valid_bounds(bounds)) {
        return {};
    }

    const auto graph = to_voronoi_topology(sites, bounds);
    const auto polygons = to_voronoi_polygons(sites, graph, bounds);
    return build_embedding_from_polygons(
        point_sites_view{ sites },
        graph,
        polygons,
        bounds,
        k_clip_epsilon
    );
}

std::vector<cz::point> cz::perform_lloyd_relaxation(
    std::span<const point> sites,
    double min_delta_thresh,
    int max_iterations,
    const rect& bounds)
{
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

std::vector<std::vector<std::size_t>> cz::to_voronoi_topology(
    std::span<const weighted_point> sites,
    const rect& bounds)
{
    if (sites.empty()) {
        return {};
    }

    if (!is_valid_bounds(bounds)) {
        return {};
    }

    return build_power_neighbor_lists(sites);
}

std::vector<cz::polygon> cz::to_voronoi_polygons(
    std::span<const weighted_point> sites,
    const std::vector<std::vector<std::size_t>>& graph,
    const rect& bounds)
{
    if (sites.empty() || graph.size() != sites.size() || !is_valid_bounds(bounds)) {
        return {};
    }

    std::vector<cz::polygon> result(sites.size());
    std::vector<std::size_t> indices(sites.size());
    std::iota(indices.begin(), indices.end(), std::size_t{ 0 });

    const weighted_sites_view site_view{ sites };

    std::for_each(
        std::execution::par,
        indices.begin(),
        indices.end(),
        [&](std::size_t i) {
            result[i] = construct_cell_polygon(
                site_view,
                i,
                graph[i],
                bounds,
                k_clip_epsilon
            );
        }
    );

    return result;
}

cz::voronoi_embedding cz::to_voronoi_embedding(
    std::span<const weighted_point> sites,
    const rect& bounds)
{
    if (sites.empty() || !is_valid_bounds(bounds)) {
        return {};
    }

    const auto graph = to_voronoi_topology(sites, bounds);
    const auto polygons = to_voronoi_polygons(sites, graph, bounds);
    return build_embedding_from_polygons(
        weighted_sites_view{ sites },
        graph,
        polygons,
        bounds,
        k_clip_epsilon
    );
}

std::vector<cz::point> cz::perform_lloyd_relaxation(
    std::span<const weighted_point> sites,
    double min_delta_thresh,
    int max_iterations,
    const rect& bounds)
{
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

    return relaxed
        | rv::transform([](const cz::weighted_point& p) -> cz::point {
        return p.pt;
            })
        | r::to<std::vector>();
}