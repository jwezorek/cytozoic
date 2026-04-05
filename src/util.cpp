#include "util.hpp"
#include <random>
#include <ranges>

namespace r = std::ranges;
namespace rv = std::ranges::views;

/*------------------------------------------------------------------------------------------------*/

namespace {

    std::mt19937_64 make_rng(std::optional<uint64_t> seed = std::nullopt) {
        if (seed.has_value()) {
            return std::mt19937_64(*seed);
        }

        std::random_device rd;
        std::seed_seq seq{
            rd(), rd(), rd(), rd(),
            rd(), rd(), rd(), rd()
        };
        return std::mt19937_64(seq);
    }

}

std::vector<cz::point> cz::random_points(
        size_t n,
        double bounds_wd,
        double bounds_hgt,
        std::optional<uint64_t> seed) {

    if (n == 0 || bounds_wd <= 0.0 || bounds_hgt <= 0.0) {
        return {};
    }

    auto rng = make_rng(seed);

    std::uniform_real_distribution<double> x_dist(0.0, bounds_wd);
    std::uniform_real_distribution<double> y_dist(0.0, bounds_hgt);

    std::vector<cz::point> points;
    points.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        points.push_back({
            .x = x_dist(rng),
            .y = y_dist(rng)
            });
    }

    return points;
}

cz::cyto_frame cz::to_cyto_frame(const cyto_state& state, const color_table& colors) {
    return rv::zip(state.cells, state.states) |
        rv::transform(
            [&](auto&& cell_state)->cell {
                const auto& [vor_cell, s] = cell_state;
                return {
                    .shape = vor_cell.cell,
                    .color = colors[s],
                    .seed = vor_cell.site
                };
            }
        ) | r::to<cyto_frame>();
}

cz::cyto_frame cz::to_cyto_frame(
        std::span<const point> pts, std::span<const polygon> polys, std::span<const color> colors) {

    return rv::zip(polys, colors, pts) | rv::transform(
            [](const auto& v)->cell {
                return {
                    .shape = std::get<0>(v),
                    .color = std::get<1>(v),
                    .seed = std::get<2>(v)
                };
            }
        ) | r::to<std::vector>();
}

cz::cyto_state cz::blank_state(const voronoi_diagram& v) {
    return {
        .cells = v,
        .states = std::vector<int8_t>(v.size(), int8_t{0})
    };
}

cz::point cz::centroid(std::span<const point> pts) {
    if (pts.empty()) {
        throw std::runtime_error(
            "attempted to find centroid of an empty point set."
        );
    }
    auto sum = r::fold_left(  pts, point{ 0.0,0.0 }, std::plus<>{} );
    const double n = static_cast<double>(pts.size());
    return {
        sum.x / n, sum.y / n
    };
}

double cz::dot(const point& u, const point& v) {
    return u.x * v.x + u.y * v.y;
}

double cz::distance(const point& u, const point& v) {
    auto diff = v - u;
    return std::sqrt(dot(diff,diff));
}

double cz::magnitude(const point& u) {
    return std::sqrt(dot(u, u));
}

cz::point cz::interpolate_point(const point& from, const point& to, double t) {
    t = std::clamp(t, 0.0, 1.0);
    return from + (to - from) * t;
}

cz::color cz::interpolate_color(const color& from, const color& to, double t) {
    t = std::clamp(t, 0.0, 1.0);

    auto lerp_channel = [t](uint8_t a, uint8_t b) -> uint8_t {
        const double value = static_cast<double>(a)
            + (static_cast<double>(b) - static_cast<double>(a)) * t;
        return static_cast<uint8_t>(std::round(value));
        };

    return {
        lerp_channel(from.r, to.r),
        lerp_channel(from.g, to.g),
        lerp_channel(from.b, to.b)
    };
}
