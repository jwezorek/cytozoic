#include "cytozoic.hpp"
#include "voronoi.hpp"
#include <ranges>
#include <cmath>
#include <algorithm>

namespace r = std::ranges;
namespace rv = std::ranges::views;

/*------------------------------------------------------------------------------------------------*/

namespace {

    cz::color interpolate_color(const cz::color& from, const cz::color& to, double t) {
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

}

cz::cyto_frame cz::to_cyto_frame(
    std::span<const point> pts, std::span<const color> colors) {

    auto polys = to_voronoi_polygons(pts);
    return rv::zip(polys, colors, pts) | rv::transform(
        [](const auto& v)->frame_cell {
            return {
                .shape = std::get<0>(v),
                .color = std::get<1>(v),
                .site = std::get<2>(v)
            };
        }
    ) | r::to<std::vector>();
}



cz::cyto_frame cz::interpolate_cyto_frames( std::span<const cz::frame_cell> from,
        std::span<const cz::frame_cell> to, double t) {

    if (from.size() != to.size() || from.empty()) {
        return {};
    }

    t = std::clamp(t, 0.0, 1.0);

    std::vector<cz::point> sites;
    std::vector<cz::color> colors;
    sites.reserve(from.size());
    colors.reserve(from.size());

    for (size_t i = 0; i < from.size(); ++i) {
        sites.push_back(cz::interpolate_point(from[i].site, to[i].site, t));
        colors.push_back(interpolate_color(from[i].color, to[i].color, t));
    }

    return cz::to_cyto_frame(sites, colors);
}