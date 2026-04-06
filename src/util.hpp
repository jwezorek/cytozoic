#pragma once

#include "types.hpp"
#include <vector>
#include <optional>
#include <span>

namespace cz {

    std::vector<point> random_points(size_t n,
        double bounds_wd = 1.0, double bounds_hgt = 1.0,
        std::optional<uint64_t> seed = {}
    );

    cyto_frame to_cyto_frame(
        std::span<const point> sites,
        std::span<const color> colors
    );

    //cyto_state blank_state(const voronoi_diagram& v);
    point centroid(std::span<const point> pts);
    double dot(const point& u, const point& v);
    double distance(const point& u, const point& v);
    double magnitude(const point& v);

    point interpolate_point(const point& from, const point& to, double t);
    color interpolate_color(const color& from, const color& to, double t);
}