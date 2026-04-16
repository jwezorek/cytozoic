#pragma once

#include <vector>
#include <optional>
#include <span>
#include "vec2.hpp"

namespace cz {

    using point = vec2<double>;
    using polygon = std::vector<point>;

    struct rect {
        point min_point;
        point max_point;
    };

    std::vector<point> random_points(size_t n,
        std::optional<uint64_t> seed = {},
        double bounds_wd = 1.0, 
        double bounds_hgt = 1.0
    );


    point centroid(std::span<const point> pts);
    double dot(const point& u, const point& v);
    double distance(const point& u, const point& v);
    double magnitude(const point& v);

    point interpolate_point(const point& from, const point& to, double t);
}