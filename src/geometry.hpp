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

    enum class center_type {
        area_centroid = 0,
        mean_vertex,
        chebyshev,
        geometric_median
    };
    constexpr size_t k_num_center_types = 4;

    std::string center_type_to_string(center_type ct);
    center_type center_type_from_string(const std::string& str);
    std::vector<center_type> center_types();

    point center(std::span<const point>, center_type which_center = center_type::area_centroid);
    double dot(const point& u, const point& v);
    double distance(const point& u, const point& v);
    double magnitude(const point& v);

    point interpolate_point(const point& from, const point& to, double t);
}