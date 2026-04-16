#include "geometry.hpp"
#include "voronoi.hpp"
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
        std::optional<uint64_t> seed,
        double bounds_wd,
        double bounds_hgt) {

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

cz::point cz::centroid(std::span<const point> pts) {
    if (pts.empty()) {
        throw std::runtime_error(
            "attempted to find centroid of an empty point set."
        );
    }
    auto sum = r::fold_left(pts, point{ 0.0,0.0 }, std::plus<>{});
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


