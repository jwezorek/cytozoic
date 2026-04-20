#include "geometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>

namespace r = std::ranges;
namespace rv = std::ranges::views;

/*------------------------------------------------------------------------------------------------*/

namespace {

    constexpr double k_epsilon = 1e-12;
    constexpr int k_weiszfeld_max_iterations = 256;
    constexpr double k_weiszfeld_tolerance = 1e-9;

    struct half_plane_constraint {
        double a;
        double b;
        double c;
        double d;
    };

    double signed_area_twice(std::span<const cz::point> poly)
    {
        if (poly.size() < 3) {
            return 0.0;
        }

        double area_twice = 0.0;

        for (std::size_t i = 0; i < poly.size(); ++i) {
            const cz::point& p = poly[i];
            const cz::point& q = poly[(i + 1) % poly.size()];

            area_twice += p.x * q.y - q.x * p.y;
        }

        return area_twice;
    }

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

    cz::point mean_vertex(std::span<const cz::point> pts) {
        if (pts.empty()) {
            throw std::runtime_error(
                "mean_vertex: attempted to average an empty point set."
            );
        }

        auto sum = r::fold_left(pts, cz::point{ 0.0, 0.0 }, std::plus<>{});
        const double n = static_cast<double>(pts.size());

        return {
            .x = sum.x / n,
            .y = sum.y / n
        };
    }

    cz::point area_centroid(std::span<const cz::point> poly) {
        if (poly.size() < 3) {
            throw std::runtime_error(
                "area_centroid: polygon must contain at least 3 vertices."
            );
        }

        double area_twice = 0.0;
        double cx_times_6a = 0.0;
        double cy_times_6a = 0.0;

        for (std::size_t i = 0; i < poly.size(); ++i) {
            const cz::point& a = poly[i];
            const cz::point& b = poly[(i + 1) % poly.size()];

            const double cross = a.x * b.y - b.x * a.y;

            area_twice += cross;
            cx_times_6a += (a.x + b.x) * cross;
            cy_times_6a += (a.y + b.y) * cross;
        }

        if (std::abs(area_twice) < k_epsilon) {
            throw std::runtime_error(
                "area_centroid: degenerate polygon has zero area."
            );
        }

        return {
            .x = cx_times_6a / (3.0 * area_twice),
            .y = cy_times_6a / (3.0 * area_twice)
        };
    }

    bool solve_3x3(
        const double m[3][3],
        const double rhs[3],
        double out[3])
    {
        double aug[3][4] = {
            { m[0][0], m[0][1], m[0][2], rhs[0] },
            { m[1][0], m[1][1], m[1][2], rhs[1] },
            { m[2][0], m[2][1], m[2][2], rhs[2] }
        };

        for (int col = 0; col < 3; ++col) {
            int pivot_row = col;

            for (int row = col + 1; row < 3; ++row) {
                if (std::abs(aug[row][col]) > std::abs(aug[pivot_row][col])) {
                    pivot_row = row;
                }
            }

            if (std::abs(aug[pivot_row][col]) < k_epsilon) {
                return false;
            }

            if (pivot_row != col) {
                for (int k = col; k < 4; ++k) {
                    std::swap(aug[col][k], aug[pivot_row][k]);
                }
            }

            const double pivot = aug[col][col];
            for (int k = col; k < 4; ++k) {
                aug[col][k] /= pivot;
            }

            for (int row = 0; row < 3; ++row) {
                if (row == col) {
                    continue;
                }

                const double factor = aug[row][col];
                for (int k = col; k < 4; ++k) {
                    aug[row][k] -= factor * aug[col][k];
                }
            }
        }

        out[0] = aug[0][3];
        out[1] = aug[1][3];
        out[2] = aug[2][3];
        return true;
    }

    std::vector<half_plane_constraint> make_chebyshev_constraints(
        std::span<const cz::point> poly)
    {
        if (poly.size() < 3) {
            throw std::runtime_error(
                "chebyshev_center: polygon must contain at least 3 vertices."
            );
        }

        const double area2 = signed_area_twice(poly);
        if (std::abs(area2) < k_epsilon) {
            throw std::runtime_error(
                "chebyshev_center: degenerate polygon has zero area."
            );
        }

        const bool ccw = area2 > 0.0;

        std::vector<half_plane_constraint> constraints;
        constraints.reserve(poly.size() + 1);

        for (std::size_t i = 0; i < poly.size(); ++i) {
            const cz::point& a = poly[i];
            const cz::point& b = poly[(i + 1) % poly.size()];

            const cz::point edge = b - a;
            const double length = std::sqrt(edge.x * edge.x + edge.y * edge.y);

            if (length < k_epsilon) {
                throw std::runtime_error(
                    "chebyshev_center: polygon contains a zero-length edge."
                );
            }

            cz::point inward_normal;
            if (ccw) {
                inward_normal = {
                    .x = -edge.y / length,
                    .y = edge.x / length
                };
            }
            else {
                inward_normal = {
                    .x = edge.y / length,
                    .y = -edge.x / length
                };
            }

            constraints.push_back(
                half_plane_constraint{
                    .a = -inward_normal.x,
                    .b = -inward_normal.y,
                    .c = 1.0,
                    .d = -(inward_normal.x * a.x + inward_normal.y * a.y)
                }
            );
        }

        constraints.push_back(
            half_plane_constraint{
                .a = 0.0,
                .b = 0.0,
                .c = -1.0,
                .d = 0.0
            }
        );

        return constraints;
    }

    bool satisfies_constraint(
        const half_plane_constraint& constraint,
        double x,
        double y,
        double radius)
    {
        const double lhs =
            constraint.a * x +
            constraint.b * y +
            constraint.c * radius;

        return lhs <= constraint.d + 1e-9;
    }

    cz::point chebyshev_center(std::span<const cz::point> poly) {
        const auto constraints = make_chebyshev_constraints(poly);

        double best_radius = -std::numeric_limits<double>::infinity();
        std::optional<cz::point> best_center;

        for (std::size_t i = 0; i < constraints.size(); ++i) {
            for (std::size_t j = i + 1; j < constraints.size(); ++j) {
                for (std::size_t k = j + 1; k < constraints.size(); ++k) {
                    const double matrix[3][3] = {
                        {
                            constraints[i].a,
                            constraints[i].b,
                            constraints[i].c
                        },
                        {
                            constraints[j].a,
                            constraints[j].b,
                            constraints[j].c
                        },
                        {
                            constraints[k].a,
                            constraints[k].b,
                            constraints[k].c
                        }
                    };

                    const double rhs[3] = {
                        constraints[i].d,
                        constraints[j].d,
                        constraints[k].d
                    };

                    double solution[3];
                    if (!solve_3x3(matrix, rhs, solution)) {
                        continue;
                    }

                    const double x = solution[0];
                    const double y = solution[1];
                    const double radius = solution[2];

                    if (radius < -1e-9) {
                        continue;
                    }

                    bool feasible = true;
                    for (const auto& constraint : constraints) {
                        if (!satisfies_constraint(constraint, x, y, radius)) {
                            feasible = false;
                            break;
                        }
                    }

                    if (!feasible) {
                        continue;
                    }

                    if (radius > best_radius) {
                        best_radius = radius;
                        best_center = cz::point{ .x = x, .y = y };
                    }
                }
            }
        }

        if (!best_center.has_value()) {
            throw std::runtime_error(
                "chebyshev_center: failed to solve for polygon center."
            );
        }

        return *best_center;
    }

    cz::point geometric_median(std::span<const cz::point> poly) {
        if (poly.empty()) {
            throw std::runtime_error(
                "geometric_median: polygon must contain at least 1 vertex."
            );
        }

        if (poly.size() == 1) {
            return poly.front();
        }

        cz::point current = mean_vertex(poly);

        for (int iter = 0; iter < k_weiszfeld_max_iterations; ++iter) {
            double weight_sum = 0.0;
            cz::point weighted_sum{ 0.0, 0.0 };

            for (const cz::point& vertex : poly) {
                const double dist = cz::distance(current, vertex);

                if (dist < k_weiszfeld_tolerance) {
                    return vertex;
                }

                const double weight = 1.0 / dist;
                weight_sum += weight;
                weighted_sum += vertex * weight;
            }

            if (weight_sum < k_epsilon) {
                break;
            }

            const cz::point next = weighted_sum / weight_sum;

            if (cz::distance(current, next) < k_weiszfeld_tolerance) {
                return next;
            }

            current = next;
        }

        return current;
    }

} // namespace

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

std::string cz::center_type_to_string(center_type ct) {
    switch (ct) {
        case center_type::area_centroid:
            return "area centroid";

        case center_type::mean_vertex:
            return "mean vertex";

        case center_type::chebyshev:
            return "chebyshev center";

        case center_type::geometric_median:
            return "geometric median";
    }

    std::unreachable();
}

cz::center_type cz::center_type_from_string(const std::string& value) {
    if (value == "chebyshev center") {
        return cz::center_type::chebyshev;
    }

    if (value == "area centroid") {
        return cz::center_type::area_centroid;
    }

    if (value == "mean vertex") {
        return cz::center_type::mean_vertex;
    }

    if (value == "geometric median") {
        return cz::center_type::geometric_median;
    }

    throw std::runtime_error("unrecognized center_type string.");
}

cz::point cz::center(std::span<const point> poly, center_type ct) {

    switch (ct) {
        case center_type::area_centroid:
            return area_centroid(poly);

        case center_type::mean_vertex:
            return mean_vertex(poly);

        case center_type::chebyshev:
            return chebyshev_center(poly);

        case center_type::geometric_median:
            return geometric_median(poly);
    }

    std::unreachable();
}

double cz::dot(const point& u, const point& v) {
    return u.x * v.x + u.y * v.y;
}

double cz::distance(const point& u, const point& v) {
    auto diff = v - u;
    return std::sqrt(dot(diff, diff));
}

double cz::magnitude(const point& u) {
    return std::sqrt(dot(u, u));
}

cz::point cz::interpolate_point(const point& from, const point& to, double t) {
    t = std::clamp(t, 0.0, 1.0);
    return from + (to - from) * t;
}