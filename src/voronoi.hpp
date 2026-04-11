#include <span>
#include <vector>

#include "geometry.hpp"

namespace cz {

    std::vector<std::vector<size_t>> to_voronoi_topology(
        std::span<const point> sites,
        const rect& bounds = { {0.0, 0.0}, {1.0, 1.0} }
    );

    std::vector<polygon> to_voronoi_polygons(
        std::span<const point> sites,
        const std::vector<std::vector<size_t>>& graph,
        const rect& bounds = { {0.0, 0.0}, {1.0, 1.0} }
    );

    std::vector<polygon> to_voronoi_polygons(
        std::span<const point> sites,
        const rect& bounds = { {0.0, 0.0}, {1.0, 1.0} }
    );

    std::vector<point> perform_lloyd_relaxation(
        std::span<const point> sites,
        double min_delta_thresh,
        int max_iterations,
        const rect& bounds = { {0.0, 0.0}, {1.0, 1.0} }
    );

    struct weighted_point {
        point pt;
        double weight;
    };

    std::vector<std::vector<size_t>> to_voronoi_topology(
        std::span<const weighted_point> sites,
        const rect& bounds = { {0.0, 0.0}, {1.0, 1.0} }
    );

    std::vector<polygon> to_voronoi_polygons(
        std::span<const weighted_point> sites,
        const std::vector<std::vector<size_t>>& graph,
        const rect& bounds = { {0.0, 0.0}, {1.0, 1.0} }
    );

    std::vector<polygon> to_voronoi_polygons(
        std::span<const weighted_point> sites,
        const rect& bounds = { {0.0, 0.0}, {1.0, 1.0} }
    );

    std::vector<point> perform_lloyd_relaxation(
        std::span<const weighted_point> sites,
        double min_delta_thresh,
        int max_iterations,
        const rect& bounds = { {0.0, 0.0}, {1.0, 1.0} }
    );

}