#pragma once

#include <span>
#include <vector>

#include "geometry.hpp"

/*------------------------------------------------------------------------------------------------*/

namespace cz {

    using polygon_vertex_indices = std::vector<size_t>;

    struct voronoi_embedding {
        std::vector<point> vertices;
        std::vector<polygon_vertex_indices> cells;
    };

    struct voronoi_diagram {
        std::vector<std::vector<size_t>> graph;
        std::vector<polygon> polygons;
        voronoi_embedding embedding;
    };

    struct weighted_point {
        point pt;
        double weight;
        double scale;
    };

    std::vector<std::vector<size_t>> to_voronoi_topology(
        std::span<const point> sites,
        const rect& bounds = { {0.0, 0.0}, {1.0, 1.0} }
    );

    voronoi_embedding to_voronoi_embedding(
        std::span<const point> sites,
        const rect& bounds = { {0.0, 0.0}, {1.0, 1.0} }
    );

    std::vector<point> perform_lloyd_relaxation(
        std::span<const point> sites,
        double min_delta_thresh,
        int max_iterations,
        const rect& bounds = { {0.0, 0.0}, {1.0, 1.0} }
    );

    std::vector<std::vector<size_t>> to_voronoi_topology(
        std::span<const weighted_point> sites,
        const rect& bounds = { {0.0, 0.0}, {1.0, 1.0} }
    );

    voronoi_embedding to_voronoi_embedding(
        std::span<const weighted_point> sites,
        const rect& bounds = { {0.0, 0.0}, {1.0, 1.0} }
    );

    std::vector<point> perform_lloyd_relaxation(
        std::span<const weighted_point> sites,
        double min_delta_thresh,
        int max_iterations,
        const rect& bounds = { {0.0, 0.0}, {1.0, 1.0} }
    );

    voronoi_diagram to_voronoi_diagram(
        std::span<const point> sites,
        const rect& bounds = { {0.0, 0.0}, {1.0, 1.0} });

    voronoi_diagram to_voronoi_diagram(
        std::span<const weighted_point> sites,
        const rect& bounds = { {0.0, 0.0}, {1.0, 1.0} }
    );


} // namespace cz