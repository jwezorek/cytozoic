#pragma once

#include <vector>
#include <span>
#include "types.hpp"

namespace cz {

    std::vector<std::vector<size_t>> to_voronoi_topology(
        std::span<const point> sites, const rect& bounds = { {0.0,0.0},{1.0,{1.0}} }
    );

    std::vector<polygon> to_voronoi_polygons(
        std::span<const point> sites, const rect& bounds = {{0.0,0.0},{1.0,{1.0}}}
    );

    std::vector<point> perform_lloyd_relaxation(
        std::span<const point> sites,
        double min_delta_thresh, int max_iterations, 
        const rect& bounds = { {0.0,0.0},{1.0,{1.0}} }
    );
}