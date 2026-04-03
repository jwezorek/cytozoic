#pragma once

#include <vector>
#include <span>
#include "types.hpp"

namespace cz {

    voronoi_diagram construct_voronoi_diagram(
        std::span<const point> sites, const rect& bounds
    );

    
}