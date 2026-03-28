#pragma once

#include <vector>
#include "geometric_types.hpp"

namespace cz {

    voronoi_diagram construct_voronoi_diagram(
        const std::vector<point>& sites,
        const rect& bounds,
        double coordinate_scale = 1000.0
    );

}