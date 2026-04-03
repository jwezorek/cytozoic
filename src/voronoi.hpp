#pragma once

#include <vector>
#include "types.hpp"

namespace cz {

    voronoi_diagram construct_voronoi_diagram(
        const std::vector<point>& sites, const rect& bounds
    );

}