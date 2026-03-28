#pragma once

#include <vector>

namespace cz {

    struct point {
        double x;
        double y;
    };

    using polygon = std::vector<point>;

    struct voronoi_cell {
        point site;
        polygon cell;
        std::vector<size_t> neighbors;
    };

    struct rect {
        point min_point;
        point max_point;
    };

    using voronoi_diagram = std::vector<voronoi_cell>;

}