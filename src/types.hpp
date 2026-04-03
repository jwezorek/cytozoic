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

    struct color {
        uint8_t r;
        uint8_t g;
        uint8_t b;
    };

    struct cell {
        polygon shape;
        color color;
        point seed;
    };

    using cyto_frame = std::vector<cell>;
    using color_table = std::vector<color>;

    struct cyto_state {
        voronoi_diagram cells;
        std::vector<int8_t> states;
    };
    
}