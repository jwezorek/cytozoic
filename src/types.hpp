#pragma once

#include <cstdint>
#include <vector>
#include "vec2.hpp"

namespace cz {

    using point = vec2<double>;
    using polygon = std::vector<point>;

    struct rect {
        point min_point;
        point max_point;
    };

    struct voronoi_region {
        point site;
        polygon region;
        std::vector<size_t> neighbors;
    };

    using voronoi_diagram = std::vector<voronoi_region>;

    struct color {
        uint8_t r;
        uint8_t g;
        uint8_t b;
    };

    using color_table = std::vector<color>;

    using cell_id = uint64_t;

    struct frame_cell {
        polygon shape;
        color color;
        point site;
    };

    using cyto_frame = std::vector<frame_cell>;

    struct cell_state {
        cell_id id;
        point site;
        int8_t state;
    };

    using cyto_state = std::vector<cell_state>;

}