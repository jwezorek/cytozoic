#pragma once

#include "geometry.hpp"

/*------------------------------------------------------------------------------------------------*/

namespace cz {

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

    cyto_frame to_cyto_frame(
        std::span<const point> sites,
        std::span<const color> colors
    );

    cz::cyto_frame interpolate_cyto_frames(
        std::span<const cz::frame_cell> from,
        std::span<const cz::frame_cell> to,
        double t
    );

}