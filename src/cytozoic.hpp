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
        cell_id id;
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

    class cell_id_source {
    public:
        cell_id_source();
        cell_id acquire();
        void release(cell_id id);
        void reset();

    private:
        cell_id next_id_;
        std::vector<cell_id> free_ids_;
    };

    cyto_state random_cyto_state(int num_cells, int num_states, cell_id_source& ids);

    cyto_frame to_cyto_frame(const cyto_state& state, const color_table& palette);

    cz::cyto_frame interpolate_cyto_frames(
        std::span<const cz::frame_cell> from,
        std::span<const cz::frame_cell> to,
        double t
    );

}