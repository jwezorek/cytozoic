#pragma once

#include "geometry.hpp"
#include "neighborhood_indexer.hpp"
#include <optional>
#include <variant>

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
        double weight;
        double scale;
    };

    using cyto_frame = std::vector<frame_cell>;

    enum class life_stage {
        new_born,
        normal,
        dying
    };

    struct cell_state {
        cell_id id;
        point site;
        int8_t state;
        life_stage phase;
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

    cyto_state random_cyto_state(
        int num_cells, int num_states, cell_id_source& ids,
        std::optional<uint64_t> seed = {}
    );

    cyto_state initial_cyto_state(
        const struct cyto_params& params,
        cell_id_source& ids
    );

    cyto_frame to_cyto_frame(
        const cyto_state& state,
        const color_table& palette,
        const std::vector<double>& scales);

    cyto_frame interpolate_cyto_frames(
        std::span<const frame_cell> from,
        std::span<const frame_cell> to,
        double t
    );

    std::vector<cell_id> deleted_cell_ids(std::span<const frame_cell> frame);

    cyto_frame remove_deleted_cells(std::span<const frame_cell> frame);

    using state_table_row = std::vector<int8_t>;
    using state_table = std::vector<state_table_row>;

    struct cyto_state_transition {
        cyto_state from;
        cyto_state to;
    };

    cyto_state_transition generate_transition(
        const cyto_state& state,
        const cyto_state& next_state,
        const std::unordered_set<cz::cell_id>& deleted_set,
        const std::vector<cell_state> add_cells
    );

    struct state_table_result {
        cyto_state next_state;
        cyto_frame anim_start;
        cyto_frame anim_end;
    };

    enum class center_type {
        incircle,
        johnson_ellipse,
        center_of_mass
    };

    struct cell_based_birth {
        center_type spawn_site;
        state_table state_table;
    };

    struct vertex_based_birth {
        neighborhood_indexer vertex_indexer;
        state_table_row vertex_table;
    };

    using birth_parameters = std::variant<vertex_based_birth, cell_based_birth>;

    struct cyto_params {

        neighborhood_indexer cell_indexer;
        state_table cell_state_table;
        int num_states;
        int num_initial_cells;
        std::vector<double> initial_state_density;
        color_table palette;
        birth_parameters birth_params;

        cyto_params();

    };

    state_table_result apply_state_tables_animated(
        cell_id_source& id_source,
        const cyto_state& state,
        const state_table& cell_tbl,
        const neighborhood_indexer& cell_indexer,
        const birth_parameters& birth,
        const color_table& palette
    );

    cyto_state apply_state_tables_quick(
        cell_id_source& id_source,
        const cyto_state& state,
        const state_table& cell_tbl,
        const neighborhood_indexer& cell_indexer,
        const birth_parameters& birth
    );

}