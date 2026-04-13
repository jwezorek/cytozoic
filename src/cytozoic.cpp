#include "cytozoic.hpp"
#include "geometry.hpp"
#include "voronoi.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <ranges>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace r = std::ranges;
namespace rv = std::ranges::views;

/*------------------------------------------------------------------------------------------------*/

namespace
{
    constexpr auto k_small_weight = 0.01;
    constexpr auto k_lloyd_min_delta = 0.001;
    constexpr auto k_max_iterations = 20;

    struct cell
    {
        cz::cell_id id;
        cz::point site;
        int8_t state;
        std::vector<cz::cell_id> neighbors;
    };

    struct vertex_neighborhood
    {
        cz::point vertex;
        std::vector<cz::cell_id> cells;
    };

    using cell_graph = std::unordered_map<cz::cell_id, cell>;

    struct topology_snapshot
    {
        cell_graph graph;
        std::vector<vertex_neighborhood> vertex_neighborhoods;
    };

    cz::color interpolate_color(const cz::color& from, const cz::color& to, double t)
    {
        t = std::clamp(t, 0.0, 1.0);

        auto lerp_channel = [t](uint8_t a, uint8_t b) -> uint8_t {
            const double value =
                static_cast<double>(a) +
                (static_cast<double>(b) - static_cast<double>(a)) * t;
            return static_cast<uint8_t>(std::round(value));
            };

        return {
            lerp_channel(from.r, to.r),
            lerp_channel(from.g, to.g),
            lerp_channel(from.b, to.b)
        };
    }

    cz::frame_cell interpolate_frame_cell(
        const cz::frame_cell& from,
        const cz::frame_cell& to,
        double t)
    {
        if (from.id != to.id) {
            throw std::runtime_error(
                "interpolate_frame_cell: cell ids do not match."
            );
        }

        t = std::clamp(t, 0.0, 1.0);

        return {
            .id = from.id,
            .shape = {},
            .color = interpolate_color(from.color, to.color, t),
            .site = cz::interpolate_point(from.site, to.site, t),
            .weight = from.weight + (to.weight - from.weight) * t
        };
    }

    std::unordered_map<cz::cell_id, const cz::frame_cell*> make_frame_cell_map(
        std::span<const cz::frame_cell> frame)
    {
        std::unordered_map<cz::cell_id, const cz::frame_cell*> map;
        map.reserve(frame.size());

        for (const cz::frame_cell& cell : frame) {
            auto [it, inserted] = map.emplace(cell.id, &cell);
            if (!inserted) {
                throw std::runtime_error("duplicate frame_cell id.");
            }
        }

        return map;
    }

    int8_t random_cell_state(int num_states)
    {
        if (num_states <= 0) {
            throw std::invalid_argument("num_states must be positive.");
        }

        static thread_local std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dist(0, num_states - 1);
        return static_cast<int8_t>(dist(rng));
    }

    std::vector<cz::cell_state> collect_live_cells_and_release_dead(
        const cz::cyto_state& state,
        cz::cell_id_source& id_source)
    {
        std::vector<cz::cell_state> live_cells;
        live_cells.reserve(state.size());

        std::unordered_set<cz::cell_id> released_ids;
        released_ids.reserve(state.size());

        for (const auto& cell : state) {
            if (cell.phase == cz::life_stage::dying) {
                if (released_ids.insert(cell.id).second) {
                    id_source.release(cell.id);
                }
            }
            else {
                live_cells.push_back(cell);
            }
        }

        return live_cells;
    }

    topology_snapshot build_topology_snapshot(const std::vector<cz::cell_state>& live_cells)
    {
        topology_snapshot snapshot;

        if (live_cells.empty()) {
            return snapshot;
        }

        const auto sites = live_cells
            | rv::transform([](const cz::cell_state& cell) -> cz::point {
            return cell.site;
                })
            | r::to<std::vector>();

        const auto diagram = cz::to_voronoi_diagram(sites);

        if (diagram.graph.size() != live_cells.size()) {
            throw std::runtime_error(
                "build_topology_snapshot: neighbor count mismatch."
            );
        }

        if (diagram.embedding.cells.size() != live_cells.size()) {
            throw std::runtime_error(
                "build_topology_snapshot: embedding cell count mismatch."
            );
        }

        snapshot.graph.reserve(live_cells.size());

        for (std::size_t i = 0; i < live_cells.size(); ++i) {
            const auto& src = live_cells[i];

            std::vector<cz::cell_id> neighbor_ids;
            neighbor_ids.reserve(diagram.graph[i].size());

            for (std::size_t adj_index : diagram.graph[i]) {
                if (adj_index >= live_cells.size()) {
                    throw std::runtime_error(
                        "build_topology_snapshot: adjacency index out of range."
                    );
                }

                neighbor_ids.push_back(live_cells[adj_index].id);
            }

            auto [it, inserted] = snapshot.graph.emplace(
                src.id,
                cell{
                    .id = src.id,
                    .site = src.site,
                    .state = src.state,
                    .neighbors = std::move(neighbor_ids)
                }
            );

            if (!inserted) {
                throw std::runtime_error(
                    "build_topology_snapshot: duplicate cell id."
                );
            }
        }

        snapshot.vertex_neighborhoods.reserve(diagram.embedding.vertices.size());

        for (std::size_t vertex_index = 0;
            vertex_index < diagram.embedding.vertices.size();
            ++vertex_index) {
            std::vector<cz::cell_id> incident_cells;

            for (std::size_t cell_index = 0;
                cell_index < diagram.embedding.cells.size();
                ++cell_index) {
                const auto& polygon_vertices = diagram.embedding.cells[cell_index];

                if (r::find(polygon_vertices, vertex_index) != polygon_vertices.end()) {
                    incident_cells.push_back(live_cells[cell_index].id);
                }
            }

            if (incident_cells.size() < 3) {
                continue;
            }

            snapshot.vertex_neighborhoods.push_back(
                vertex_neighborhood{
                    .vertex = diagram.embedding.vertices[vertex_index],
                    .cells = std::move(incident_cells)
                }
            );
        }

        return snapshot;
    }

    double relaxation_weight(cz::life_stage phase)
    {
        switch (phase) {
        case cz::life_stage::normal:
            return 1.0;

        case cz::life_stage::new_born:
        case cz::life_stage::dying:
            return k_small_weight;
        }

        return 1.0;
    }

    void relax_cyto_state(cz::cyto_state& state)
    {
        if (state.empty()) {
            return;
        }

        const auto weighted_sites = state
            | rv::transform([](const cz::cell_state& cell) -> cz::weighted_point {
            return {
                .pt = cell.site,
                .weight = relaxation_weight(cell.phase)
            };
                })
            | r::to<std::vector>();

        const std::vector<cz::point> relaxed_sites = cz::perform_lloyd_relaxation(
            weighted_sites,
            k_lloyd_min_delta,
            k_max_iterations
        );

        if (relaxed_sites.size() != state.size()) {
            throw std::runtime_error(
                "relax_cyto_state: relaxed site count mismatch."
            );
        }

        for (auto&& [cell, site] : rv::zip(state, relaxed_sites)) {
            cell.site = site;
        }
    }

    void apply_cell_table(
        cz::cyto_state& current,
        std::vector<cz::cell_id>& delete_list,
        const cell_graph& graph,
        const cz::state_table& cell_tbl,
        const cz::neighborhood_indexer& cell_indexer)
    {
        for (const auto& [id, cell] : graph) {
            (void)id;

            if (cell.state < 0) {
                throw std::runtime_error(
                    "apply_cell_table: cell state cannot be negative."
                );
            }

            const auto state_index = static_cast<std::size_t>(cell.state);
            if (state_index >= cell_tbl.size()) {
                throw std::runtime_error(
                    "apply_cell_table: cell state is out of range for cell_tbl."
                );
            }

            auto neighborhood = cell.neighbors
                | rv::transform([&](cz::cell_id neighbor_id) -> int8_t {
                return graph.at(neighbor_id).state;
                    })
                | r::to<std::vector>();

            const std::size_t column =
                cell_indexer->column_index(neighborhood, cell_tbl.size());

            if (column >= cell_tbl[state_index].size()) {
                throw std::runtime_error(
                    "apply_cell_table: column index is out of range."
                );
            }

            const int8_t new_state = cell_tbl[state_index][column];

            current.push_back(
                cz::cell_state{
                    .id = cell.id,
                    .site = cell.site,
                    .state = new_state >= 0 ? new_state : cell.state,
                    .phase = cz::life_stage::normal
                }
            );

            if (new_state < 0) {
                delete_list.push_back(cell.id);
            }
        }
    }

    std::vector<cz::cell_state> apply_vertex_table(
        cz::cell_id_source& id_source,
        const cell_graph& graph,
        const std::vector<vertex_neighborhood>& neighborhoods,
        const cz::state_table_row& vert_tbl,
        const cz::neighborhood_indexer& vert_indexer,
        std::size_t num_states)
    {
        std::vector<cz::cell_state> add_list;

        for (const auto& neighborhood : neighborhoods) {
            const auto neighborhood_states = neighborhood.cells
                | rv::transform([&](cz::cell_id id) -> int8_t {
                return graph.at(id).state;
                    })
                | r::to<std::vector>();

            const std::size_t column =
                vert_indexer->column_index(neighborhood_states, num_states);

            if (column >= vert_tbl.size()) {
                throw std::runtime_error(
                    "apply_vertex_table: column index is out of range."
                );
            }

            const int8_t spawn_state = vert_tbl.at(column);

            if (spawn_state >= 0) {
                add_list.push_back(
                    cz::cell_state{
                        .id = id_source.acquire(),
                        .site = neighborhood.vertex,
                        .state = spawn_state,
                        .phase = cz::life_stage::new_born
                    }
                );
            }
        }

        return add_list;
    }

} // namespace

/*------------------------------------------------------------------------------------------------*/

cz::cyto_state cz::random_cyto_state(int num_cells, int num_states, cell_id_source& ids)
{
    if (num_cells <= 0) {
        return {};
    }

    const auto sites = perform_lloyd_relaxation(
        random_points(static_cast<std::size_t>(num_cells)),
        k_lloyd_min_delta,
        k_max_iterations
    );

    return sites
        | rv::transform([&ids, num_states](const auto& pt) -> cell_state {
        return {
            .id = ids.acquire(),
            .site = pt,
            .state = random_cell_state(num_states),
            .phase = life_stage::normal
        };
            })
        | r::to<std::vector>();
}

/*------------------------------------------------------------------------------------------------*/

cz::cell_id_source::cell_id_source()
    : next_id_(0)
{}

cz::cell_id cz::cell_id_source::acquire()
{
    if (!free_ids_.empty()) {
        const cell_id id = free_ids_.back();
        free_ids_.pop_back();
        return id;
    }

    return next_id_++;
}

void cz::cell_id_source::release(cell_id id)
{
    free_ids_.push_back(id);
}

void cz::cell_id_source::reset()
{
    next_id_ = 0;
    free_ids_.clear();
}

/*------------------------------------------------------------------------------------------------*/

cz::cyto_frame cz::to_cyto_frame(
    const cyto_state& state,
    const color_table& palette)
{
    if (state.empty()) {
        return {};
    }

    const auto sites = state
        | rv::transform([](const cell_state& cell) -> point {
        return cell.site;
            })
        | r::to<std::vector>();

    const auto diagram = to_voronoi_diagram(sites);

    if (diagram.polygons.size() != state.size()) {
        throw std::runtime_error(
            "to_cyto_frame: voronoi polygon count did not match cyto_state size."
        );
    }

    return rv::zip(state, diagram.polygons)
        | rv::transform([&palette](const auto& v) -> frame_cell {
        const auto& [cell, poly] = v;

        if (cell.state < 0) {
            throw std::out_of_range("cell state cannot be negative.");
        }

        const auto palette_index = static_cast<std::size_t>(cell.state);
        if (palette_index >= palette.size()) {
            throw std::out_of_range(
                "cell state is out of range for palette."
            );
        }

        return {
            .id = cell.id,
            .shape = poly,
            .color = palette[palette_index],
            .site = cell.site,
            .weight = (cell.phase == life_stage::normal) ? 1.0 : k_small_weight
        };
            })
        | r::to<std::vector>();
}

cz::cyto_frame cz::interpolate_cyto_frames(
    std::span<const cz::frame_cell> from,
    std::span<const cz::frame_cell> to,
    double t)
{
    if (from.size() != to.size()) {
        throw std::runtime_error(
            "interpolate_cyto_frames: frame cell counts differ."
        );
    }

    if (from.empty()) {
        return {};
    }

    t = std::clamp(t, 0.0, 1.0);

    const auto to_by_id = make_frame_cell_map(to);

    cz::cyto_frame result;
    result.reserve(from.size());

    std::vector<cz::weighted_point> sites;
    sites.reserve(from.size());

    for (const cz::frame_cell& from_cell : from) {
        auto it = to_by_id.find(from_cell.id);
        if (it == to_by_id.end()) {
            throw std::runtime_error(
                "interpolate_cyto_frames: frame cell ids do not match."
            );
        }

        cz::frame_cell cell = interpolate_frame_cell(from_cell, *it->second, t);

        sites.push_back(
            cz::weighted_point{
                .pt = cell.site,
                .weight = cell.weight
            }
        );

        result.push_back(std::move(cell));
    }

    const auto diagram = cz::to_voronoi_diagram(sites);

    if (diagram.polygons.size() != result.size()) {
        throw std::runtime_error(
            "interpolate_cyto_frames: polygon count mismatch."
        );
    }

    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i].shape = diagram.polygons[i];
    }

    return result;
}

cz::cyto_state_transition cz::generate_transition(
    const cyto_state& state,
    const std::vector<cell_id>& delete_cells,
    const std::vector<cell_state> add_cells)
{
    cyto_state from;
    cyto_state to;

    const auto deletion_set = delete_cells | r::to<std::unordered_set>();

    for (const auto& cell : state) {
        if (cell.phase == life_stage::dying) {
            continue;
        }

        auto new_cell = cell;
        new_cell.phase = life_stage::normal;
        from.push_back(new_cell);

        if (deletion_set.contains(cell.id)) {
            auto deleted = cell;
            deleted.phase = life_stage::dying;
            to.push_back(deleted);
        }
        else {
            to.push_back(new_cell);
        }
    }

    for (const auto& new_cell : add_cells) {
        auto addee = new_cell;
        addee.phase = life_stage::new_born;
        from.push_back(addee);

        addee.phase = life_stage::normal;
        to.push_back(addee);
    }

    relax_cyto_state(from);
    relax_cyto_state(to);

    return { from, to };
}

cz::state_table_result cz::apply_state_tables(
    cell_id_source& id_source,
    const cyto_state& state,
    const state_table& cell_tbl,
    const neighborhood_indexer& cell_indexer,
    const state_table_row& vert_tbl,
    const neighborhood_indexer& vert_indexer,
    const color_table& palette)
{
    const auto live_cells = collect_live_cells_and_release_dead(state, id_source);
    const auto snapshot = build_topology_snapshot(live_cells);

    cyto_state current;
    current.reserve(live_cells.size());

    std::vector<cell_id> delete_list;
    apply_cell_table(current, delete_list, snapshot.graph, cell_tbl, cell_indexer);

    const auto add_list = apply_vertex_table(
        id_source,
        snapshot.graph,
        snapshot.vertex_neighborhoods,
        vert_tbl,
        vert_indexer,
        cell_tbl.size()
    );

    auto trans = generate_transition(current, delete_list, add_list);

    return {
        trans.to,
        to_cyto_frame(trans.from, palette),
        to_cyto_frame(trans.to, palette)
    };
}