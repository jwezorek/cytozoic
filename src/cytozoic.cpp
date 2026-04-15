#include "cytozoic.hpp"
#include "geometry.hpp"
#include "voronoi.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <random>
#include <ranges>
#include <stdexcept>
#include <tuple>
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

    enum class normal_scale_source
    {
        current,
        next
    };

    double polygon_area(const cz::polygon& poly)
    {
        if (poly.size() < 3) {
            return 0.0;
        }

        double area2 = 0.0;

        for (std::size_t i = 0; i < poly.size(); ++i) {
            const cz::point& a = poly[i];
            const cz::point& b = poly[(i + 1) % poly.size()];
            area2 += a.x * b.y - b.x * a.y;
        }

        return std::abs(area2) * 0.5;
    }

    double polygon_scale(const cz::polygon& poly)
    {
        constexpr double k_min_scale = 1e-9;
        return std::max(
            polygon_area(poly) / std::numbers::pi_v<double>,
            k_min_scale
        );
    }

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

        const bool weight_is_growing = from.weight < to.weight;

        const double weight_t = weight_is_growing
            ? 1.0 - (1.0 - t) * (1.0 - t) * (1.0 - t) * (1.0 - t) * (1.0 - t)
            : t;

        const double interpolated_weight =
            from.weight + (to.weight - from.weight) * weight_t;

        const double scale =
            weight_is_growing ? to.scale : from.scale;

        return {
            .id = from.id,
            .shape = {},
            .color = interpolate_color(from.color, to.color, t),
            .site = cz::interpolate_point(from.site, to.site, t),
            .weight = interpolated_weight,
            .scale = scale
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

    std::unordered_map<cz::cell_id, double> make_scale_map(const cz::cyto_state& state)
    {
        std::unordered_map<cz::cell_id, double> scale_by_id;
        scale_by_id.reserve(state.size());

        if (state.empty()) {
            return scale_by_id;
        }

        const auto sites = state
            | rv::transform([](const cz::cell_state& cell) -> cz::point {
            return cell.site;
                })
            | r::to<std::vector>();

        const auto polygons = cz::to_voronoi_polygons(sites);

        if (polygons.size() != state.size()) {
            throw std::runtime_error(
                "make_scale_map: polygon count did not match state size."
            );
        }

        for (std::size_t i = 0; i < state.size(); ++i) {
            auto [it, inserted] = scale_by_id.emplace(
                state[i].id,
                polygon_scale(polygons[i])
            );

            if (!inserted) {
                throw std::runtime_error(
                    "make_scale_map: duplicate cell id."
                );
            }
        }

        return scale_by_id;
    }

    std::vector<double> make_transition_scales(
        const cz::cyto_state& state,
        const std::unordered_map<cz::cell_id, double>& current_scale_by_id,
        const std::unordered_map<cz::cell_id, double>& next_scale_by_id,
        normal_scale_source normal_source)
    {
        std::vector<double> scales;
        scales.reserve(state.size());

        for (const auto& cell : state) {
            if (cell.phase == cz::life_stage::new_born) {
                auto it = next_scale_by_id.find(cell.id);
                if (it == next_scale_by_id.end()) {
                    throw std::runtime_error(
                        "make_transition_scales: missing next scale for newborn cell."
                    );
                }

                scales.push_back(it->second);
                continue;
            }

            if (cell.phase == cz::life_stage::dying) {
                auto it = current_scale_by_id.find(cell.id);
                if (it == current_scale_by_id.end()) {
                    throw std::runtime_error(
                        "make_transition_scales: missing current scale for dying cell."
                    );
                }

                scales.push_back(it->second);
                continue;
            }

            if (normal_source == normal_scale_source::current) {
                auto it = current_scale_by_id.find(cell.id);
                if (it == current_scale_by_id.end()) {
                    throw std::runtime_error(
                        "make_transition_scales: missing current scale for normal cell."
                    );
                }

                scales.push_back(it->second);
            }
            else {
                auto it = next_scale_by_id.find(cell.id);
                if (it == next_scale_by_id.end()) {
                    throw std::runtime_error(
                        "make_transition_scales: missing next scale for normal cell."
                    );
                }

                scales.push_back(it->second);
            }
        }

        return scales;
    }

    void relax_cyto_state(cz::cyto_state& state, const std::vector<double>& scales)
    {
        if (state.empty()) {
            return;
        }

        if (!scales.empty() && scales.size() != state.size()) {
            throw std::runtime_error(
                "relax_cyto_state: scale count did not match state size."
            );
        }

        std::vector<cz::point> relaxed_sites;

        if (scales.empty()) {
            const auto sites = state
                | rv::transform([](const cz::cell_state& cell) -> cz::point {
                return cell.site;
                    })
                | r::to<std::vector>();

            relaxed_sites = cz::perform_lloyd_relaxation(
                sites,
                k_lloyd_min_delta,
                k_max_iterations
            );
        }
        else {
            std::vector<cz::weighted_point> weighted_sites;
            weighted_sites.reserve(state.size());

            for (std::size_t i = 0; i < state.size(); ++i) {
                weighted_sites.push_back(
                    cz::weighted_point{
                        .pt = state[i].site,
                        .weight = relaxation_weight(state[i].phase),
                        .scale = scales[i]
                    }
                );
            }

            relaxed_sites = cz::perform_lloyd_relaxation(
                weighted_sites,
                k_lloyd_min_delta,
                k_max_iterations
            );
        }

        if (relaxed_sites.size() != state.size()) {
            throw std::runtime_error(
                "relax_cyto_state: relaxed site count mismatch."
            );
        }

        for (auto&& [cell, site] : rv::zip(state, relaxed_sites)) {
            cell.site = site;
        }
    }

    std::tuple<cz::cyto_state, std::vector<cz::cell_id>> apply_cell_table(
        const cell_graph& graph,
        const cz::state_table& cell_tbl,
        const cz::neighborhood_indexer& cell_indexer)
    {
        cz::cyto_state next;
        std::vector<cz::cell_id> delete_list;

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

            next.push_back(
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

        return { next, delete_list };
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
    const color_table& palette,
    const std::vector<double>& scales)
{
    if (state.empty()) {
        return {};
    }

    if (!scales.empty() && scales.size() != state.size()) {
        throw std::runtime_error(
            "to_cyto_frame: scale count did not match cyto_state size."
        );
    }

    std::vector<cz::polygon> polygons;
    std::vector<double> frame_scales;
    frame_scales.reserve(state.size());

    if (scales.empty()) {
        const auto sites = state
            | rv::transform([](const cz::cell_state& cell) -> cz::point {
            return cell.site;
                })
            | r::to<std::vector>();

        polygons = cz::to_voronoi_polygons(sites);

        if (polygons.size() != state.size()) {
            throw std::runtime_error(
                "to_cyto_frame: voronoi polygon count did not match cyto_state size."
            );
        }

        for (const auto& poly : polygons) {
            frame_scales.push_back(polygon_scale(poly));
        }
    }
    else {
        std::vector<cz::weighted_point> weighted_sites;
        weighted_sites.reserve(state.size());

        for (std::size_t i = 0; i < state.size(); ++i) {
            weighted_sites.push_back(
                cz::weighted_point{
                    .pt = state[i].site,
                    .weight = relaxation_weight(state[i].phase),
                    .scale = scales[i]
                }
            );
        }

        polygons = cz::to_voronoi_polygons(weighted_sites);

        if (polygons.size() != state.size()) {
            throw std::runtime_error(
                "to_cyto_frame: weighted voronoi polygon count did not match cyto_state size."
            );
        }

        frame_scales = scales;
    }

    return rv::zip(state, polygons, frame_scales)
        | rv::transform([&palette](const auto& v) -> cz::frame_cell {
        const auto& [cell, poly, scale] = v;

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
            .weight = relaxation_weight(cell.phase),
            .scale = scale
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
                .weight = cell.weight,
                .scale = cell.scale
            }
        );

        result.push_back(std::move(cell));
    }

    const auto polygons = cz::to_voronoi_polygons(sites);

    if (polygons.size() != result.size()) {
        throw std::runtime_error(
            "interpolate_cyto_frames: polygon count mismatch."
        );
    }

    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i].shape = polygons[i];
    }

    return result;
}

cz::cyto_state_transition cz::generate_transition(
    const cyto_state& state,
    const cyto_state& next_state,
    const std::vector<cell_id>& delete_cells,
    const std::vector<cell_state> add_cells)
{
    cyto_state from;
    cyto_state to;

    const auto deletion_set = delete_cells | r::to<std::unordered_set>();

    for (const auto& cell : state) {
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

    const auto current_scale_by_id = make_scale_map(state);
    const auto next_scale_by_id = make_scale_map(next_state);

    const auto to_scales = make_transition_scales(
        to,
        current_scale_by_id,
        next_scale_by_id,
        normal_scale_source::next
    );

    relax_cyto_state(to, to_scales);

    return { from, to };
}

cz::state_table_result cz::apply_state_tables(
    cell_id_source& id_source,
    const cyto_state& current_state,
    const state_table& cell_tbl,
    const neighborhood_indexer& cell_indexer,
    const state_table_row& vert_tbl,
    const neighborhood_indexer& vert_indexer,
    const color_table& palette)
{
    const auto snapshot = build_topology_snapshot(current_state);

    auto [next_state, delete_list] = apply_cell_table(
        snapshot.graph,
        cell_tbl,
        cell_indexer
    );

    const auto add_list = apply_vertex_table(
        id_source,
        snapshot.graph,
        snapshot.vertex_neighborhoods,
        vert_tbl,
        vert_indexer,
        cell_tbl.size()
    );

    for (const auto& new_cell : add_list) {
        auto canonical_cell = new_cell;
        canonical_cell.phase = life_stage::normal;
        next_state.push_back(canonical_cell);
    }

    relax_cyto_state(next_state, {});

    auto trans = generate_transition(current_state, next_state, delete_list, add_list);

    const auto current_scale_by_id = make_scale_map(current_state);
    const auto next_scale_by_id = make_scale_map(next_state);

    const auto from_scales = make_transition_scales(
        trans.from,
        current_scale_by_id,
        next_scale_by_id,
        normal_scale_source::current
    );

    const auto to_scales = make_transition_scales(
        trans.to,
        current_scale_by_id,
        next_scale_by_id,
        normal_scale_source::next
    );

    return {
        next_state,
        to_cyto_frame(trans.from, palette, from_scales),
        to_cyto_frame(trans.to, palette, to_scales)
    };
}