#include "cytozoic.hpp"
#include "voronoi.hpp"
#include "geometry.hpp"
#include <boost/container_hash/hash.hpp>
#include <algorithm>
#include <cmath>
#include <ranges>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <random>

namespace r = std::ranges;
namespace rv = std::ranges::views;

/*------------------------------------------------------------------------------------------------*/

namespace {

    constexpr auto k_small_weight = 0.00001;
    constexpr auto k_lloyd_min_delta = 0.001;
    constexpr auto k_max_iterations = 20;

    struct point_key {
        static constexpr double k_scale = 1000000.0;

        std::int64_t x;
        std::int64_t y;

        point_key() = default;

        explicit point_key(const cz::point& pt)
            : x(static_cast<std::int64_t>(std::llround(pt.x* k_scale))),
            y(static_cast<std::int64_t>(std::llround(pt.y* k_scale)))
        {}

        friend bool operator==(const point_key&, const point_key&) = default;
    };

    struct point_key_hasher {
        std::size_t operator()(const point_key& key) const noexcept
        {
            std::size_t seed = 0;
            boost::hash_combine(seed, key.x);
            boost::hash_combine(seed, key.y);
            return seed;
        }
    };

    template<typename V>
    using point_map = std::unordered_map<point_key, V, point_key_hasher>;

    cz::color interpolate_color(const cz::color& from, const cz::color& to, double t) {

        t = std::clamp(t, 0.0, 1.0);

        auto lerp_channel = [t](uint8_t a, uint8_t b) -> uint8_t {
            const double value = static_cast<double>(a)
                + (static_cast<double>(b) - static_cast<double>(a)) * t;
            return static_cast<uint8_t>(std::round(value));
            };

        return {
            lerp_channel(from.r, to.r),
            lerp_channel(from.g, to.g),
            lerp_channel(from.b, to.b)
        };
    }

    cz::frame_cell interpolate_frame_cell( const cz::frame_cell& from, const cz::frame_cell& to,
            double t) {
        if (from.id != to.id) {
            throw std::runtime_error("interpolate_frame_cell: cell ids do not match.");
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
            std::span<const cz::frame_cell> frame)  {
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

    int8_t random_cell_state(int num_states) {
        if (num_states <= 0) {
            throw std::invalid_argument("num_states must be positive.");
        }

        static thread_local std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dist(0, num_states - 1);
        return static_cast<int8_t>(dist(rng));
    }

    struct cell {
        cz::cell_id id;
        cz::point site;
        int8_t state;
        cz::polygon shape;
        std::vector<cz::cell_id> neighbors;
    };

    using cell_graph = std::unordered_map<cz::cell_id, cell>;

    cell_graph build_live_cell_graph(const cz::cyto_state& state) {
        cell_graph graph;
        auto live_cells = state | rv::filter(
            [](const auto& c) {
                // negative ids indicate dead cells...
                return c.state >= 0;
            }
        ) | r::to<std::vector>();
        auto sites = live_cells | rv::transform(
            [](const auto& c) {
                return c.site;
            }
        ) | r::to<std::vector>();
        auto neighbors = to_voronoi_topology(sites);
        auto polygons = to_voronoi_polygons(sites, neighbors);
        
        return rv::zip(live_cells, neighbors, polygons) | rv::transform(
            [&](const auto& triple)->cell_graph::value_type {
                const auto& [cell, adj, shape] = triple;
                auto adj_list_as_ids = adj | rv::transform(
                    [&](size_t adj_index)->cz::cell_id {
                        return live_cells.at(adj_index).id;
                    }
                ) | r::to<std::vector>();
                return {
                    cell.id,
                    ::cell{
                        .id = cell.id,
                        .site = cell.site,
                        .state = cell.state,
                        .shape = shape,
                        .neighbors = adj_list_as_ids
                    }
                };
            }
        ) | r::to< cell_graph>();
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

    void relax_cyto_state(cz::cyto_state& state) {
        if (state.empty()) {
            return;
        }

        const auto weighted_sites = state
            | rv::transform(
                [](const cz::cell_state& cell) -> cz::weighted_point {
                    return {
                        .pt = cell.site,
                        .weight = relaxation_weight(cell.phase)
                    };
                }
            )
            | r::to<std::vector>();

        const std::vector<cz::point> relaxed_sites = cz::perform_lloyd_relaxation(
            weighted_sites,
            k_lloyd_min_delta,
            k_max_iterations
        );  

        for (const auto& [cell, site] : rv::zip(state, relaxed_sites)) {
            cell.site = site;
        }
    }
} 

cz::cyto_state cz::random_cyto_state(int num_cells, int num_states, cell_id_source& ids) {
    auto sites = perform_lloyd_relaxation(
        random_points(num_cells),
        k_lloyd_min_delta, k_max_iterations
    );
    return sites | rv::transform(
        [&ids, num_states](const auto& pt)->cell_state {
            return {
                .id = ids.acquire(),
                .site = pt,
                .state = random_cell_state(num_states)
            };
        }
    ) | r::to<std::vector>();
}

/*------------------------------------------------------------------------------------------------*/

cz::cell_id_source::cell_id_source() : next_id_(0) {

}

cz::cell_id cz::cell_id_source::acquire() {
    if (!free_ids_.empty()) {
        cell_id id = free_ids_.back();
        free_ids_.pop_back();
        return id;
    }

    return next_id_++;
}

void cz::cell_id_source::release(cell_id id) {
    free_ids_.push_back(id);
}

void cz::cell_id_source::reset() {
    next_id_ = 0;
    free_ids_.clear();
}

/*------------------------------------------------------------------------------------------------*/


cz::cyto_frame cz::to_cyto_frame(const cyto_state& state, const color_table& palette) {
    if (state.empty()) {
        return {};
    }

    auto sites = state | rv::transform(
            [](const cell_state& cell) -> point {
                return cell.site;
            }
        ) | r::to<std::vector>();

    auto polygons = to_voronoi_polygons(sites);

    if (polygons.size() != state.size()) {
        throw std::runtime_error("voronoi polygon count did not match cyto_state size.");
    }

    return rv::zip(state, polygons) | rv::transform(
        [&palette](const auto& v) -> frame_cell {
            const auto& [cell, poly] = v;

            if (cell.state < 0) {
                throw std::out_of_range("cell state cannot be negative.");
            }

            const auto palette_index = static_cast<size_t>(cell.state);
            if (palette_index >= palette.size()) {
                throw std::out_of_range("cell state is out of range for palette.");
            }

            return {
                .id = cell.id,
                .shape = poly,
                .color = palette[palette_index],
                .site = cell.site,
                .weight = (cell.phase == life_stage::normal) ? 1.0 : k_small_weight
            };
        }
    ) | r::to<std::vector>();
}

cz::cyto_frame cz::interpolate_cyto_frames( std::span<const cz::frame_cell> from,
        std::span<const cz::frame_cell> to, double t) 
{
    if (from.size() != to.size()) {
        throw std::runtime_error("interpolate_cyto_frames: frame cell counts differ.");
    }

    if (from.empty()) {
        return {};
    }

    t = std::clamp(t, 0.0, 1.0);

    const auto to_by_id = make_frame_cell_map(to);

    cz::cyto_frame result;
    result.reserve(from.size());

    std::vector<cz::point> sites;
    sites.reserve(from.size());

    for (const cz::frame_cell& from_cell : from) {
        auto it = to_by_id.find(from_cell.id);
        if (it == to_by_id.end()) {
            throw std::runtime_error("interpolate_cyto_frames: frame cell ids do not match.");
        }

        cz::frame_cell cell = interpolate_frame_cell(from_cell, *it->second, t);
        sites.push_back(cell.site);
        result.push_back(std::move(cell));
    }

    const std::vector<cz::polygon> polygons = cz::to_voronoi_polygons(sites);

    if (polygons.size() != result.size()) {
        throw std::runtime_error("interpolate_cyto_frames: polygon count mismatch.");
    }

    for (size_t i = 0; i < result.size(); ++i) {
        result[i].shape = polygons[i];
    }

    return result;
}

cz::cyto_state_transition cz::generate_transition(
        const cyto_state& state, const std::vector<cell_id>& delete_cells, 
        const std::vector<cell_state> add_cells) {

    cyto_state from;
    cyto_state to;

    auto deletion_set = delete_cells | r::to<std::unordered_set>();

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
        } else {
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
        const cyto_state& state, const state_table& cell_tbl, 
        const neighborhood_indexer& cell_indexer, const state_table_row& vert_tbl, 
        const neighborhood_indexer& vert_indexer ) {

    auto cell_graph = build_live_cell_graph(state);

    cyto_state from_state;
    cyto_state to_state;

    for (const auto& cell : cell_graph | rv::values) {
        auto neighborhood = cell.neighbors | rv::transform(
                [&](auto neighbor) {
                    return cell_graph.at(neighbor).state;
                }
            ) | r::to<std::vector>();
        auto column = cell_indexer->column_index(neighborhood, cell_tbl.size());
    }

    return {};
}
