#include "main_window.hpp"
#include "cytozoic_widget.hpp"
#include "cytozoic.hpp"
#include "voronoi.hpp"

#include <QEvent>
#include <QKeyEvent>
#include <QMenu>
#include <QMenuBar>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <ranges>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace r = std::ranges;
namespace rv = std::ranges::views;

namespace
{
    constexpr double k_lloyd_min_delta = 0.001;
    constexpr int k_max_iterations = 20;

    struct pending_animation
    {
        bool initialized = false;
        bool started = false;
        cz::cyto_frame from;
        cz::cyto_frame to;
    };

    std::unordered_map<cz::main_window*, pending_animation> g_pending_animations;

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
                throw std::runtime_error("make_scale_map: duplicate cell id.");
            }
        }

        return scale_by_id;
    }

    std::vector<double> make_from_scales(
        const cz::cyto_state& from_state,
        const std::unordered_map<cz::cell_id, double>& current_scale_by_id,
        const std::unordered_map<cz::cell_id, double>& next_scale_by_id)
    {
        std::vector<double> scales;
        scales.reserve(from_state.size());

        for (const auto& cell : from_state) {
            if (cell.phase == cz::life_stage::new_born) {
                auto it = next_scale_by_id.find(cell.id);
                if (it == next_scale_by_id.end()) {
                    throw std::runtime_error(
                        "make_from_scales: missing next scale for newborn cell."
                    );
                }

                scales.push_back(it->second);
            }
            else {
                auto it = current_scale_by_id.find(cell.id);
                if (it == current_scale_by_id.end()) {
                    throw std::runtime_error(
                        "make_from_scales: missing current scale for existing cell."
                    );
                }

                scales.push_back(it->second);
            }
        }

        return scales;
    }

    std::vector<double> make_to_scales(
        const cz::cyto_state& to_state,
        const std::unordered_map<cz::cell_id, double>& current_scale_by_id,
        const std::unordered_map<cz::cell_id, double>& next_scale_by_id)
    {
        std::vector<double> scales;
        scales.reserve(to_state.size());

        for (const auto& cell : to_state) {
            if (cell.phase == cz::life_stage::dying) {
                auto it = current_scale_by_id.find(cell.id);
                if (it == current_scale_by_id.end()) {
                    throw std::runtime_error(
                        "make_to_scales: missing current scale for dying cell."
                    );
                }

                scales.push_back(it->second);
            }
            else {
                auto it = next_scale_by_id.find(cell.id);
                if (it == next_scale_by_id.end()) {
                    throw std::runtime_error(
                        "make_to_scales: missing next scale for non-dying cell."
                    );
                }

                scales.push_back(it->second);
            }
        }

        return scales;
    }

    void relax_state_unweighted(cz::cyto_state& state)
    {
        if (state.empty()) {
            return;
        }

        const auto sites = state
            | rv::transform([](const cz::cell_state& cell) -> cz::point {
            return cell.site;
                })
            | r::to<std::vector>();

        const auto relaxed_sites = cz::perform_lloyd_relaxation(
            sites,
            k_lloyd_min_delta,
            k_max_iterations
        );

        if (relaxed_sites.size() != state.size()) {
            throw std::runtime_error(
                "relax_state_unweighted: relaxed site count mismatch."
            );
        }

        for (auto&& [cell, site] : rv::zip(state, relaxed_sites)) {
            cell.site = site;
        }
    }

    class key_start_filter : public QObject
    {
    public:
        explicit key_start_filter(cz::main_window* window)
            : QObject(window),
            window_(window)
        {}

    protected:
        bool eventFilter(QObject* watched, QEvent* event) override
        {
            Q_UNUSED(watched);

            if (event->type() != QEvent::KeyPress) {
                return QObject::eventFilter(watched, event);
            }

            auto it = g_pending_animations.find(window_);
            if (it == g_pending_animations.end()) {
                return QObject::eventFilter(watched, event);
            }

            pending_animation& anim = it->second;
            if (!anim.initialized || anim.started) {
                return QObject::eventFilter(watched, event);
            }

            anim.started = true;
            window_->centralWidget()->removeEventFilter(this);

            auto* canvas = qobject_cast<cz::cytozoic_widget*>(window_->centralWidget());
            if (canvas != nullptr) {
                canvas->start_transition(anim.from, anim.to);
            }

            return true;
        }

    private:
        cz::main_window* window_;
    };

} // namespace

cz::main_window::main_window(QWidget* parent)
    : QMainWindow(parent)
{
    setCentralWidget(canvas_ = new cytozoic_widget(this));
    create_menus();
    setWindowTitle(tr("cytozoic"));
    resize(1200, 1200);

    setFocusPolicy(Qt::StrongFocus);
    centralWidget()->setFocusPolicy(Qt::StrongFocus);
    centralWidget()->installEventFilter(new key_start_filter(this));
}

cz::main_window::~main_window()
{
    g_pending_animations.erase(this);
}

void cz::main_window::showEvent(QShowEvent* event)
{
    QMainWindow::showEvent(event);

    pending_animation& anim = g_pending_animations[this];
    if (anim.initialized) {
        return;
    }

    anim.initialized = true;

    cz::cell_id_source ids;
    auto current_state = cz::random_cyto_state(500, 4, ids);

    for (auto& cell : current_state) {
        cell.state = 0;
        cell.phase = cz::life_stage::normal;
    }

    // Pick some deletions by index into the current canonical state.
    const std::vector<std::size_t> delete_indices = { 5, 17, 42, 80, 111, 150, 222, 301 };
    std::vector<cz::cell_id> delete_list;
    delete_list.reserve(delete_indices.size());

    for (std::size_t index : delete_indices) {
        if (index >= current_state.size()) {
            throw std::runtime_error("showEvent: delete index out of range.");
        }

        delete_list.push_back(current_state[index].id);
        current_state[index].state = 2; // mark them visually in the initial frame
    }

    const auto current_sites = current_state
        | rv::transform([](const cz::cell_state& cell) -> cz::point {
        return cell.site;
            })
        | r::to<std::vector>();

    const auto diagram = cz::to_voronoi_diagram(current_sites);

    std::vector<cz::cell_state> add_list;
    add_list.reserve(8);

    constexpr std::size_t num_to_add = 8;

    if (!diagram.embedding.vertices.empty()) {
        const std::size_t step = std::max<std::size_t>(
            1,
            diagram.embedding.vertices.size() / num_to_add
        );

        for (std::size_t i = 0;
            i < diagram.embedding.vertices.size() && add_list.size() < num_to_add;
            i += step) {
            add_list.push_back(
                cz::cell_state{
                    .id = ids.acquire(),
                    .site = diagram.embedding.vertices[i],
                    .state = 1,
                    .phase = cz::life_stage::new_born
                }
            );
        }
    }

    // Build canonical next_state:
    // survivors from current_state (excluding deletions) + births normalized to normal,
    // then unweighted relaxation.
    const auto deletion_set = delete_list | r::to<std::unordered_set>();

    cz::cyto_state next_state;
    next_state.reserve(current_state.size() - delete_list.size() + add_list.size());

    for (const auto& cell : current_state) {
        if (deletion_set.contains(cell.id)) {
            continue;
        }

        auto survivor = cell;
        survivor.phase = cz::life_stage::normal;
        survivor.state = 0;
        next_state.push_back(survivor);
    }

    for (const auto& new_cell : add_list) {
        auto canonical_birth = new_cell;
        canonical_birth.phase = cz::life_stage::normal;
        next_state.push_back(canonical_birth);
    }

    relax_state_unweighted(next_state);

    const auto trans = cz::generate_transition(
        current_state,
        next_state,
        delete_list,
        add_list
    );

    const auto current_scale_by_id = make_scale_map(current_state);
    const auto next_scale_by_id = make_scale_map(next_state);

    const auto from_scales = make_from_scales(
        trans.from,
        current_scale_by_id,
        next_scale_by_id
    );

    const auto to_scales = make_to_scales(
        trans.to,
        current_scale_by_id,
        next_scale_by_id
    );

    cz::color_table palette = {
        {0, 255, 0},     // state 0
        {255, 0, 0},     // state 1 (births)
        {255, 255, 130}, // state 2 (cells chosen for deletion marker in initial frame)
        {0, 60, 200}
    };

    canvas_->set_show_cell_nuceli(false);

    const auto initial_frame = cz::to_cyto_frame(current_state, palette, {});
    anim.from = cz::to_cyto_frame(trans.from, palette, from_scales);
    anim.to = cz::to_cyto_frame(trans.to, palette, to_scales);

    canvas_->set(initial_frame);
    centralWidget()->setFocus();
}

void cz::main_window::create_menus()
{
    QMenu* file_menu = menuBar()->addMenu(tr("&File"));
    Q_UNUSED(file_menu);
}