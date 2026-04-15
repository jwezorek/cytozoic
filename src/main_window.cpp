#include "main_window.hpp"

#include "cytozoic.hpp"
#include "voronoi.hpp"

#include <QDialog>
#include <QDialogButtonBox>
#include <QFile>
#include <QFileDialog>
#include <QJsonDocument>
#include <QJsonObject>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QVBoxLayout>

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

    std::unordered_map<cz::cell_id, double> make_scale_map(
        const cz::cyto_state& state)
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

    void show_empty_modal_dialog(QWidget* parent, const QString& title)
    {
        QDialog dialog(parent);
        dialog.setWindowTitle(title);
        dialog.setModal(true);

        auto* layout = new QVBoxLayout(&dialog);
        auto* buttons = new QDialogButtonBox(QDialogButtonBox::Ok, &dialog);

        QObject::connect(
            buttons,
            &QDialogButtonBox::accepted,
            &dialog,
            &QDialog::accept
        );

        layout->addWidget(buttons);

        dialog.exec();
    }

} // namespace

cz::main_window::main_window(QWidget* parent)
    : QMainWindow(parent)
{
    setCentralWidget(canvas_ = new cytozoic_widget(this));
    create_menus();
    setWindowTitle(tr("cytozoic"));
    resize(1200, 1200);
}

cz::main_window::~main_window() = default;

void cz::main_window::set_params(const cyto_params& params) {
    params_ = params;
}

cz::cyto_params cz::main_window::get_params() const {
    return params_;
}

void cz::main_window::create_menus()
{
    QMenu* file_menu = menuBar()->addMenu(tr("&File"));

    QAction* open_ruleset_action =
        file_menu->addAction(tr("&Open Ruleset..."));
    QObject::connect(
        open_ruleset_action,
        &QAction::triggered,
        this,
        &cz::main_window::open_ruleset
    );

    QAction* save_ruleset_action =
        file_menu->addAction(tr("&Save Ruleset"));
    QObject::connect(
        save_ruleset_action,
        &QAction::triggered,
        this,
        &cz::main_window::save_ruleset
    );

    QAction* save_ruleset_as_action =
        file_menu->addAction(tr("Save Ruleset &As..."));
    QObject::connect(
        save_ruleset_as_action,
        &QAction::triggered,
        this,
        &cz::main_window::save_ruleset_as
    );

    file_menu->addSeparator();

    QAction* exit_action = file_menu->addAction(tr("E&xit"));
    QObject::connect(
        exit_action,
        &QAction::triggered,
        this,
        &QWidget::close
    );

    QMenu* cytozoic_menu = menuBar()->addMenu(tr("&Cytozoic"));

    QAction* rules_action =
        cytozoic_menu->addAction(tr("View/Edit Current &Rules..."));
    QObject::connect(
        rules_action,
        &QAction::triggered,
        this,
        &cz::main_window::view_edit_current_rules
    );

    cytozoic_menu->addSeparator();

    QAction* run_simulation_action =
        cytozoic_menu->addAction(tr("&Run Simulation"));
    QObject::connect(
        run_simulation_action,
        &QAction::triggered,
        this,
        &cz::main_window::run_simulation
    );

    QAction* debug_action =
        cytozoic_menu->addAction(tr("&Debug"));
    QObject::connect(
        debug_action,
        &QAction::triggered,
        this,
        &cz::main_window::run_debug_demo
    );
}

void cz::main_window::open_ruleset()
{
    const QString file_path = QFileDialog::getOpenFileName(
        this,
        tr("Open Ruleset"),
        current_ruleset_path_,
        tr("JSON Files (*.json);;All Files (*.*)")
    );

    if (file_path.isEmpty()) {
        return;
    }

    load_ruleset_from_file(file_path);
}

void cz::main_window::save_ruleset()
{
    if (current_ruleset_path_.isEmpty()) {
        save_ruleset_as();
        return;
    }

    save_ruleset_to_file(current_ruleset_path_);
}

void cz::main_window::save_ruleset_as()
{
    const QString file_path = QFileDialog::getSaveFileName(
        this,
        tr("Save Ruleset As"),
        current_ruleset_path_.isEmpty() ? QStringLiteral("ruleset.json")
        : current_ruleset_path_,
        tr("JSON Files (*.json);;All Files (*.*)")
    );

    if (file_path.isEmpty()) {
        return;
    }

    save_ruleset_to_file(file_path);
}

void cz::main_window::view_edit_current_rules()
{
    show_empty_modal_dialog(this, tr("Current Rules"));
}

void cz::main_window::view_edit_current_start_conditions()
{
    show_empty_modal_dialog(this, tr("Current Start Conditions"));
}

void cz::main_window::run_simulation()
{
    QMessageBox::information(
        this,
        tr("Run Simulation"),
        tr("run_simulation() is currently stubbed in.")
    );

    // TODO:
    // Implement the actual simulation run here.
}

void cz::main_window::run_debug_demo()
{
    try {
        cz::cell_id_source ids;
        auto current_state = cz::random_cyto_state(500, 4, ids);

        for (auto& cell : current_state) {
            cell.state = 0;
            cell.phase = cz::life_stage::normal;
        }

        const std::vector<std::size_t> delete_indices = {
            5, 17, 42, 80, 111, 150, 222, 301
        };

        std::vector<cz::cell_id> delete_list;
        delete_list.reserve(delete_indices.size());

        for (std::size_t index : delete_indices) {
            if (index >= current_state.size()) {
                throw std::runtime_error(
                    "run_debug_demo: delete index out of range."
                );
            }

            delete_list.push_back(current_state[index].id);
            current_state[index].state = 2;
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
                i < diagram.embedding.vertices.size() &&
                add_list.size() < num_to_add;
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

        const auto deletion_set = delete_list | r::to<std::unordered_set>();

        cz::cyto_state next_state;
        next_state.reserve(
            current_state.size() - delete_list.size() + add_list.size()
        );

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
            {0, 255, 0},
            {255, 0, 0},
            {255, 255, 130},
            {0, 60, 200}
        };

        canvas_->set_show_cell_nuceli(false);

        const auto initial_frame = cz::to_cyto_frame(current_state, palette, {});
        const auto from_frame = cz::to_cyto_frame(trans.from, palette, from_scales);
        const auto to_frame = cz::to_cyto_frame(trans.to, palette, to_scales);

        canvas_->set(initial_frame);
        canvas_->start_transition(from_frame, to_frame);
    }
    catch (const std::exception& ex) {
        QMessageBox::critical(
            this,
            tr("Debug"),
            tr("Debug run failed:\n%1").arg(QString::fromUtf8(ex.what()))
        );
    }
}

bool cz::main_window::load_ruleset_from_file(const QString& file_path)
{
    QFile file(file_path);

    if (!file.open(QIODevice::ReadOnly)) {
        QMessageBox::critical(
            this,
            tr("Open Ruleset"),
            tr("Could not open file:\n%1").arg(file_path)
        );
        return false;
    }

    const QByteArray bytes = file.readAll();
    file.close();

    QJsonParseError parse_error;
    const QJsonDocument doc = QJsonDocument::fromJson(bytes, &parse_error);

    if (parse_error.error != QJsonParseError::NoError) {
        QMessageBox::critical(
            this,
            tr("Open Ruleset"),
            tr("Invalid JSON in file:\n%1\n\n%2")
            .arg(file_path, parse_error.errorString())
        );
        return false;
    }

    if (!doc.isObject()) {
        QMessageBox::critical(
            this,
            tr("Open Ruleset"),
            tr("Ruleset file must contain a top-level JSON object.")
        );
        return false;
    }

    current_ruleset_path_ = file_path;

    // TODO:
    // Read the JSON object and populate the current in-memory ruleset.
    [[maybe_unused]] const QJsonObject obj = doc.object();

    QMessageBox::information(
        this,
        tr("Open Ruleset"),
        tr("Ruleset loading is currently stubbed in.\n\nLoaded JSON from:\n%1")
        .arg(file_path)
    );

    return true;
}

bool cz::main_window::save_ruleset_to_file(const QString& file_path)
{
    QJsonObject obj;

    // TODO:
    // Populate obj from the current in-memory ruleset before writing.
    obj.insert("stub", true);

    const QJsonDocument doc(obj);

    QFile file(file_path);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
        QMessageBox::critical(
            this,
            tr("Save Ruleset"),
            tr("Could not write file:\n%1").arg(file_path)
        );
        return false;
    }

    file.write(doc.toJson(QJsonDocument::Indented));
    file.close();

    current_ruleset_path_ = file_path;

    QMessageBox::information(
        this,
        tr("Save Ruleset"),
        tr("Ruleset saving is currently stubbed in.\n\nWrote JSON to:\n%1")
        .arg(file_path)
    );

    return true;
}