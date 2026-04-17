#include "main_window.hpp"
#include "rules_dialog.hpp"
#include "third-party/json.hpp"
#include "cytozoic.hpp"
#include "voronoi.hpp"

#include <QDialog>
#include <QDialogButtonBox>
#include <QFile>
#include <QFileDialog>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QVBoxLayout>

#include <QColorDialog>
#include <QComboBox>
#include <QGridLayout>
#include <QHeaderView>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QSignalBlocker>
#include <QSlider>
#include <QSpinBox>
#include <QTabWidget>
#include <QTableWidget>
#include <QWidget>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <ranges>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <fstream>

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



    using json = nlohmann::json;

    json color_to_json(const cz::color& color)
    {
        return json{
            { "r", color.r },
            { "g", color.g },
            { "b", color.b }
        };
    }

    cz::color color_from_json(const json& j)
    {
        return cz::color{
            static_cast<uint8_t>(j.at("r").get<int>()),
            static_cast<uint8_t>(j.at("g").get<int>()),
            static_cast<uint8_t>(j.at("b").get<int>())
        };
    }

    json state_table_row_to_json(const cz::state_table_row& row)
    {
        json j = json::array();

        for (int8_t value : row) {
            j.push_back(static_cast<int>(value));
        }

        return j;
    }

    cz::state_table_row state_table_row_from_json(const json& j)
    {
        if (!j.is_array()) {
            throw std::runtime_error("state table row must be a JSON array.");
        }

        cz::state_table_row row;
        row.reserve(j.size());

        for (const auto& value : j) {
            row.push_back(static_cast<int8_t>(value.get<int>()));
        }

        return row;
    }

    json state_table_to_json(const cz::state_table& table)
    {
        json j = json::array();

        for (const auto& row : table) {
            j.push_back(state_table_row_to_json(row));
        }

        return j;
    }

    cz::state_table state_table_from_json(const json& j)
    {
        if (!j.is_array()) {
            throw std::runtime_error("state table must be a JSON array.");
        }

        cz::state_table table;
        table.reserve(j.size());

        for (const auto& row : j) {
            table.push_back(state_table_row_from_json(row));
        }

        return table;
    }

    std::vector<double> double_vector_from_json(const json& j)
    {
        if (!j.is_array()) {
            throw std::runtime_error("expected a JSON array of numbers.");
        }

        std::vector<double> values;
        values.reserve(j.size());

        for (const auto& value : j) {
            values.push_back(value.get<double>());
        }

        return values;
    }

    json double_vector_to_json(const std::vector<double>& values)
    {
        json j = json::array();

        for (double value : values) {
            j.push_back(value);
        }

        return j;
    }

    cz::color_table color_table_from_json(const json& j)
    {
        if (!j.is_array()) {
            throw std::runtime_error("palette must be a JSON array.");
        }

        cz::color_table palette;
        palette.reserve(j.size());

        for (const auto& value : j) {
            palette.push_back(color_from_json(value));
        }

        return palette;
    }

    json color_table_to_json(const cz::color_table& palette)
    {
        json j = json::array();

        for (const auto& color : palette) {
            j.push_back(color_to_json(color));
        }

        return j;
    }

    void validate_loaded_params(const cz::cyto_params& params)
    {
        if (params.num_states <= 0) {
            throw std::runtime_error("num_states must be positive.");
        }

        if (!params.cell_indexer) {
            throw std::runtime_error("missing cell_indexer.");
        }

        if (!params.vertex_indexer) {
            throw std::runtime_error("missing vertex_indexer.");
        }

        if (params.cell_state_table.size() != static_cast<std::size_t>(params.num_states)) {
            throw std::runtime_error(
                "cell_state_table row count must equal num_states."
            );
        }

        const std::size_t expected_cell_columns =
            params.cell_indexer->num_columns(static_cast<std::size_t>(params.num_states));

        for (const auto& row : params.cell_state_table) {
            if (row.size() != expected_cell_columns) {
                throw std::runtime_error(
                    "cell_state_table width does not match cell_indexer."
                );
            }
        }

        const std::size_t expected_vertex_columns =
            params.vertex_indexer->num_columns(static_cast<std::size_t>(params.num_states));

        if (params.vertex_table.size() != expected_vertex_columns) {
            throw std::runtime_error(
                "vertex_table width does not match vertex_indexer."
            );
        }

        if (params.palette.size() < static_cast<std::size_t>(params.num_states)) {
            throw std::runtime_error(
                "palette must contain at least num_states colors."
            );
        }
    }

} // namespace

cz::main_window::main_window(QWidget* parent)
    : QMainWindow(parent)
{
    setCentralWidget(canvas_ = new cytozoic_widget(this));
    create_menus();
    setWindowTitle(tr("cytozoic"));
    resize(1200, 1200);

    connect(
        canvas_,
        &cz::cytozoic_widget::transition_finished,
        this,
        &cz::main_window::on_transition_finished
    );
}

cz::main_window::~main_window() = default;

void cz::main_window::set_params(const cyto_params& params) {
    params_ = params;
    reset_simulation_session();
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

    run_simulation_action_ =
        cytozoic_menu->addAction(tr("&Run Simulation"));
    QObject::connect(
        run_simulation_action_,
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

void cz::main_window::view_edit_current_rules() {

    rules_dialog dialog(this, get_params());

    if (dialog.exec() != QDialog::Accepted) {
        return;
    }

    auto edited = dialog.get();
    set_params(edited);

}

void cz::main_window::view_edit_current_start_conditions() {

    show_empty_modal_dialog(this, tr("Current Start Conditions"));

}

void cz::main_window::run_simulation()
{
    try {
        if (simulation_running_) {
            simulation_running_ = false;
            update_run_action_text();
            return;
        }

        if (transition_in_flight_) {
            return;
        }

        id_source_.reset();
        pending_next_state_.clear();
        current_state_ = cz::initial_cyto_state(params_, id_source_);
        simulation_initialized_ = true;

        if (current_state_.empty()) {
            canvas_->clear(Qt::black);
            simulation_running_ = false;
            update_run_action_text();
            return;
        }

        canvas_->set_show_cell_nuceli(false);
        canvas_->set(cz::to_cyto_frame(current_state_, params_.palette, {}));

        simulation_running_ = true;
        update_run_action_text();

        advance_simulation();
    }
    catch (const std::exception& ex) {
        reset_simulation_session(false);

        QMessageBox::critical(
            this,
            tr("Run Simulation"),
            tr("Simulation failed:\n%1").arg(QString::fromUtf8(ex.what()))
        );
    }
}

void cz::main_window::advance_simulation()
{
    if (!simulation_running_ || transition_in_flight_) {
        return;
    }

    if (!simulation_initialized_) {
        return;
    }

    if (!params_.cell_indexer) {
        throw std::runtime_error("advance_simulation: missing cell indexer.");
    }

    if (!params_.vertex_indexer) {
        throw std::runtime_error("advance_simulation: missing vertex indexer.");
    }

    if (params_.num_states <= 0) {
        throw std::runtime_error("advance_simulation: num_states must be positive.");
    }

    if (params_.cell_state_table.size() != static_cast<std::size_t>(params_.num_states)) {
        throw std::runtime_error(
            "advance_simulation: cell_state_table row count must equal num_states."
        );
    }

    if (params_.vertex_table.size() !=
        params_.vertex_indexer->num_columns(static_cast<std::size_t>(params_.num_states))) {
        throw std::runtime_error(
            "advance_simulation: vertex_table width does not match vertex indexer."
        );
    }

    if (params_.palette.size() < static_cast<std::size_t>(params_.num_states)) {
        throw std::runtime_error("advance_simulation: palette is too small.");
    }

    const auto result = cz::apply_state_tables(
        id_source_,
        current_state_,
        params_.cell_state_table,
        params_.cell_indexer,
        params_.vertex_table,
        params_.vertex_indexer,
        params_.palette
    );

    if (!simulation_running_) {
        return;
    }

    pending_next_state_ = result.next_state;
    transition_in_flight_ = true;

    canvas_->start_transition(result.anim_start, result.anim_end);
}

void cz::main_window::on_transition_finished()
{
    for (const cz::cell_id id : canvas_->take_reclaimable_ids()) {
        id_source_.release(id);
    }

    current_state_ = std::move(pending_next_state_);
    pending_next_state_.clear();
    transition_in_flight_ = false;

    if (!simulation_running_) {
        update_run_action_text();
        return;
    }

    QTimer::singleShot(
        0,
        this,
        [this]() {
            advance_simulation();
        }
    );
}

void cz::main_window::run_debug_demo()
{   /*
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
    } */
}

bool cz::main_window::load_ruleset_from_file(const QString& file_path)
{
    try {
        std::ifstream file(file_path.toStdString(), std::ios::binary);

        if (!file) {
            QMessageBox::critical(
                this,
                tr("Open Ruleset"),
                tr("Could not open file:\n%1").arg(file_path)
            );
            return false;
        }

        const json j = json::parse(file);

        if (!j.is_object()) {
            throw std::runtime_error(
                "ruleset file must contain a top-level JSON object."
            );
        }

        const std::string format = j.value("format", std::string{});
        const int version = j.value("version", 0);

        if (format != "cytozoic_ruleset") {
            throw std::runtime_error("unrecognized ruleset format.");
        }

        if (version != 1) {
            throw std::runtime_error("unsupported ruleset version.");
        }

        cz::cyto_params loaded;

        loaded.num_states = j.at("num_states").get<int>();
        loaded.num_initial_cells = j.at("num_initial_cells").get<int>();
        loaded.initial_state_density =
            double_vector_from_json(j.at("initial_state_density"));
        loaded.palette =
            color_table_from_json(j.at("palette"));

        loaded.cell_indexer = cz::indexer_from_name(
            j.at("cell_indexer").get<std::string>()
        );

        loaded.vertex_indexer = cz::indexer_from_name(
            j.at("vertex_indexer").get<std::string>()
        );

        loaded.cell_state_table =
            state_table_from_json(j.at("cell_state_table"));

        loaded.vertex_table =
            state_table_row_from_json(j.at("vertex_table"));

        validate_loaded_params(loaded);

        set_params(loaded);
        current_ruleset_path_ = file_path;

        return true;
    }
    catch (const std::exception& ex) {
        QMessageBox::critical(
            this,
            tr("Open Ruleset"),
            tr("Failed to load ruleset from:\n%1\n\n%2")
            .arg(file_path, QString::fromUtf8(ex.what()))
        );
        return false;
    }
}

bool cz::main_window::save_ruleset_to_file(const QString& file_path)
{
    try {
        const cz::cyto_params params = get_params();

        validate_loaded_params(params);

        json j = {
            { "format", "cytozoic_ruleset" },
            { "version", 1 },
            { "num_states", params.num_states },
            { "num_initial_cells", params.num_initial_cells },
            { "initial_state_density", double_vector_to_json(params.initial_state_density) },
            { "palette", color_table_to_json(params.palette) },
            { "cell_indexer", params.cell_indexer->name() },
            { "cell_state_table", state_table_to_json(params.cell_state_table) },
            { "vertex_indexer", params.vertex_indexer->name() },
            { "vertex_table", state_table_row_to_json(params.vertex_table) }
        };

        std::ofstream file(file_path.toStdString(), std::ios::binary | std::ios::trunc);

        if (!file) {
            QMessageBox::critical(
                this,
                tr("Save Ruleset"),
                tr("Could not write file:\n%1").arg(file_path)
            );
            return false;
        }

        file << j.dump(4);
        file.close();

        current_ruleset_path_ = file_path;
        return true;
    }
    catch (const std::exception& ex) {
        QMessageBox::critical(
            this,
            tr("Save Ruleset"),
            tr("Failed to save ruleset to:\n%1\n\n%2")
            .arg(file_path, QString::fromUtf8(ex.what()))
        );
        return false;
    }
}

void cz::main_window::reset_simulation_session(bool clear_canvas)
{
    simulation_initialized_ = false;
    simulation_running_ = false;
    transition_in_flight_ = false;

    current_state_.clear();
    pending_next_state_.clear();
    id_source_.reset();

    if (clear_canvas && canvas_ != nullptr) {
        canvas_->clear(Qt::black);
    }

    update_run_action_text();
}

void cz::main_window::update_run_action_text()
{
    if (run_simulation_action_ == nullptr) {
        return;
    }

    run_simulation_action_->setText(
        simulation_running_ ? tr("S&top") : tr("&Run Simulation")
    );
}