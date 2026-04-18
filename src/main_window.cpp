#include "main_window.hpp"
#include "rules_dialog.hpp"
#include "cytozoic.hpp"
#include "serialize.hpp"
#include <stdexcept>
#include <utility>
#include <vector>
#include <ranges>

namespace r = std::ranges;
namespace rv = std::ranges::views;

namespace
{
    constexpr double k_lloyd_min_delta = 0.001;
    constexpr int k_max_iterations = 20;

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

bool cz::main_window::load_ruleset_from_file(const QString& file_path) {

    auto loaded = cz::load_ruleset_from_file( file_path.toStdString() );
    if (!loaded) {
        QMessageBox::critical(
            this,
            tr("Open Ruleset"),
            tr("Could not open file:\n%1").arg(file_path)
        );
        return false;
    }
    set_params( *loaded );
    return true;

}

bool cz::main_window::save_ruleset_to_file(const QString& file_path) {

    auto result = cz::save_ruleset_to_file(file_path.toStdString(), get_params());
    if (!result) {
        QMessageBox::critical(
            this,
            tr("Save Ruleset"),
            tr("Could not write file:\n%1").arg(file_path)
        );
        return false;
    }
    return true;
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