#include "main_window.hpp"

#include "rules_dialog.hpp"
#include "serialize.hpp"
#include "simulation.hpp"

#include <QAction>
#include <QFileDialog>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>

#include <stdexcept>
#include <utility>

namespace
{
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

        if (params.cell_state_table.size() !=
            static_cast<std::size_t>(params.num_states)) {
            throw std::runtime_error(
                "cell_state_table row count must equal num_states."
            );
        }

        const std::size_t expected_cell_columns =
            params.cell_indexer->num_columns(
                static_cast<std::size_t>(params.num_states)
            );

        for (const auto& row : params.cell_state_table) {
            if (row.size() != expected_cell_columns) {
                throw std::runtime_error(
                    "cell_state_table width does not match cell_indexer."
                );
            }
        }

        const std::size_t expected_vertex_columns =
            params.vertex_indexer->num_columns(
                static_cast<std::size_t>(params.num_states)
            );

        if (params.vertex_table.size() != expected_vertex_columns) {
            throw std::runtime_error(
                "vertex_table width does not match vertex_indexer."
            );
        }

        if (params.palette.size() <
            static_cast<std::size_t>(params.num_states)) {
            throw std::runtime_error(
                "palette must contain at least num_states colors."
            );
        }
    }
} // namespace

cz::main_window::main_window(QWidget* parent) : 
        QMainWindow(parent),
        simulation_thread_(new cz::simulation_thread(this)) {

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

    connect(
        simulation_thread_,
        &cz::simulation_thread::initial_frame_ready,
        this,
        &cz::main_window::on_initial_frame_ready,
        Qt::QueuedConnection
    );

    connect(
        simulation_thread_,
        &cz::simulation_thread::transition_ready,
        this,
        &cz::main_window::on_transition_ready,
        Qt::QueuedConnection
    );

    connect(
        simulation_thread_,
        &cz::simulation_thread::simulation_failed,
        this,
        &cz::main_window::on_simulation_failed,
        Qt::QueuedConnection
    );

    connect(
        simulation_thread_,
        &cz::simulation_thread::simulation_stopped,
        this,
        &cz::main_window::on_simulation_stopped,
        Qt::QueuedConnection
    );
}

cz::main_window::~main_window()
{
    stop_simulation_thread(true);
}

void cz::main_window::set_params(const cyto_params& params)
{
    params_ = params;
    reset_simulation_session();
}

cz::cyto_params cz::main_window::get_params() const
{
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
        current_ruleset_path_.isEmpty()
        ? QStringLiteral("ruleset.json")
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
    rules_dialog dialog(this, get_params());

    if (dialog.exec() != QDialog::Accepted) {
        return;
    }

    set_params(dialog.get());
}

void cz::main_window::run_simulation()
{
    try {
        if (simulation_running_) {
            simulation_running_ = false;
            update_run_action_text();
            stop_simulation_thread(false);
            return;
        }

        if (transition_in_flight_) {
            return;
        }

        if (simulation_thread_ == nullptr || simulation_thread_->isRunning()) {
            return;
        }

        canvas_->set_show_cell_nuceli(false);

        simulation_running_ = true;
        transition_in_flight_ = false;
        update_run_action_text();

        simulation_thread_->start_simulation(params_);
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

void cz::main_window::on_initial_frame_ready(const cyto_frame& frame)
{
    if (!simulation_running_) {
        return;
    }

    if (frame.empty()) {
        canvas_->clear(Qt::black);
    }
    else {
        canvas_->set(frame);
    }
}

void cz::main_window::on_transition_ready(
    const cyto_frame& anim_start,
    const cyto_frame& anim_end)
{
    if (!simulation_running_) {
        return;
    }

    transition_in_flight_ = true;
    canvas_->start_transition(anim_start, anim_end);
}

void cz::main_window::on_transition_finished()
{
    if (!transition_in_flight_) {
        return;
    }

    transition_in_flight_ = false;

    if (simulation_thread_ != nullptr) {
        simulation_thread_->notify_transition_finished();
    }
}

void cz::main_window::on_simulation_failed(const QString& message)
{
    simulation_running_ = false;
    update_run_action_text();

    QMessageBox::critical(
        this,
        tr("Run Simulation"),
        tr("Simulation failed:\n%1").arg(message)
    );
}

void cz::main_window::on_simulation_stopped()
{
    simulation_running_ = false;
    update_run_action_text();
}

bool cz::main_window::load_ruleset_from_file(const QString& file_path)
{
    auto loaded = cz::load_ruleset_from_file(file_path.toStdString());

    if (!loaded) {
        QMessageBox::critical(
            this,
            tr("Open Ruleset"),
            tr("Could not open file:\n%1").arg(file_path)
        );
        return false;
    }

    try {
        validate_loaded_params(*loaded);
    }
    catch (const std::exception& ex) {
        QMessageBox::critical(
            this,
            tr("Open Ruleset"),
            tr("Invalid ruleset:\n%1").arg(QString::fromUtf8(ex.what()))
        );
        return false;
    }

    set_params(*loaded);
    current_ruleset_path_ = file_path;
    return true;
}

bool cz::main_window::save_ruleset_to_file(const QString& file_path)
{
    const bool result =
        cz::save_ruleset_to_file(file_path.toStdString(), get_params());

    if (!result) {
        QMessageBox::critical(
            this,
            tr("Save Ruleset"),
            tr("Could not write file:\n%1").arg(file_path)
        );
        return false;
    }

    current_ruleset_path_ = file_path;
    return true;
}

void cz::main_window::stop_simulation_thread(bool wait_for_finish)
{
    if (simulation_thread_ == nullptr) {
        return;
    }

    simulation_thread_->request_stop();

    if (transition_in_flight_ && canvas_ != nullptr) {
        canvas_->cancel_transition();
    }

    if (wait_for_finish && simulation_thread_->isRunning()) {
        simulation_thread_->wait();
    }
}

void cz::main_window::reset_simulation_session(bool clear_canvas)
{
    simulation_running_ = false;
    update_run_action_text();

    stop_simulation_thread(true);

    transition_in_flight_ = false;

    if (clear_canvas && canvas_ != nullptr) {
        canvas_->clear(Qt::black);
    }
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