#include "simulation.hpp"

#include <QMetaType>
#include <QMutexLocker>

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

cz::simulation_thread::simulation_thread(QObject* parent)
    : QThread(parent)
{
    qRegisterMetaType<cz::cyto_frame>("cz::cyto_frame");
}

cz::simulation_thread::~simulation_thread()
{
    request_stop();
    wait();
}

void cz::simulation_thread::start_simulation(const cyto_params& params)
{
    {
        QMutexLocker locker(&mutex_);
        params_ = params;
        stop_requested_ = false;
        transition_finished_ = false;
        pending_reclaim_ids_.clear();
    }

    if (!isRunning()) {
        start();
    }
}

void cz::simulation_thread::request_stop()
{
    QMutexLocker locker(&mutex_);
    stop_requested_ = true;
    transition_finished_ = true;
    condition_.wakeAll();
}

void cz::simulation_thread::notify_transition_finished()
{
    QMutexLocker locker(&mutex_);
    transition_finished_ = true;
    condition_.wakeAll();
}

void cz::simulation_thread::run()
{
    cyto_params params;

    {
        QMutexLocker locker(&mutex_);
        params = params_;
        stop_requested_ = false;
        transition_finished_ = false;
        pending_reclaim_ids_.clear();
    }

    try {
        validate_params(params);

        id_source_.reset();
        current_state_.clear();
        current_state_ = cz::initial_cyto_state(params, id_source_);

        if (is_stop_requested()) {
            emit simulation_stopped();
            return;
        }

        emit initial_frame_ready(
            cz::to_cyto_frame(current_state_, params.palette, {})
        );

        if (current_state_.empty()) {
            emit simulation_stopped();
            return;
        }

        while (!is_stop_requested()) {
            auto result = cz::apply_state_tables(
                id_source_,
                current_state_,
                params.cell_state_table,
                params.cell_indexer,
                params.vertex_table,
                params.vertex_indexer,
                params.palette
            );

            if (is_stop_requested()) {
                break;
            }

            {
                QMutexLocker locker(&mutex_);
                pending_reclaim_ids_ = cz::deleted_cell_ids(result.anim_end);
                transition_finished_ = false;
            }

            emit transition_ready(
                std::move(result.anim_start),
                std::move(result.anim_end)
            );

            if (!wait_for_transition_finished()) {
                break;
            }

            release_pending_ids();
            current_state_ = std::move(result.next_state);

            if (current_state_.empty()) {
                break;
            }
        }
    } catch (const std::exception& ex) {
        emit simulation_failed(QString::fromUtf8(ex.what()));
    } catch (...) {
        emit simulation_failed(QStringLiteral("Unknown simulation error."));
    }

    emit simulation_stopped();
}

void cz::simulation_thread::validate_params(const cyto_params& params) const
{
    validate_loaded_params(params);
}

bool cz::simulation_thread::is_stop_requested() const
{
    QMutexLocker locker(&mutex_);
    return stop_requested_;
}

bool cz::simulation_thread::wait_for_transition_finished()
{
    QMutexLocker locker(&mutex_);

    while (!stop_requested_ && !transition_finished_) {
        condition_.wait(&mutex_);
    }

    const bool keep_running = !stop_requested_;
    transition_finished_ = false;
    return keep_running;
}

void cz::simulation_thread::release_pending_ids()
{
    std::vector<cell_id> reclaim_ids;

    {
        QMutexLocker locker(&mutex_);
        reclaim_ids = std::move(pending_reclaim_ids_);
        pending_reclaim_ids_.clear();
    }

    for (cz::cell_id id : reclaim_ids) {
        id_source_.release(id);
    }
}