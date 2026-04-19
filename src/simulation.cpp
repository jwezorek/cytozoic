#include "simulation.hpp"

#include <QMetaType>
#include <QMutexLocker>

#include <stdexcept>
#include <utility>

namespace {

    constexpr int k_min_frame_spacing_ms = 16;


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

void cz::simulation_thread::start_simulation(
    const cyto_params& params,
    simulation_presentation_mode mode)
{
    {
        QMutexLocker locker(&mutex_);
        params_ = params;
        mode_ = mode;
        stop_requested_ = false;
        transition_finished_ = false;
        pending_reclaim_ids_.clear();
        frame_spacing_started_ = false;
        frame_spacing_timer_.invalidate();
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
    simulation_presentation_mode mode;

    {
        QMutexLocker locker(&mutex_);
        params = params_;
        mode = mode_;
        stop_requested_ = false;
        transition_finished_ = false;
        pending_reclaim_ids_.clear();
        frame_spacing_started_ = false;
        frame_spacing_timer_.invalidate();
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

        if (!throttle_frame_rate(k_min_frame_spacing_ms)) {
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
            if (mode == simulation_presentation_mode::animated) {
                auto result = cz::apply_state_tables_animated(
                    id_source_,
                    current_state_,
                    params.cell_state_table,
                    params.cell_indexer,
                    params.birth_params,
                    params.palette
                );

                if (is_stop_requested()) {
                    break;
                }

                if (!throttle_frame_rate(k_min_frame_spacing_ms)) {
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
            }
            else {
                current_state_ = cz::apply_state_tables_quick(
                    id_source_,
                    current_state_,
                    params.cell_state_table,
                    params.cell_indexer,
                    params.birth_params
                );

                if (is_stop_requested()) {
                    break;
                }

                if (!throttle_frame_rate(k_min_frame_spacing_ms)) {
                    break;
                }

                emit canonical_frame_ready(
                    cz::to_cyto_frame(current_state_, params.palette, {})
                );
            }

            if (current_state_.empty()) {
                break;
            }
        }
    }
    catch (const std::exception& ex) {
        emit simulation_failed(QString::fromUtf8(ex.what()));
    }
    catch (...) {
        emit simulation_failed(QStringLiteral("Unknown simulation error."));
    }

    emit simulation_stopped();
}

void cz::simulation_thread::validate_params(const cyto_params& params) const
{
    //TODO: call a function in cytozoic.cpp to do this...
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

bool cz::simulation_thread::throttle_frame_rate(int min_frame_spacing_ms)
{
    if (min_frame_spacing_ms <= 0) {
        QMutexLocker locker(&mutex_);
        frame_spacing_timer_.restart();
        frame_spacing_started_ = true;
        return !stop_requested_;
    }

    QMutexLocker locker(&mutex_);

    if (!frame_spacing_started_) {
        frame_spacing_timer_.start();
        frame_spacing_started_ = true;
        return !stop_requested_;
    }

    const qint64 elapsed_ms = frame_spacing_timer_.elapsed();

    if (elapsed_ms < min_frame_spacing_ms) {
        const unsigned long remaining_ms =
            static_cast<unsigned long>(min_frame_spacing_ms - elapsed_ms);

        condition_.wait(&mutex_, remaining_ms);
    }

    if (stop_requested_) {
        return false;
    }

    frame_spacing_timer_.restart();
    return true;
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