#pragma once

#include <QMutex>
#include <QThread>
#include <QWaitCondition>

#include "cytozoic.hpp"

namespace cz {

    enum class simulation_presentation_mode {
        animated,
        quick
    };

    class simulation_thread : public QThread {
        Q_OBJECT

    public:
        explicit simulation_thread(QObject* parent = nullptr);
        ~simulation_thread() override;

        void start_simulation(
            const cyto_params& params,
            simulation_presentation_mode mode
        );

        void request_stop();
        void notify_transition_finished();

    signals:
        void initial_frame_ready(cyto_frame frame);
        void canonical_frame_ready(cyto_frame frame);
        void transition_ready(cyto_frame anim_start, cyto_frame anim_end);
        void simulation_failed(QString message);
        void simulation_stopped();

    protected:
        void run() override;

    private:
        void validate_params(const cyto_params& params) const;
        bool is_stop_requested() const;
        bool wait_for_transition_finished();
        void release_pending_ids();

        mutable QMutex mutex_;
        QWaitCondition condition_;

        cyto_params params_;
        simulation_presentation_mode mode_ =
            simulation_presentation_mode::animated;

        cell_id_source id_source_;
        cyto_state current_state_;
        std::vector<cell_id> pending_reclaim_ids_;

        bool stop_requested_ = false;
        bool transition_finished_ = false;
    };

} // namespace cz

Q_DECLARE_METATYPE(cz::cyto_frame)