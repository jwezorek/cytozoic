#pragma once

#include <QMainWindow>
#include <QString>

#include "cytozoic.hpp"
#include "cytozoic_widget.hpp"

class QAction;

namespace cz {

    class simulation_thread;

    class main_window : public QMainWindow {
        Q_OBJECT

    public:
        explicit main_window(QWidget* parent = nullptr);
        ~main_window() override;

        void set_params(const cyto_params& params);
        cyto_params get_params() const;

    private:
        void create_menus();

        void open_ruleset();
        void save_ruleset();
        void save_ruleset_as();
        void view_edit_current_rules();
        void run_simulation();

        void on_initial_frame_ready(const cyto_frame& frame);
        void on_transition_ready(
            const cyto_frame& anim_start,
            const cyto_frame& anim_end
        );
        void on_transition_finished();
        void on_simulation_failed(const QString& message);
        void on_simulation_stopped();

        bool load_ruleset_from_file(const QString& file_path);
        bool save_ruleset_to_file(const QString& file_path);

        void stop_simulation_thread(bool wait_for_finish);
        void reset_simulation_session(bool clear_canvas = true);
        void update_run_action_text();

        cytozoic_widget* canvas_ = nullptr;
        simulation_thread* simulation_thread_ = nullptr;
        QAction* run_simulation_action_ = nullptr;

        QString current_ruleset_path_;
        cyto_params params_;

        bool simulation_running_ = false;
        bool transition_in_flight_ = false;
    };

} // namespace cz