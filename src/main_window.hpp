#pragma once

#include <QMainWindow>
#include <QString>

#include "cytozoic_widget.hpp"
#include "cytozoic.hpp"

class QAction;

namespace cz {

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

        void advance_simulation();
        void on_transition_finished();

        bool load_ruleset_from_file(const QString& file_path);
        bool save_ruleset_to_file(const QString& file_path);

        void reset_simulation_session(bool clear_canvas = true);
        void update_run_action_text();

        cytozoic_widget* canvas_ = nullptr;
        QAction* run_simulation_action_ = nullptr;

        QString current_ruleset_path_;
        cyto_params params_;

        cell_id_source id_source_;
        cyto_state current_state_;
        cyto_state pending_next_state_;

        bool simulation_initialized_ = false;
        bool simulation_running_ = false;
        bool transition_in_flight_ = false;
    };

} // namespace cz