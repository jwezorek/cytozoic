#pragma once

#include <QMainWindow>
#include <QString>

#include "cytozoic_widget.hpp"

class QAction;

namespace cz {

    class main_window : public QMainWindow {
        Q_OBJECT

    public:
        explicit main_window(QWidget* parent = nullptr);
        ~main_window() override;

    private:
        void create_menus();

        void open_ruleset();
        void save_ruleset();
        void save_ruleset_as();
        void view_edit_current_rules();
        void view_edit_current_start_conditions();
        void run_simulation();
        void run_debug_demo();

        bool load_ruleset_from_file(const QString& file_path);
        bool save_ruleset_to_file(const QString& file_path);

        cytozoic_widget* canvas_ = nullptr;
        QString current_ruleset_path_;
    };

} // namespace cz