#pragma once

#include <QMainWindow>
#include "cytozoic_widget.hpp"

namespace cz{

    class main_window : public QMainWindow {
        Q_OBJECT

    public:
        main_window(QWidget* parent = nullptr);
        ~main_window();

    protected:
        void showEvent(QShowEvent* event) override;

    private:
        void create_menus();
        cytozoic_widget* canvas_;
    };
}