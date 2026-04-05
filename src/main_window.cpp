#include "main_window.hpp"
#include "cytozoic_widget.hpp"
#include "util.hpp"
#include "voronoi.hpp"
#include <QMenuBar>
#include <QMenu>

namespace {

}

cz::main_window::main_window(QWidget* parent)
    : QMainWindow(parent)
{
    setCentralWidget(canvas_ = new cytozoic_widget(this));
    create_menus();
    setWindowTitle(tr("cytozoic"));
    resize(1200, 1200);

}

cz::main_window::~main_window()
{
}

void cz::main_window::showEvent(QShowEvent* event) {
    static bool initialized = false;
    QMainWindow::showEvent(event);

    if (initialized) {
        return;
    }

    initialized = true;

    auto seeds = cz::random_points(1000, 1.0, 1.0);
    auto unrelaxed = cz::construct_voronoi_diagram(
        seeds,
        { {0,0},{1.0,1.0} }
    );

    seeds = perform_lloyd_relaxation(seeds, { {0,0},{1.0,1.0} }, 0.001, 20);

    auto relaxed = cz::construct_voronoi_diagram(
        seeds,
        { {0,0},{1.0,1.0} }
    );
    auto from = to_cyto_frame(blank_state(unrelaxed),
        std::vector<color>{1, color{ 255,128,55 }}
    );
    auto to = to_cyto_frame(blank_state(relaxed),
        std::vector<color>{1, color{ 55,128,255 }}
    );

    canvas_->set_show_cell_nuceli(true);
    canvas_->start_transition(from, to);
}

void cz::main_window::create_menus() {
    QMenu* file_menu = menuBar()->addMenu(tr("&File"));
	/*
    QAction* open_act = new QAction(tr("&Open..."), this);
    open_act->setShortcut(QKeySequence::Open);
    connect(open_act, &QAction::triggered, this, &main_window::open_file);
    file_menu->addAction(open_act);

    file_menu->addSeparator();

    QAction* exit_act = new QAction(tr("E&xit"), this);
    exit_act->setShortcut(QKeySequence::Quit);
    connect(exit_act, &QAction::triggered, this, &QWidget::close);
    file_menu->addAction(exit_act);

    // Optional: View menu to toggle docks
    QMenu* view_menu = menuBar()->addMenu(tr("&View"));
    view_menu->addAction(tr("Toggle Source Palette"), [this](bool) {
        if (auto* dw = findChild<QDockWidget*>("Source Dock"))
            dw->setVisible(!dw->isVisible());
        });
    view_menu->addAction(tr("Toggle Target Palette"), [this](bool) {
        if (auto* dw = findChild<QDockWidget*>("Target Dock"))
            dw->setVisible(!dw->isVisible());
        });
	*/
}
