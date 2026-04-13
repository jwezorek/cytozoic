#include "main_window.hpp"
#include "cytozoic_widget.hpp"
#include "voronoi.hpp"

#include <QEvent>
#include <QKeyEvent>
#include <QMenu>
#include <QMenuBar>
#include <ranges>
#include <unordered_map>
#include <utility>

namespace r = std::ranges;
namespace rv = std::ranges::views;

namespace {

    struct pending_animation {
        bool initialized = false;
        bool started = false;
        cz::cyto_frame from;
        cz::cyto_frame to;
    };

    std::unordered_map<cz::main_window*, pending_animation> g_pending_animations;

    class key_start_filter : public QObject {
    public:
        explicit key_start_filter(cz::main_window* window)
            : QObject(window),
            window_(window) {}

    protected:
        bool eventFilter(QObject* watched, QEvent* event) override {
            Q_UNUSED(watched);

            if (event->type() != QEvent::KeyPress) {
                return QObject::eventFilter(watched, event);
            }

            auto it = g_pending_animations.find(window_);
            if (it == g_pending_animations.end()) {
                return QObject::eventFilter(watched, event);
            }

            pending_animation& anim = it->second;
            if (!anim.initialized || anim.started) {
                return QObject::eventFilter(watched, event);
            }

            anim.started = true;
            window_->centralWidget()->removeEventFilter(this);

            auto* canvas = qobject_cast<cz::cytozoic_widget*>(window_->centralWidget());
            if (canvas != nullptr) {
                canvas->start_transition(anim.from, anim.to);
            }

            return true;
        }

    private:
        cz::main_window* window_;
    };

} // namespace

cz::main_window::main_window(QWidget* parent)
    : QMainWindow(parent)
{
    setCentralWidget(canvas_ = new cytozoic_widget(this));
    create_menus();
    setWindowTitle(tr("cytozoic"));
    resize(1200, 1200);

    setFocusPolicy(Qt::StrongFocus);
    centralWidget()->setFocusPolicy(Qt::StrongFocus);
    centralWidget()->installEventFilter(new key_start_filter(this));
}

cz::main_window::~main_window()
{
    g_pending_animations.erase(this);
}

void cz::main_window::showEvent(QShowEvent* event) {
    QMainWindow::showEvent(event);

    pending_animation& anim = g_pending_animations[this];
    if (anim.initialized) {
        return;
    }

    anim.initialized = true;

    cz::cell_id_source ids;
    auto state = random_cyto_state(500, 4, ids);

    for (auto& c : state) {
        c.state = 0;
    }

    const auto sites = state
        | rv::transform([](const cz::cell_state& cell) -> cz::point {
        return cell.site;
            })
        | r::to<std::vector>();

    const auto diagram = cz::to_voronoi_diagram(sites);

    std::vector<cz::cell_state> add_list;
    add_list.reserve(8);

    // Pick a few reasonably separated Voronoi vertices.
    // Stepping through the embedding vertices is a simple way to avoid
    // clustering all births in one tiny region.
    constexpr std::size_t num_to_add = 8;

    if (!diagram.embedding.vertices.empty()) {
        const std::size_t step = std::max<std::size_t>(
            1,
            diagram.embedding.vertices.size() / num_to_add
        );

        for (std::size_t i = 0;
            i < diagram.embedding.vertices.size() && add_list.size() < num_to_add;
            i += step) {
            add_list.push_back(
                cz::cell_state{
                    .id = ids.acquire(),
                    .site = diagram.embedding.vertices[i],
                    .state = 1,
                    .phase = cz::life_stage::new_born
                }
            );
        }
    }

    auto trans = generate_transition(
        state,
        {},
        add_list
    );

    cz::color_table palette = {
        {0, 255, 0},
        {255, 0, 0},
        {255, 255, 130},
        {0, 60, 200}
    };

    canvas_->set_show_cell_nuceli(false);

    const auto initial_frame = to_cyto_frame(state, palette);

    anim.from = to_cyto_frame(trans.from, palette);
    anim.to = to_cyto_frame(trans.to, palette);

    canvas_->set(initial_frame);
    centralWidget()->setFocus();
}

void cz::main_window::create_menus() {
    QMenu* file_menu = menuBar()->addMenu(tr("&File"));
    Q_UNUSED(file_menu);
}