#include "main_window.hpp"

#include "cytozoic.hpp"
#include "voronoi.hpp"

#include <QDialog>
#include <QDialogButtonBox>
#include <QFile>
#include <QFileDialog>
#include <QJsonDocument>
#include <QJsonObject>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QVBoxLayout>

#include <QColorDialog>
#include <QComboBox>
#include <QGridLayout>
#include <QHeaderView>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QSignalBlocker>
#include <QSlider>
#include <QSpinBox>
#include <QTabWidget>
#include <QTableWidget>
#include <QWidget>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <ranges>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace r = std::ranges;
namespace rv = std::ranges::views;

namespace
{
    constexpr double k_lloyd_min_delta = 0.001;
    constexpr int k_max_iterations = 20;

    double polygon_area(const cz::polygon& poly)
    {
        if (poly.size() < 3) {
            return 0.0;
        }

        double area2 = 0.0;

        for (std::size_t i = 0; i < poly.size(); ++i) {
            const cz::point& a = poly[i];
            const cz::point& b = poly[(i + 1) % poly.size()];
            area2 += a.x * b.y - b.x * a.y;
        }

        return std::abs(area2) * 0.5;
    }

    double polygon_scale(const cz::polygon& poly)
    {
        constexpr double k_min_scale = 1e-9;

        return std::max(
            polygon_area(poly) / std::numbers::pi_v<double>,
            k_min_scale
        );
    }

    std::unordered_map<cz::cell_id, double> make_scale_map(
        const cz::cyto_state& state)
    {
        std::unordered_map<cz::cell_id, double> scale_by_id;
        scale_by_id.reserve(state.size());

        if (state.empty()) {
            return scale_by_id;
        }

        const auto sites = state
            | rv::transform([](const cz::cell_state& cell) -> cz::point {
            return cell.site;
                })
            | r::to<std::vector>();

        const auto polygons = cz::to_voronoi_polygons(sites);

        if (polygons.size() != state.size()) {
            throw std::runtime_error(
                "make_scale_map: polygon count did not match state size."
            );
        }

        for (std::size_t i = 0; i < state.size(); ++i) {
            auto [it, inserted] = scale_by_id.emplace(
                state[i].id,
                polygon_scale(polygons[i])
            );

            if (!inserted) {
                throw std::runtime_error("make_scale_map: duplicate cell id.");
            }
        }

        return scale_by_id;
    }

    std::vector<double> make_from_scales(
        const cz::cyto_state& from_state,
        const std::unordered_map<cz::cell_id, double>& current_scale_by_id,
        const std::unordered_map<cz::cell_id, double>& next_scale_by_id)
    {
        std::vector<double> scales;
        scales.reserve(from_state.size());

        for (const auto& cell : from_state) {
            if (cell.phase == cz::life_stage::new_born) {
                auto it = next_scale_by_id.find(cell.id);
                if (it == next_scale_by_id.end()) {
                    throw std::runtime_error(
                        "make_from_scales: missing next scale for newborn cell."
                    );
                }

                scales.push_back(it->second);
            }
            else {
                auto it = current_scale_by_id.find(cell.id);
                if (it == current_scale_by_id.end()) {
                    throw std::runtime_error(
                        "make_from_scales: missing current scale for existing cell."
                    );
                }

                scales.push_back(it->second);
            }
        }

        return scales;
    }

    std::vector<double> make_to_scales(
        const cz::cyto_state& to_state,
        const std::unordered_map<cz::cell_id, double>& current_scale_by_id,
        const std::unordered_map<cz::cell_id, double>& next_scale_by_id)
    {
        std::vector<double> scales;
        scales.reserve(to_state.size());

        for (const auto& cell : to_state) {
            if (cell.phase == cz::life_stage::dying) {
                auto it = current_scale_by_id.find(cell.id);
                if (it == current_scale_by_id.end()) {
                    throw std::runtime_error(
                        "make_to_scales: missing current scale for dying cell."
                    );
                }

                scales.push_back(it->second);
            }
            else {
                auto it = next_scale_by_id.find(cell.id);
                if (it == next_scale_by_id.end()) {
                    throw std::runtime_error(
                        "make_to_scales: missing next scale for non-dying cell."
                    );
                }

                scales.push_back(it->second);
            }
        }

        return scales;
    }

    void relax_state_unweighted(cz::cyto_state& state)
    {
        if (state.empty()) {
            return;
        }

        const auto sites = state
            | rv::transform([](const cz::cell_state& cell) -> cz::point {
            return cell.site;
                })
            | r::to<std::vector>();

        const auto relaxed_sites = cz::perform_lloyd_relaxation(
            sites,
            k_lloyd_min_delta,
            k_max_iterations
        );

        if (relaxed_sites.size() != state.size()) {
            throw std::runtime_error(
                "relax_state_unweighted: relaxed site count mismatch."
            );
        }

        for (auto&& [cell, site] : rv::zip(state, relaxed_sites)) {
            cell.site = site;
        }
    }

    void show_empty_modal_dialog(QWidget* parent, const QString& title)
    {
        QDialog dialog(parent);
        dialog.setWindowTitle(title);
        dialog.setModal(true);

        auto* layout = new QVBoxLayout(&dialog);
        auto* buttons = new QDialogButtonBox(QDialogButtonBox::Ok, &dialog);

        QObject::connect(
            buttons,
            &QDialogButtonBox::accepted,
            &dialog,
            &QDialog::accept
        );

        layout->addWidget(buttons);

        dialog.exec();
    }

    class state_density_editor : public QWidget
    {
    public:
        explicit state_density_editor(QWidget* parent = nullptr)
            : QWidget(parent)
        {
            auto* root = new QVBoxLayout(this);
            root->setContentsMargins(0, 0, 0, 0);
            layout_ = new QGridLayout();
            root->addLayout(layout_);
            root->addStretch();
        }

        void set_num_states(int num_states)
        {
            if (num_states < 1) {
                num_states = 1;
            }

            while (sliders_.size() > static_cast<std::size_t>(num_states)) {
                remove_last_row();
            }

            while (sliders_.size() < static_cast<std::size_t>(num_states)) {
                add_row(static_cast<int>(sliders_.size()));
            }

            refresh_labels();
        }

        void set_density(const std::vector<double>& density)
        {
            if (sliders_.empty()) {
                return;
            }

            std::vector<int> discrete(sliders_.size(), 0);

            if (!density.empty()) {
                const std::size_t count = std::min(density.size(), sliders_.size());

                for (std::size_t i = 0; i < count; ++i) {
                    discrete[i] = std::clamp(
                        static_cast<int>(std::lround(density[i] * 100.0)),
                        0,
                        100
                    );
                }
            }

            if (std::ranges::all_of(discrete, [](int v) { return v == 0; })) {
                discrete.front() = 100;
            }

            for (std::size_t i = 0; i < sliders_.size(); ++i) {
                QSignalBlocker blocker(sliders_[i]);
                sliders_[i]->setValue(discrete[i]);
            }

            refresh_labels();
        }

        std::vector<double> density() const
        {
            std::vector<double> result;
            result.reserve(sliders_.size());

            int total = 0;
            for (const auto* slider : sliders_) {
                total += slider->value();
            }

            if (total <= 0) {
                result.resize(sliders_.size(), 0.0);
                if (!result.empty()) {
                    result[0] = 1.0;
                }
                return result;
            }

            for (const auto* slider : sliders_) {
                result.push_back(
                    static_cast<double>(slider->value()) /
                    static_cast<double>(total)
                );
            }

            return result;
        }

    private:
        void add_row(int state_index)
        {
            auto* state_label = new QLabel(
                QStringLiteral("State %1").arg(state_index),
                this
            );

            auto* slider = new QSlider(Qt::Horizontal, this);
            slider->setRange(0, 100);
            slider->setValue(state_index == 0 ? 100 : 0);

            auto* value_label = new QLabel(this);
            value_label->setMinimumWidth(120);

            const int row = static_cast<int>(sliders_.size());
            layout_->addWidget(state_label, row, 0);
            layout_->addWidget(slider, row, 1);
            layout_->addWidget(value_label, row, 2);

            labels_.push_back(state_label);
            sliders_.push_back(slider);
            value_labels_.push_back(value_label);

            connect(
                slider,
                &QSlider::valueChanged,
                this,
                [this]() {
                    refresh_labels();
                }
            );
        }

        void remove_last_row()
        {
            if (sliders_.empty()) {
                return;
            }

            delete labels_.back();
            delete sliders_.back();
            delete value_labels_.back();

            labels_.pop_back();
            sliders_.pop_back();
            value_labels_.pop_back();
        }

        void refresh_labels()
        {
            int total = 0;
            for (const auto* slider : sliders_) {
                total += slider->value();
            }

            for (std::size_t i = 0; i < sliders_.size(); ++i) {
                labels_[i]->setText(QStringLiteral("State %1").arg(i));

                const int raw = sliders_[i]->value();
                const double normalized = total > 0
                    ? static_cast<double>(raw) / static_cast<double>(total)
                    : 0.0;

                value_labels_[i]->setText(
                    QStringLiteral("%1  (%2)")
                    .arg(raw)
                    .arg(normalized, 0, 'f', 3)
                );
            }
        }

        QGridLayout* layout_ = nullptr;
        std::vector<QLabel*> labels_;
        std::vector<QSlider*> sliders_;
        std::vector<QLabel*> value_labels_;
    };

    class palette_editor : public QWidget
    {
    public:
        explicit palette_editor(QWidget* parent = nullptr)
            : QWidget(parent)
        {
            auto* root = new QGridLayout(this);
            root->setContentsMargins(0, 0, 0, 0);
            layout_ = root;
        }

        void set_num_states(int num_states)
        {
            if (num_states < 1) {
                num_states = 1;
            }

            while (buttons_.size() > static_cast<std::size_t>(num_states)) {
                remove_last_row();
            }

            while (buttons_.size() < static_cast<std::size_t>(num_states)) {
                add_row(static_cast<int>(buttons_.size()));
            }

            refresh_buttons();
        }

        void set_palette(const cz::color_table& palette)
        {
            colors_.clear();
            colors_.reserve(buttons_.size());

            for (std::size_t i = 0; i < buttons_.size(); ++i) {
                if (i < palette.size()) {
                    colors_.push_back(QColor(
                        palette[i].r,
                        palette[i].g,
                        palette[i].b
                    ));
                }
                else {
                    const int hue = static_cast<int>((360 * i) / std::max<std::size_t>(1, buttons_.size()));
                    colors_.push_back(QColor::fromHsv(hue, 255, 220));
                }
            }

            refresh_buttons();
        }

        cz::color_table palette() const
        {
            cz::color_table result;
            result.reserve(colors_.size());

            for (const auto& color : colors_) {
                result.push_back(cz::color{
                    static_cast<uint8_t>(color.red()),
                    static_cast<uint8_t>(color.green()),
                    static_cast<uint8_t>(color.blue())
                    });
            }

            return result;
        }

    private:
        void add_row(int state_index)
        {
            auto* state_label = new QLabel(
                QStringLiteral("State %1").arg(state_index),
                this
            );

            auto* button = new QPushButton(this);
            button->setMinimumWidth(80);
            button->setText(QStringLiteral("Pick..."));

            const int row = static_cast<int>(buttons_.size());
            layout_->addWidget(state_label, row, 0);
            layout_->addWidget(button, row, 1);

            labels_.push_back(state_label);
            buttons_.push_back(button);
            colors_.push_back(QColor::fromHsv((row * 36) % 360, 255, 220));

            connect(
                button,
                &QPushButton::clicked,
                this,
                [this, row]() {
                    const QColor chosen = QColorDialog::getColor(
                        colors_[row],
                        this,
                        QStringLiteral("Choose color for state %1").arg(row)
                    );

                    if (!chosen.isValid()) {
                        return;
                    }

                    colors_[row] = chosen;
                    refresh_buttons();
                }
            );
        }

        void remove_last_row()
        {
            if (buttons_.empty()) {
                return;
            }

            delete labels_.back();
            delete buttons_.back();

            labels_.pop_back();
            buttons_.pop_back();
            colors_.pop_back();
        }

        void refresh_buttons()
        {
            for (std::size_t i = 0; i < buttons_.size(); ++i) {
                labels_[i]->setText(QStringLiteral("State %1").arg(i));

                const QColor& color = colors_[i];
                buttons_[i]->setStyleSheet(
                    QStringLiteral(
                        "background-color: rgb(%1, %2, %3);"
                        "border: 1px solid black;"
                        "min-height: 24px;"
                    )
                    .arg(color.red())
                    .arg(color.green())
                    .arg(color.blue())
                );

                buttons_[i]->setText(
                    QStringLiteral("(%1, %2, %3)")
                    .arg(color.red())
                    .arg(color.green())
                    .arg(color.blue())
                );
            }
        }

        QGridLayout* layout_ = nullptr;
        std::vector<QLabel*> labels_;
        std::vector<QPushButton*> buttons_;
        std::vector<QColor> colors_;
    };

    cz::neighborhood_indexer make_indexer_from_dialog_name(const QString& name)
    {
        const std::string value = name.toStdString();

        if (value == "sum of states") {
            return pro::make_proxy<cz::neighborhood_indexer_facade>(
                cz::sum_of_states_indexer{ 10 }
            );
        }

        if (value == "max state") {
            return pro::make_proxy<cz::neighborhood_indexer_facade>(
                cz::max_state_indexer{}
            );
        }

        throw std::runtime_error("unknown neighborhood indexer name.");
    }

    QString indexer_name_or_default(
        const cz::neighborhood_indexer& indexer,
        const QString& fallback)
    {
        if (indexer) {
            return QString::fromStdString(indexer->name());
        }

        return fallback;
    }

    int safe_cell_value(const QTableWidget* table, int row, int col, int fallback = -1)
    {
        if (auto* widget = qobject_cast<QSpinBox*>(table->cellWidget(row, col))) {
            return widget->value();
        }

        return fallback;
    }

    void rebuild_state_table(
        QTableWidget* table,
        int rows,
        int columns,
        int min_value,
        int max_value)
    {
        std::vector<std::vector<int>> old_values(
            static_cast<std::size_t>(table->rowCount()),
            std::vector<int>(static_cast<std::size_t>(table->columnCount()), -1)
        );

        for (int r = 0; r < table->rowCount(); ++r) {
            for (int c = 0; c < table->columnCount(); ++c) {
                old_values[r][c] = safe_cell_value(table, r, c, -1);
            }
        }

        table->clear();
        table->setRowCount(rows);
        table->setColumnCount(columns);

        QStringList headers;
        headers.reserve(columns);
        for (int c = 0; c < columns; ++c) {
            headers.push_back(QString::number(c));
        }

        table->setHorizontalHeaderLabels(headers);

        QStringList row_headers;
        row_headers.reserve(rows);
        for (int r = 0; r < rows; ++r) {
            row_headers.push_back(QString::number(r));
        }

        table->setVerticalHeaderLabels(row_headers);

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < columns; ++c) {
                auto* spin = new QSpinBox(table);
                spin->setRange(min_value, max_value);

                const int preserved =
                    (r < static_cast<int>(old_values.size()) &&
                        c < static_cast<int>(old_values[r].size()))
                    ? old_values[r][c]
                    : -1;

                spin->setValue(std::clamp(preserved, min_value, max_value));
                table->setCellWidget(r, c, spin);
            }
        }

        table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
        table->verticalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    }

    void set_state_table_values(
        QTableWidget* table,
        const cz::state_table& values)
    {
        const int row_count = std::min(table->rowCount(), static_cast<int>(values.size()));

        for (int r = 0; r < row_count; ++r) {
            const int col_count = std::min(
                table->columnCount(),
                static_cast<int>(values[r].size())
            );

            for (int c = 0; c < col_count; ++c) {
                if (auto* spin = qobject_cast<QSpinBox*>(table->cellWidget(r, c))) {
                    spin->setValue(values[r][c]);
                }
            }
        }
    }

    cz::state_table read_state_table(const QTableWidget* table)
    {
        cz::state_table result(
            static_cast<std::size_t>(table->rowCount()),
            cz::state_table_row(static_cast<std::size_t>(table->columnCount()), -1)
        );

        for (int r = 0; r < table->rowCount(); ++r) {
            for (int c = 0; c < table->columnCount(); ++c) {
                result[r][c] = static_cast<int8_t>(safe_cell_value(table, r, c, -1));
            }
        }

        return result;
    }

    void set_state_table_row_values(
        QTableWidget* table,
        const cz::state_table_row& values)
    {
        if (table->rowCount() < 1) {
            return;
        }

        const int col_count = std::min(
            table->columnCount(),
            static_cast<int>(values.size())
        );

        for (int c = 0; c < col_count; ++c) {
            if (auto* spin = qobject_cast<QSpinBox*>(table->cellWidget(0, c))) {
                spin->setValue(values[c]);
            }
        }
    }

    cz::state_table_row read_state_table_row(const QTableWidget* table)
    {
        cz::state_table_row result(
            static_cast<std::size_t>(table->columnCount()),
            -1
        );

        if (table->rowCount() < 1) {
            return result;
        }

        for (int c = 0; c < table->columnCount(); ++c) {
            result[c] = static_cast<int8_t>(safe_cell_value(table, 0, c, -1));
        }

        return result;
    }

} // namespace

cz::main_window::main_window(QWidget* parent)
    : QMainWindow(parent)
{
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
}

cz::main_window::~main_window() = default;

void cz::main_window::set_params(const cyto_params& params) {
    params_ = params;
}

cz::cyto_params cz::main_window::get_params() const {
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

    QAction* run_simulation_action =
        cytozoic_menu->addAction(tr("&Run Simulation"));
    QObject::connect(
        run_simulation_action,
        &QAction::triggered,
        this,
        &cz::main_window::run_simulation
    );

    QAction* debug_action =
        cytozoic_menu->addAction(tr("&Debug"));
    QObject::connect(
        debug_action,
        &QAction::triggered,
        this,
        &cz::main_window::run_debug_demo
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
        current_ruleset_path_.isEmpty() ? QStringLiteral("ruleset.json")
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
    cyto_params edited = get_params();

    if (edited.num_states < 1) {
        edited.num_states = 1;
    }

    if (edited.num_states > 10) {
        edited.num_states = 10;
    }

    if (edited.num_initial_cells < 10) {
        edited.num_initial_cells = 10;
    }

    if (edited.num_initial_cells > 5000) {
        edited.num_initial_cells = 5000;
    }

    if (!edited.cell_indexer) {
        edited.cell_indexer = pro::make_proxy<cz::neighborhood_indexer_facade>(
            cz::sum_of_states_indexer{ 10 }
        );
    }

    if (!edited.vertex_indexer) {
        edited.vertex_indexer = pro::make_proxy<cz::neighborhood_indexer_facade>(
            cz::sum_of_states_indexer{ 10 }
        );
    }

    QDialog dialog(this);
    dialog.setWindowTitle(tr("Current Rules"));
    dialog.setModal(true);
    dialog.resize(1100, 700);

    auto* root_layout = new QVBoxLayout(&dialog);
    auto* tabs = new QTabWidget(&dialog);
    root_layout->addWidget(tabs);

    auto* buttons = new QDialogButtonBox(
        QDialogButtonBox::Ok | QDialogButtonBox::Cancel,
        &dialog
    );
    root_layout->addWidget(buttons);

    connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);

    //--------------------------------------------------------------------------
    // Initial configuration tab
    //--------------------------------------------------------------------------

    auto* initial_tab = new QWidget(tabs);
    auto* initial_layout = new QVBoxLayout(initial_tab);
    auto* initial_grid = new QGridLayout();
    initial_layout->addLayout(initial_grid);

    auto* num_states_label = new QLabel(tr("Number of states"), initial_tab);
    auto* num_states_combo = new QComboBox(initial_tab);
    for (int i = 1; i <= 10; ++i) {
        num_states_combo->addItem(QString::number(i), i);
    }

    auto* num_initial_cells_label = new QLabel(tr("Initial cell count"), initial_tab);
    auto* num_initial_cells_spin = new QSpinBox(initial_tab);
    num_initial_cells_spin->setRange(10, 5000);

    auto* state_density_label = new QLabel(tr("State density"), initial_tab);
    auto* density_editor = new state_density_editor(initial_tab);

    auto* palette_label = new QLabel(tr("Palette"), initial_tab);
    auto* palette_editor_widget = new palette_editor(initial_tab);

    initial_grid->addWidget(num_states_label, 0, 0);
    initial_grid->addWidget(num_states_combo, 0, 1);
    initial_grid->addWidget(num_initial_cells_label, 1, 0);
    initial_grid->addWidget(num_initial_cells_spin, 1, 1);
    initial_grid->addWidget(state_density_label, 2, 0, Qt::AlignTop);
    initial_grid->addWidget(density_editor, 2, 1);
    initial_grid->addWidget(palette_label, 3, 0, Qt::AlignTop);
    initial_grid->addWidget(palette_editor_widget, 3, 1);

    initial_layout->addStretch();
    tabs->addTab(initial_tab, tr("Initial Configuration"));

    //--------------------------------------------------------------------------
    // Cell table tab
    //--------------------------------------------------------------------------

    auto* cell_tab = new QWidget(tabs);
    auto* cell_layout = new QVBoxLayout(cell_tab);
    auto* cell_top_row = new QHBoxLayout();

    auto* cell_indexer_label = new QLabel(tr("Indexer"), cell_tab);
    auto* cell_indexer_combo = new QComboBox(cell_tab);
    for (const auto& name : cz::named_indexers()) {
        cell_indexer_combo->addItem(QString::fromStdString(name));
    }

    cell_top_row->addWidget(cell_indexer_label);
    cell_top_row->addWidget(cell_indexer_combo);
    cell_top_row->addStretch();

    auto* cell_table = new QTableWidget(cell_tab);
    cell_table->setAlternatingRowColors(true);

    cell_layout->addLayout(cell_top_row);
    cell_layout->addWidget(cell_table);
    tabs->addTab(cell_tab, tr("Cell Table"));

    //--------------------------------------------------------------------------
    // Vertex table tab
    //--------------------------------------------------------------------------

    auto* vertex_tab = new QWidget(tabs);
    auto* vertex_layout = new QVBoxLayout(vertex_tab);
    auto* vertex_top_row = new QHBoxLayout();

    auto* vertex_indexer_label = new QLabel(tr("Indexer"), vertex_tab);
    auto* vertex_indexer_combo = new QComboBox(vertex_tab);
    for (const auto& name : cz::named_indexers()) {
        vertex_indexer_combo->addItem(QString::fromStdString(name));
    }

    vertex_top_row->addWidget(vertex_indexer_label);
    vertex_top_row->addWidget(vertex_indexer_combo);
    vertex_top_row->addStretch();

    auto* vertex_table = new QTableWidget(vertex_tab);
    vertex_table->setAlternatingRowColors(true);

    vertex_layout->addLayout(vertex_top_row);
    vertex_layout->addWidget(vertex_table);
    tabs->addTab(vertex_tab, tr("Vertex Table"));

    //--------------------------------------------------------------------------
    // Initialization from current params
    //--------------------------------------------------------------------------

    num_states_combo->setCurrentIndex(
        std::max(0, num_states_combo->findData(edited.num_states))
    );
    num_initial_cells_spin->setValue(edited.num_initial_cells);

    density_editor->set_num_states(edited.num_states);
    density_editor->set_density(edited.initial_state_density);

    palette_editor_widget->set_num_states(edited.num_states);
    palette_editor_widget->set_palette(edited.palette);

    cell_indexer_combo->setCurrentText(
        indexer_name_or_default(edited.cell_indexer, QStringLiteral("sum of states"))
    );

    vertex_indexer_combo->setCurrentText(
        indexer_name_or_default(edited.vertex_indexer, QStringLiteral("sum of states"))
    );

    auto rebuild_tables = [&]() {
        const int num_states = num_states_combo->currentData().toInt();
        const int min_value = -1;
        const int max_value = num_states - 1;

        density_editor->set_num_states(num_states);
        palette_editor_widget->set_num_states(num_states);

        const auto cell_indexer =
            make_indexer_from_dialog_name(cell_indexer_combo->currentText());
        const auto vertex_indexer =
            make_indexer_from_dialog_name(vertex_indexer_combo->currentText());

        const int cell_columns =
            static_cast<int>(cell_indexer->num_columns(static_cast<std::size_t>(num_states)));

        const int vertex_columns =
            static_cast<int>(vertex_indexer->num_columns(static_cast<std::size_t>(num_states)));

        rebuild_state_table(
            cell_table,
            num_states,
            cell_columns,
            min_value,
            max_value
        );

        rebuild_state_table(
            vertex_table,
            1,
            vertex_columns,
            min_value,
            max_value
        );
        };

    rebuild_tables();
    set_state_table_values(cell_table, edited.cell_state_table);
    set_state_table_row_values(vertex_table, edited.vertex_table);

    connect(
        num_states_combo,
        &QComboBox::currentIndexChanged,
        &dialog,
        [&](int) {
            rebuild_tables();
        }
    );

    connect(
        cell_indexer_combo,
        &QComboBox::currentIndexChanged,
        &dialog,
        [&](int) {
            rebuild_tables();
        }
    );

    connect(
        vertex_indexer_combo,
        &QComboBox::currentIndexChanged,
        &dialog,
        [&](int) {
            rebuild_tables();
        }
    );

    if (dialog.exec() != QDialog::Accepted) {
        return;
    }

    edited.num_states = num_states_combo->currentData().toInt();
    edited.num_initial_cells = num_initial_cells_spin->value();
    edited.initial_state_density = density_editor->density();
    edited.palette = palette_editor_widget->palette();

    edited.cell_indexer =
        make_indexer_from_dialog_name(cell_indexer_combo->currentText());
    edited.vertex_indexer =
        make_indexer_from_dialog_name(vertex_indexer_combo->currentText());

    edited.cell_state_table = read_state_table(cell_table);
    edited.vertex_table = read_state_table_row(vertex_table);

    set_params(edited);
}

void cz::main_window::view_edit_current_start_conditions()
{
    show_empty_modal_dialog(this, tr("Current Start Conditions"));
}

void cz::main_window::run_simulation()
{
    try {
        if (transition_in_flight_) {
            return;
        }

        if (!simulation_initialized_) {
            id_source_.reset();
            current_state_ = cz::initial_cyto_state(params_, id_source_);
            pending_next_state_.clear();
            simulation_initialized_ = true;

            if (current_state_.empty()) {
                canvas_->clear(Qt::black);
                return;
            }

            canvas_->set_show_cell_nuceli(false);
            canvas_->set(cz::to_cyto_frame(current_state_, params_.palette, {}));
        }

        simulation_running_ = !simulation_running_;

        if (simulation_running_) {
            advance_simulation();
        }
    }
    catch (const std::exception& ex) {
        simulation_running_ = false;
        transition_in_flight_ = false;

        QMessageBox::critical(
            this,
            tr("Run Simulation"),
            tr("Simulation failed:\n%1").arg(QString::fromUtf8(ex.what()))
        );
    }
}

void cz::main_window::advance_simulation()
{
    if (!simulation_running_ || transition_in_flight_) {
        return;
    }

    if (!params_.cell_indexer) {
        throw std::runtime_error("advance_simulation: missing cell indexer.");
    }

    if (!params_.vertex_indexer) {
        throw std::runtime_error("advance_simulation: missing vertex indexer.");
    }

    if (params_.num_states <= 0) {
        throw std::runtime_error("advance_simulation: num_states must be positive.");
    }

    if (params_.cell_state_table.size() != static_cast<std::size_t>(params_.num_states)) {
        throw std::runtime_error(
            "advance_simulation: cell_state_table row count must equal num_states."
        );
    }

    if (params_.vertex_table.size() !=
        params_.vertex_indexer->num_columns(static_cast<std::size_t>(params_.num_states))) {
        throw std::runtime_error(
            "advance_simulation: vertex_table width does not match vertex indexer."
        );
    }

    if (params_.palette.size() < static_cast<std::size_t>(params_.num_states)) {
        throw std::runtime_error("advance_simulation: palette is too small.");
    }

    const auto result = cz::apply_state_tables(
        id_source_,
        current_state_,
        params_.cell_state_table,
        params_.cell_indexer,
        params_.vertex_table,
        params_.vertex_indexer,
        params_.palette
    );

    pending_next_state_ = result.next_state;
    transition_in_flight_ = true;

    canvas_->start_transition(result.anim_start, result.anim_end);
}

void cz::main_window::on_transition_finished()
{
    for (const cz::cell_id id : canvas_->take_reclaimable_ids()) {
        id_source_.release(id);
    }

    current_state_ = std::move(pending_next_state_);
    pending_next_state_.clear();
    transition_in_flight_ = false;

    if (simulation_running_) {
        QTimer::singleShot(
            0,
            this,
            [this]() {
                advance_simulation();
            }
        );
    }
}

void cz::main_window::run_debug_demo()
{
    try {
        cz::cell_id_source ids;
        auto current_state = cz::random_cyto_state(500, 4, ids);

        for (auto& cell : current_state) {
            cell.state = 0;
            cell.phase = cz::life_stage::normal;
        }

        const std::vector<std::size_t> delete_indices = {
            5, 17, 42, 80, 111, 150, 222, 301
        };

        std::vector<cz::cell_id> delete_list;
        delete_list.reserve(delete_indices.size());

        for (std::size_t index : delete_indices) {
            if (index >= current_state.size()) {
                throw std::runtime_error(
                    "run_debug_demo: delete index out of range."
                );
            }

            delete_list.push_back(current_state[index].id);
            current_state[index].state = 2;
        }

        const auto current_sites = current_state
            | rv::transform([](const cz::cell_state& cell) -> cz::point {
            return cell.site;
                })
            | r::to<std::vector>();

        const auto diagram = cz::to_voronoi_diagram(current_sites);

        std::vector<cz::cell_state> add_list;
        add_list.reserve(8);

        constexpr std::size_t num_to_add = 8;

        if (!diagram.embedding.vertices.empty()) {
            const std::size_t step = std::max<std::size_t>(
                1,
                diagram.embedding.vertices.size() / num_to_add
            );

            for (std::size_t i = 0;
                i < diagram.embedding.vertices.size() &&
                add_list.size() < num_to_add;
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

        const auto deletion_set = delete_list | r::to<std::unordered_set>();

        cz::cyto_state next_state;
        next_state.reserve(
            current_state.size() - delete_list.size() + add_list.size()
        );

        for (const auto& cell : current_state) {
            if (deletion_set.contains(cell.id)) {
                continue;
            }

            auto survivor = cell;
            survivor.phase = cz::life_stage::normal;
            survivor.state = 0;
            next_state.push_back(survivor);
        }

        for (const auto& new_cell : add_list) {
            auto canonical_birth = new_cell;
            canonical_birth.phase = cz::life_stage::normal;
            next_state.push_back(canonical_birth);
        }

        relax_state_unweighted(next_state);

        const auto trans = cz::generate_transition(
            current_state,
            next_state,
            delete_list,
            add_list
        );

        const auto current_scale_by_id = make_scale_map(current_state);
        const auto next_scale_by_id = make_scale_map(next_state);

        const auto from_scales = make_from_scales(
            trans.from,
            current_scale_by_id,
            next_scale_by_id
        );

        const auto to_scales = make_to_scales(
            trans.to,
            current_scale_by_id,
            next_scale_by_id
        );

        cz::color_table palette = {
            {0, 255, 0},
            {255, 0, 0},
            {255, 255, 130},
            {0, 60, 200}
        };

        canvas_->set_show_cell_nuceli(false);

        const auto initial_frame = cz::to_cyto_frame(current_state, palette, {});
        const auto from_frame = cz::to_cyto_frame(trans.from, palette, from_scales);
        const auto to_frame = cz::to_cyto_frame(trans.to, palette, to_scales);

        canvas_->set(initial_frame);
        canvas_->start_transition(from_frame, to_frame);
    }
    catch (const std::exception& ex) {
        QMessageBox::critical(
            this,
            tr("Debug"),
            tr("Debug run failed:\n%1").arg(QString::fromUtf8(ex.what()))
        );
    }
}

bool cz::main_window::load_ruleset_from_file(const QString& file_path)
{
    QFile file(file_path);

    if (!file.open(QIODevice::ReadOnly)) {
        QMessageBox::critical(
            this,
            tr("Open Ruleset"),
            tr("Could not open file:\n%1").arg(file_path)
        );
        return false;
    }

    const QByteArray bytes = file.readAll();
    file.close();

    QJsonParseError parse_error;
    const QJsonDocument doc = QJsonDocument::fromJson(bytes, &parse_error);

    if (parse_error.error != QJsonParseError::NoError) {
        QMessageBox::critical(
            this,
            tr("Open Ruleset"),
            tr("Invalid JSON in file:\n%1\n\n%2")
            .arg(file_path, parse_error.errorString())
        );
        return false;
    }

    if (!doc.isObject()) {
        QMessageBox::critical(
            this,
            tr("Open Ruleset"),
            tr("Ruleset file must contain a top-level JSON object.")
        );
        return false;
    }

    current_ruleset_path_ = file_path;

    // TODO:
    // Read the JSON object and populate the current in-memory ruleset.
    [[maybe_unused]] const QJsonObject obj = doc.object();

    QMessageBox::information(
        this,
        tr("Open Ruleset"),
        tr("Ruleset loading is currently stubbed in.\n\nLoaded JSON from:\n%1")
        .arg(file_path)
    );

    return true;
}

bool cz::main_window::save_ruleset_to_file(const QString& file_path)
{
    QJsonObject obj;

    // TODO:
    // Populate obj from the current in-memory ruleset before writing.
    obj.insert("stub", true);

    const QJsonDocument doc(obj);

    QFile file(file_path);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
        QMessageBox::critical(
            this,
            tr("Save Ruleset"),
            tr("Could not write file:\n%1").arg(file_path)
        );
        return false;
    }

    file.write(doc.toJson(QJsonDocument::Indented));
    file.close();

    current_ruleset_path_ = file_path;

    QMessageBox::information(
        this,
        tr("Save Ruleset"),
        tr("Ruleset saving is currently stubbed in.\n\nWrote JSON to:\n%1")
        .arg(file_path)
    );

    return true;
}