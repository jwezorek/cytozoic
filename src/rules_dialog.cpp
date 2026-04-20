#include "rules_dialog.hpp"
#include "cytozoic.hpp"

#include <QButtonGroup>
#include <QColorDialog>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QGridLayout>
#include <QGroupBox>
#include <QHeaderView>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QRandomGenerator>
#include <QRadioButton>
#include <QSignalBlocker>
#include <QSlider>
#include <QSpinBox>
#include <QStackedWidget>
#include <QTabWidget>
#include <QTableWidget>
#include <QVBoxLayout>
#include <QWidget>

#include <algorithm>
#include <cmath>
#include <ranges>
#include <vector>

namespace {

    constexpr int k_min_cells = 20;
    constexpr int k_max_cells = 10000;
    int g_cell_death_probability_percent = 25;

    class palette_editor : public QWidget {
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
                    const int hue = static_cast<int>(
                        (360 * i) / std::max<std::size_t>(1, buttons_.size())
                        );
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
                result.push_back(
                    cz::color{
                        static_cast<uint8_t>(color.red()),
                        static_cast<uint8_t>(color.green()),
                        static_cast<uint8_t>(color.blue())
                    }
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

    class state_density_editor : public QWidget {
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

    int safe_cell_value(
        const QTableWidget* table,
        int row,
        int col,
        int fallback = -1)
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

        table->horizontalHeader()->setSectionResizeMode(
            QHeaderView::ResizeToContents
        );
        table->verticalHeader()->setSectionResizeMode(
            QHeaderView::ResizeToContents
        );
    }

    void set_state_table_values(
        QTableWidget* table,
        const cz::state_table& values)
    {
        const int row_count = std::min(
            table->rowCount(),
            static_cast<int>(values.size())
        );

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
            cz::state_table_row(
                static_cast<std::size_t>(table->columnCount()),
                -1
            )
        );

        for (int r = 0; r < table->rowCount(); ++r) {
            for (int c = 0; c < table->columnCount(); ++c) {
                result[r][c] = static_cast<int8_t>(
                    safe_cell_value(table, r, c, -1)
                    );
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

    cz::neighborhood_indexer make_indexer_from_dialog_name(const QString& name)
    {
        const std::string value = name.toStdString();
        return cz::indexer_from_name(value);
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

    void randomize_cell_table(
        QTableWidget* table,
        int num_states,
        int death_probability_percent)
    {
        if (table == nullptr) {
            return;
        }

        death_probability_percent = std::clamp(
            death_probability_percent,
            0,
            100
        );

        for (int r = 0; r < table->rowCount(); ++r) {
            for (int c = 0; c < table->columnCount(); ++c) {
                auto* spin = qobject_cast<QSpinBox*>(table->cellWidget(r, c));
                if (spin == nullptr) {
                    continue;
                }

                const bool make_dead =
                    QRandomGenerator::global()->bounded(100) <
                    death_probability_percent;

                if (make_dead) {
                    spin->setValue(-1);
                }
                else {
                    spin->setValue(
                        QRandomGenerator::global()->bounded(
                            std::max(1, num_states)
                        )
                    );
                }
            }
        }
    }

    void refresh_vertex_birth_count_combo(QComboBox* combo, int column_count)
    {
        if (combo == nullptr) {
            return;
        }

        const int clamped_old_value = std::clamp(
            combo->currentData().toInt(),
            0,
            std::max(0, column_count)
        );

        QSignalBlocker blocker(combo);
        combo->clear();

        for (int i = 0; i <= column_count; ++i) {
            combo->addItem(QString::number(i), i);
        }

        combo->setCurrentIndex(combo->findData(clamped_old_value));
    }

    void randomize_vertex_table_impl(
        QTableWidget* table,
        int num_states,
        int num_live_entries)
    {
        if (table == nullptr || table->rowCount() < 1) {
            return;
        }

        const int column_count = table->columnCount();
        num_live_entries = std::clamp(num_live_entries, 0, column_count);

        for (int c = 0; c < column_count; ++c) {
            if (auto* spin = qobject_cast<QSpinBox*>(table->cellWidget(0, c))) {
                spin->setValue(-1);
            }
        }

        if (num_states <= 1 || num_live_entries <= 0) {
            return;
        }

        std::vector<int> columns;
        columns.reserve(static_cast<std::size_t>(column_count));
        for (int c = 0; c < column_count; ++c) {
            columns.push_back(c);
        }

        for (int i = column_count - 1; i > 0; --i) {
            const int j = QRandomGenerator::global()->bounded(i + 1);
            std::swap(columns[i], columns[j]);
        }

        for (int i = 0; i < num_live_entries; ++i) {
            const int c = columns[static_cast<std::size_t>(i)];
            if (auto* spin = qobject_cast<QSpinBox*>(table->cellWidget(0, c))) {
                spin->setValue(
                    1 + QRandomGenerator::global()->bounded(num_states - 1)
                );
            }
        }
    }

    QString center_type_to_string(cz::center_type type)
    {
        return { cz::center_type_to_string(type).c_str() };
    }

    cz::center_type center_type_from_combo(const QComboBox* combo)
    {
        const QString value = combo->currentData().toString();
        return cz::center_type_from_string(value.toStdString());
    }

} // namespace

void cz::rules_dialog::refresh_birth_mode_ui()
{
    if (birth_mode_stack_ == nullptr) {
        return;
    }

    if (vertex_birth_radio_ != nullptr && vertex_birth_radio_->isChecked()) {
        birth_mode_stack_->setCurrentWidget(vertex_birth_page_);
        if (tabs_ != nullptr && birth_tab_ != nullptr) {
            tabs_->setTabText(tabs_->indexOf(birth_tab_), tr("Vertex Birth"));
        }
    }
    else {
        birth_mode_stack_->setCurrentWidget(cell_birth_page_);
        if (tabs_ != nullptr && birth_tab_ != nullptr) {
            tabs_->setTabText(tabs_->indexOf(birth_tab_), tr("Cell Birth"));
        }
    }
}

void cz::rules_dialog::rebuild_tables()
{
    const int num_states = num_states_combo_->currentData().toInt();
    const int min_value = -1;
    const int max_value = num_states - 1;

    static_cast<state_density_editor*>(density_editor_)->set_num_states(num_states);
    static_cast<palette_editor*>(palette_editor_widget_)->set_num_states(num_states);

    const auto cell_indexer =
        make_indexer_from_dialog_name(cell_indexer_combo_->currentText());

    const int cell_columns =
        static_cast<int>(
            cell_indexer->num_columns(static_cast<std::size_t>(num_states))
            );

    rebuild_state_table(
        cell_table_,
        num_states,
        cell_columns,
        min_value,
        max_value
    );

    if (vert_indexer_combo_ != nullptr) {
        const auto vertex_indexer =
            make_indexer_from_dialog_name(vert_indexer_combo_->currentText());

        const int vertex_columns =
            static_cast<int>(
                vertex_indexer->num_columns(static_cast<std::size_t>(num_states))
                );

        rebuild_state_table(
            vertex_table_,
            1,
            vertex_columns,
            min_value,
            max_value
        );

        refresh_vertex_birth_count_combo(
            vertex_birth_probability_combo_,
            vertex_columns
        );
    }

    if (cell_birth_table_ != nullptr) {
        rebuild_state_table(
            cell_birth_table_,
            num_states,
            1,
            min_value,
            max_value
        );
    }

    refresh_birth_mode_ui();
}

cz::rules_dialog::rules_dialog(QWidget* parent, const cyto_params& params)
    : QDialog(parent)
{
    setWindowTitle(tr("Current Rules"));
    setModal(true);
    resize(1100, 700);

    auto* root_layout = new QVBoxLayout(this);
    tabs_ = new QTabWidget(this);
    root_layout->addWidget(tabs_);

    auto* buttons = new QDialogButtonBox(
        QDialogButtonBox::Ok | QDialogButtonBox::Cancel,
        this
    );
    root_layout->addWidget(buttons);

    connect(buttons, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);

    //----------------------------------------------------------------------
    // Initial configuration tab
    //----------------------------------------------------------------------

    auto* initial_tab = new QWidget(tabs_);
    auto* initial_layout = new QVBoxLayout(initial_tab);
    auto* initial_grid = new QGridLayout();
    initial_layout->addLayout(initial_grid);

    auto* num_states_label = new QLabel(tr("Number of states"), initial_tab);
    num_states_combo_ = new QComboBox(initial_tab);
    for (int i = 1; i <= 10; ++i) {
        num_states_combo_->addItem(QString::number(i), i);
    }

    auto* num_initial_cells_label = new QLabel(tr("Initial cell count"), initial_tab);
    num_initial_cells_spin_ = new QSpinBox(initial_tab);
    num_initial_cells_spin_->setRange(k_min_cells, k_max_cells);

    auto* state_density_label = new QLabel(tr("State density"), initial_tab);
    density_editor_ = new state_density_editor(initial_tab);
    auto* density_editor = static_cast<state_density_editor*>(density_editor_);

    auto* palette_label = new QLabel(tr("Palette"), initial_tab);
    palette_editor_widget_ = new palette_editor(initial_tab);
    auto* palette_editor =
        static_cast<::palette_editor*>(palette_editor_widget_);

    auto* birth_mode_label = new QLabel(tr("Birth mode"), initial_tab);
    auto* birth_mode_widget = new QWidget(initial_tab);
    auto* birth_mode_layout = new QVBoxLayout(birth_mode_widget);
    birth_mode_layout->setContentsMargins(0, 0, 0, 0);

    vertex_birth_radio_ =
        new QRadioButton(tr("Vertex-based birth"), birth_mode_widget);
    cell_birth_radio_ =
        new QRadioButton(tr("Cell-based birth"), birth_mode_widget);

    auto* birth_mode_group = new QButtonGroup(this);
    birth_mode_group->addButton(vertex_birth_radio_);
    birth_mode_group->addButton(cell_birth_radio_);

    birth_mode_layout->addWidget(vertex_birth_radio_);
    birth_mode_layout->addWidget(cell_birth_radio_);

    initial_grid->addWidget(num_states_label, 0, 0);
    initial_grid->addWidget(num_states_combo_, 0, 1);
    initial_grid->addWidget(num_initial_cells_label, 1, 0);
    initial_grid->addWidget(num_initial_cells_spin_, 1, 1);
    initial_grid->addWidget(state_density_label, 2, 0, Qt::AlignTop);
    initial_grid->addWidget(density_editor, 2, 1);
    initial_grid->addWidget(palette_label, 3, 0, Qt::AlignTop);
    initial_grid->addWidget(palette_editor, 3, 1);
    initial_grid->addWidget(birth_mode_label, 4, 0, Qt::AlignTop);
    initial_grid->addWidget(birth_mode_widget, 4, 1);

    initial_layout->addStretch();
    tabs_->addTab(initial_tab, tr("Initial Configuration"));

    //----------------------------------------------------------------------
    // Cell table tab
    //----------------------------------------------------------------------

    auto* cell_tab = new QWidget(tabs_);
    auto* cell_layout = new QVBoxLayout(cell_tab);
    auto* cell_top_row = new QHBoxLayout();

    auto* cell_indexer_label = new QLabel(tr("Indexer"), cell_tab);
    cell_indexer_combo_ = new QComboBox(cell_tab);
    for (const auto& name : cz::named_indexers()) {
        cell_indexer_combo_->addItem(QString::fromStdString(name));
    }

    cell_top_row->addWidget(cell_indexer_label);
    cell_top_row->addWidget(cell_indexer_combo_);
    cell_top_row->addStretch();

    cell_table_ = new QTableWidget(cell_tab);
    cell_table_->setAlternatingRowColors(true);

    auto* cell_random_group = new QGroupBox(tr("Random Generation"), cell_tab);
    auto* cell_random_layout = new QHBoxLayout(cell_random_group);

    auto* cell_death_probability_label =
        new QLabel(tr("Death probability"), cell_random_group);
    auto* cell_death_probability_spin = new QSpinBox(cell_random_group);
    cell_death_probability_spin->setRange(0, 100);
    cell_death_probability_spin->setSuffix(tr("%"));
    cell_death_probability_spin->setValue(g_cell_death_probability_percent);

    auto* cell_randomize_button =
        new QPushButton(tr("Generate Random States"), cell_random_group);

    cell_random_layout->addWidget(cell_death_probability_label);
    cell_random_layout->addWidget(cell_death_probability_spin);
    cell_random_layout->addSpacing(16);
    cell_random_layout->addWidget(cell_randomize_button);
    cell_random_layout->addStretch();

    cell_layout->addLayout(cell_top_row);
    cell_layout->addWidget(cell_table_);
    cell_layout->addWidget(cell_random_group);
    tabs_->addTab(cell_tab, tr("Cell Table"));

    //----------------------------------------------------------------------
    // Birth tab
    //----------------------------------------------------------------------

    birth_tab_ = new QWidget(tabs_);
    auto* birth_layout = new QVBoxLayout(birth_tab_);

    birth_mode_stack_ = new QStackedWidget(birth_tab_);

    //----------------------------------------------------------------------
    // Vertex-based birth page
    //----------------------------------------------------------------------

    vertex_birth_page_ = new QWidget(birth_mode_stack_);
    auto* vertex_layout = new QVBoxLayout(vertex_birth_page_);
    auto* vertex_top_row = new QHBoxLayout();

    auto* vertex_indexer_label = new QLabel(tr("Indexer"), vertex_birth_page_);
    vert_indexer_combo_ = new QComboBox(vertex_birth_page_);
    for (const auto& name : cz::named_indexers()) {
        vert_indexer_combo_->addItem(QString::fromStdString(name));
    }

    vertex_top_row->addWidget(vertex_indexer_label);
    vertex_top_row->addWidget(vert_indexer_combo_);
    vertex_top_row->addStretch();

    vertex_table_ = new QTableWidget(vertex_birth_page_);
    vertex_table_->setAlternatingRowColors(true);

    auto* vertex_random_group =
        new QGroupBox(tr("Random Generation"), vertex_birth_page_);
    auto* vertex_random_layout = new QHBoxLayout(vertex_random_group);

    auto* vertex_birth_probability_label =
        new QLabel(tr("Birth probability"), vertex_random_group);
    vertex_birth_probability_combo_ = new QComboBox(vertex_random_group);

    auto* vertex_randomize_button =
        new QPushButton(tr("Generate Random States"), vertex_random_group);

    vertex_random_layout->addWidget(vertex_birth_probability_label);
    vertex_random_layout->addWidget(vertex_birth_probability_combo_);
    vertex_random_layout->addSpacing(16);
    vertex_random_layout->addWidget(vertex_randomize_button);
    vertex_random_layout->addStretch();

    vertex_layout->addLayout(vertex_top_row);
    vertex_layout->addWidget(vertex_table_);
    vertex_layout->addWidget(vertex_random_group);

    //----------------------------------------------------------------------
    // Cell-based birth page
    //----------------------------------------------------------------------

    cell_birth_page_ = new QWidget(birth_mode_stack_);
    auto* cell_birth_layout = new QVBoxLayout(cell_birth_page_);
    auto* cell_birth_top_grid = new QGridLayout();

    auto* cell_spawn_site_label =
        new QLabel(tr("Spawn site"), cell_birth_page_);
    cell_spawn_site_combo_ = new QComboBox(cell_birth_page_);
    cell_spawn_site_combo_->addItem(tr("Incircle"), QStringLiteral("incircle"));
    cell_spawn_site_combo_->addItem(
        tr("Johnson ellipse"),
        QStringLiteral("johnson_ellipse")
    );
    cell_spawn_site_combo_->addItem(
        tr("Center of mass"),
        QStringLiteral("center_of_mass")
    );

    auto* cell_birth_info = new QLabel(
        tr("Cell-based birth is not yet implemented in the simulation core. "
            "These controls are present so the dialog matches the new params "
            "layout and can round-trip the new variant type."),
        cell_birth_page_
    );
    cell_birth_info->setWordWrap(true);

    cell_birth_table_ = new QTableWidget(cell_birth_page_);
    cell_birth_table_->setAlternatingRowColors(true);

    cell_birth_top_grid->addWidget(cell_spawn_site_label, 0, 0);
    cell_birth_top_grid->addWidget(cell_spawn_site_combo_, 0, 1);

    cell_birth_layout->addLayout(cell_birth_top_grid);
    cell_birth_layout->addWidget(cell_birth_info);
    cell_birth_layout->addWidget(cell_birth_table_);

    birth_mode_stack_->addWidget(vertex_birth_page_);
    birth_mode_stack_->addWidget(cell_birth_page_);
    birth_layout->addWidget(birth_mode_stack_);

    tabs_->addTab(birth_tab_, tr("Birth"));

    //----------------------------------------------------------------------
    // Initialization from current params
    //----------------------------------------------------------------------

    num_states_combo_->setCurrentIndex(
        std::max(0, num_states_combo_->findData(params.num_states))
    );
    num_initial_cells_spin_->setValue(params.num_initial_cells);

    density_editor->set_num_states(params.num_states);
    density_editor->set_density(params.initial_state_density);

    palette_editor->set_num_states(params.num_states);
    palette_editor->set_palette(params.palette);

    cell_indexer_combo_->setCurrentText(
        indexer_name_or_default(
            params.cell_indexer,
            QStringLiteral("sum of states")
        )
    );

    if (const auto* vertex_birth =
        std::get_if<cz::vertex_based_birth>(&params.birth_params)) {
        vertex_birth_radio_->setChecked(true);

        vert_indexer_combo_->setCurrentText(
            indexer_name_or_default(
                vertex_birth->vertex_indexer,
                QStringLiteral("sum of states")
            )
        );
    }
    else if (const auto* cell_birth =
        std::get_if<cz::cell_based_birth>(&params.birth_params)) {
        cell_birth_radio_->setChecked(true);
        cell_spawn_site_combo_->setCurrentText(
            cz::center_type_to_string(cell_birth->spawn_site).c_str()
        );
    }
    else {
        vertex_birth_radio_->setChecked(true);
    }

    rebuild_tables();
    set_state_table_values(cell_table_, params.cell_state_table);

    if (const auto* vertex_birth =
        std::get_if<cz::vertex_based_birth>(&params.birth_params)) {
        set_state_table_row_values(
            vertex_table_,
            vertex_birth->vertex_table
        );
    }

    if (const auto* cell_birth =
        std::get_if<cz::cell_based_birth>(&params.birth_params)) {
        set_state_table_values(
            cell_birth_table_,
            cell_birth->state_table
        );
    }

    refresh_birth_mode_ui();

    connect(
        num_states_combo_,
        &QComboBox::currentIndexChanged,
        this,
        [&](int) {
            rebuild_tables();
        }
    );

    connect(
        cell_indexer_combo_,
        &QComboBox::currentIndexChanged,
        this,
        [&](int) {
            rebuild_tables();
        }
    );

    connect(
        vert_indexer_combo_,
        &QComboBox::currentIndexChanged,
        this,
        [&](int) {
            rebuild_tables();
        }
    );

    connect(
        vertex_birth_radio_,
        &QRadioButton::toggled,
        this,
        [&](bool) {
            refresh_birth_mode_ui();
        }
    );

    connect(
        cell_birth_radio_,
        &QRadioButton::toggled,
        this,
        [&](bool) {
            refresh_birth_mode_ui();
        }
    );

    connect(
        cell_death_probability_spin,
        qOverload<int>(&QSpinBox::valueChanged),
        this,
        [](int value) {
            g_cell_death_probability_percent = value;
        }
    );

    connect(
        cell_randomize_button,
        &QPushButton::clicked,
        this,
        [this, cell_death_probability_spin]() {
            g_cell_death_probability_percent =
                cell_death_probability_spin->value();

            randomize_cell_table(
                cell_table_,
                num_states_combo_->currentData().toInt(),
                g_cell_death_probability_percent
            );
        }
    );

    connect(
        vertex_randomize_button,
        &QPushButton::clicked,
        this,
        [this]() {
            randomize_vertex_table_impl(
                vertex_table_,
                num_states_combo_->currentData().toInt(),
                vertex_birth_probability_combo_->currentData().toInt()
            );
        }
    );
}

cz::cyto_params cz::rules_dialog::get() const
{
    cz::cyto_params edited;

    edited.num_states = num_states_combo_->currentData().toInt();
    edited.num_initial_cells = num_initial_cells_spin_->value();
    edited.initial_state_density =
        static_cast<state_density_editor*>(density_editor_)->density();
    edited.palette =
        static_cast<palette_editor*>(palette_editor_widget_)->palette();

    edited.cell_indexer =
        make_indexer_from_dialog_name(cell_indexer_combo_->currentText());

    edited.cell_state_table = read_state_table(cell_table_);

    if (vertex_birth_radio_ != nullptr && vertex_birth_radio_->isChecked()) {
        edited.birth_params = cz::vertex_based_birth{
            .vertex_indexer =
                make_indexer_from_dialog_name(vert_indexer_combo_->currentText()),
            .vertex_table = read_state_table_row(vertex_table_)
        };
    }
    else {
        edited.birth_params = cz::cell_based_birth{
            .spawn_site = center_type_from_combo(cell_spawn_site_combo_),
            .state_table = read_state_table(cell_birth_table_)
        };
    }

    return edited;
}