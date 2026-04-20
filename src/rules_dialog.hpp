#pragma once

#include <QDialog>
#include <QtWidgets>

#include "cytozoic.hpp"

namespace cz {

    class rules_dialog : public QDialog {
        Q_OBJECT

        QComboBox* num_states_combo_ = nullptr;
        QSpinBox* num_initial_cells_spin_ = nullptr;
        QWidget* density_editor_ = nullptr;
        QWidget* palette_editor_widget_ = nullptr;

        QRadioButton* vertex_birth_radio_ = nullptr;
        QRadioButton* cell_birth_radio_ = nullptr;

        QComboBox* cell_indexer_combo_ = nullptr;
        QTableWidget* cell_table_ = nullptr;

        QTabWidget* tabs_ = nullptr;
        QWidget* birth_tab_ = nullptr;
        QStackedWidget* birth_mode_stack_ = nullptr;

        QWidget* vertex_birth_page_ = nullptr;
        QComboBox* vert_indexer_combo_ = nullptr;
        QTableWidget* vertex_table_ = nullptr;
        QComboBox* vertex_birth_probability_combo_ = nullptr;

        QWidget* cell_birth_page_ = nullptr;
        QComboBox* cell_spawn_site_combo_ = nullptr;
        QTableWidget* cell_birth_table_ = nullptr;

        void build_initial_tab();
        void build_cell_tab();
        void build_birth_tab();

        void initialize_from_params(const cyto_params& params);
        void connect_rebuild_signals();
        void connect_birth_mode_signals();

        void rebuild_tables();
        void refresh_birth_mode_ui();

    public:
        rules_dialog(QWidget* parent = nullptr, const cyto_params& params = {});
        cyto_params get() const;
    };

} // namespace cz