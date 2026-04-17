#include <QDialog>
#include <QtWidgets>
#include "cytozoic.hpp"

namespace cz {

    class rules_dialog : public QDialog {

        QComboBox* num_states_combo_;
        QSpinBox* num_initial_cells_spin_;
        QWidget* density_editor_;
        QWidget* palette_editor_widget_;
        QComboBox* cell_indexer_combo_;
        QComboBox* vert_indexer_combo_;
        QTableWidget* cell_table_;
        QTableWidget* vertex_table_;

        void rebuild_tables();

    public:

        rules_dialog(QWidget* parent = nullptr, const cyto_params& params = {});
        cyto_params get() const;
    };


}