#pragma once

#include <QImage>
#include <QPoint>
#include <QSize>
#include <QTimer>
#include <QWidget>

#include "cytozoic.hpp"

namespace cz {

    class cytozoic_widget : public QWidget {
        Q_OBJECT

    public:
        explicit cytozoic_widget(
            QWidget* parent = nullptr,
            int animation_duration_ms = 500,
            int animation_frame_interval_ms = 16
        );

        QSize minimumSizeHint() const override;
        QSize sizeHint() const override;

        void clear(const QColor& color = Qt::black);
        void resize_framebuffer(int width, int height);

        const QImage& framebuffer() const noexcept;

        void set(const cyto_frame& cs);
        bool show_cell_nuclei() const;
        void set_show_cell_nuceli(bool show);

        void start_transition(
            const cyto_frame& from,
            const cyto_frame& to
        );
        void cancel_transition();
        std::vector<cell_id> take_reclaimable_ids();

    signals:
        void transition_finished();

    protected:
        void paintEvent(QPaintEvent* event) override;
        void resizeEvent(QResizeEvent* event) override;

    private:
        void ensure_framebuffer_matches_widget_size();
        void advance_animation();
        void finish_transition();

        QImage framebuffer_;
        int animation_duration_ms_;
        int animation_frame_interval_ms_;
        int animation_elapsed_ms_;
        cyto_frame anim_start_;
        cyto_frame anim_end_;
        bool show_cell_nuclei_;
        QTimer animation_timer_;
        std::vector<cell_id> reclaimable_ids_;
    };

} // namespace cz