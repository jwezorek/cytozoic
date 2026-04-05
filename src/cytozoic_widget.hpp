#pragma once

#include <QImage>
#include <QPoint>
#include <QSize>
#include <QTimer>
#include <QWidget>
#include "types.hpp"

namespace cz {

    class cytozoic_widget : public QWidget {
        Q_OBJECT

    public:
        explicit cytozoic_widget(
            QWidget* parent = nullptr,
            int animation_duration_ms = 400,
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

    protected:
        void paintEvent(QPaintEvent* event) override;
        void resizeEvent(QResizeEvent* event) override;

    private:
        void ensure_framebuffer_matches_widget_size();

        QImage framebuffer_;
        int animation_duration_ms_;
        int animation_frame_interval_ms_;
        int animation_elapsed_ms_;
        cyto_frame anim_start_;
        cyto_frame anim_end_;
        bool show_cell_nuclei_;
        QTimer animation_timer_;
    };

}