#pragma once

#include <QImage>
#include <QPoint>
#include <QSize>
#include <QWidget>
#include "geometric_types.hpp"

namespace cz {

    class cytozoic_widget : public QWidget {
        Q_OBJECT

    public:
        explicit cytozoic_widget(QWidget* parent = nullptr,
            double log_wd = 1.0, double log_hgt = 1.0);

        QSize minimumSizeHint() const override;
        QSize sizeHint() const override;

        void clear(const QColor& color = Qt::black);
        void resize_framebuffer(int width, int height);

        const QImage& framebuffer() const noexcept;

        void set(const voronoi_diagram& v);

    protected:
        void paintEvent(QPaintEvent* event) override;
        void resizeEvent(QResizeEvent* event) override;

    private:
        void ensure_framebuffer_matches_widget_size();

        QImage framebuffer_;
        double logical_wd_;
        double logical_hgt_;
    };

}