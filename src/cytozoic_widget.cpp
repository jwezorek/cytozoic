#include "cytozoic_widget.hpp"
#include "voronoi.hpp"

#include <QColor>
#include <QPaintEvent>
#include <QPainter>
#include <QPen>
#include <QResizeEvent>

cz::cytozoic_widget::cytozoic_widget(QWidget* parent, double log_wd, double log_hgt) : 
        QWidget(parent), 
        logical_wd_(log_wd),
        logical_hgt_(log_hgt) {
    // These help Qt treat this as a fully self-painted widget and reduce flicker.
    setAttribute(Qt::WA_OpaquePaintEvent);
    setAttribute(Qt::WA_NoSystemBackground);
    setAutoFillBackground(false);

    ensure_framebuffer_matches_widget_size();
    clear(Qt::black);
}

QSize cz::cytozoic_widget::minimumSizeHint() const {
    return { 64, 64 };
}

QSize cz::cytozoic_widget::sizeHint() const {
    return { 640, 480 };
}

void cz::cytozoic_widget::clear(const QColor& color) {
    if (framebuffer_.isNull()) {
        return;
    }

    framebuffer_.fill(color);
    update();
}

void cz::cytozoic_widget::resize_framebuffer(int width, int height) {
    if (width <= 0 || height <= 0) {
        framebuffer_ = QImage();
        update();
        return;
    }

    QImage new_buffer(width, height, QImage::Format_ARGB32_Premultiplied);
    new_buffer.fill(Qt::black);

    if (!framebuffer_.isNull()) {
        QPainter painter(&new_buffer);
        painter.drawImage(0, 0, framebuffer_);
    }

    framebuffer_ = std::move(new_buffer);
    update();
}

const QImage& cz::cytozoic_widget::framebuffer() const noexcept {
    return framebuffer_;
}

void cz::cytozoic_widget::set(const voronoi_diagram& v) {
    ensure_framebuffer_matches_widget_size();

    if (framebuffer_.isNull()) {
        return;
    }

    framebuffer_.fill(Qt::black);

    QPainter painter(&framebuffer_);
    painter.setRenderHint(QPainter::Antialiasing, true);

    const double sx = logical_wd_ != 0.0
        ? static_cast<double>(framebuffer_.width()) / logical_wd_
        : 1.0;

    const double sy = logical_hgt_ != 0.0
        ? static_cast<double>(framebuffer_.height()) / logical_hgt_
        : 1.0;

    auto to_qpointf = [&](const cz::point& p) -> QPointF {
        return {
            p.x * sx,
            p.y * sy
        };
        };

    QPen edge_pen(QColor(200, 200, 200));
    edge_pen.setWidthF(1.0);
    painter.setPen(edge_pen);
    painter.setBrush(Qt::NoBrush);

    for (const auto& cell : v) {
        if (cell.cell.empty()) {
            continue;
        }

        QPolygonF poly;
        poly.reserve(static_cast<int>(cell.cell.size()));

        for (const auto& p : cell.cell) {
            poly << to_qpointf(p);
        }

        painter.drawPolygon(poly);
    }

    painter.setPen(Qt::NoPen);
    painter.setBrush(QColor(255, 80, 80));

    constexpr double site_radius = 2.5;

    for (const auto& cell : v) {
        const auto q = to_qpointf(cell.site);
        painter.drawEllipse(q, site_radius, site_radius);
    }

    update();
}

void cz::cytozoic_widget::paintEvent(QPaintEvent* event) {
    Q_UNUSED(event);

    QPainter painter(this);
    painter.setRenderHint(QPainter::SmoothPixmapTransform, false);

    if (framebuffer_.isNull()) {
        painter.fillRect(rect(), Qt::black);
        return;
    }

    // Draw the backing image directly to the widget.
    // Right now this stretches to fit the widget exactly.
    // If later you want 1:1 pixels, integer scaling, camera transforms, etc.,
    // this is the place to change it.
    painter.drawImage(rect(), framebuffer_, framebuffer_.rect());
}

void cz::cytozoic_widget::resizeEvent(QResizeEvent* event) {
    QWidget::resizeEvent(event);
    ensure_framebuffer_matches_widget_size();
}

void cz::cytozoic_widget::ensure_framebuffer_matches_widget_size() {
    const QSize s = size();
    if (s.width() <= 0 || s.height() <= 0) {
        return;
    }

    if (framebuffer_.size() == s) {
        return;
    }

    resize_framebuffer(s.width(), s.height());
}