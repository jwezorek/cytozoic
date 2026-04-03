#include "cytozoic_widget.hpp"
#include "voronoi.hpp"

#include <QColor>
#include <QPaintEvent>
#include <QPainter>
#include <QPen>
#include <QResizeEvent>

#include <ranges>

namespace r = std::ranges;
namespace rv = std::ranges::views;

/*------------------------------------------------------------------------------------------------*/

cz::cytozoic_widget::cytozoic_widget(QWidget* parent, double log_wd, double log_hgt) : 
        QWidget(parent), 
        logical_wd_(log_wd),
        logical_hgt_(log_hgt),
        show_cell_nuclei_(false) {
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

void cz::cytozoic_widget::set(const cyto_frame& v) {
    ensure_framebuffer_matches_widget_size();

    if (framebuffer_.isNull()) {
        return;
    }

    framebuffer_.fill(Qt::black);

    auto to_pixel = [this](const cz::point& p) -> QPointF {
        const double sx = framebuffer_.width() / logical_wd_;
        const double sy = framebuffer_.height() / logical_hgt_;
        return QPointF(p.x * sx, p.y * sy);
        };

    QPainter painter(&framebuffer_);
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setPen(QPen(Qt::black, 1.0));

    
    for (const auto& c : v ) {

        if (c.shape.empty()) {
            continue;
        }

        QPolygonF poly;
        poly.reserve(static_cast<int>(c.shape.size()));

        for (const auto& p : c.shape) {
            poly.push_back(to_pixel(p));
        }

        painter.setBrush(QColor(c.color.r, c.color.g, c.color.b));
        painter.drawPolygon(poly);
    }

    if (show_cell_nuclei_) {
        painter.setPen(Qt::NoPen);
        painter.setBrush(Qt::white);

        constexpr double nucleus_radius = 3.0;

        for (const auto& c : v) {
            const QPointF center = to_pixel(c.seed);
            painter.drawEllipse(center, nucleus_radius, nucleus_radius);
        }
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

    painter.drawImage(rect(), framebuffer_, framebuffer_.rect());
}

void cz::cytozoic_widget::resizeEvent(QResizeEvent* event) {
    QWidget::resizeEvent(event);
    ensure_framebuffer_matches_widget_size();
}

bool cz::cytozoic_widget::show_cell_nuclei() const {
    return show_cell_nuclei_;
}

void cz::cytozoic_widget::set_show_cell_nuceli(bool show) {
    show_cell_nuclei_ = show;
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