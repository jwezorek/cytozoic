#include "cytozoic_widget.hpp"
#include "voronoi.hpp"
#include "util.hpp"

#include <QColor>
#include <QPaintEvent>
#include <QPainter>
#include <QPen>
#include <QResizeEvent>

#include <algorithm>
#include <ranges>
#include <span>
#include <vector>

namespace r = std::ranges;
namespace rv = std::ranges::views;

/*------------------------------------------------------------------------------------------------*/

cz::cytozoic_widget::cytozoic_widget(QWidget* parent, int duration_ms, int interval_ms) :
    QWidget(parent),
    animation_duration_ms_(duration_ms),
    animation_frame_interval_ms_(interval_ms),
    animation_elapsed_ms_(0),
    show_cell_nuclei_(false) {

    setAttribute(Qt::WA_OpaquePaintEvent);
    setAttribute(Qt::WA_NoSystemBackground);
    setAutoFillBackground(false);

    animation_timer_.setSingleShot(false);
    animation_timer_.setInterval(animation_frame_interval_ms_);

    connect(&animation_timer_, &QTimer::timeout, this, [this]() {
        if (anim_start_.size() != anim_end_.size() || anim_start_.empty()) {
            animation_timer_.stop();
            if (!anim_end_.empty()) {
                set(anim_end_);
            }
            return;
        }

        animation_elapsed_ms_ = std::min(
            animation_elapsed_ms_ + animation_frame_interval_ms_,
            animation_duration_ms_
        );

        double t = 1.0;
        if (animation_duration_ms_ > 0) {
            t = static_cast<double>(animation_elapsed_ms_) /
                static_cast<double>(animation_duration_ms_);
        }

        t = std::clamp(t, 0.0, 1.0);

        std::vector<point> sites;
        std::vector<color> colors;
        sites.reserve(anim_start_.size());
        colors.reserve(anim_start_.size());

        for (size_t i = 0; i < anim_start_.size(); ++i) {
            sites.push_back(interpolate_point(
                anim_start_[i].seed,
                anim_end_[i].seed,
                t
            ));

            colors.push_back(interpolate_color(
                anim_start_[i].color,
                anim_end_[i].color,
                t
            ));
        }

        const cz::rect bounds{ {0.0, 0.0}, {1.0, 1.0} };
        auto polys = cz::to_voronoi_polygons(sites, bounds);
        auto frame = cz::to_cyto_frame(sites, polys, colors);
        set(frame);

        if (animation_elapsed_ms_ >= animation_duration_ms_) {
            animation_timer_.stop();
            set(anim_end_);
        }
        });

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
        const double sx = framebuffer_.width();
        const double sy = framebuffer_.height();
        return QPointF(p.x * sx, p.y * sy);
        };

    QPainter painter(&framebuffer_);
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setPen(QPen(Qt::black, 1.0));

    for (const auto& c : v) {
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
    update();
}

void cz::cytozoic_widget::start_transition(const cyto_frame& from, const cyto_frame& to) {
    animation_timer_.stop();

    anim_start_ = from;
    anim_end_ = to;
    animation_elapsed_ms_ = 0;

    if (anim_start_.size() != anim_end_.size()) {
        set(anim_end_);
        return;
    }

    if (anim_end_.empty()) {
        clear(Qt::black);
        return;
    }

    if (animation_duration_ms_ <= 0 || animation_frame_interval_ms_ <= 0) {
        set(anim_end_);
        return;
    }

    animation_timer_.setInterval(animation_frame_interval_ms_);

    set(anim_start_);
    animation_timer_.start();
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