#include <cassert>
#include <cmath>

#include "geom/Circle.h"
#include "geom/CurveVisitor.h"
#include "geom/Line.h"

namespace geom {

// ---- Line ----

Line::Line(const Vec3& p0, const Vec3& p1) : p0_(p0), p1_(p1) {
    assert(!nearlyEqual(p0, p1) && "Line endpoints must be distinct");
}

Result<Line> Line::create(const Vec3& p0, const Vec3& p1) {
    if (nearlyEqual(p0, p1)) return Error{ErrorCode::Degenerate, "line endpoints coincide"};
    return Line(p0, p1);
}

Vec3 Line::derivative(double, int order) const {
    assert(order == 1 || order == 2);
    return order == 1 ? p1_ - p0_ : Vec3();
}

void Line::accept(CurveVisitor& visitor) const { visitor.visit(*this); }

// ---- Circle ----

Circle::Circle(const Vec3& center, const Vec3& normal, const Vec3& xAxis, double radius,
               double start, double end)
    : c_(center), r_(radius), dom_{start, end} {
    const auto n = normal.normalized();
    assert(n && "Circle normal must be non-zero");
    n_ = *n;
    const auto x = (xAxis - n_ * xAxis.dot(n_)).normalized();
    assert(x && "Circle xAxis must not be parallel to the normal");
    x_ = *x;
    y_ = n_.cross(x_);
    assert(radius > 0 && end > start);
}

Result<Circle> Circle::create(const Vec3& center, const Vec3& normal, const Vec3& xAxis,
                              double radius, double start, double end) {
    const auto n = normal.normalized();
    if (!n) return Error{ErrorCode::Degenerate, "circle normal is zero"};
    if (!(xAxis - *n * xAxis.dot(*n)).normalized())
        return Error{ErrorCode::Degenerate, "circle xAxis is parallel to the normal"};
    if (radius <= kLinear) return invalidInput("circle radius must be positive");
    if (end <= start) return invalidInput("arc end angle must exceed start angle");
    if (end - start > 2 * kPi + kAngular) return invalidInput("arc spans more than 2 pi");
    return Circle(center, normal, xAxis, radius, start, end);
}

Vec3 Circle::point(double t) const {
    return c_ + (x_ * std::cos(t) + y_ * std::sin(t)) * r_;
}

Vec3 Circle::derivative(double t, int order) const {
    assert(order == 1 || order == 2);
    const double c = std::cos(t), s = std::sin(t);
    if (order == 1) return (y_ * c - x_ * s) * r_;
    return (x_ * c + y_ * s) * -r_;
}

void Circle::accept(CurveVisitor& visitor) const { visitor.visit(*this); }

}  // namespace geom
