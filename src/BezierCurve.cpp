#include "geom/BezierCurve.h"

#include <algorithm>
#include <cassert>

#include "geom/CurveVisitor.h"

namespace geom {

BezierCurve::BezierCurve(std::vector<Vec3> poles) : poles_(std::move(poles)) {
    assert(!poles_.empty());
}

Result<BezierCurve> BezierCurve::create(std::vector<Vec3> poles) {
    if (poles.size() < 2) return invalidInput("Bezier curve needs at least 2 poles");
    return BezierCurve(std::move(poles));
}

Vec3 BezierCurve::deCasteljau(std::vector<Vec3> pts, double t) {
    for (std::size_t k = pts.size(); k > 1; --k)
        for (std::size_t i = 0; i + 1 < k; ++i) pts[i] = pts[i] * (1 - t) + pts[i + 1] * t;
    return pts[0];
}

Vec3 BezierCurve::point(double t) const { return deCasteljau(poles_, t); }

Vec3 BezierCurve::derivative(double t, int order) const {
    assert(order == 1 || order == 2);
    if (order > degree()) return {};
    BezierCurve d = derivativeCurve();
    return order == 1 ? d.point(t) : d.derivative(t, 1);
}

BezierCurve BezierCurve::derivativeCurve() const {
    if (poles_.size() == 1) return BezierCurve({Vec3()});
    const double n = degree();
    std::vector<Vec3> d;
    d.reserve(poles_.size() - 1);
    for (std::size_t i = 0; i + 1 < poles_.size(); ++i) d.push_back((poles_[i + 1] - poles_[i]) * n);
    return BezierCurve(std::move(d));
}

std::pair<BezierCurve, BezierCurve> BezierCurve::subdivide(double t) const {
    // The left and right edges of the de Casteljau triangle are the poles of
    // the two halves.
    std::vector<Vec3> pts = poles_;
    const std::size_t n = pts.size();
    std::vector<Vec3> left, right;
    left.reserve(n);
    right.reserve(n);
    left.push_back(pts.front());
    right.push_back(pts.back());
    for (std::size_t k = 1; k < n; ++k) {
        for (std::size_t i = 0; i + k < n; ++i) pts[i] = pts[i] * (1 - t) + pts[i + 1] * t;
        left.push_back(pts[0]);
        right.push_back(pts[n - 1 - k]);
    }
    std::reverse(right.begin(), right.end());
    return {BezierCurve(std::move(left)), BezierCurve(std::move(right))};
}

void BezierCurve::accept(CurveVisitor& visitor) const { visitor.visit(*this); }

}  // namespace geom
