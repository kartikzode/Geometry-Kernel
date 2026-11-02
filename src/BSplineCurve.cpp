#include "geom/BSplineCurve.h"

#include <algorithm>
#include <cassert>

#include "geom/BSplineBasis.h"
#include "geom/CurveVisitor.h"

namespace geom {

BSplineCurve::BSplineCurve(int degree, std::vector<double> knots, std::vector<Vec3> poles)
    : degree_(degree), knots_(std::move(knots)), poles_(std::move(poles)) {
    assert(!bspline::validate(degree_, knots_, poles_.size()));
}

Result<BSplineCurve> BSplineCurve::create(int degree, std::vector<double> knots,
                                          std::vector<Vec3> poles) {
    if (auto err = bspline::validate(degree, knots, poles.size())) return invalidInput(*err);
    return BSplineCurve(degree, std::move(knots), std::move(poles));
}

Vec3 BSplineCurve::evaluate(double t, int order) const {
    const int span = bspline::findSpan(degree_, knots_, poles_.size(), t);
    const auto N = bspline::basisDerivatives(span, t, degree_, order, knots_);
    Vec3 sum;
    for (int j = 0; j <= degree_; ++j) sum += poles_[span - degree_ + j] * N(order, j);
    return sum;
}

Vec3 BSplineCurve::point(double t) const { return evaluate(domain().clamp(t), 0); }

Vec3 BSplineCurve::derivative(double t, int order) const {
    assert(order == 1 || order == 2);
    return evaluate(domain().clamp(t), order);
}

Interval BSplineCurve::domain() const { return bspline::domain(degree_, knots_); }

std::vector<double> BSplineCurve::breakpoints() const {
    return bspline::breakpoints(degree_, knots_);
}

std::optional<BSplineCurve> BSplineCurve::insertKnot(double u) const {
    if (!domain().contains(u)) return std::nullopt;
    const int p = degree_;
    const int k = bspline::findSpan(p, knots_, poles_.size(), u);
    const int s = static_cast<int>(std::count(knots_.begin(), knots_.end(), u));
    if (s >= p) return std::nullopt;

    const int n = static_cast<int>(poles_.size()) - 1;
    std::vector<Vec3> q(poles_.size() + 1);
    for (int i = 0; i <= k - p; ++i) q[i] = poles_[i];
    for (int i = k - s; i <= n; ++i) q[i + 1] = poles_[i];
    for (int i = k - p + 1; i <= k - s; ++i) {
        const double alpha = (u - knots_[i]) / (knots_[i + p] - knots_[i]);
        q[i] = poles_[i] * alpha + poles_[i - 1] * (1 - alpha);
    }
    std::vector<double> knots = knots_;
    knots.insert(knots.begin() + k + 1, u);
    return BSplineCurve(p, std::move(knots), std::move(q));
}

void BSplineCurve::accept(CurveVisitor& visitor) const { visitor.visit(*this); }

}  // namespace geom
