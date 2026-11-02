#include "geom/NurbsCurve.h"

#include <cassert>
#include <cmath>

#include "geom/BSplineBasis.h"
#include "geom/CurveVisitor.h"

namespace geom {

NurbsCurve::NurbsCurve(int degree, std::vector<double> knots, std::vector<Vec3> poles,
                       std::vector<double> weights)
    : degree_(degree),
      knots_(std::move(knots)),
      poles_(std::move(poles)),
      weights_(std::move(weights)) {
    assert(!bspline::validate(degree_, knots_, poles_.size()));
    assert(weights_.size() == poles_.size());
}

Result<NurbsCurve> NurbsCurve::create(int degree, std::vector<double> knots,
                                      std::vector<Vec3> poles, std::vector<double> weights) {
    if (auto err = bspline::validate(degree, knots, poles.size())) return invalidInput(*err);
    if (weights.size() != poles.size()) return invalidInput("need one weight per pole");
    for (double w : weights)
        if (!(w > 0)) return invalidInput("NURBS weights must be positive");
    return NurbsCurve(degree, std::move(knots), std::move(poles), std::move(weights));
}

NurbsCurve NurbsCurve::fullCircle(const Vec3& center, const Vec3& normal, const Vec3& xAxis,
                                  double radius) {
    const Vec3 n = *normal.normalized();
    const Vec3 x = *(xAxis - n * xAxis.dot(n)).normalized();
    const Vec3 y = n.cross(x);
    // Square of side 2r around the circle; corners carry weight sqrt(2)/2.
    const double w = std::sqrt(2.0) / 2;
    const Vec3 X = x * radius, Y = y * radius;
    std::vector<Vec3> poles = {center + X,     center + X + Y, center + Y,
                               center - X + Y, center - X,     center - X - Y,
                               center - Y,     center + X - Y, center + X};
    std::vector<double> weights = {1, w, 1, w, 1, w, 1, w, 1};
    std::vector<double> knots = {0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1};
    return NurbsCurve(2, std::move(knots), std::move(poles), std::move(weights));
}

void NurbsCurve::evaluate(double t, int maxOrder, Vec3 out[3]) const {
    t = domain().clamp(t);
    const int span = bspline::findSpan(degree_, knots_, poles_.size(), t);
    const auto N = bspline::basisDerivatives(span, t, degree_, maxOrder, knots_);
    Vec3 A[3];
    double w[3] = {0, 0, 0};
    for (int k = 0; k <= maxOrder; ++k)
        for (int j = 0; j <= degree_; ++j) {
            const int i = span - degree_ + j;
            A[k] += poles_[i] * (N(k, j) * weights_[i]);
            w[k] += N(k, j) * weights_[i];
        }
    out[0] = A[0] / w[0];
    if (maxOrder >= 1) out[1] = (A[1] - out[0] * w[1]) / w[0];
    if (maxOrder >= 2) out[2] = (A[2] - out[1] * (2 * w[1]) - out[0] * w[2]) / w[0];
}

Vec3 NurbsCurve::point(double t) const {
    Vec3 out[3];
    evaluate(t, 0, out);
    return out[0];
}

Vec3 NurbsCurve::derivative(double t, int order) const {
    assert(order == 1 || order == 2);
    Vec3 out[3];
    evaluate(t, order, out);
    return out[order];
}

Interval NurbsCurve::domain() const { return bspline::domain(degree_, knots_); }

std::vector<double> NurbsCurve::breakpoints() const {
    return bspline::breakpoints(degree_, knots_);
}

void NurbsCurve::accept(CurveVisitor& visitor) const { visitor.visit(*this); }

}  // namespace geom
