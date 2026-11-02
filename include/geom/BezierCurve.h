#pragma once
#include <utility>
#include <vector>

#include "geom/Curve.h"
#include "geom/Result.h"

namespace geom {

// Bezier curve of degree n = poles.size() - 1 on t in [0, 1], evaluated with
// de Casteljau's algorithm (stable: only convex combinations of points).
class BezierCurve final : public Curve {
public:
    // Precondition: at least one pole. A single pole is the degree-0 curve
    // that the hodograph of a line produces. create() requires two.
    explicit BezierCurve(std::vector<Vec3> poles);
    static Result<BezierCurve> create(std::vector<Vec3> poles);

    Vec3 point(double t) const override;
    Vec3 derivative(double t, int order = 1) const override;
    Interval domain() const override { return {0, 1}; }
    std::unique_ptr<Curve> clone() const override { return std::make_unique<BezierCurve>(*this); }
    void accept(CurveVisitor& visitor) const override;
    std::string_view typeName() const override { return "Bezier"; }

    int degree() const { return static_cast<int>(poles_.size()) - 1; }
    const std::vector<Vec3>& poles() const { return poles_; }

    // Hodograph: the derivative is a Bezier curve of degree n - 1 with poles
    // n (P[i+1] - P[i]).
    BezierCurve derivativeCurve() const;
    // Split at t into two curves covering [0, t] and [t, 1].
    std::pair<BezierCurve, BezierCurve> subdivide(double t) const;

private:
    static Vec3 deCasteljau(std::vector<Vec3> pts, double t);
    std::vector<Vec3> poles_;
};

}  // namespace geom
