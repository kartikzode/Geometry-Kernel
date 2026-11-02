#pragma once
#include <optional>
#include <vector>

#include "geom/Curve.h"
#include "geom/Result.h"

namespace geom {

// Non-rational B-spline curve of any degree up to bspline::kMaxDegree.
// Domain is [knots[p], knots[m - p]] with m = knots.size() - 1.
class BSplineCurve final : public Curve {
public:
    // Precondition: bspline::validate(degree, knots, poles.size()) passes.
    BSplineCurve(int degree, std::vector<double> knots, std::vector<Vec3> poles);
    static Result<BSplineCurve> create(int degree, std::vector<double> knots,
                                       std::vector<Vec3> poles);

    Vec3 point(double t) const override;
    Vec3 derivative(double t, int order = 1) const override;
    Interval domain() const override;
    std::vector<double> breakpoints() const override;
    std::unique_ptr<Curve> clone() const override { return std::make_unique<BSplineCurve>(*this); }
    void accept(CurveVisitor& visitor) const override;
    std::string_view typeName() const override { return "BSpline"; }

    int degree() const { return degree_; }
    const std::vector<double>& knots() const { return knots_; }
    const std::vector<Vec3>& poles() const { return poles_; }

    // Boehm knot insertion: same curve, one more knot and pole.
    // nullopt if u is outside the domain or already has multiplicity = degree.
    std::optional<BSplineCurve> insertKnot(double u) const;

private:
    Vec3 evaluate(double t, int order) const;
    int degree_;
    std::vector<double> knots_;
    std::vector<Vec3> poles_;
};

}  // namespace geom
