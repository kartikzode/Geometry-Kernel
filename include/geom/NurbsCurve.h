#pragma once
#include <vector>

#include "geom/Curve.h"
#include "geom/Result.h"

namespace geom {

// Rational B-spline: C(t) = sum(N_i w_i P_i) / sum(N_i w_i).
// Weights must be positive, which keeps the curve inside the convex hull of
// its poles. Represents conics (circles, ellipses) exactly.
class NurbsCurve final : public Curve {
public:
    // Precondition: knots valid for the degree and pole count, weights.size()
    // == poles.size(), all weights > 0.
    NurbsCurve(int degree, std::vector<double> knots, std::vector<Vec3> poles,
               std::vector<double> weights);
    static Result<NurbsCurve> create(int degree, std::vector<double> knots,
                                     std::vector<Vec3> poles, std::vector<double> weights);
    // Exact full circle as a 9-pole quadratic NURBS.
    static NurbsCurve fullCircle(const Vec3& center, const Vec3& normal, const Vec3& xAxis,
                                 double radius);

    Vec3 point(double t) const override;
    Vec3 derivative(double t, int order = 1) const override;
    Interval domain() const override;
    std::vector<double> breakpoints() const override;
    std::unique_ptr<Curve> clone() const override { return std::make_unique<NurbsCurve>(*this); }
    void accept(CurveVisitor& visitor) const override;
    std::string_view typeName() const override { return "NURBS"; }

    int degree() const { return degree_; }
    const std::vector<double>& knots() const { return knots_; }
    const std::vector<Vec3>& poles() const { return poles_; }
    const std::vector<double>& weights() const { return weights_; }

private:
    // C, C', C'' (up to maxOrder) via the quotient rule on A(t) / w(t).
    void evaluate(double t, int maxOrder, Vec3 out[3]) const;
    int degree_;
    std::vector<double> knots_;
    std::vector<Vec3> poles_;
    std::vector<double> weights_;
};

}  // namespace geom
