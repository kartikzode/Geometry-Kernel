#pragma once
#include "geom/Curve.h"
#include "geom/Result.h"

namespace geom {

// Circle or arc: P(t) = c + r (cos t X + sin t Y), t in [start, end] radians,
// with Y = N x X, so t increases counter-clockwise seen from +N.
class Circle final : public Curve {
public:
    // Preconditions: normal non-zero, xAxis not parallel to normal, radius > 0,
    // 0 < end - start <= 2 pi. xAxis is projected into the plane.
    Circle(const Vec3& center, const Vec3& normal, const Vec3& xAxis, double radius,
           double start = 0, double end = 2 * kPi);
    static Result<Circle> create(const Vec3& center, const Vec3& normal, const Vec3& xAxis,
                                 double radius, double start = 0, double end = 2 * kPi);
    // Circle in the XY plane, X axis along +x.
    static Circle xy(const Vec3& center, double radius, double start = 0, double end = 2 * kPi) {
        return Circle(center, {0, 0, 1}, {1, 0, 0}, radius, start, end);
    }

    Vec3 point(double t) const override;
    Vec3 derivative(double t, int order = 1) const override;
    Interval domain() const override { return dom_; }
    std::unique_ptr<Curve> clone() const override { return std::make_unique<Circle>(*this); }
    void accept(CurveVisitor& visitor) const override;
    std::string_view typeName() const override { return "Circle"; }
    double length() const override { return r_ * dom_.length(); }

    const Vec3& center() const { return c_; }
    const Vec3& normal() const { return n_; }
    const Vec3& xAxis() const { return x_; }
    const Vec3& yAxis() const { return y_; }
    double radius() const { return r_; }
    bool isFull() const { return dom_.length() >= 2 * kPi - kAngular; }

private:
    Vec3 c_, n_, x_, y_;
    double r_;
    Interval dom_;
};

}  // namespace geom
