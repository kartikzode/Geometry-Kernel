#pragma once
#include "geom/Curve.h"
#include "geom/Result.h"

namespace geom {

// Segment P(t) = p0 + t (p1 - p0), t in [0, 1].
class Line final : public Curve {
public:
    // Precondition: p0 and p1 are distinct. Use create() for unchecked input.
    Line(const Vec3& p0, const Vec3& p1);
    static Result<Line> create(const Vec3& p0, const Vec3& p1);

    Vec3 point(double t) const override { return p0_ + (p1_ - p0_) * t; }
    Vec3 derivative(double t, int order = 1) const override;
    Interval domain() const override { return {0, 1}; }
    std::unique_ptr<Curve> clone() const override { return std::make_unique<Line>(*this); }
    void accept(CurveVisitor& visitor) const override;
    std::string_view typeName() const override { return "Line"; }
    double length() const override { return p0_.distanceTo(p1_); }

    const Vec3& p0() const { return p0_; }
    const Vec3& p1() const { return p1_; }
    Vec3 direction() const { return *(p1_ - p0_).normalized(0); }

private:
    Vec3 p0_, p1_;
};

}  // namespace geom
