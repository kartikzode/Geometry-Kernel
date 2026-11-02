#pragma once
#include <memory>
#include <optional>
#include <string_view>
#include <vector>

#include "geom/Vec3.h"

namespace geom {

class CurveVisitor;

struct Interval {
    double lo = 0, hi = 1;
    double length() const { return hi - lo; }
    double clamp(double t) const { return t < lo ? lo : (t > hi ? hi : t); }
    bool contains(double t, double tol = 0) const { return t >= lo - tol && t <= hi + tol; }
};

// Abstract parametric curve C(t), t in domain().
//
// Concrete curves are immutable after construction, so every const member
// function is safe to call from several threads at once (see Parallel.h).
// Polymorphic copies go through clone(); the base copy operations are
// protected so a Curve can't be sliced by accident.
class Curve {
public:
    virtual ~Curve() = default;

    virtual Vec3 point(double t) const = 0;
    // order must be 1 or 2.
    virtual Vec3 derivative(double t, int order = 1) const = 0;
    virtual Interval domain() const = 0;
    virtual std::unique_ptr<Curve> clone() const = 0;
    virtual void accept(CurveVisitor& visitor) const = 0;
    virtual std::string_view typeName() const = 0;

    // Parameters where the curve may lose smoothness (distinct knots),
    // always including both domain ends. Integration and sampling work
    // span by span between breakpoints.
    virtual std::vector<double> breakpoints() const;

    // Total arc length. Overridden with exact formulas where they exist.
    virtual double length() const;
    // Arc length over [t0, t1] by adaptive Gauss-Legendre quadrature.
    double arcLength(double t0, double t1) const;

    Vec3 startPoint() const { return point(domain().lo); }
    Vec3 endPoint() const { return point(domain().hi); }

    // Unit tangent; nullopt where the first derivative vanishes.
    std::optional<Vec3> tangent(double t) const;
    // |C' x C''| / |C'|^3; 0 where the first derivative vanishes.
    double curvature(double t) const;

protected:
    Curve() = default;
    Curve(const Curve&) = default;
    Curve(Curve&&) = default;
    Curve& operator=(const Curve&) = default;
    Curve& operator=(Curve&&) = default;
};

}  // namespace geom
