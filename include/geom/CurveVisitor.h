#pragma once

namespace geom {

class Line;
class Circle;
class BezierCurve;
class BSplineCurve;
class NurbsCurve;
class OffsetCurve;

// Visitor over the closed set of curve types. Adding an operation means
// writing a new visitor, not touching every curve class. The price is that
// adding a new curve type means touching every visitor.
class CurveVisitor {
public:
    virtual ~CurveVisitor() = default;
    virtual void visit(const Line&) = 0;
    virtual void visit(const Circle&) = 0;
    virtual void visit(const BezierCurve&) = 0;
    virtual void visit(const BSplineCurve&) = 0;
    virtual void visit(const NurbsCurve&) = 0;
    virtual void visit(const OffsetCurve&) = 0;
};

}  // namespace geom
