# PES-02: Curves

## Scope
An abstract `Curve` interface and five implementations: `Line`, `Circle` (arcs too),
`BezierCurve`, `BSplineCurve`, `NurbsCurve`. Derivatives, tangent, curvature and arc
length are covered, plus a `CurveFactory` and a `CurveVisitor`. `OffsetCurve` is
specified in PES-03.

## Interface
```cpp
class Curve {
    virtual Vec3 point(double t) const = 0;
    virtual Vec3 derivative(double t, int order = 1) const = 0;  // order 1 or 2
    virtual Interval domain() const = 0;
    virtual std::unique_ptr<Curve> clone() const = 0;
    virtual void accept(CurveVisitor&) const = 0;
    virtual std::vector<double> breakpoints() const;  // knots; default = domain ends
    virtual double length() const;                    // exact override where possible
    double arcLength(double t0, double t1) const;     // adaptive Gauss-Legendre
    std::optional<Vec3> tangent(double t) const;
    double curvature(double t) const;                 // |C' x C''| / |C'|^3
};
```

## Algorithms chosen
| Need | Choice | Why |
|---|---|---|
| Bezier evaluation | de Casteljau | Only convex combinations, so numerically stable. It also yields subdivision. |
| Bezier derivative | Hodograph (derivative poles `n (P[i+1] - P[i])`) | Exact, and reuses de Casteljau. |
| B-spline basis | Piegl & Tiller A2.1 (span search) and A2.3 (basis + derivatives) | Standard reference algorithms that compute all derivatives in one pass. |
| NURBS derivatives | Quotient rule on homogeneous sums `A(t)/w(t)` | Exact; avoids differentiating the rational form numerically. |
| Circle as NURBS | 9-pole quadratic, corner weights `sqrt(2)/2` | Exact conic, and a strong test oracle. |
| Knot insertion | Boehm's algorithm | Leaves the shape unchanged, which is easy to test. |
| Arc length | Adaptive 5-point Gauss-Legendre, span by span | Converges quickly on smooth spans. Splitting at knots avoids the derivative kinks there. |

## Ownership and copying
- Curves are immutable after construction. They have no caches and no mutable state, so
  every const method is thread-safe (see PES-05 / `Parallel.h`).
- Polymorphic curves live in `std::unique_ptr<Curve>`. `clone()` gives a deep copy.
- The base class copy and move operations are `protected`, so `Curve c = someCircle;`
  (slicing) doesn't compile. Concrete classes are `final` value types with the rule of zero.

## Construction and validation
Each concrete class has two ways in:
- A **public constructor** with documented preconditions, checked by `assert`. It's for
  trusted code such as tests and internal construction.
- A **`static create(...)` that returns `Result<T>`**, for untrusted input. It reports why
  it failed, for example "knot count must equal control points + degree + 1 (expected 10,
  got 9)".

`CurveFactory::create(CurveSpec)` takes a `std::variant` of plain spec structs, dispatches
with `std::visit`, and returns `Result<std::unique_ptr<Curve>>`. Callers then see only the
interface and an error message.

## Visitor vs virtual functions
Core evaluation (`point`, `derivative`) is virtual: every curve must support it. Operations
that are type-specific *and* external to the curve go in visitors:
- `BoundingBoxVisitor` uses each type's guarantee: exact arc extremes, and the convex-hull
  property of poles.
- The SVG `PathWriter` emits exact SVG arcs for circles.

The tradeoff: adding an operation means writing a new visitor without touching the
curves, but adding a curve type means updating every visitor.

## Alternatives considered
- **`std::variant<Line, Circle, ...>` instead of inheritance.** It avoids heap allocation
  and virtual calls, but closes the set of types at compile time and makes user-defined
  curves impossible. CAD kernels (Open CASCADE's `Geom_Curve`, CATIA's geometry) use
  class hierarchies, so the interface version is closer to what the job involves.
- **Templates over curve type.** Fastest, but every algorithm becomes a template and a
  heterogeneous container of curves needs type erasure anyway.

## Limits
- Derivatives above order 2 are not part of the interface.
- No degree elevation, knot removal or curve fitting yet.
