# PES-01: Core math layer

## Scope
`Vec3`, `Matrix4`/`Transform`, `Quaternion`, `BoundingBox`, and the tolerance policy.
Everything later (curves, intersections, unfolding) builds on this.

## Tolerance strategy
- `kLinear = 1e-6` (model units, mm): points closer than this are coincident.
- `kAngular = 1e-9` rad: directions closer than this are parallel.
- `nearlyEqual(a, b)` accepts if the absolute difference <= `kLinear` **or** the relative
  difference <= `kRelative`. Absolute-only breaks for large coordinates; relative-only breaks
  near zero.
- `Vec3::operator==` is `= delete`d, so exact comparison of geometry is a compile error.
- Why `==` fails: `0.1 + 0.2 != 0.3` in binary floating point. Every computed coordinate
  carries rounding error, so equality must mean "within tolerance". Tolerance is also not
  transitive (a~b and b~c does not imply a~c), so avoid using it as a map key.
- Angles use `atan2(|a x b|, a.b)` rather than `acos(a.b/(|a||b|))`; `acos` loses precision
  near 0 and pi (test: `AngleBetweenAccurateForTinyAngles`).

## Error handling strategy
Operations that can fail on valid-looking input return `std::optional` (normalizing a zero
vector, a rotation about a zero axis, inverting a singular matrix). No exceptions in the math
layer: failure is an expected outcome, the caller must handle it, and it keeps hot paths cheap.
Programmer errors (e.g. out-of-range matrix index) are not checked at runtime in this layer.
Operations that need to say *why* they failed (factories, sheet-metal unfolding) return
`Result<T>` (`Result.h`): a value or an `Error{code, message}`.

## Ownership
All types here are small value types (copyable, no heap, no pointers). There is nothing to own
or free. Ownership questions start in Phase 2 with `std::unique_ptr<Curve>`.

## Conventions
- Column vectors, row-major storage: `p' = M * p`; `A * B` applies `B` first.
- `transformPoint` uses w = 1 (translation applies), `transformDirection` uses w = 0.
- Rotations are right-handed, angles in radians.
- Matrix inverse is general Gauss-Jordan with partial pivoting; singular = pivot < 1e-12.
  A faster rigid-body inverse can be added if profiling asks for it.

## Alternatives considered
- Eigen: excellent, but the point of this project is to show the implementation.
- Separate `Point3`/`Vector3` types: safer, more verbose. Deferred; `Vec3` plus the
  point/direction split on `Transform` is enough at this size.
