#pragma once
#include <algorithm>
#include <cmath>

namespace geom {

// Global tolerance policy. Lengths are in model units (mm by convention).
// Two points closer than kLinear are "the same point"; two directions whose
// angle is below kAngular are "parallel". See docs/PES-01-core-math.md.
inline constexpr double kLinear = 1e-6;
inline constexpr double kAngular = 1e-9;  // radians
inline constexpr double kRelative = 1e-12;
inline constexpr double kPi = 3.14159265358979323846;  // M_PI is not portable to MSVC

// Mixed absolute/relative comparison. Absolute alone fails for large values,
// relative alone fails near zero; accept if either test passes.
inline bool nearlyEqual(double a, double b, double absTol = kLinear, double relTol = kRelative) {
    const double diff = std::abs(a - b);
    if (diff <= absTol) return true;
    return diff <= relTol * std::max(std::abs(a), std::abs(b));
}

inline bool isZero(double a, double absTol = kLinear) { return std::abs(a) <= absTol; }

}  // namespace geom
