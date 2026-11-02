#include "geom/Curve.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>

namespace geom {

namespace {

// 5-point Gauss-Legendre on [a, b]; exact for polynomials up to degree 9.
constexpr std::array<double, 5> kNodes = {0.0, -0.5384693101056831, 0.5384693101056831,
                                          -0.9061798459386640, 0.9061798459386640};
constexpr std::array<double, 5> kWeights = {0.5688888888888889, 0.4786286704993665,
                                            0.4786286704993665, 0.2369268850561891,
                                            0.2369268850561891};

double gauss5(const std::function<double(double)>& f, double a, double b) {
    const double mid = 0.5 * (a + b), half = 0.5 * (b - a);
    double sum = 0;
    for (size_t i = 0; i < kNodes.size(); ++i) sum += kWeights[i] * f(mid + half * kNodes[i]);
    return sum * half;
}

double adaptive(const std::function<double(double)>& f, double a, double b, double whole,
                double tol, int depth) {
    const double m = 0.5 * (a + b);
    const double left = gauss5(f, a, m), right = gauss5(f, m, b);
    if (depth <= 0 || std::abs(left + right - whole) <= tol) return left + right;
    return adaptive(f, a, m, left, 0.5 * tol, depth - 1) +
           adaptive(f, m, b, right, 0.5 * tol, depth - 1);
}

}  // namespace

std::vector<double> Curve::breakpoints() const {
    const Interval d = domain();
    return {d.lo, d.hi};
}

double Curve::length() const {
    const Interval d = domain();
    return arcLength(d.lo, d.hi);
}

double Curve::arcLength(double t0, double t1) const {
    if (t1 < t0) return -arcLength(t1, t0);
    const auto speed = [this](double t) { return derivative(t, 1).norm(); };
    // Integrate span by span so kinks at knots don't slow convergence.
    double total = 0;
    double a = t0;
    for (double bp : breakpoints()) {
        if (bp <= a) continue;
        const double b = std::min(bp, t1);
        total += adaptive(speed, a, b, gauss5(speed, a, b), 1e-11, 24);
        a = b;
        if (a >= t1) break;
    }
    return total;
}

std::optional<Vec3> Curve::tangent(double t) const {
    return derivative(t, 1).normalized(1e-14);
}

double Curve::curvature(double t) const {
    const Vec3 d1 = derivative(t, 1), d2 = derivative(t, 2);
    const double speed = d1.norm();
    if (speed < 1e-14) return 0;
    return d1.cross(d2).norm() / (speed * speed * speed);
}

}  // namespace geom
