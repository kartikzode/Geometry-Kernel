#pragma once
#include <cmath>
#include <optional>

#include "geom/Tolerance.h"

namespace geom {

// Value type: 3 doubles, trivially copyable. Used for both points and
// directions; Transform distinguishes them (w = 1 vs w = 0).
struct Vec3 {
    double x = 0, y = 0, z = 0;

    constexpr Vec3() = default;
    constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    constexpr Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    constexpr Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    constexpr Vec3 operator-() const { return {-x, -y, -z}; }
    constexpr Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    constexpr Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }

    // Exact comparison is deliberately not provided; use nearlyEqual().
    bool operator==(const Vec3&) const = delete;

    constexpr double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
    constexpr Vec3 cross(const Vec3& o) const {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }
    constexpr double normSquared() const { return dot(*this); }
    double norm() const { return std::sqrt(normSquared()); }

    // Returns nullopt for (near-)zero vectors instead of dividing by ~0.
    std::optional<Vec3> normalized(double tol = kLinear) const {
        const double n = norm();
        if (n <= tol) return std::nullopt;
        return *this / n;
    }

    double distanceTo(const Vec3& o) const { return (*this - o).norm(); }
};

constexpr Vec3 operator*(double s, const Vec3& v) { return v * s; }

inline bool nearlyEqual(const Vec3& a, const Vec3& b, double tol = kLinear) {
    return a.distanceTo(b) <= tol;
}

// Angle between two vectors in [0, pi]. Uses atan2(|a x b|, a.b), which stays
// accurate near 0 and pi where acos(dot) loses precision. Zero vectors -> nullopt.
inline std::optional<double> angleBetween(const Vec3& a, const Vec3& b) {
    if (a.norm() <= kLinear || b.norm() <= kLinear) return std::nullopt;
    return std::atan2(a.cross(b).norm(), a.dot(b));
}

// Parallel or anti-parallel within kAngular. False if either vector is zero.
inline bool isParallel(const Vec3& a, const Vec3& b, double angTol = kAngular) {
    const auto ang = angleBetween(a, b);
    if (!ang) return false;
    return *ang <= angTol || (kPi - *ang) <= angTol;
}

inline bool isPerpendicular(const Vec3& a, const Vec3& b, double angTol = kAngular) {
    const auto ang = angleBetween(a, b);
    return ang && std::abs(*ang - kPi / 2) <= angTol;
}

}  // namespace geom
