#pragma once
#include <cmath>
#include <optional>

#include "geom/Matrix4.h"
#include "geom/Vec3.h"

namespace geom {

// Unit quaternion (w, x, y, z) representing a rotation. Factories normalise.
struct Quaternion {
    double w = 1, x = 0, y = 0, z = 0;

    static std::optional<Quaternion> fromAxisAngle(const Vec3& axis, double angle) {
        const auto u = axis.normalized();
        if (!u) return std::nullopt;
        const double h = angle / 2, s = std::sin(h);
        return Quaternion{std::cos(h), u->x * s, u->y * s, u->z * s};
    }

    Quaternion conjugate() const { return {w, -x, -y, -z}; }

    // Hamilton product: (a*b) applies b first, then a.
    Quaternion operator*(const Quaternion& o) const {
        return {w * o.w - x * o.x - y * o.y - z * o.z, w * o.x + x * o.w + y * o.z - z * o.y,
                w * o.y - x * o.z + y * o.w + z * o.x, w * o.z + x * o.y - y * o.x + z * o.w};
    }

    Vec3 rotate(const Vec3& v) const {
        const Quaternion p{0, v.x, v.y, v.z};
        const Quaternion r = (*this) * p * conjugate();
        return {r.x, r.y, r.z};
    }

    Matrix4 toMatrix() const {
        Matrix4 m;
        m(0, 0) = 1 - 2 * (y * y + z * z);
        m(0, 1) = 2 * (x * y - z * w);
        m(0, 2) = 2 * (x * z + y * w);
        m(1, 0) = 2 * (x * y + z * w);
        m(1, 1) = 1 - 2 * (x * x + z * z);
        m(1, 2) = 2 * (y * z - x * w);
        m(2, 0) = 2 * (x * z - y * w);
        m(2, 1) = 2 * (y * z + x * w);
        m(2, 2) = 1 - 2 * (x * x + y * y);
        return m;
    }
};

}  // namespace geom
