#pragma once
#include <algorithm>

#include "geom/Vec3.h"

namespace geom {

// Axis-aligned box. Default-constructed is empty (min > max) so that
// expand() on the first point works without a special case.
class BoundingBox {
public:
    BoundingBox() = default;
    BoundingBox(const Vec3& a, const Vec3& b)
        : min_(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)),
          max_(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)) {}

    bool isEmpty() const { return min_.x > max_.x; }
    const Vec3& min() const { return min_; }
    const Vec3& max() const { return max_; }
    Vec3 center() const { return (min_ + max_) * 0.5; }
    Vec3 size() const { return max_ - min_; }

    void expand(const Vec3& p) {
        min_ = {std::min(min_.x, p.x), std::min(min_.y, p.y), std::min(min_.z, p.z)};
        max_ = {std::max(max_.x, p.x), std::max(max_.y, p.y), std::max(max_.z, p.z)};
    }
    void expand(const BoundingBox& b) {
        if (b.isEmpty()) return;
        expand(b.min_);
        expand(b.max_);
    }

    bool contains(const Vec3& p, double tol = kLinear) const {
        return !isEmpty() && p.x >= min_.x - tol && p.x <= max_.x + tol && p.y >= min_.y - tol &&
               p.y <= max_.y + tol && p.z >= min_.z - tol && p.z <= max_.z + tol;
    }

    // Touching boxes count as intersecting (within tol).
    bool intersects(const BoundingBox& o, double tol = kLinear) const {
        if (isEmpty() || o.isEmpty()) return false;
        return min_.x <= o.max_.x + tol && max_.x >= o.min_.x - tol && min_.y <= o.max_.y + tol &&
               max_.y >= o.min_.y - tol && min_.z <= o.max_.z + tol && max_.z >= o.min_.z - tol;
    }

    BoundingBox inflated(double d) const {
        if (isEmpty()) return *this;
        return BoundingBox(min_ - Vec3(d, d, d), max_ + Vec3(d, d, d));
    }

private:
    static constexpr double kInf = 1e300;
    Vec3 min_{kInf, kInf, kInf};
    Vec3 max_{-kInf, -kInf, -kInf};
};

}  // namespace geom
