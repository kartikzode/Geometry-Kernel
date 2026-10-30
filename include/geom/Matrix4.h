#pragma once
#include <array>
#include <optional>

#include "geom/Vec3.h"

namespace geom {

// 4x4 row-major matrix for affine transforms; column-vector convention
// (p' = M * p), so composition A*B applies B first.
class Matrix4 {
public:
    Matrix4();  // identity
    explicit Matrix4(const std::array<double, 16>& rowMajor) : m_(rowMajor) {}

    static Matrix4 identity() { return Matrix4(); }

    double operator()(int r, int c) const { return m_[static_cast<size_t>(r * 4 + c)]; }
    double& operator()(int r, int c) { return m_[static_cast<size_t>(r * 4 + c)]; }

    Matrix4 operator*(const Matrix4& o) const;
    Matrix4 transposed() const;
    // Gauss-Jordan with partial pivoting; nullopt if singular.
    std::optional<Matrix4> inverse() const;

    Vec3 transformPoint(const Vec3& p) const;      // w = 1, translation applies
    Vec3 transformDirection(const Vec3& d) const;  // w = 0, translation ignored

private:
    std::array<double, 16> m_;
};

bool nearlyEqual(const Matrix4& a, const Matrix4& b, double tol = kLinear);

// Factories for common transforms.
namespace Transform {
Matrix4 translation(const Vec3& t);
Matrix4 scale(double sx, double sy, double sz);
inline Matrix4 scale(double s) { return scale(s, s, s); }
// Rotation by `angle` radians about `axis` (need not be unit). nullopt for a zero axis.
std::optional<Matrix4> rotation(const Vec3& axis, double angle);
}  // namespace Transform

}  // namespace geom
