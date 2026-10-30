#include "geom/Matrix4.h"

#include <cmath>
#include <utility>

namespace geom {

Matrix4::Matrix4() : m_{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1} {}

Matrix4 Matrix4::operator*(const Matrix4& o) const {
    Matrix4 r(std::array<double, 16>{});
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            double s = 0;
            for (int k = 0; k < 4; ++k) s += (*this)(i, k) * o(k, j);
            r(i, j) = s;
        }
    return r;
}

Matrix4 Matrix4::transposed() const {
    Matrix4 r;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) r(i, j) = (*this)(j, i);
    return r;
}

std::optional<Matrix4> Matrix4::inverse() const {
    Matrix4 a = *this;
    Matrix4 inv;
    for (int col = 0; col < 4; ++col) {
        int piv = col;
        for (int r = col + 1; r < 4; ++r)
            if (std::abs(a(r, col)) > std::abs(a(piv, col))) piv = r;
        if (std::abs(a(piv, col)) < 1e-12) return std::nullopt;
        if (piv != col)
            for (int c = 0; c < 4; ++c) {
                std::swap(a(piv, c), a(col, c));
                std::swap(inv(piv, c), inv(col, c));
            }
        const double d = a(col, col);
        for (int c = 0; c < 4; ++c) {
            a(col, c) /= d;
            inv(col, c) /= d;
        }
        for (int r = 0; r < 4; ++r) {
            if (r == col) continue;
            const double f = a(r, col);
            if (f == 0) continue;
            for (int c = 0; c < 4; ++c) {
                a(r, c) -= f * a(col, c);
                inv(r, c) -= f * inv(col, c);
            }
        }
    }
    return inv;
}

Vec3 Matrix4::transformPoint(const Vec3& p) const {
    const auto& m = *this;
    return {m(0, 0) * p.x + m(0, 1) * p.y + m(0, 2) * p.z + m(0, 3),
            m(1, 0) * p.x + m(1, 1) * p.y + m(1, 2) * p.z + m(1, 3),
            m(2, 0) * p.x + m(2, 1) * p.y + m(2, 2) * p.z + m(2, 3)};
}

Vec3 Matrix4::transformDirection(const Vec3& d) const {
    const auto& m = *this;
    return {m(0, 0) * d.x + m(0, 1) * d.y + m(0, 2) * d.z,
            m(1, 0) * d.x + m(1, 1) * d.y + m(1, 2) * d.z,
            m(2, 0) * d.x + m(2, 1) * d.y + m(2, 2) * d.z};
}

bool nearlyEqual(const Matrix4& a, const Matrix4& b, double tol) {
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c)
            if (!nearlyEqual(a(r, c), b(r, c), tol)) return false;
    return true;
}

namespace Transform {

Matrix4 translation(const Vec3& t) {
    Matrix4 m;
    m(0, 3) = t.x;
    m(1, 3) = t.y;
    m(2, 3) = t.z;
    return m;
}

Matrix4 scale(double sx, double sy, double sz) {
    Matrix4 m;
    m(0, 0) = sx;
    m(1, 1) = sy;
    m(2, 2) = sz;
    return m;
}

// Rodrigues' rotation formula.
std::optional<Matrix4> rotation(const Vec3& axis, double angle) {
    const auto u = axis.normalized();
    if (!u) return std::nullopt;
    const double c = std::cos(angle), s = std::sin(angle), t = 1 - c;
    Matrix4 m;
    m(0, 0) = t * u->x * u->x + c;
    m(0, 1) = t * u->x * u->y - s * u->z;
    m(0, 2) = t * u->x * u->z + s * u->y;
    m(1, 0) = t * u->x * u->y + s * u->z;
    m(1, 1) = t * u->y * u->y + c;
    m(1, 2) = t * u->y * u->z - s * u->x;
    m(2, 0) = t * u->x * u->z - s * u->y;
    m(2, 1) = t * u->y * u->z + s * u->x;
    m(2, 2) = t * u->z * u->z + c;
    return m;
}

}  // namespace Transform
}  // namespace geom
