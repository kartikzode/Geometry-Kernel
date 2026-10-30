#include <gtest/gtest.h>

#include "geom/Quaternion.h"

using namespace geom;

TEST(Quaternion, RotatesLikeAxisAngleMatrix) {
    const Vec3 axis(1, -2, 0.5), p(3, 1, -4);
    const double ang = 1.1;
    const auto q = Quaternion::fromAxisAngle(axis, ang);
    const auto m = Transform::rotation(axis, ang);
    ASSERT_TRUE(q && m);
    EXPECT_TRUE(nearlyEqual(q->rotate(p), m->transformPoint(p), 1e-12));
    EXPECT_TRUE(nearlyEqual(q->toMatrix(), *m, 1e-12));
}

TEST(Quaternion, ZeroAxisRejected) {
    EXPECT_FALSE(Quaternion::fromAxisAngle({0, 0, 0}, 1).has_value());
}

TEST(Quaternion, CompositionMatchesMatrixProduct) {
    const auto a = *Quaternion::fromAxisAngle({1, 0, 0}, 0.5);
    const auto b = *Quaternion::fromAxisAngle({0, 1, 0}, 0.8);
    const Matrix4 ma = *Transform::rotation({1, 0, 0}, 0.5);
    const Matrix4 mb = *Transform::rotation({0, 1, 0}, 0.8);
    EXPECT_TRUE(nearlyEqual((a * b).toMatrix(), ma * mb, 1e-12));
}

TEST(Quaternion, ConjugateInverts) {
    const auto q = *Quaternion::fromAxisAngle({1, 1, 1}, 2.0);
    const Vec3 p(1, 2, 3);
    EXPECT_TRUE(nearlyEqual(q.conjugate().rotate(q.rotate(p)), p, 1e-12));
}
