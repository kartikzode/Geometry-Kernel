#include <gtest/gtest.h>

#include "geom/Matrix4.h"

using namespace geom;

TEST(Matrix4, IdentityLeavesPointsAlone) {
    const Vec3 p(1, 2, 3);
    EXPECT_TRUE(nearlyEqual(Matrix4::identity().transformPoint(p), p));
}

TEST(Matrix4, TranslationMovesPointsNotDirections) {
    const auto t = Transform::translation({10, 20, 30});
    EXPECT_TRUE(nearlyEqual(t.transformPoint({1, 1, 1}), Vec3(11, 21, 31)));
    EXPECT_TRUE(nearlyEqual(t.transformDirection({1, 1, 1}), Vec3(1, 1, 1)));
}

TEST(Matrix4, Scale) {
    const auto s = Transform::scale(2, 3, 4);
    EXPECT_TRUE(nearlyEqual(s.transformPoint({1, 1, 1}), Vec3(2, 3, 4)));
    EXPECT_TRUE(nearlyEqual(Transform::scale(2).transformPoint({1, 2, 3}), Vec3(2, 4, 6)));
}

TEST(Matrix4, RotationAboutZ90) {
    const auto r = Transform::rotation({0, 0, 1}, kPi / 2);
    ASSERT_TRUE(r.has_value());
    EXPECT_TRUE(nearlyEqual(r->transformPoint({1, 0, 0}), Vec3(0, 1, 0)));
    EXPECT_TRUE(nearlyEqual(r->transformPoint({0, 1, 0}), Vec3(-1, 0, 0)));
}

TEST(Matrix4, RotationAxisNeedNotBeUnit) {
    const auto r = Transform::rotation({0, 0, 5}, kPi / 2);
    ASSERT_TRUE(r.has_value());
    EXPECT_TRUE(nearlyEqual(r->transformPoint({1, 0, 0}), Vec3(0, 1, 0)));
}

TEST(Matrix4, RotationAboutZeroAxisIsRejected) {
    EXPECT_FALSE(Transform::rotation({0, 0, 0}, 1.0).has_value());
}

TEST(Matrix4, RotationPreservesLengthAndAxis) {
    const Vec3 axis(1, 2, 3);
    const auto r = Transform::rotation(axis, 0.7);
    ASSERT_TRUE(r.has_value());
    const Vec3 p(4, -1, 2);
    EXPECT_NEAR(r->transformPoint(p).norm(), p.norm(), 1e-12);
    EXPECT_TRUE(nearlyEqual(r->transformPoint(axis), axis));  // axis is fixed
}

TEST(Matrix4, CompositionOrder) {
    // A*B applies B first. Rotate then translate != translate then rotate.
    const auto rot = *Transform::rotation({0, 0, 1}, kPi / 2);
    const auto tr = Transform::translation({1, 0, 0});
    EXPECT_TRUE(nearlyEqual((tr * rot).transformPoint({1, 0, 0}), Vec3(1, 1, 0)));
    EXPECT_TRUE(nearlyEqual((rot * tr).transformPoint({1, 0, 0}), Vec3(0, 2, 0)));
}

TEST(Matrix4, InverseTimesSelfIsIdentity) {
    const auto m = Transform::translation({3, -2, 5}) * *Transform::rotation({1, 1, 0}, 0.9) *
                   Transform::scale(2, 3, 0.5);
    const auto inv = m.inverse();
    ASSERT_TRUE(inv.has_value());
    EXPECT_TRUE(nearlyEqual(m * *inv, Matrix4::identity(), 1e-12));
    EXPECT_TRUE(nearlyEqual(*inv * m, Matrix4::identity(), 1e-12));
}

TEST(Matrix4, InverseUndoesTransform) {
    const auto m = Transform::translation({1, 2, 3}) * *Transform::rotation({0, 1, 0}, 1.2);
    const Vec3 p(7, 8, 9);
    EXPECT_TRUE(nearlyEqual(m.inverse()->transformPoint(m.transformPoint(p)), p, 1e-12));
}

TEST(Matrix4, SingularMatrixHasNoInverse) {
    EXPECT_FALSE(Transform::scale(1, 0, 1).inverse().has_value());
}

TEST(Matrix4, Transpose) {
    const auto r = *Transform::rotation({0, 0, 1}, 0.4);
    // For a pure rotation, transpose == inverse.
    EXPECT_TRUE(nearlyEqual(r.transposed(), *r.inverse(), 1e-12));
}
