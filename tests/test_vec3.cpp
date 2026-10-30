#include <gtest/gtest.h>

#include "geom/Vec3.h"

using namespace geom;

TEST(Vec3, Arithmetic) {
    const Vec3 a(1, 2, 3), b(4, 5, 6);
    EXPECT_TRUE(nearlyEqual(a + b, Vec3(5, 7, 9)));
    EXPECT_TRUE(nearlyEqual(b - a, Vec3(3, 3, 3)));
    EXPECT_TRUE(nearlyEqual(a * 2, Vec3(2, 4, 6)));
    EXPECT_TRUE(nearlyEqual(2 * a, Vec3(2, 4, 6)));
    EXPECT_TRUE(nearlyEqual(-a, Vec3(-1, -2, -3)));
    EXPECT_TRUE(nearlyEqual(b / 2, Vec3(2, 2.5, 3)));
}

TEST(Vec3, DotAndCross) {
    const Vec3 x(1, 0, 0), y(0, 1, 0), z(0, 0, 1);
    EXPECT_DOUBLE_EQ(x.dot(y), 0);
    EXPECT_DOUBLE_EQ(Vec3(1, 2, 3).dot(Vec3(4, 5, 6)), 32);
    EXPECT_TRUE(nearlyEqual(x.cross(y), z));
    EXPECT_TRUE(nearlyEqual(y.cross(x), -z));
    EXPECT_TRUE(nearlyEqual(x.cross(x), Vec3()));
}

TEST(Vec3, CrossIsPerpendicularToBoth) {
    const Vec3 a(1.3, -2.1, 0.7), b(-0.4, 3.3, 2.2);
    const Vec3 c = a.cross(b);
    EXPECT_NEAR(c.dot(a), 0, 1e-12);
    EXPECT_NEAR(c.dot(b), 0, 1e-12);
}

TEST(Vec3, Norm) {
    EXPECT_DOUBLE_EQ(Vec3(3, 4, 0).norm(), 5);
    EXPECT_DOUBLE_EQ(Vec3(3, 4, 12).normSquared(), 169);
}

TEST(Vec3, NormalizedUnitLength) {
    const auto n = Vec3(3, 4, 12).normalized();
    ASSERT_TRUE(n.has_value());
    EXPECT_NEAR(n->norm(), 1.0, 1e-15);
}

TEST(Vec3, NormalizeZeroVectorFailsSafely) {
    EXPECT_FALSE(Vec3().normalized().has_value());
    EXPECT_FALSE(Vec3(1e-9, 0, 0).normalized().has_value());
}

TEST(Vec3, Distance) { EXPECT_DOUBLE_EQ(Vec3(1, 1, 1).distanceTo(Vec3(4, 5, 1)), 5); }

TEST(Vec3, AngleBetween) {
    EXPECT_NEAR(*angleBetween({1, 0, 0}, {0, 1, 0}), kPi / 2, 1e-15);
    EXPECT_NEAR(*angleBetween({1, 0, 0}, {-1, 0, 0}), kPi, 1e-15);
    EXPECT_NEAR(*angleBetween({1, 0, 0}, {1, 1, 0}), kPi / 4, 1e-15);
    EXPECT_FALSE(angleBetween({0, 0, 0}, {1, 0, 0}).has_value());
}

TEST(Vec3, AngleBetweenAccurateForTinyAngles) {
    // acos(dot) would return 0 or ~1e-8 noise here; atan2 form is exact.
    const double eps = 1e-10;
    EXPECT_NEAR(*angleBetween({1, 0, 0}, {1, eps, 0}), eps, 1e-20);
}

TEST(Vec3, Parallel) {
    EXPECT_TRUE(isParallel({1, 2, 3}, {2, 4, 6}));
    EXPECT_TRUE(isParallel({1, 2, 3}, {-2, -4, -6}));  // anti-parallel counts
    EXPECT_FALSE(isParallel({1, 0, 0}, {0, 1, 0}));
    EXPECT_FALSE(isParallel({0, 0, 0}, {1, 0, 0}));
}

TEST(Vec3, Perpendicular) {
    EXPECT_TRUE(isPerpendicular({1, 0, 0}, {0, 3, 0}));
    EXPECT_FALSE(isPerpendicular({1, 0, 0}, {1, 1, 0}));
}

TEST(Vec3, CompoundAssign) {
    Vec3 a(1, 1, 1);
    a += Vec3(1, 2, 3);
    a -= Vec3(0, 1, 0);
    a *= 2;
    EXPECT_TRUE(nearlyEqual(a, Vec3(4, 4, 8)));
}
