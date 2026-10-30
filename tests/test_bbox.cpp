#include <gtest/gtest.h>

#include "geom/BoundingBox.h"

using namespace geom;

TEST(BoundingBox, DefaultIsEmpty) {
    BoundingBox b;
    EXPECT_TRUE(b.isEmpty());
    EXPECT_FALSE(b.contains({0, 0, 0}));
}

TEST(BoundingBox, ExpandFromEmpty) {
    BoundingBox b;
    b.expand({1, 2, 3});
    EXPECT_FALSE(b.isEmpty());
    b.expand({-1, 5, 0});
    EXPECT_TRUE(nearlyEqual(b.min(), Vec3(-1, 2, 0)));
    EXPECT_TRUE(nearlyEqual(b.max(), Vec3(1, 5, 3)));
    EXPECT_TRUE(nearlyEqual(b.center(), Vec3(0, 3.5, 1.5)));
    EXPECT_TRUE(nearlyEqual(b.size(), Vec3(2, 3, 3)));
}

TEST(BoundingBox, ConstructorOrdersCorners) {
    const BoundingBox b({5, 5, 5}, {0, 0, 0});
    EXPECT_TRUE(nearlyEqual(b.min(), Vec3(0, 0, 0)));
}

TEST(BoundingBox, Contains) {
    const BoundingBox b({0, 0, 0}, {1, 1, 1});
    EXPECT_TRUE(b.contains({0.5, 0.5, 0.5}));
    EXPECT_TRUE(b.contains({1, 1, 1}));  // boundary
    EXPECT_FALSE(b.contains({1.1, 0.5, 0.5}));
}

TEST(BoundingBox, Intersects) {
    const BoundingBox a({0, 0, 0}, {1, 1, 1});
    EXPECT_TRUE(a.intersects(BoundingBox({0.5, 0.5, 0.5}, {2, 2, 2})));
    EXPECT_TRUE(a.intersects(BoundingBox({1, 0, 0}, {2, 1, 1})));  // touching face
    EXPECT_FALSE(a.intersects(BoundingBox({1.5, 0, 0}, {2, 1, 1})));
    EXPECT_FALSE(a.intersects(BoundingBox()));
}

TEST(BoundingBox, MergeAndInflate) {
    BoundingBox a({0, 0, 0}, {1, 1, 1});
    a.expand(BoundingBox({2, 2, 2}, {3, 3, 3}));
    EXPECT_TRUE(nearlyEqual(a.max(), Vec3(3, 3, 3)));
    a.expand(BoundingBox());  // merging empty is a no-op
    EXPECT_TRUE(nearlyEqual(a.min(), Vec3(0, 0, 0)));
    EXPECT_TRUE(nearlyEqual(a.inflated(1).min(), Vec3(-1, -1, -1)));
}
