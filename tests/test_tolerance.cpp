#include <gtest/gtest.h>

#include "geom/Tolerance.h"

using namespace geom;

TEST(Tolerance, DoublesDoNotCompareExactly) {
    // The canonical reason raw == fails on geometry.
    EXPECT_NE(0.1 + 0.2, 0.3);
    EXPECT_TRUE(nearlyEqual(0.1 + 0.2, 0.3));
}

TEST(Tolerance, AbsoluteNearZero) {
    EXPECT_TRUE(nearlyEqual(0.0, 1e-9));
    EXPECT_FALSE(nearlyEqual(0.0, 1e-3));
}

TEST(Tolerance, RelativeForLargeValues) {
    // 1e12 + 1e-3 differs by more than kLinear in absolute terms but is
    // within the relative tolerance.
    EXPECT_TRUE(nearlyEqual(1e12, 1e12 + 0.5));
    EXPECT_FALSE(nearlyEqual(1e12, 1e12 + 100.0));
}

TEST(Tolerance, Symmetric) {
    EXPECT_EQ(nearlyEqual(1.0, 1.0 + 5e-7), nearlyEqual(1.0 + 5e-7, 1.0));
}

TEST(Tolerance, IsZero) {
    EXPECT_TRUE(isZero(-1e-8));
    EXPECT_FALSE(isZero(1e-3));
}
