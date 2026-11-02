#include <gtest/gtest.h>

#include "geom/Circle.h"
#include "geom/Line.h"

using namespace geom;

TEST(Line, EvaluatesAndDifferentiates) {
    const Line l({0, 0, 0}, {3, 4, 0});
    EXPECT_TRUE(nearlyEqual(l.point(0.5), Vec3(1.5, 2, 0)));
    EXPECT_TRUE(nearlyEqual(l.derivative(0.3), Vec3(3, 4, 0)));
    EXPECT_TRUE(nearlyEqual(l.derivative(0.3, 2), Vec3()));
    EXPECT_DOUBLE_EQ(l.length(), 5);
    EXPECT_DOUBLE_EQ(l.curvature(0.5), 0);
}

TEST(Line, NumericLengthMatchesExact) {
    const Line l({1, 2, 3}, {4, -2, 15});
    EXPECT_NEAR(l.arcLength(0, 1), l.length(), 1e-12);
    EXPECT_NEAR(l.arcLength(0.25, 0.75), 0.5 * l.length(), 1e-12);
}

TEST(Line, CreateRejectsDegenerate) {
    const auto r = Line::create({1, 1, 1}, {1, 1, 1 + 1e-9});
    ASSERT_FALSE(r.ok());
    EXPECT_EQ(r.error().code, ErrorCode::Degenerate);
}

TEST(Circle, PointsLieOnCircle) {
    const Circle c({1, 2, 3}, {1, 1, 0}, {0, 0, 1}, 2.5);
    for (int i = 0; i <= 20; ++i) {
        const Vec3 p = c.point(i * 2 * kPi / 20);
        EXPECT_NEAR(p.distanceTo(c.center()), 2.5, 1e-12);
        EXPECT_NEAR((p - c.center()).dot(c.normal()), 0, 1e-12);  // stays in plane
    }
}

TEST(Circle, CounterClockwiseAboutNormal) {
    const Circle c = Circle::xy({0, 0, 0}, 1);
    EXPECT_TRUE(nearlyEqual(c.point(0), Vec3(1, 0, 0)));
    EXPECT_TRUE(nearlyEqual(c.point(kPi / 2), Vec3(0, 1, 0)));
}

TEST(Circle, LengthAndCurvature) {
    const Circle full = Circle::xy({0, 0, 0}, 3);
    EXPECT_NEAR(full.length(), 6 * kPi, 1e-12);
    EXPECT_NEAR(full.arcLength(0, 2 * kPi), 6 * kPi, 1e-9);  // quadrature agrees
    EXPECT_NEAR(full.curvature(1.234), 1.0 / 3, 1e-12);
    const Circle arc = Circle::xy({0, 0, 0}, 2, 0, kPi / 2);
    EXPECT_NEAR(arc.length(), kPi, 1e-12);
    EXPECT_FALSE(arc.isFull());
}

TEST(Circle, TangentIsPerpendicularToRadius) {
    const Circle c({0, 0, 0}, {0, 1, 1}, {1, 0, 0}, 4);
    for (double t : {0.0, 0.7, 2.0, 5.5}) {
        const auto tan = c.tangent(t);
        ASSERT_TRUE(tan);
        EXPECT_NEAR(tan->dot(c.point(t) - c.center()), 0, 1e-12);
        EXPECT_NEAR(tan->norm(), 1, 1e-15);
    }
}

TEST(Circle, SecondDerivativePointsToCenter) {
    const Circle c = Circle::xy({5, 5, 0}, 2);
    const Vec3 d2 = c.derivative(0.3, 2);
    const Vec3 toCenter = c.center() - c.point(0.3);
    EXPECT_TRUE(isParallel(d2, toCenter));
    EXPECT_GT(d2.dot(toCenter), 0);
}

TEST(Circle, CreateValidates) {
    EXPECT_FALSE(Circle::create({0, 0, 0}, {0, 0, 0}, {1, 0, 0}, 1).ok());
    EXPECT_FALSE(Circle::create({0, 0, 0}, {0, 0, 1}, {0, 0, 2}, 1).ok());
    EXPECT_FALSE(Circle::create({0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 0).ok());
    EXPECT_FALSE(Circle::create({0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 1, 1, 1).ok());
    EXPECT_FALSE(Circle::create({0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 1, 0, 7).ok());
    EXPECT_TRUE(Circle::create({0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 1, 0, kPi).ok());
}

TEST(Circle, XAxisIsProjectedIntoPlane) {
    const Circle c({0, 0, 0}, {0, 0, 1}, {1, 0, 1}, 1);  // tilted x axis
    EXPECT_TRUE(nearlyEqual(c.xAxis(), Vec3(1, 0, 0)));
}

TEST(Curve, CloneIsIndependentDeepCopy) {
    const Circle c = Circle::xy({1, 0, 0}, 2);
    std::unique_ptr<Curve> copy = c.clone();
    EXPECT_EQ(copy->typeName(), "Circle");
    EXPECT_TRUE(nearlyEqual(copy->point(1.0), c.point(1.0)));
}
