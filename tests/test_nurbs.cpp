#include <gtest/gtest.h>

#include <cmath>

#include "geom/BSplineCurve.h"
#include "geom/NurbsCurve.h"

using namespace geom;

TEST(Nurbs, QuarterCircleIsExact) {
    const double w = std::sqrt(2.0) / 2;
    const NurbsCurve q(2, {0, 0, 0, 1, 1, 1}, {{1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, {1, w, 1});
    for (double t = 0; t <= 1.0; t += 0.01) EXPECT_NEAR(q.point(t).norm(), 1.0, 1e-14);
    EXPECT_NEAR(q.length(), kPi / 2, 1e-10);
}

TEST(Nurbs, FullCircleRadiusCurvatureLength) {
    const Vec3 c(1, -2, 3);
    const auto circle = NurbsCurve::fullCircle(c, {0, 1, 1}, {1, 0, 0}, 5);
    for (double t = 0; t <= 1.0; t += 0.01) {
        EXPECT_NEAR(circle.point(t).distanceTo(c), 5, 1e-12);
        EXPECT_NEAR(circle.curvature(t), 0.2, 1e-10);
    }
    EXPECT_NEAR(circle.length(), 10 * kPi, 1e-9);
    EXPECT_TRUE(nearlyEqual(circle.startPoint(), circle.endPoint()));
}

TEST(Nurbs, UnitWeightsEqualBSpline) {
    const std::vector<double> knots = {0, 0, 0, 0.5, 1, 1, 1};
    const std::vector<Vec3> poles = {{0, 0, 0}, {1, 2, 0}, {3, 1, 1}, {4, 0, 0}};
    const NurbsCurve n(2, knots, poles, {1, 1, 1, 1});
    const BSplineCurve b(2, knots, poles);
    for (double t = 0; t <= 1.0; t += 0.05) {
        EXPECT_TRUE(nearlyEqual(n.point(t), b.point(t), 1e-13));
        EXPECT_TRUE(nearlyEqual(n.derivative(t), b.derivative(t), 1e-12));
        EXPECT_TRUE(nearlyEqual(n.derivative(t, 2), b.derivative(t, 2), 1e-11));
    }
}

TEST(Nurbs, ScalingAllWeightsChangesNothing) {
    const std::vector<double> knots = {0, 0, 0, 1, 1, 1};
    const std::vector<Vec3> poles = {{0, 0, 0}, {1, 2, 0}, {3, 1, 1}};
    const NurbsCurve a(2, knots, poles, {1, 2, 0.5});
    const NurbsCurve b(2, knots, poles, {3, 6, 1.5});
    for (double t = 0; t <= 1.0; t += 0.1) EXPECT_TRUE(nearlyEqual(a.point(t), b.point(t), 1e-13));
}

TEST(Nurbs, RationalDerivativesMatchFiniteDifference) {
    const auto c = NurbsCurve::fullCircle({0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 2);
    const double h = 1e-6;
    for (double t : {0.1, 0.3, 0.6, 0.9}) {
        const Vec3 fd = (c.point(t + h) - c.point(t - h)) / (2 * h);
        EXPECT_TRUE(nearlyEqual(c.derivative(t), fd, 1e-6));
        const Vec3 fd2 = (c.derivative(t + h) - c.derivative(t - h)) / (2 * h);
        EXPECT_TRUE(nearlyEqual(c.derivative(t, 2), fd2, 1e-4));
    }
}

TEST(Nurbs, CreateRejectsBadWeights) {
    const std::vector<double> knots = {0, 0, 0, 1, 1, 1};
    const std::vector<Vec3> poles = {{0, 0, 0}, {1, 1, 0}, {2, 0, 0}};
    EXPECT_FALSE(NurbsCurve::create(2, knots, poles, {1, 0, 1}).ok());
    EXPECT_FALSE(NurbsCurve::create(2, knots, poles, {1, -1, 1}).ok());
    EXPECT_FALSE(NurbsCurve::create(2, knots, poles, {1, 1}).ok());
    EXPECT_TRUE(NurbsCurve::create(2, knots, poles, {1, 3, 1}).ok());
}

TEST(Nurbs, HigherWeightPullsTowardsPole) {
    const std::vector<double> knots = {0, 0, 0, 1, 1, 1};
    const std::vector<Vec3> poles = {{0, 0, 0}, {1, 1, 0}, {2, 0, 0}};
    const NurbsCurve lo(2, knots, poles, {1, 1, 1});
    const NurbsCurve hi(2, knots, poles, {1, 5, 1});
    EXPECT_LT(hi.point(0.5).distanceTo(poles[1]), lo.point(0.5).distanceTo(poles[1]));
}
