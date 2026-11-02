#include <gtest/gtest.h>

#include <numeric>

#include "geom/BSplineBasis.h"
#include "geom/BSplineCurve.h"
#include "geom/BezierCurve.h"

using namespace geom;

namespace {
BSplineCurve sampleCubic() {
    std::vector<Vec3> poles = {{0, 0, 0}, {1, 3, 0}, {3, 4, 1}, {5, 1, 0}, {7, 2, -1}, {9, 0, 0}};
    return BSplineCurve(3, {0, 0, 0, 0, 0.3, 0.6, 1, 1, 1, 1}, poles);
}
}  // namespace

TEST(BSplineBasis, PartitionOfUnity) {
    const std::vector<double> knots = {0, 0, 0, 0, 1, 2, 2, 3, 4, 4, 4, 4};
    for (double u = 0; u <= 4.0; u += 0.1) {
        const int span = bspline::findSpan(3, knots, 8, u);
        const auto N = bspline::basisDerivatives(span, u, 3, 2, knots);
        double sum = 0, dsum = 0, d2sum = 0;
        for (int j = 0; j <= 3; ++j) {
            EXPECT_GE(N(0, j), -1e-15);  // non-negative
            sum += N(0, j);
            dsum += N(1, j);
            d2sum += N(2, j);
        }
        EXPECT_NEAR(sum, 1, 1e-14);
        EXPECT_NEAR(dsum, 0, 1e-12);  // derivatives of a constant sum
        EXPECT_NEAR(d2sum, 0, 1e-10);
    }
}

TEST(BSplineBasis, FindSpanHandlesEnds) {
    const std::vector<double> knots = {0, 0, 0, 0.5, 1, 1, 1};
    EXPECT_EQ(bspline::findSpan(2, knots, 4, 0.0), 2);
    EXPECT_EQ(bspline::findSpan(2, knots, 4, 0.5), 3);
    EXPECT_EQ(bspline::findSpan(2, knots, 4, 1.0), 3);  // end maps to last span
    EXPECT_EQ(bspline::findSpan(2, knots, 4, -1.0), 2);
}

TEST(BSplineBasis, Validation) {
    EXPECT_FALSE(bspline::validate(2, {0, 0, 0, 1, 1, 1}, 3));
    EXPECT_TRUE(bspline::validate(0, {0, 1}, 1));              // degree 0
    EXPECT_TRUE(bspline::validate(2, {0, 0, 0, 1, 1}, 3));     // wrong count
    EXPECT_TRUE(bspline::validate(2, {0, 0, 1, 0, 1, 1}, 3));  // decreasing
    EXPECT_TRUE(bspline::validate(1, {0, 0, 0.5, 0.5, 1, 1}, 4));  // multiplicity 2 > degree 1
    EXPECT_TRUE(bspline::validate(2, {0, 0, 0, 0, 0, 0}, 3));  // empty domain
    EXPECT_TRUE(bspline::validate(3, {0, 0, 0, 0, 1, 1, 1}, 2));  // too few poles
}

TEST(BSpline, BezierKnotsGiveBezierCurve) {
    const std::vector<Vec3> poles = {{0, 0, 0}, {1, 2, 0}, {3, 2, 1}, {4, 0, 0}};
    const BSplineCurve bs(3, {0, 0, 0, 0, 1, 1, 1, 1}, poles);
    const BezierCurve bz(poles);
    for (double t = 0; t <= 1.0; t += 0.05) {
        EXPECT_TRUE(nearlyEqual(bs.point(t), bz.point(t), 1e-13));
        EXPECT_TRUE(nearlyEqual(bs.derivative(t), bz.derivative(t), 1e-12));
        EXPECT_TRUE(nearlyEqual(bs.derivative(t, 2), bz.derivative(t, 2), 1e-11));
    }
}

TEST(BSpline, ClampedEndsInterpolate) {
    const auto c = sampleCubic();
    EXPECT_TRUE(nearlyEqual(c.startPoint(), c.poles().front()));
    EXPECT_TRUE(nearlyEqual(c.endPoint(), c.poles().back()));
}

TEST(BSpline, DerivativesMatchFiniteDifference) {
    const auto c = sampleCubic();
    const double h = 1e-6;
    for (double t : {0.1, 0.45, 0.8}) {
        const Vec3 fd = (c.point(t + h) - c.point(t - h)) / (2 * h);
        EXPECT_TRUE(nearlyEqual(c.derivative(t), fd, 1e-6));
        const Vec3 fd2 = (c.derivative(t + h) - c.derivative(t - h)) / (2 * h);
        EXPECT_TRUE(nearlyEqual(c.derivative(t, 2), fd2, 1e-5));
    }
}

TEST(BSpline, KnotInsertionPreservesShape) {
    const auto c = sampleCubic();
    const auto refined = c.insertKnot(0.45);
    ASSERT_TRUE(refined);
    EXPECT_EQ(refined->poles().size(), c.poles().size() + 1);
    EXPECT_EQ(refined->knots().size(), c.knots().size() + 1);
    for (double t = 0; t <= 1.0; t += 0.02)
        EXPECT_TRUE(nearlyEqual(refined->point(t), c.point(t), 1e-12));
}

TEST(BSpline, KnotInsertionAtExistingKnot) {
    const auto c = sampleCubic();
    auto once = c.insertKnot(0.3);
    ASSERT_TRUE(once);
    auto twice = once->insertKnot(0.3);
    ASSERT_TRUE(twice);
    EXPECT_FALSE(twice->insertKnot(0.3));  // multiplicity would exceed degree
    for (double t = 0; t <= 1.0; t += 0.05)
        EXPECT_TRUE(nearlyEqual(twice->point(t), c.point(t), 1e-12));
    // At full multiplicity the curve passes through a pole.
    bool hitsPole = false;
    for (const Vec3& p : twice->poles()) hitsPole |= nearlyEqual(p, c.point(0.3), 1e-12);
    EXPECT_TRUE(hitsPole);
}

TEST(BSpline, KnotInsertionOutsideDomainFails) { EXPECT_FALSE(sampleCubic().insertKnot(1.5)); }

TEST(BSpline, LengthMatchesPolylineLimit) {
    const auto c = sampleCubic();
    double poly = 0;
    const int n = 20000;
    for (int i = 0; i < n; ++i)
        poly += c.point(static_cast<double>(i) / n).distanceTo(c.point(static_cast<double>(i + 1) / n));
    EXPECT_NEAR(c.length(), poly, 1e-6);
}

TEST(BSpline, Breakpoints) {
    const auto bps = sampleCubic().breakpoints();
    ASSERT_EQ(bps.size(), 4u);
    EXPECT_DOUBLE_EQ(bps[1], 0.3);
}

TEST(BSpline, LinearSplineIsPolyline) {
    const BSplineCurve c(1, {0, 0, 1, 2, 2}, {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}});
    EXPECT_TRUE(nearlyEqual(c.point(0.5), Vec3(0.5, 0, 0)));
    EXPECT_TRUE(nearlyEqual(c.point(1.5), Vec3(1, 0.5, 0)));
    EXPECT_NEAR(c.length(), 2, 1e-12);
}
