#include <gtest/gtest.h>

#include "geom/BezierCurve.h"
#include "geom/Line.h"

using namespace geom;

namespace {
BezierCurve sampleCubic() { return BezierCurve({{0, 0, 0}, {1, 2, 0}, {3, 2, 1}, {4, 0, 0}}); }

// Explicit Bernstein form for checking de Casteljau.
Vec3 bernsteinCubic(const BezierCurve& c, double t) {
    const auto& P = c.poles();
    const double u = 1 - t;
    return P[0] * (u * u * u) + P[1] * (3 * u * u * t) + P[2] * (3 * u * t * t) + P[3] * (t * t * t);
}
}  // namespace

TEST(Bezier, InterpolatesEndPoles) {
    const auto c = sampleCubic();
    EXPECT_TRUE(nearlyEqual(c.point(0), c.poles().front()));
    EXPECT_TRUE(nearlyEqual(c.point(1), c.poles().back()));
    EXPECT_EQ(c.degree(), 3);
}

TEST(Bezier, MatchesBernsteinForm) {
    const auto c = sampleCubic();
    for (double t = 0; t <= 1.0; t += 0.125)
        EXPECT_TRUE(nearlyEqual(c.point(t), bernsteinCubic(c, t), 1e-13));
}

TEST(Bezier, EndTangentsFollowControlPolygon) {
    const auto c = sampleCubic();
    const auto& P = c.poles();
    EXPECT_TRUE(nearlyEqual(c.derivative(0), (P[1] - P[0]) * 3));
    EXPECT_TRUE(nearlyEqual(c.derivative(1), (P[3] - P[2]) * 3));
}

TEST(Bezier, DerivativeMatchesFiniteDifference) {
    const auto c = sampleCubic();
    const double t = 0.37, h = 1e-6;
    const Vec3 fd = (c.point(t + h) - c.point(t - h)) / (2 * h);
    EXPECT_TRUE(nearlyEqual(c.derivative(t), fd, 1e-7));
    const Vec3 fd2 = (c.derivative(t + h) - c.derivative(t - h)) / (2 * h);
    EXPECT_TRUE(nearlyEqual(c.derivative(t, 2), fd2, 1e-6));
}

TEST(Bezier, LinearBezierIsALine) {
    const BezierCurve b({{0, 0, 0}, {2, 2, 2}});
    const Line l({0, 0, 0}, {2, 2, 2});
    EXPECT_TRUE(nearlyEqual(b.point(0.3), l.point(0.3)));
    EXPECT_NEAR(b.length(), l.length(), 1e-12);
    EXPECT_TRUE(nearlyEqual(b.derivative(0.5, 2), Vec3()));
}

TEST(Bezier, SubdivisionReproducesCurve) {
    const auto c = sampleCubic();
    const double s = 0.4;
    const auto [left, right] = c.subdivide(s);
    for (double u = 0; u <= 1.0; u += 0.1) {
        EXPECT_TRUE(nearlyEqual(left.point(u), c.point(s * u), 1e-13));
        EXPECT_TRUE(nearlyEqual(right.point(u), c.point(s + (1 - s) * u), 1e-13));
    }
    EXPECT_NEAR(left.length() + right.length(), c.length(), 1e-10);
}

TEST(Bezier, CreateNeedsTwoPoles) {
    EXPECT_FALSE(BezierCurve::create({{1, 1, 1}}).ok());
    EXPECT_TRUE(BezierCurve::create({{0, 0, 0}, {1, 1, 1}}).ok());
}
