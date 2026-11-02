#pragma once
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "geom/Curve.h"

// B-spline basis functions shared by BSplineCurve, NurbsCurve and NurbsSurface.
// Algorithms follow Piegl & Tiller, "The NURBS Book", A2.1 and A2.3.
namespace geom::bspline {

inline constexpr int kMaxDegree = 10;
inline constexpr int kMaxDerivative = 2;

// Index of the knot span containing u: knots[span] <= u < knots[span + 1].
// u is clamped into the domain; u at the domain end maps to the last span.
int findSpan(int degree, const std::vector<double>& knots, std::size_t numPoles, double u);

// Values and derivatives of the degree+1 non-zero basis functions on a span:
// (k, j) -> k-th derivative of N_{span-degree+j}(u), for k <= nDers.
class BasisDerivs {
public:
    double operator()(int k, int j) const { return d_[static_cast<std::size_t>(k * kStride + j)]; }
    double& operator()(int k, int j) { return d_[static_cast<std::size_t>(k * kStride + j)]; }

private:
    static constexpr int kStride = kMaxDegree + 1;
    double d_[(kMaxDerivative + 1) * kStride] = {};
};

BasisDerivs basisDerivatives(int span, double u, int degree, int nDers,
                             const std::vector<double>& knots);

// nullopt if valid, otherwise a description of the problem.
std::optional<std::string> validate(int degree, const std::vector<double>& knots,
                                    std::size_t numPoles);

Interval domain(int degree, const std::vector<double>& knots);
// Distinct knot values inside the domain, ends included.
std::vector<double> breakpoints(int degree, const std::vector<double>& knots);
// Clamped knot vector with uniformly spaced interior knots on [0, 1].
std::vector<double> clampedUniformKnots(int degree, std::size_t numPoles);

}  // namespace geom::bspline
