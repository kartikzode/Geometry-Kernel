#include "geom/BSplineBasis.h"

#include <algorithm>
#include <utility>

namespace geom::bspline {

int findSpan(int degree, const std::vector<double>& knots, std::size_t numPoles, double u) {
    const int n = static_cast<int>(numPoles) - 1;
    const int p = degree;
    if (u >= knots[n + 1]) {
        int k = n;
        while (k > p && knots[k] >= knots[n + 1]) --k;
        return k;
    }
    if (u < knots[p]) u = knots[p];
    int low = p, high = n + 1;
    int mid = (low + high) / 2;
    while (u < knots[mid] || u >= knots[mid + 1]) {
        if (u < knots[mid])
            high = mid;
        else
            low = mid;
        mid = (low + high) / 2;
    }
    return mid;
}

BasisDerivs basisDerivatives(int span, double u, int degree, int nDers,
                             const std::vector<double>& knots) {
    const int p = degree;
    const int i = span;
    BasisDerivs out;

    // Fixed-size scratch on the stack: degree is capped at kMaxDegree, and heap
    // vectors here cost more than the arithmetic (see docs/performance.md).
    constexpr int M = kMaxDegree + 1;
    double ndu[M][M] = {};
    double left[M] = {}, right[M] = {};
    double a[2][M] = {};

    // ndu: basis values (upper triangle) and knot differences (lower triangle).
    ndu[0][0] = 1.0;
    for (int j = 1; j <= p; ++j) {
        left[j] = u - knots[i + 1 - j];
        right[j] = knots[i + j] - u;
        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            ndu[j][r] = right[r + 1] + left[j - r];
            const double temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    for (int j = 0; j <= p; ++j) out(0, j) = ndu[j][p];

    const int n = std::min(nDers, p);  // derivatives above the degree are zero
    for (int r = 0; r <= p; ++r) {
        int s1 = 0, s2 = 1;
        a[0][0] = 1.0;
        for (int k = 1; k <= n; ++k) {
            double d = 0.0;
            const int rk = r - k, pk = p - k;
            if (r >= k) {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }
            const int j1 = rk >= -1 ? 1 : -rk;
            const int j2 = (r - 1 <= pk) ? k - 1 : p - r;
            for (int j = j1; j <= j2; ++j) {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                d += a[s2][j] * ndu[rk + j][pk];
            }
            if (r <= pk) {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }
            out(k, r) = d;
            std::swap(s1, s2);
        }
    }
    double factor = p;
    for (int k = 1; k <= n; ++k) {
        for (int j = 0; j <= p; ++j) out(k, j) *= factor;
        factor *= (p - k);
    }
    return out;
}

std::optional<std::string> validate(int degree, const std::vector<double>& knots,
                                    std::size_t numPoles) {
    if (degree < 1 || degree > kMaxDegree)
        return "degree must be between 1 and " + std::to_string(kMaxDegree);
    const auto p = static_cast<std::size_t>(degree);
    if (numPoles < p + 1) return "need at least degree + 1 control points";
    if (knots.size() != numPoles + p + 1)
        return "knot count must equal control points + degree + 1 (expected " +
               std::to_string(numPoles + p + 1) + ", got " + std::to_string(knots.size()) + ")";
    if (!std::is_sorted(knots.begin(), knots.end())) return "knots must be non-decreasing";
    // Interior multiplicity above the degree would make the curve discontinuous.
    for (std::size_t i = p + 1; i < numPoles;) {
        std::size_t j = i;
        while (j < numPoles && knots[j] == knots[i]) ++j;
        if (j - i > p) return "interior knot multiplicity exceeds degree";
        i = j;
    }
    if (!(knots[p] < knots[numPoles])) return "knot vector has an empty domain";
    return std::nullopt;
}

Interval domain(int degree, const std::vector<double>& knots) {
    const auto p = static_cast<std::size_t>(degree);
    return {knots[p], knots[knots.size() - p - 1]};
}

std::vector<double> breakpoints(int degree, const std::vector<double>& knots) {
    const Interval d = domain(degree, knots);
    std::vector<double> out;
    for (double k : knots)
        if (k >= d.lo && k <= d.hi && (out.empty() || k > out.back())) out.push_back(k);
    return out;
}

std::vector<double> clampedUniformKnots(int degree, std::size_t numPoles) {
    const auto p = static_cast<std::size_t>(degree);
    const std::size_t interior = numPoles - p - 1;
    std::vector<double> k(p + 1, 0.0);
    for (std::size_t i = 1; i <= interior; ++i)
        k.push_back(static_cast<double>(i) / static_cast<double>(interior + 1));
    k.insert(k.end(), p + 1, 1.0);
    return k;
}

}  // namespace geom::bspline
