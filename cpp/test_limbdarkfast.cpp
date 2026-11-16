// Check the C++ version against the shared test vectors.
//
//     ./test_limbdarkfast [path/to/tests/vectors]
//
// Exits non-zero if any point differs from the reference by more than TOL.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "limbdarkfast.hpp"

static const double TOL = 1e-12;

static std::vector<std::vector<double> > read_rows(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) throw std::runtime_error("cannot open " + path);
    std::vector<std::vector<double> > rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::vector<double> row;
        double v;
        while (ss >> v) row.push_back(v);
        rows.push_back(row);
    }
    return rows;
}

// ---------------------------------------------------------------------------
// brute-force integration, as in Agol, Luger & Foreman-Mackey (2020)
// (test/test_transit_poly.jl of Limbdark.jl): sum the intensity over NS
// annuli of the star, each weighted by the angle Phi(s) it spends inside
// the planet, with the midpoint rule
// ---------------------------------------------------------------------------

static const int NS = 5000;          // annuli
static const int M_ARC = 48;         // samples in angle before bisecting for crossings
static const double NUM_TOL = 1e-6;
static const double PI = 3.141592653589793;
static const double NUM_F = 0.2, NUM_BETA = 0.6, NUM_THETA = 0.4;

static double num_total(const double* w, limbdarkfast::Law law) {
    double t = 0.0;
    for (int j = 0; j < NS; ++j) {
        const double s = (j + 0.5) / NS;
        t += 2.0 * PI * s / NS * limbdarkfast::intensity(std::sqrt(1.0 - s * s), w, law);
    }
    return t;
}

// spherical planet, by brute force
static double numeric_flux(double z, double r, const double* w, limbdarkfast::Law law) {
    if (z >= 1.0 + r) return 1.0;
    const double s1 = std::max(0.0, z - r), s2 = std::min(1.0, z + r);
    const double ds = (s2 - s1) / NS;
    double occ = 0.0;
    for (int j = 0; j < NS; ++j) {
        const double s = s1 + (j + 0.5) * ds;
        double phi = 2.0 * PI;
        if (!(s < r - z)) {
            const double c = (s * s + z * z - r * r) / (2.0 * z * s + 1e-300);
            phi = 2.0 * std::acos(std::max(-1.0, std::min(1.0, c)));
        }
        occ += s * ds * phi * limbdarkfast::intensity(std::sqrt(1.0 - s * s), w, law);
    }
    return 1.0 - occ / num_total(w, law);
}

// The angle of the circle of radius s inside the ellipse
// ((x - x0)/A)^2 + ((y - y0)/B)^2 <= 1: sample M_ARC + 1 angles across the
// window the ellipse can occupy, bisect each sign change, and add up the
// parts of the intervals inside.
static double arc_inside(double s, double x0, double y0, double A, double B) {
    const double d = std::hypot(x0, y0);
    double lo = 0.0, hi = 2.0 * PI;
    if (d > A) {
        const double half = std::asin(A / d), c = std::atan2(y0, x0);
        lo = c - half;
        hi = c + half;
    }
    const double h = (hi - lo) / M_ARC;
    struct In {
        double s, x0, y0, A, B;
        bool operator()(double p) const {
            const double u = (s * std::cos(p) - x0) / A, v = (s * std::sin(p) - y0) / B;
            return u * u + v * v < 1.0;
        }
    } inside = {s, x0, y0, A, B};
    double total = 0.0;
    bool in_a = inside(lo);
    for (int k = 0; k < M_ARC; ++k) {
        const double pa = lo + h * k, pb = lo + h * (k + 1);
        const bool in_b = inside(pb);
        if (in_a && in_b) {
            total += h;
        } else if (in_a != in_b) {
            double a = pa, b = pb;
            for (int it = 0; it < 50; ++it) {
                const double m = 0.5 * (a + b);
                if (inside(m) == in_a) a = m; else b = m;
            }
            const double cross = 0.5 * (a + b);
            total += in_a ? cross - pa : pb - cross;
        }
        in_a = in_b;
    }
    return total;
}

// oblate planet, by brute force
static double numeric_flux_oblate(double X, double Y, double r, double f, double theta,
                                  const double* w, limbdarkfast::Law law) {
    const double A = r / std::sqrt(1.0 - f), B = A * (1.0 - f);
    const double x0 = X * std::cos(theta) + Y * std::sin(theta);
    const double y0 = -X * std::sin(theta) + Y * std::cos(theta);
    const double d = std::hypot(x0, y0);
    const double s1 = std::max(0.0, d - A), s2 = std::min(1.0, d + A);
    if (s1 >= s2) return 1.0;
    const double ds = (s2 - s1) / NS;
    double occ = 0.0;
    for (int j = 0; j < NS; ++j) {
        const double s = s1 + (j + 0.5) * ds;
        occ += s * ds * arc_inside(s, x0, y0, A, B)
            * limbdarkfast::intensity(std::sqrt(1.0 - s * s), w, law);
    }
    return 1.0 - occ / num_total(w, law);
}

// Largest difference from brute force: their table of (r, b) with every
// contact (r <= 2), its oblate twin, then a light curve for r = 0.1
// (1000 points spherical, 201 oblate).
static double numeric_checks() {
    static const double W[3][5] = {
        {0.12123702682335066, 0.22967167533748584},
        {0.3709047608728925, 0.00711438176464726, 0.3162293789763985, 0.01888876666845794},
        {0.136725258105043, 0.25114827117024674, 0.1088217352833381,
         0.030866335871229982, 0.32786772243688317}};
    static const double RB[16][2] = {
        {0.01, 0.0}, {0.01, 0.01}, {0.01, 0.99}, {0.01, 1.0}, {0.01, 1.01},
        {0.1, 0.0}, {0.1, 0.1}, {0.1, 0.9}, {0.1, 1.0}, {0.1, 1.1},
        {1.0, 0.0}, {1.0, 1.0}, {1.0, 2.0}, {2.0, 1.0}, {2.0, 2.0}, {2.0, 3.0}};
    static const double RB_OBL[4][2] = {{0.5, 0.0}, {0.5, 0.5}, {0.5, 1.0}, {0.5, 1.5}};
    const limbdarkfast::Law laws[3] = {limbdarkfast::Law::hybrid2, limbdarkfast::Law::hybrid4,
                                 limbdarkfast::Law::hybrid5};
    double worst = 0.0;
    for (int il = 0; il < 3; ++il) {
        const double* w = W[il];
        const limbdarkfast::Law law = laws[il];
        for (int i = 0; i < 16; ++i)
            worst = std::max(worst, std::fabs(limbdarkfast::hybrid_flux(RB[i][1], RB[i][0], w, law)
                                              - numeric_flux(RB[i][1], RB[i][0], w, law)));
        for (int i = 0; i < 14; ++i) {
            const double r = i < 10 ? RB[i][0] : RB_OBL[i - 10][0];
            const double b = i < 10 ? RB[i][1] : RB_OBL[i - 10][1];
            const double X = b * std::cos(NUM_BETA), Y = b * std::sin(NUM_BETA);
            worst = std::max(worst, std::fabs(
                limbdarkfast::hybrid_flux_oblate(X, Y, r, NUM_F, NUM_THETA, w, law)
                - numeric_flux_oblate(X, Y, r, NUM_F, NUM_THETA, w, law)));
        }
    }
    const double* w = W[1];
    const limbdarkfast::Law law = limbdarkfast::Law::hybrid4;
    for (int i = 1; i <= 1000; ++i) {
        const double b = std::fabs(3.0 * (i / 1000.0 - 0.5));
        worst = std::max(worst, std::fabs(limbdarkfast::hybrid_flux(b, 0.1, w, law)
                                          - numeric_flux(b, 0.1, w, law)));
    }
    for (int i = 0; i <= 200; ++i) {
        const double X = -1.5 + 3.0 * i / 200.0;
        worst = std::max(worst, std::fabs(
            limbdarkfast::hybrid_flux_oblate(X, 0.3, 0.1, NUM_F, NUM_THETA, w, law)
            - numeric_flux_oblate(X, 0.3, 0.1, NUM_F, NUM_THETA, w, law)));
    }
    return worst;
}

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? argv[1] : "../tests/vectors";
    const char* laws[3] = {"hybrid2", "hybrid4", "hybrid5"};
    const char* kinds[2] = {"sphere", "oblate"};
    int failures = 0;
    for (int ik = 0; ik < 2; ++ik) {
        for (int il = 0; il < 3; ++il) {
            const std::string kind = kinds[ik];
            const limbdarkfast::Law law = limbdarkfast::law_from_name(laws[il]);
            const int n_w = limbdarkfast::n_weights(law);
            const int n_geo = kind == "sphere" ? 2 : 5;
            const std::vector<std::vector<double> > rows =
                read_rows(dir + "/" + kind + "_" + laws[il] + ".txt");
            double worst = 0.0;
            for (std::size_t i = 0; i < rows.size(); ++i) {
                const std::vector<double>& d = rows[i];
                const double* w = &d[n_geo];
                const double F = kind == "sphere"
                    ? limbdarkfast::hybrid_flux(d[0], d[1], w, law)
                    : limbdarkfast::hybrid_flux_oblate(d[0], d[1], d[2], d[3], d[4], w, law);
                const double err = std::fabs(F - d[n_geo + n_w]);
                if (!(err <= worst)) worst = err;      // NaN counts as worst
            }
            const bool ok = worst <= TOL;
            failures += !ok;
            std::printf("%-6s %s: %5zu points, max |dF| = %.2e  %s\n", kind.c_str(),
                        laws[il], rows.size(), worst, ok ? "ok" : "FAIL");
        }
    }
    // every shape is 1 at the limb and 0 at the centre
    for (int il = 0; il < 3; ++il) {
        const limbdarkfast::Law law = limbdarkfast::law_from_name(laws[il]);
        for (int j = 0; j < limbdarkfast::n_weights(law); ++j) {
            double w[5] = {0, 0, 0, 0, 0};
            w[j] = 1.0;
            if (std::fabs(limbdarkfast::intensity(0.0, w, law)) > 1e-12 ||
                std::fabs(limbdarkfast::intensity(1.0, w, law) - 1.0) > 1e-12) {
                std::printf("shape %d of %s is not 1 at the limb and 0 at the centre\n",
                            j + 1, laws[il]);
                ++failures;
            }
        }
    }
    std::printf("law shapes: %s\n", failures ? "see above" : "ok");
    const double worst = numeric_checks();
    failures += !(worst <= NUM_TOL);
    std::printf("brute-force integration: max |dF| = %.2e  %s\n", worst,
                worst <= NUM_TOL ? "ok" : "FAIL");
    return failures ? 1 : 0;
}
