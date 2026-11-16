// Example: a spherical and an oblate planet on a circular orbit.
//
// limbdarkfast computes the flux at given sky positions; the orbit is up to the
// caller. On a circular orbit with orbital phase phi = 2 pi (t - t0) / P,
// semi-major axis a and impact parameter b (both in stellar radii), the
// planet's sky position, X along its direction of motion, is
//
//     X = a sin(phi),   Y = -b cos(phi),   in front of the star when cos(phi) > 0.
//
// Writes t, F_sphere, F_oblate as three columns to stdout.

#include <cmath>
#include <cstdio>
#include <vector>

#include "limbdarkfast.hpp"

static const double PI = 3.141592653589793;

int main() {
    const double P = 3.0, t0 = 0.0, a = 10.0, b = 0.3;   // days, days, R*, R*
    const double r = 0.1;                                // area-equivalent if oblate
    const double f = 0.2, theta = 30.0 * PI / 180.0;   // flattening; long axis from X
    const limbdarkfast::Law law = limbdarkfast::Law::hybrid4;
    const std::vector<double> w = {0.30, 0.20, 0.10, 0.05};

    const int n = 2001;
    for (int i = 0; i < n; ++i) {
        const double t = -0.1 + 0.2 * i / (n - 1);
        const double phi = 2.0 * PI * (t - t0) / P;
        const double X = a * std::sin(phi), Y = -b * std::cos(phi);
        double Fs = 1.0, Fo = 1.0;
        if (std::cos(phi) > 0.0) {
            Fs = limbdarkfast::hybrid_flux(std::hypot(X, Y), r, w, law);
            Fo = limbdarkfast::hybrid_flux_oblate(X, Y, r, f, theta, w, law);
        }
        std::printf("%.10f %.17g %.17g\n", t, Fs, Fo);
    }
    return 0;
}
