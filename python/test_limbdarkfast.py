"""Check the Python version against the shared test vectors.

    python test_limbdarkfast.py          (or: pytest test_limbdarkfast.py)

Exits non-zero if any point differs from the reference by more than TOL.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import limbdarkfast  # noqa: E402

VECTORS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "tests", "vectors")
TOL = 1e-12
LAWS = ("hybrid2", "hybrid4", "hybrid5")


def _check(kind, law):
    data = np.loadtxt(os.path.join(VECTORS, f"{kind}_{law}.txt"))
    n_w = limbdarkfast.get_law(law).n_w
    F_ref = data[:, -1]
    W = data[:, -1 - n_w:-1]
    F = np.empty_like(F_ref)
    # rows come in blocks sharing one weight vector
    for w in np.unique(W, axis=0):
        sel = np.all(W == w, axis=1)
        d = data[sel]
        if kind == "sphere":
            F[sel] = limbdarkfast.hybrid_flux(d[:, 0], d[:, 1], w, law)
        else:
            F[sel] = limbdarkfast.hybrid_flux_oblate(d[:, 0], d[:, 1], d[:, 2],
                                               d[:, 3], d[:, 4], w, law)
    err = np.abs(F - F_ref)
    return len(F), float(err.max())


def test_vectors():
    for kind in ("sphere", "oblate"):
        for law in LAWS:
            n, err = _check(kind, law)
            assert err <= TOL, f"{kind} {law}: max error {err:.2e}"


def test_law_shapes():
    """Every shape is 1 at the limb and 0 at the centre, and the
    normalisation is the integral of the intensity over the disc."""
    for law in LAWS:
        L = limbdarkfast.get_law(law)
        for j in range(L.n_w):
            w = np.zeros(L.n_w)
            w[j] = 1.0
            assert abs(limbdarkfast.intensity(0.0, w, law)) < 1e-12
            assert abs(limbdarkfast.intensity(1.0, w, law) - 1.0) < 1e-12
            # int_disc S_j dA = 2 pi int_0^1 S_j(mu) mu dmu
            mu = np.linspace(0.0, 1.0, 400001)
            S = 1.0 - limbdarkfast.intensity(mu, w, law)
            y = S * mu
            integral = 2.0 * np.pi * np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(mu))
            assert abs(integral - L.norms[j]) < 1e-8 * max(1.0, L.norms[j])


# ---------------------------------------------------------------------------
# brute-force integration, as in Agol, Luger & Foreman-Mackey (2020)
# (test/test_transit_poly.jl of Limbdark.jl): sum the intensity over NS
# annuli of the star, each weighted by the angle Phi(s) it spends inside
# the planet, with the midpoint rule
# ---------------------------------------------------------------------------

NS = 5000                 # annuli
M_ARC = 48                # samples in angle before bisecting for crossings
NUM_TOL = 1e-6
NUM_W = {"hybrid2": [0.12123702682335066, 0.22967167533748584],
         "hybrid4": [0.3709047608728925, 0.00711438176464726,
                     0.3162293789763985, 0.01888876666845794],
         "hybrid5": [0.136725258105043, 0.25114827117024674, 0.1088217352833381,
                     0.030866335871229982, 0.32786772243688317]}
# their table of (r, b), every contact included, with r <= 2
NUM_RB = [(0.01, 0.0), (0.01, 0.01), (0.01, 0.99), (0.01, 1.0), (0.01, 1.01),
          (0.1, 0.0), (0.1, 0.1), (0.1, 0.9), (0.1, 1.0), (0.1, 1.1),
          (1.0, 0.0), (1.0, 1.0), (1.0, 2.0), (2.0, 1.0), (2.0, 2.0), (2.0, 3.0)]
# the oblate table: the same contacts for r <= 0.5, at f = 0.2, the planet
# displaced along the direction 0.6 rad and tilted by 0.4 rad
NUM_OBL = [(r, b) for r, b in NUM_RB if r < 1.0] + \
    [(0.5, 0.0), (0.5, 0.5), (0.5, 1.0), (0.5, 1.5)]
NUM_F, NUM_BETA, NUM_THETA = 0.2, 0.6, 0.4


def _num_total(w, law):
    s = (np.arange(NS) + 0.5) / NS
    return np.sum(2.0 * np.pi * s / NS * limbdarkfast.intensity(np.sqrt(1.0 - s * s), w, law))


def numeric_flux(z, r, w, law):
    """Spherical planet, by brute force."""
    if z >= 1.0 + r:
        return 1.0
    s1, s2 = max(0.0, z - r), min(1.0, z + r)
    ds = (s2 - s1) / NS
    s = s1 + (np.arange(NS) + 0.5) * ds
    c = np.clip((s * s + z * z - r * r) / (2.0 * z * s + 1e-300), -1.0, 1.0)
    phi = np.where(s < r - z, 2.0 * np.pi, 2.0 * np.arccos(c))
    occ = np.sum(s * ds * phi * limbdarkfast.intensity(np.sqrt(1.0 - s * s), w, law))
    return 1.0 - occ / _num_total(w, law)


def _arc_inside(s, x0, y0, A, B):
    """The angle of the circle of radius s (array) inside the ellipse
    ((x - x0)/A)^2 + ((y - y0)/B)^2 <= 1: sample M_ARC + 1 angles across the
    window the ellipse can occupy, bisect each sign change of the ellipse
    function, and add up the parts of the intervals inside."""
    d = np.hypot(x0, y0)
    if d > A:
        half = np.arcsin(A / d)
        lo, hi = np.arctan2(y0, x0) - half, np.arctan2(y0, x0) + half
    else:
        lo, hi = 0.0, 2.0 * np.pi
    p = lo + (hi - lo) * np.arange(M_ARC + 1) / M_ARC
    gp = ((s[:, None] * np.cos(p) - x0) / A) ** 2 + \
        ((s[:, None] * np.sin(p) - y0) / B) ** 2 - 1.0
    inside = gp < 0.0
    part = np.where(inside[:, :-1] & inside[:, 1:], (hi - lo) / M_ARC, 0.0)
    # bisect the intervals with a sign change
    i, k = np.nonzero(inside[:, :-1] != inside[:, 1:])
    si, a, b = s[i], p[k], p[k + 1]
    ina = inside[i, k]
    for _ in range(50):
        m = 0.5 * (a + b)
        inm = ((si * np.cos(m) - x0) / A) ** 2 + ((si * np.sin(m) - y0) / B) ** 2 < 1.0
        a, b = np.where(inm == ina, m, a), np.where(inm == ina, b, m)
    cross = 0.5 * (a + b)
    np.add.at(part, (i, k), np.where(ina, cross - p[k], p[k + 1] - cross))
    return part.sum(axis=1)


def numeric_flux_oblate(X, Y, r, f, theta, w, law):
    """Oblate planet, by brute force."""
    A = r / np.sqrt(1.0 - f)
    B = A * (1.0 - f)
    x0 = X * np.cos(theta) + Y * np.sin(theta)
    y0 = -X * np.sin(theta) + Y * np.cos(theta)
    d = np.hypot(x0, y0)
    s1, s2 = max(0.0, d - A), min(1.0, d + A)
    if s1 >= s2:
        return 1.0
    ds = (s2 - s1) / NS
    s = s1 + (np.arange(NS) + 0.5) * ds
    phi = _arc_inside(s, x0, y0, A, B)
    occ = np.sum(s * ds * phi * limbdarkfast.intensity(np.sqrt(1.0 - s * s), w, law))
    return 1.0 - occ / _num_total(w, law)


def _numeric_checks():
    """Largest difference from brute force: the table of contacts, then a
    light curve (1000 points spherical, 201 oblate) for r = 0.1."""
    worst = 0.0
    for law in LAWS:
        w = NUM_W[law]
        for r, b in NUM_RB:
            worst = max(worst, abs(limbdarkfast.hybrid_flux(b, r, w, law)
                                   - numeric_flux(b, r, w, law)))
        for r, b in NUM_OBL:
            X, Y = b * np.cos(NUM_BETA), b * np.sin(NUM_BETA)
            worst = max(worst, abs(
                limbdarkfast.hybrid_flux_oblate(X, Y, r, NUM_F, NUM_THETA, w, law)
                - numeric_flux_oblate(X, Y, r, NUM_F, NUM_THETA, w, law)))
    w, law = NUM_W["hybrid4"], "hybrid4"
    for b in 3.0 * (np.arange(1, 1001) / 1000.0 - 0.5):
        worst = max(worst, abs(limbdarkfast.hybrid_flux(abs(b), 0.1, w, law)
                               - numeric_flux(abs(b), 0.1, w, law)))
    for X in np.linspace(-1.5, 1.5, 201):
        worst = max(worst, abs(
            limbdarkfast.hybrid_flux_oblate(X, 0.3, 0.1, NUM_F, NUM_THETA, w, law)
            - numeric_flux_oblate(X, 0.3, 0.1, NUM_F, NUM_THETA, w, law)))
    return worst


def test_numeric():
    worst = _numeric_checks()
    assert worst <= NUM_TOL, f"brute-force integration differs by {worst:.2e}"


def main():
    bad = False
    for kind in ("sphere", "oblate"):
        for law in LAWS:
            n, err = _check(kind, law)
            flag = "ok" if err <= TOL else "FAIL"
            bad |= err > TOL
            print(f"{kind:6s} {law}: {n:5d} points, max |dF| = {err:.2e}  {flag}")
    test_law_shapes()
    print("law shapes and norms: ok")
    worst = _numeric_checks()
    bad |= worst > NUM_TOL
    print(f"brute-force integration: max |dF| = {worst:.2e}  "
          f"{'ok' if worst <= NUM_TOL else 'FAIL'}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
