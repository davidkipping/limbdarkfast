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
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
