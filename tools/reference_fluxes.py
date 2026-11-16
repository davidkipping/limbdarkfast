"""Write the shared test vectors from high-precision quadrature.

    python tools/reference_fluxes.py [--workers N]

Needs NumPy and mpmath. Each reference flux is a direct numerical integral
of the intensity over the occulted part of the star, in 40-digit
arithmetic, so it shares nothing with the closed forms it checks; this
mirrors the BigFloat and numerical-integration checks of Agol, Luger &
Foreman-Mackey (2020). Only the laws' constants (the poles eps and the
matrix G, as stored doubles) come from python/limbdarkfast.py; the
unocculted fluxes are exact.

The occulted flux of a generator g(mu^2) is the integral over annuli,

    O_g = int g(1 - s^2) Phi(s) s ds,

where Phi(s) is the angle of the circle of radius s (about the star's
centre) that lies inside the planet. For a circle Phi is an arccosine; for
an ellipse its ends are the circle's crossings with the outline, the
unit-modulus roots of a quartic in exp(i phi), found with mpmath.polyroots.
Phi is smooth between the radii at which the circle touches the outline
(the critical distances of the outline from the star's centre, again roots
of a quartic), so the integral is split there and done by tanh-sinh
quadrature, which converges exponentially despite the square-root ends.

It writes, for each law L in hybrid2, hybrid4, hybrid5,

    tests/vectors/sphere_L.txt   columns: z r w_1..w_n F
    tests/vectors/oblate_L.txt   columns: X Y r f theta w_1..w_n F

as plain text: '#' header lines, then one point per row, 17 significant
digits (F is the correctly rounded double of the reference). Angles are in
radians. The rows come in three blocks, one per weight vector.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp_pool
import os
import sys

import mpmath as mpm
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python"))
import limbdarkfast  # noqa: E402  (law constants only)

OUT = os.path.join(HERE, "..", "tests", "vectors")
LAWS = ("hybrid2", "hybrid4", "hybrid5")
DPS = 40                     # working precision, decimal digits
ROOT_DPS = 60                # for the quartics (near-double roots)
QUAD_TOL = mpm.mpf(10) ** -30
OFFSETS = (0.0, 1e-13, 1e-10, 1e-7, 1e-4, 1e-2)
rng = np.random.default_rng(20261116)

# three weight vectors per law, inside each law's physical region
WEIGHTS = {
    "hybrid2": [[0.12123702682335066, 0.22967167533748584],
                [0.38587171994566877, 0.5837433480177472],
                [0.15298147401582932, 0.256902579205786]],
    "hybrid4": [[0.3709047608728925, 0.00711438176464726, 0.3162293789763985,
                 0.01888876666845794],
                [0.0044589268157508785, 0.13609141729031174, 0.20138397018412915,
                 0.34602619662222606],
                [0.04782893284628414, 0.6777281164409152, 0.06043401764564811,
                 0.06618308428644647]],
    "hybrid5": [[0.136725258105043, 0.25114827117024674, 0.1088217352833381,
                 0.030866335871229982, 0.32786772243688317],
                [0.18267605012057253, 0.04871020763852243, 0.06971269006071566,
                 0.289247163088596, 0.1936049284750606],
                [0.06602840959945644, 0.15891045247993513, 0.06325441258349657,
                 0.34624553813206715, 0.3110643395711629]],
}


def eps_all():
    return sorted({e for name in LAWS for e in limbdarkfast.get_law(name).eps})


EPS_ALL = eps_all()          # every pole of every law: generators 3, 4, ...


# ---------------------------------------------------------------------------
# geometry: the angle Phi(s) inside the planet, and where to split
# ---------------------------------------------------------------------------

def sphere_setup(z, r):
    """Pieces [(s_a, s_b)] and Phi(s) for a circle of radius r at distance z."""
    z, r = mpm.mpf(z), mpm.mpf(r)
    lo, hi = abs(z - r), min(z + r, mpm.mpf(1))
    pieces = []
    if r > z:                                     # star's centre covered
        pieces.append((mpm.mpf(0), min(r - z, mpm.mpf(1)), "full"))
    if lo < hi:
        pieces.append((lo, hi, "arc"))

    def phi(s, kind):
        if kind == "full":
            return 2 * mpm.pi
        c = (s * s + z * z - r * r) / (2 * z * s)
        return 2 * mpm.acos(max(mpm.mpf(-1), min(mpm.mpf(1), c)))
    return pieces, phi


def _unit_roots(coeffs, near=False):
    """Arguments of the roots of a polynomial (highest power first) that lie
    on the unit circle (with near=True, also of those off it, as a second
    list)."""
    with mpm.workdps(ROOT_DPS):
        c = [mpm.mpc(v) for v in coeffs]
        while abs(c[0]) < mpm.mpf(10) ** (-ROOT_DPS + 10) * max(abs(v) for v in c):
            c = c[1:]
        if len(c) < 2:
            return []
        try:
            roots = mpm.polyroots(c, maxsteps=200, extraprec=2 * ROOT_DPS)
        except mpm.libmp.NoConvergence:
            # roots of very different sizes (a nearly round planet): the
            # eigenvalues of the companion matrix
            n = len(c) - 1
            M = mpm.matrix(n, n)
            for j in range(n):
                M[0, j] = -c[j + 1] / c[0]
            for j in range(1, n):
                M[j, j - 1] = 1
            roots = mpm.eig(M, left=False, right=False)
        tol = mpm.mpf(10) ** -24
        on = [mpm.arg(w) for w in roots if abs(abs(w) - 1) < tol]
        if near:
            return on, [mpm.arg(w) for w in roots if not abs(abs(w) - 1) < tol]
        return on


def ellipse_setup(x0, y0, A, B):
    """Pieces and Phi(s) for the ellipse (x - x0)^2/A^2 + (y - y0)^2/B^2 <= 1."""
    x0, y0, A, B = (mpm.mpf(v) for v in (x0, y0, A, B))

    def inside(x, y):
        return (x - x0) ** 2 / A ** 2 + (y - y0) ** 2 / B ** 2 - 1

    # critical distances of the outline (x0 + A cos t, y0 + B sin t):
    # d(rho^2)/dt = 0, times i w^2 with w = exp(i t), is a quartic in w
    a2 = (B * B - A * A)
    crit, off = _unit_roots([a2, -2 * A * x0 + 2j * B * y0, 0,
                             2 * A * x0 + 2j * B * y0, -a2], near=True)

    def rho(t):
        return mpm.sqrt((x0 + A * mpm.cos(t)) ** 2 + (y0 + B * mpm.sin(t)) ** 2)
    rhos = sorted(rho(t) for t in crit)
    covered = inside(0, 0) < 0
    lo = mpm.mpf(0) if covered else rhos[0]
    hi = min(rhos[-1], mpm.mpf(1))
    # also split at the distances of the off-circle critical points: near an
    # inflection of the distance along the outline, Phi falls steeply there
    # (a singularity just off the real axis)
    extra = [rho(t) for t in off]
    cuts = [lo] + sorted(p for p in rhos + extra if lo < p < hi) + [hi]
    pieces = [(a, b, "ellipse") for a, b in zip(cuts[:-1], cuts[1:]) if b > a]

    def phi(s, kind):
        # 4 A^2 B^2 w^2 * inside(s cos p, s sin p) with w = exp(i p)
        BB, AA = B * B, A * A
        c4 = (BB - AA) * s * s
        c3 = -4 * BB * s * x0 + 4j * AA * s * y0
        c2 = 2 * (BB + AA) * s * s + 4 * BB * x0 * x0 + 4 * AA * y0 * y0 - 4 * AA * BB
        c1 = -4 * BB * s * x0 - 4j * AA * s * y0
        ang = sorted(_unit_roots([c4, c3, c2, c1, c4]))
        if not ang:
            return 2 * mpm.pi if inside(s, 0) < 0 else mpm.mpf(0)
        total = mpm.mpf(0)
        ends = ang + [ang[0] + 2 * mpm.pi]
        for a, b in zip(ends[:-1], ends[1:]):
            m = (a + b) / 2
            if inside(s * mpm.cos(m), s * mpm.sin(m)) < 0:
                total += b - a
        return total
    return pieces, phi


# ---------------------------------------------------------------------------
# tanh-sinh quadrature of a vector integrand
# ---------------------------------------------------------------------------

def tanh_sinh(fvec, a, b, n_out):
    """int_a^b fvec(s) ds for a vector-valued fvec, refining the step until
    successive estimates agree to QUAD_TOL (relative to the largest)."""
    half = (b - a) / 2
    if half == 0:
        return [mpm.mpf(0)] * n_out, mpm.mpf(0)
    T = mpm.mpf(4)                       # weights ~ exp(-pi/2 e^T): far below 1e-40
    cache = {}

    def node(k, h):
        t = k * h
        if t in cache:
            return cache[t]
        u = mpm.pi / 2 * mpm.sinh(t)
        w = mpm.pi / 2 * mpm.cosh(t) / mpm.cosh(u) ** 2
        one_minus = 2 / (1 + mpm.exp(2 * u))           # 1 - tanh(u), accurately
        if t <= 0:
            s = a + half * (2 / (1 + mpm.exp(-2 * u)))  # (1 + x) = 2/(1 + e^{-2u})
        else:
            s = b - half * one_minus
        val = fvec(s)
        cache[t] = [w * v for v in val]
        return cache[t]

    h = mpm.mpf(1) / 2
    n = int(T / h)
    total = [mpm.mpf(0)] * n_out
    for k in range(-n, n + 1):
        total = [x + y for x, y in zip(total, node(k, h))]
    est = [half * h * x for x in total]
    err = mpm.inf
    for _ in range(9):
        h /= 2
        n = int(T / h)
        for k in range(1 - n, n, 2):         # the new, odd nodes
            total = [x + y for x, y in zip(total, node(k, h))]
        new = [half * h * x for x in total]
        scale = max(max(abs(x) for x in new), mpm.mpf(10) ** -30)
        err = max(abs(x - y) for x, y in zip(new, est)) / scale
        est = new
        if err < QUAD_TOL:
            break
    return est, err


def occulted(geom):
    """Occulted flux of every generator [1, mu^2, mu^4, (mu^2 + eps)^-2 ...]
    for one geometry: ('sphere', z, r) or ('oblate', X, Y, r, f, theta).
    Returns the values as strings, and the quadrature's error estimate."""
    mpm.mp.dps = DPS
    if geom[0] == "sphere" or geom[4] == 0.0:
        if geom[0] == "sphere":
            z, r = geom[1], geom[2]
        else:
            z, r = mpm.sqrt(mpm.mpf(geom[1]) ** 2 + mpm.mpf(geom[2]) ** 2), geom[3]
        pieces, phi = sphere_setup(z, r)
    else:
        X, Y, r, f, th = (mpm.mpf(v) for v in geom[1:])
        A = r / mpm.sqrt(1 - f)
        B = A * (1 - f)
        x0 = X * mpm.cos(th) + Y * mpm.sin(th)
        y0 = -X * mpm.sin(th) + Y * mpm.cos(th)
        pieces, phi = ellipse_setup(x0, y0, A, B)
    eps = [mpm.mpf(e) for e in EPS_ALL]
    n_out = 3 + len(eps)
    total = [mpm.mpf(0)] * n_out
    worst = mpm.mpf(0)
    for a, b, kind in pieces:
        def fvec(s, kind=kind):
            u = 1 - s * s                                  # mu^2
            base = phi(s, kind) * s
            return [base, base * u, base * u * u] + [base / (u + e) ** 2 for e in eps]
        vals, err = tanh_sinh(fvec, a, b, n_out)
        total = [x + y for x, y in zip(total, vals)]
        worst = max(worst, err)
    return [mpm.nstr(v, 36) for v in total], float(worst)


def unocculted():
    """The star's flux in each generator: pi int_0^1 g(u) du."""
    mpm.mp.dps = DPS
    return [mpm.pi, mpm.pi / 2, mpm.pi / 3] + \
        [mpm.pi * (1 / mpm.mpf(e) - 1 / (1 + mpm.mpf(e))) for e in EPS_ALL]


def flux(occ, star, name, w):
    """F = 1 - occulted / total for a law and weights, in high precision."""
    mpm.mp.dps = DPS
    L = limbdarkfast.get_law(name)
    cols = [0, 1, 2] + [3 + EPS_ALL.index(e) for e in L.eps]
    occ = [mpm.mpf(v) for v in occ]
    o, n = occ[0], star[0]
    for j in range(L.n_w):
        for c, k in enumerate(cols):
            g = mpm.mpf(float(L.G[j, c]))
            o -= mpm.mpf(w[j]) * g * occ[k]
            n -= mpm.mpf(w[j]) * g * star[k]
    return 1 - o / n


# ---------------------------------------------------------------------------
# the points
# ---------------------------------------------------------------------------

def _segments(r, nb=10, eps=1e-14, delta=1e-3):
    """Agol et al.'s test grid in z for one radius: each contact approached
    linearly and logarithmically from 1e-3 to 1e-14 (src/IDL/test_exofast.pro
    of Limbdark.jl)."""
    lin, log = np.linspace, lambda a, b, n: np.logspace(np.log10(a), np.log10(b), n)
    if r < 1.0:
        parts = [lin(1e-15, eps, nb), lin(eps, delta, nb), lin(delta, r - delta, nb),
                 r - log(delta, eps, nb), lin(r - eps, r + eps, nb), r + log(eps, delta, nb),
                 lin(r + delta, 1 - r - delta, nb), (1 - r) - log(delta, eps, nb),
                 lin(1 - r - eps, 1 - r + eps, nb), (1 - r) + log(eps, delta, nb),
                 lin(1 - r + delta, 1 + r - delta, nb), (1 + r) - log(delta, eps, nb),
                 lin(1 + r - eps, 1 + r - 1e-15, nb)]
    else:
        parts = [lin(r - 1 + 1e-15, r - 1 + eps, nb), (r - 1) + log(eps, delta, nb),
                 lin(r - 1 + delta, r - delta, nb), r - log(delta, eps, nb),
                 lin(r - eps, r + eps, nb), r + log(eps, delta, nb),
                 lin(r + delta, r + 1 - delta, nb), (r + 1) - log(delta, eps, nb),
                 lin(r + 1 - eps, r + 1 - 1e-15, nb)]
    z = np.abs(np.concatenate(parts))
    return z[z < 1.0 + r]


def sphere_points():
    pts = []
    for r in (0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95, 1.0, 1.5, 2.0):
        # a grid across the whole event
        pts += [(z, r) for z in np.linspace(0.0, 1.0 + r + 0.01, 120)]
        # the special lines, approached from both sides: z = r (centre on
        # the limb), the contacts |1 - r| and 1 + r, and for every pole the
        # line z + r = sqrt(1 + eps) where the closed forms switch branch
        special = [r, abs(1.0 - r), 1.0 + r] + [np.sqrt(1.0 + e) - r for e in EPS_ALL]
        for z0 in special:
            for d in OFFSETS:
                for z in (z0 - d, z0 + d):
                    if 0.0 <= z <= 1.0 + r + 0.01:
                        pts.append((z, r))
    # Agol et al.'s radii and contact grids
    for r in (0.01, 0.1, 0.5, 0.9999, 1.0, 1.0001, 2.0):
        pts += [(z, r) for z in _segments(r)]
    return np.array(pts)


def _boundary_extremes(x0, y0, A, B, n=4001):
    """Largest and smallest distance from the star's centre to the outline
    of the ellipse centred at (x0, y0) (principal frame), on a dense grid."""
    t = np.linspace(0.0, 2.0 * np.pi, n)
    rho = np.hypot(x0 + A * np.cos(t), y0 + B * np.sin(t))
    return rho.max(), rho.min()


def _distance_for(rho, A, B, psi, kind):
    """The centre distance d, along principal-frame direction psi, at which
    the outline's largest (kind = 'max') or smallest ('min') distance from
    the star's centre equals rho; by bisection."""
    idx = 0 if kind == "max" else 1
    lo, hi = 0.0, rho + A + 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        v = _boundary_extremes(mid * np.cos(psi), mid * np.sin(psi), A, B)[idx]
        lo, hi = (mid, hi) if v < rho else (lo, mid)
    return 0.5 * (lo + hi)


def oblate_points():
    pts = []
    # 1. chords across the star at random geometries
    fs = (0.0, 1e-12, 3e-10, 1e-8, 1e-5, 1e-3, 0.02, 0.1, 0.2, 0.3, 0.5, 0.7,
          1e-11, 1e-13, 1e-14, 1e-15)            # either side of the switch to a circle
    for f in fs:
        for _ in range(3):
            r = rng.uniform(0.005, 1.0) * (1.0 - f) ** 1.5
            theta = rng.uniform(0.0, np.pi)
            b = rng.uniform(0.0, 1.0 + r)
            A = r / np.sqrt(1.0 - f)
            for X in np.linspace(-(1.0 + A + 0.01), 1.0 + A + 0.01, 41):
                pts.append((X, b, r, f, theta))
    # 2. the edge of the domain, r = (1 - f)^(3/2)
    for f in (0.1, 0.3, 0.5):
        r = (1.0 - f) ** 1.5
        for theta in (0.0, 0.7, np.pi / 2):
            for X in np.linspace(-1.2, 1.2, 25):
                pts.append((X, 0.3, r, f, theta))
    # 3. tangencies: the outline internally tangent to the limb (rho = 1,
    #    'max') or externally ('min'), and internally tangent to every pole
    #    circle rho = sqrt(1 + eps) ('max'), approached from both sides
    for r, f in ((0.01, 0.1), (0.1, 0.1), (0.1, 0.4), (0.2, 0.05), (0.3, 0.3),
                 (0.05, 1e-6)):
        A = r / np.sqrt(1.0 - f)
        B = A * (1.0 - f)
        for theta in (0.0, 0.4, 1.3):
            for psi in (0.3, 2.0):                    # principal-frame direction
                targets = [(1.0, "max"), (1.0, "min")] \
                    + [(np.sqrt(1.0 + e), "max") for e in EPS_ALL]
                for rho, kind in targets:
                    d0 = _distance_for(rho, A, B, psi, kind)
                    for dd in OFFSETS:
                        for d in {d0 - dd, d0 + dd}:
                            x0, y0 = d * np.cos(psi), d * np.sin(psi)
                            c, s = np.cos(theta), np.sin(theta)
                            X, Y = x0 * c - y0 * s, x0 * s + y0 * c
                            pts.append((X, Y, r, f, theta))
    # 4. wholly on the disc: centred and nearly centred, nearly round planets
    for f in (1e-9, 1e-6, 1e-3, 0.1, 0.5):
        for r in (0.01, 0.1, 0.3):
            if r > (1.0 - f) ** 1.5:
                continue
            for d in (0.0, 1e-8, 1e-4, 0.01, 0.3, 0.6):
                pts.append((d, 0.0, r, f, 0.9))
                pts.append((d * 0.6, d * 0.8, r, f, 2.1))
    return np.array(pts)


# ---------------------------------------------------------------------------

def write(path, header, rows):
    with open(path, "w") as fh:
        for line in header:
            fh.write("# " + line + "\n")
        for row in rows:
            fh.write(" ".join(f"{v:.17e}" for v in row) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    sp = sphere_points()
    op = oblate_points()
    star = unocculted()
    head = "reference: tanh-sinh quadrature over annuli, 40 digits (tools/reference_fluxes.py)"
    with mp_pool.Pool(args.workers) as pool:
        for kind, pts in (("sphere", sp), ("oblate", op)):
            geoms = [(kind, *map(float, p)) for p in pts]
            res = pool.map(occulted, geoms, chunksize=4)
            print(f"{kind}: {len(pts)} geometries, worst quadrature estimate "
                  f"{max(e for _, e in res):.1e}")
            for g, (_, e) in zip(geoms, res):
                if e > 1e-24:
                    print(f"  not converged to 1e-24 ({e:.1e}): {g}")
            for name in LAWS:
                W = WEIGHTS[name]
                n_w = len(W[0])
                wcols = " ".join(f"w{j + 1}" for j in range(n_w))
                rows = []
                for k, w in enumerate(W):
                    for i in range(k, len(pts), len(W)):
                        F = float(flux(res[i][0], star, name, w))
                        rows.append((*pts[i], *w, F))
                cols = "z r" if kind == "sphere" else "X Y r f theta"
                write(os.path.join(OUT, f"{kind}_{name}.txt"),
                      [f"limbdarkfast test vectors: {name}, "
                       f"{'spherical' if kind == 'sphere' else 'oblate'} planet",
                       head,
                       f"columns: {cols} {wcols} F"
                       + ("  (theta in radians)" if kind == "oblate" else "")], rows)
                print(f"  {kind}_{name}: {len(rows)} points")


if __name__ == "__main__":
    sys.exit(main())
