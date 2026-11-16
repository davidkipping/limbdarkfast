"""
    LimbDarkFast

Transit light curves for the hybrid limb-darkening laws, the Julia version
of limbdarkfast. It gives the flux of a star with hybrid2, hybrid4 or
hybrid5 limb darkening, occulted by a spherical or an oblate planet, in
closed form and in double precision. There are no orbits, exposures or
derivatives: give it the planet's position on the sky and it returns the
normalised flux.

If you use this code, please cite

    Kipping, D. 2027, "Limb Darkened Transits of Oblate Planets as
    Elementary Functions", AJ, submitted.

Copyright (c) 2026 David Kipping. Released under the MIT License: see
the LICENSE file, or https://opensource.org/licenses/MIT.

    F = hybrid_flux(z, r, w, :hybrid4)                        # spherical
    F = hybrid_flux_oblate(X, Y, r, f, theta, w, :hybrid4)    # oblate

Lengths are in stellar radii. `z` is the sky separation of the planet's
centre from the star's; `(X, Y)` is its sky position, X along the direction
of motion. `r` is the planet's radius, for an oblate planet its
area-equivalent radius √(AB); `f = 1 - B/A` is the projected flattening and
`theta` (radians) the angle of the long axis from X. `z` (or `X` and `Y`)
may be numbers or arrays. `w` are the law's shape weights:

    I(μ) / I(1) = 1 - Σⱼ wⱼ Sⱼ(μ),

each shape falling from 1 at the limb to 0 at the centre. The shapes
combine μ⁰, μ², μ⁴ with double poles (μ² + ε)⁻², so that the occulted flux
of every term is elementary for any elliptical planet. The oblate planet
must satisfy r ≤ (1 - f)^(3/2), so that its outline crosses the stellar
limb at most twice; otherwise a `DomainError` is thrown.

Method (python/limbdarkfast.py has the same structure, with more comments).
Green's theorem with the radial potential h(s) = (1/2s) ∫₀ˢ I, s = x² + y²,
turns the occulted flux into h(1) times the limb arc inside the planet plus
the integral of h(s) along the planet's outline inside the star. For a
circle that reduces to the lens geometry; for an ellipse the even terms are
finite Fourier sums and each pole a partial-fraction sum over the roots of
a quartic, near-double roots being integrated as one quadratic factor.
"""
module LimbDarkFast

export hybrid_flux, hybrid_flux_oblate, intensity, nweights, HYBRID2, HYBRID4, HYBRID5

const C64 = ComplexF64

const TWO_PI = 2π
const EPS = eps(Float64)
const TINY = 1e-280
const Y_SW = 1e-3        # spherical pole: series half-width in y
const F_SW = 1e-10       # below this flattening, the spherical forms
const Z_TOL = 1e-4       # |z| ~ 1: a candidate limb crossing
const G_TOL = 1e-12      # accepted crossing residual, relative
const PAIR_TOL = 1e-6    # crossings closer than this are a tangency
const PAIR_RATIO = 0.1   # near-double pole roots
const N_SERIES = 18
const N_ABERTH = 8
const N_BAIRSTOW = 2
const QR_SWEEPS4, QR_SWEEPS3, QR_BALANCE = 6, 5, 4
const PAIRS = ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))

"""
A hybrid law: `nw` weights and `K` poles `eps`; `G[j]` gives shape Sⱼ as a
combination of the generators [μ⁰, μ², μ⁴, (μ² + ε_k)⁻² ...]; `norms[j]`
is the unocculted flux of Sⱼ. From tools/law_constants.py.
"""
struct Law
    name::Symbol
    nw::Int
    K::Int
    eps::NTuple{3,Float64}
    G::NTuple{5,NTuple{6,Float64}}
    norms::NTuple{5,Float64}
end

Base.show(io::IO, L::Law) = print(io, "Law(:", L.name, ", eps = ", L.eps[1:L.K], ")")

const HYBRID2 = Law(:hybrid2, 2, 1,
    (0.35999999999999999, 0.0, 0.0),
    ((1.0, -1.0, 0.0, 0.0, 0.0, 0.0),
     (-0.075348837209302341, 0.0, 0.0, 0.13936520930232557, 0.0, 0.0),
     (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
     (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
     (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)),
    (1.5707963267948966, 0.6575426484257707, 0.0, 0.0, 0.0))
const HYBRID4 = Law(:hybrid4, 4, 2,
    (0.028601207249486824, 0.38582417584735224, 0.0),
    ((1.0, -1.0, 0.0, 0.0, 0.0, 0.0),
     (1.0, -2.0, 1.0, 0.0, 0.0, 0.0),
     (0.0031720957997626959, -0.0028080275108228465, 0.0007779295064964302,
      0.0008358238567779617, -0.003710396155608333, 0.0),
     (-0.45048769921539672, 0.51370715445091153, -0.17569337890652659,
      -4.8355348357818505e-07, 0.21600802064223382, 0.0),
     (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)),
    (1.5707963267948966, 1.0471975511965974, 0.073823572890099914,
     0.47682011132930302, 0.0))
const HYBRID5 = Law(:hybrid5, 5, 3,
    (0.0015965451792942786, 0.041275044772743019, 0.38743471584956696),
    ((1.0, 0.0, -1.0, 0.0, 0.0, 0.0),
     (1.0, -2.0, 1.0, 0.0, 0.0, 0.0),
     (-2.5408433670697262e-06, 0.0, 0.0, 2.5489629860270405e-06, 0.0, 0.0),
     (-0.0015737187012817282, 0.0, 0.0, 0.0, 0.0017063103543144378, 0.0),
     (-0.084572789617394117, 0.0, 0.0, 0.0, 0.0, 0.16280051336826454)),
    (2.0943951023931957, 1.0471975511965974, 0.0049997300166832875,
     0.11978141121311216, 0.68577554794744777))

getlaw(L::Law) = L
function getlaw(name::Union{Symbol,AbstractString})
    s = Symbol(name)
    s === :hybrid2 && return HYBRID2
    s === :hybrid4 && return HYBRID4
    s === :hybrid5 && return HYBRID5
    throw(ArgumentError("unknown law $(repr(name)); expected :hybrid2, :hybrid4 or :hybrid5"))
end

"""
    nweights(law)

Number of shape weights `law` takes: 2, 4 or 5.
"""
nweights(law) = getlaw(law).nw

function checkweights(w, L::Law)
    length(w) == L.nw ||
        throw(ArgumentError("$(L.name) takes $(L.nw) weights; got $(length(w))"))
end

"""
    intensity(mu, w, law=:hybrid4)

I(μ) / I(1) = 1 - Σⱼ wⱼ Sⱼ(μ).
"""
function intensity(mu::Real, w, law=HYBRID4)
    L = getlaw(law)
    checkweights(w, L)
    m2 = Float64(mu)^2
    gen = (1.0, m2, m2 * m2, ntuple(k -> 1.0 / (m2 + L.eps[k])^2, 3)...)
    I = 1.0
    for j in 1:L.nw
        I -= w[j] * sum(L.G[j][k] * gen[k] for k in 1:3+L.K)
    end
    return I
end

# F = 1 - (O_0 - Σⱼ wⱼ O_Sj) / (π - Σⱼ wⱼ Nⱼ), from the occulted flux of each
# generator
function flux_from_occ(occ::NTuple{6,Float64}, w, L::Law)
    blocked, norm = occ[1], Float64(π)
    for j in 1:L.nw
        blocked -= w[j] * sum(L.G[j][k] * occ[k] for k in 1:3+L.K)
        norm -= w[j] * L.norms[j]
    end
    return 1.0 - blocked / norm
end

# ===========================================================================
# spherical planet
# ===========================================================================

# 16 x (area of the triangle with sides 1, z, r)², sorted (Kahan) form;
# negative when no such triangle exists
function kite_sqarea(z, r)
    a, c = max(z, r), min(z, r)
    b = min(a, 1.0)
    a = max(a, 1.0)
    b, c = max(b, c), min(b, c)
    return (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c))
end

# occulted flux of [μ⁰, μ², μ⁴, P_ε...] for a circle of radius r at separation z
function sphere_occ(z::Float64, r::Float64, L::Law)
    z = abs(z)
    z >= 1.0 + r && return ntuple(_ -> 0.0, 6)
    r2, z2 = r * r, z * z
    alpha = 1.0 - r2 - z2
    sqarea = kite_sqarea(z, r)
    A = max((r + 1.0 - z) * (1.0 - r + z), TINY)    # 1 - (z - r)²
    Bp = (1.0 - z - r) * (1.0 + z + r)               # 1 - (z + r)²

    if z <= 1.0 - r                                 # wholly on the disc
        s0 = -π * r2
        s2 = TWO_PI * r2 * (r2 + 2.0 * z2 - 1.0)
        M0, M2 = Float64(π), π * alpha
        M4 = (6.0 * alpha * M2 + 2.0 * sqarea * M0) / 4.0
        s4 = -(2.0 * r2 * M4 - (2.0 / 3.0) * (alpha * M4 + sqarea * M2))
        poles = ntuple(3) do k
            k > L.K && return 0.0
            e = L.eps[k]
            p = 1.0 + e
            Kp = p + r2 - z2
            sqQ = sqrt((e + A) * (e + Bp))
            TWO_PI * r2 / (sqQ * (Kp + sqQ))
        end
        return (-s0, -(0.5 * s0 + 0.25 * s2), -(s0 / 3.0 + s2 / 6.0 + s4 / 6.0), poles...)
    end

    # partial overlap (or, for r > 1, total occultation: kite = 0 there, and
    # every term lands on its full-disc value)
    tot = z <= r - 1.0
    kite = tot ? 0.0 : sqrt(max(sqarea, (10.0 * EPS)^2))
    kap0 = atan(kite, r2 + z2 - 1.0)    # half-angle at the planet
    kap1 = atan(kite, 1.0 + z2 - r2)    # half-angle at the star
    s0 = -(kap1 + r2 * kap0 - 0.5 * kite)
    eta2 = r2 * (r2 + 2.0 * z2)
    s2 = 2.0 * s0 + 2.0 * (kap1 + eta2 * kap0 - 0.25 * kite * (1.0 + 5.0 * r2 + z2))
    M0, M2 = kap0, kap0 * alpha + kite
    M4 = (6.0 * alpha * M2 + 2.0 * sqarea * M0) / 4.0
    s4 = -(2.0 * r2 * M4 - (2.0 / 3.0) * (alpha * M4 + sqarea * M2))

    # the smaller of A and -Bp is taken from kite² = A (-Bp), exactly
    nBp = (!tot && A >= -Bp) ? kite * kite / A : (-Bp > 0.0 ? -Bp : kite * kite / A)
    poles = ntuple(3) do k
        k > L.K && return 0.0
        e = L.eps[k]
        p = 1.0 + e
        Kp = p + r2 - z2
        Q = (e + A) * (e + Bp)
        # J = ∫_{-kap0}^{kap0} dψ / (p - z² - r² + 2 z r cos ψ)
        apb = e + A
        U = apb * nBp
        ratio = (e + Bp) / nBp
        y = ratio * kite * kite / U
        J = if y > Y_SW
            4.0 / sqrt(Q) * atan(kite * sqrt(ratio), sqrt(U))
        elseif y < -Y_SW
            4.0 / sqrt(-Q) * log((sqrt(U) + kite * sqrt(-ratio)) /
                                 (2.0 * sqrt(max(z * r * e, TINY))))
        else                            # across Q = 0: a series
            T = kite / nBp              # tan(kap0 / 2)
            y2 = y * y
            G = 1.0 - y / 3.0 + y2 / 5.0 - y * y2 / 7.0 + y2 * y2 / 9.0 - y * y2 * y2 / 11.0
            4.0 * T * G / apb
        end
        kap1 / (p * e) + (Kp * J - 2.0 * kap0) / (4.0 * p)
    end
    return (-s0, -(0.5 * s0 + 0.25 * s2), -(s0 / 3.0 + s2 / 6.0 + s4 / 6.0), poles...)
end

# ===========================================================================
# oblate planet: the ellipse x = x0 + a cos φ, y = y0 + b sin φ. With
# z = e^{iφ}, s(φ) = x² + y² = S0 + 2 Re(S1 z) + 2 S2 cos 2φ and
# (x dy - y dx)/dφ = K0 + 2 Re(K1 z); the points where s = p are the roots
# of the self-inversive quartic D(z) = S2 z⁴ + S1 z³ + (S0 - p) z² + S1* z + S2.
# ===========================================================================

struct Ellipse
    S0::Float64
    S1::C64
    S2::Float64
    K0::Float64
    K1::C64
end

Ellipse(x0, y0, a, b) = Ellipse(x0 * x0 + y0 * y0 + 0.5 * (a * a + b * b),
                                C64(a * x0, -b * y0), 0.25 * (a * a - b * b),
                                a * b, C64(0.5 * b * x0, -0.5 * a * y0))

Dz(E::Ellipse, z, p) = (((E.S2 * z + E.S1) * z + (E.S0 - p)) * z + conj(E.S1)) * z + E.S2
dDz(E::Ellipse, z, p) = ((4.0 * E.S2 * z + 3.0 * E.S1) * z + 2.0 * (E.S0 - p)) * z + conj(E.S1)
Nz(E::Ellipse, z) = conj(E.K1) + E.K0 * z + E.K1 * z * z

# ---- the quartic's roots: a balanced, shifted 4x4 QR, then Aberth --------

# (c, s), c real: [[c, s], [-conj(s), c]] [x, y] = [ρ, 0]
function givens(x::C64, y::C64)
    xa = abs(x)
    rho = sqrt(xa * xa + abs2(y))
    ph = xa > 0.0 ? x / xa : one(C64)
    return rho > 0.0 ? (xa / rho, ph * conj(y) / rho) : (1.0, zero(C64))
end

# the eigenvalue of [[a, b], [c, d]] nearer d
function wilkinson(a, b, c, d)
    h = 0.5 * (a - d)
    disc = sqrt(h * h + b * c)
    den = abs(h + disc) >= abs(h - disc) ? h + disc : h - disc
    return abs(den) > 0.0 ? d - b * c / den : d
end

# one shifted QR step on the leading m x m block of H
function qr_sweep!(H::Matrix{C64}, m::Int)
    sig = wilkinson(H[m-1, m-1], H[m-1, m], H[m, m-1], H[m, m])
    for i in 1:m
        H[i, i] -= sig
    end
    cs = (1.0, 1.0, 1.0)
    sn = (zero(C64), zero(C64), zero(C64))
    for k in 1:m-1
        c, s = givens(H[k, k], H[k+1, k])
        cs = Base.setindex(cs, c, k)
        sn = Base.setindex(sn, s, k)
        for j in k:m
            x, y = H[k, j], H[k+1, j]
            H[k, j] = cs[k] * x + sn[k] * y
            H[k+1, j] = cs[k] * y - conj(sn[k]) * x
        end
    end
    for k in 1:m-1, i in 1:min(k + 2, m)
        x, y = H[i, k], H[i, k+1]
        H[i, k] = cs[k] * x + conj(sn[k]) * y
        H[i, k+1] = cs[k] * y - sn[k] * x
    end
    for i in 1:m
        H[i, i] += sig
    end
    return H
end

# eigenvalues of D's companion matrix
function qr_seeds(E::Ellipse, p)
    H = zeros(C64, 4, 4)
    inv = 1.0 / E.S2
    H[1, 1] = -E.S1 * inv
    H[1, 2] = -(E.S0 - p) * inv
    H[1, 3] = -conj(E.S1) * inv
    H[1, 4] = -1.0
    H[2, 1] = H[3, 2] = H[4, 3] = 1.0
    for _ in 1:QR_BALANCE                  # diagonal similarity
        for i in 1:4
            col = row = 0.0
            for j in 1:4
                j == i && continue
                col += abs(H[j, i])
                row += abs(H[i, j])
            end
            f = (col > 0.0 && row > 0.0) ? sqrt(row / col) : 1.0
            for j in 1:4
                H[j, i] *= f
                H[i, j] /= f
            end
        end
    end
    for _ in 1:QR_SWEEPS4
        qr_sweep!(H, 4)
    end
    for _ in 1:QR_SWEEPS3
        qr_sweep!(H, 3)
    end
    a, b, c, d = H[1, 1], H[1, 2], H[2, 1], H[2, 2]
    h = 0.5 * (a - d)
    disc = sqrt(h * h + b * c)
    mid = 0.5 * (a + d)
    return (mid + disc, mid - disc, H[3, 3], H[4, 4])
end

# the four roots of D at level p: eigenvalue seeds, then Aberth iterations
# (at most N_ABERTH, stopping once every correction is at rounding) and two
# Newton polishes; a correction that is not finite is skipped
function roots(E::Ellipse, p)
    z = qr_seeds(E, p)
    for _ in 1:N_ABERTH
        rat = map(zi -> Dz(E, zi, p) / dDz(E, zi, p), z)
        done = true
        nz = z
        for i in 1:4
            s = zero(C64)
            for j in 1:4
                if j != i && z[i] != z[j]
                    s += 1.0 / (z[i] - z[j])
                end
            end
            step = rat[i] / (1.0 - rat[i] * s)
            isfinite(step) || continue
            abs(step) > 4.0 * EPS * abs(z[i]) && (done = false)
            nz = Base.setindex(nz, z[i] - step, i)
        end
        z = nz
        done && break
    end
    for _ in 1:2
        z = map(z) do zi
            nw = zi - Dz(E, zi, p) / dDz(E, zi, p)
            isfinite(nw) ? nw : zi
        end
    end
    return z
end

# ---- limb crossings and the arcs ------------------------------------------

wrap2pi(t) = t - TWO_PI * floor(t / TWO_PI)

function g_dg(phi, x0, y0, a, b)
    s, c = sincos(phi)
    x, y = x0 + a * c, y0 + b * s
    return x * x + y * y - 1.0, 2.0 * (-x * a * s + y * b * c)
end

# number of limb crossings, and the first two (phi_a <= phi_b) on the ellipse
# parameter; crossings closer than PAIR_TOL are a tangency and are dropped
# together
function crossings(x0, y0, a, b)
    z = roots(Ellipse(x0, y0, a, b), 1.0)
    sc = 1.0 + (abs(x0) + a)^2 + (abs(y0) + b)^2
    phi = zeros(4)
    valid = falses(4)
    for i in 1:4
        cand = abs(abs(z[i]) - 1.0) < Z_TOL
        ph = cand ? angle(z[i]) : 0.0
        g, dg = g_dg(ph, x0, y0, a, b)
        for _ in 1:3
            trial = ph - (abs(dg) > 0.0 ? g / dg : 0.0)
            gt, dgt = g_dg(trial, x0, y0, a, b)
            if abs(gt) < abs(g)
                ph, g, dg = trial, gt, dgt
            end
        end
        valid[i] = cand && abs(g) <= G_TOL * sc
        phi[i] = wrap2pi(ph)
    end
    drop = falses(4)
    for i in 1:4, j in i+1:4
        d = abs(wrap2pi(phi[i] - phi[j] + π) - π)
        if valid[i] && valid[j] && d < PAIR_TOL
            drop[i] = drop[j] = true
        end
    end
    keep = sort([phi[i] for i in 1:4 if valid[i] && !drop[i]])
    n = length(keep)
    return n, (n > 0 ? keep[1] : 0.0), (n > 1 ? keep[2] : 0.0)
end

# the ellipse arc [lo, hi] inside the star and the length dth of the limb arc
# inside the ellipse, for two crossings. The ellipse arc is decided at the
# midpoint of the longer of its two arcs; the limb arc inside the ellipse is
# always the shorter one (the planet lies within radius A < 1).
function arcs(x0, y0, a, b, phi_a, phi_b)
    long_first = (phi_b - phi_a) >= π
    mid = long_first ? 0.5 * (phi_a + phi_b) : 0.5 * (phi_a + phi_b) + π
    xm, ym = x0 + a * cos(mid), y0 + b * sin(mid)
    mid_in = xm * xm + ym * ym < 1.0
    first = long_first ? mid_in : !mid_in       # is [phi_a, phi_b] inside?
    lo, hi = first ? (phi_a, phi_b) : (phi_b, phi_a + TWO_PI)
    ta = atan(y0 + b * sin(phi_a), x0 + a * cos(phi_a))
    tb = atan(y0 + b * sin(phi_b), x0 + a * cos(phi_b))
    dt = abs(wrap2pi(ta) - wrap2pi(tb))
    return lo, hi, min(dt, TWO_PI - dt)
end

# ---- even terms: Fourier coefficients of h_n(s) K ---------------------------

# q[n+1][m+1], n = 0..2, m = 0..5: h_n(s) K = Re Σ_m q_m e^{imφ}, with
# h_n(s) = Σ_p α[n, p] s^p, α[n, p] = (-1)^p C(n, p) / (2 (p + 1))
function q_coeffs(E::Ellipse)
    alpha = ((0.5, 0.0, 0.0), (0.5, -0.25, 0.0), (0.5, -0.5, 1.0 / 6.0))
    T = zeros(C64, 11)                  # harmonics -5..5 at index m + 6
    T[5], T[6], T[7] = conj(E.K1), E.K0, E.K1
    Sk = (C64(E.S2), conj(E.S1), C64(E.S0), E.S1, C64(E.S2))   # j = -2..2
    q = [zeros(C64, 6) for _ in 1:3]
    for p in 0:2
        if p > 0                        # T <- S * T
            Tn = zeros(C64, 11)
            for i in 1:11, j in -2:2
                1 <= i - j <= 11 && (Tn[i] += Sk[j+3] * T[i-j])
            end
            T = Tn
        end
        for n in p:2, m in 0:5
            q[n+1][m+1] += alpha[n+1][p+1] * T[m+6]
        end
    end
    return q
end

# Ψ_n(hi) - Ψ_n(lo), Ψ_n the exact antiderivative, in difference form
function even_arc(qn, lo, hi)
    mid, half = 0.5 * (lo + hi), 0.5 * (hi - lo)
    acc = real(qn[1]) * (hi - lo)
    for m in 1:5
        acc += (4.0 / m) * real(qn[m+1] * cis(m * mid)) * sin(m * half)
    end
    return acc
end

# ---- poles: partial fractions over the roots of D ---------------------------

# Newton steps on (σ, π) for the quadratic factor z² - σ z + π
function bairstow(E::Ellipse, p, sig::C64, pq::C64)
    dc = (C64(E.S2), conj(E.S1), C64(E.S0 - p), E.S1, C64(E.S2))   # low -> high
    for _ in 1:N_BAIRSTOW
        al, be = zero(C64), one(C64)    # z^k mod q = al z + be
        das = dbs = dap = dbp = zero(C64)
        r1, r0 = dc[1] * al, dc[1] * be
        J11 = J12 = J21 = J22 = zero(C64)
        for k in 2:5
            al, be, das, dbs, dap, dbp = (al * sig + be, -al * pq,
                                          das * sig + al + dbs, -das * pq,
                                          dap * sig + dbp, -dap * pq - al)
            r1 += dc[k] * al
            r0 += dc[k] * be
            J11 += dc[k] * das
            J12 += dc[k] * dap
            J21 += dc[k] * dbs
            J22 += dc[k] * dbp
        end
        det = J11 * J22 - J12 * J21
        abs(det) > 0.0 || break
        sig -= (J22 * r1 - J12 * r0) / det
        pq -= (-J21 * r1 + J11 * r0) / det
    end
    return sig, pq
end

# for the quadratic factor q of D, with N = aN z + bN and Q = D/q = aQ z + bQ
# (mod q): the sum of N/D' over its two roots (A_sum) and the C of
# ∫ (A z + B)/q dz (C_num, before the factor i); symmetric in the roots
function pair_terms(E::Ellipse, p, sig, pq)
    Q2 = C64(E.S2)
    Q1 = E.S1 + sig * Q2
    Q0 = sig * Q1 + (E.S0 - p) - pq * Q2
    aN, bN = E.K1 * sig + E.K0, conj(E.K1) - E.K1 * pq
    aQ, bQ = Q2 * sig + Q1, Q0 - Q2 * pq
    den = aQ * aQ * pq + aQ * bQ * sig + bQ * bQ
    abs(den) > 0.0 || (den = one(C64))
    A_sum = (aN * bQ - bN * aQ) / den
    C_num = 0.5 * (2.0 * aN * aQ * pq + (aN * bQ + bN * aQ) * sig + 2.0 * bN * bQ) / den
    return A_sum, C_num
end

# a whole ellipse on the disc: -(π/p) Re Σ_{|z_r| < 1} N/D'(z_r), the two
# inner roots summed as one quadratic factor
function pole_inside(E::Ellipse, p, z)
    sig, pq = zero(C64), one(C64)
    for zi in z
        if abs(zi) < 1.0
            sig += zi
            pq *= zi
        end
    end
    sig, pq = bairstow(E, p, sig, pq)
    A_sum, _ = pair_terms(E, p, sig, pq)
    return -π * real(A_sum) / p
end

# change of log(e^{iφ} - z) for φ from lo to hi, continuous on the arc
function dlog(z, lo, hi, w1, w2)
    dm = log(abs(w2 - z)) - log(abs(w1 - z))
    da = abs(z) < 1.0 ? (hi - lo) + angle(1.0 - z / w2) - angle(1.0 - z / w1) :
                        angle(1.0 - w2 / z) - angle(1.0 - w1 / z)
    return C64(dm, da)
end

# ∫ du / (u² - Δ) = -(1/u) Σ_k (Δ/u²)^k / (2k + 1)
function Jser(u, Delta)
    x = Delta / (u * u)
    acc = zero(C64)
    for k in N_SERIES-1:-1:0
        acc = acc * x + 1.0 / (2k + 1)
    end
    return -acc / u
end

# (1/2p) ∫_lo^hi K / (p - s) dφ along one ellipse arc inside the star:
# (1/2p) Re Σ_r R_r [log(e^{iφ} - z_r)], R_r = -N(z_r) / (i D'(z_r))
function pole_arc(E::Ellipse, p, z, lo, hi)
    w1, w2 = cis(lo), cis(hi)
    # near-double pairs: up to two disjoint, scored by |Δ| / d², d the
    # distance from the pair's centre to the arc
    used = falses(4)
    pairs = Tuple{Int,Int}[]
    for _ in 1:2
        best, bi = Inf, 1
        for (m, (j, l)) in enumerate(PAIRS)
            (used[j] || used[l]) && continue
            Dl = 0.25 * (z[j] - z[l])^2
            c = 0.5 * (z[j] + z[l])
            ang = wrap2pi(angle(c))
            on = (lo <= ang <= hi) || (lo <= ang + TWO_PI <= hi)
            d = on ? abs(abs(c) - 1.0) : min(abs(w1 - c), abs(w2 - c))
            score = abs(Dl) / max(d * d, 1e-300)
            if score < best
                best, bi = score, m
            end
        end
        if best <= PAIR_RATIO
            j, l = PAIRS[bi]
            used[j] = used[l] = true
            push!(pairs, (j, l))
        end
    end

    ell = zero(C64)
    # isolated roots; one beyond |z| = 2 is taken through its partner
    # ζ = 1/conj(z) inside the disc, where R(z) = conj(R(ζ)) and the log's
    # change is that of log(1 - e^{iφ} conj(ζ))
    for k in 1:4
        used[k] && continue
        a2 = abs2(z[k])
        far = a2 > 4.0
        zeta = far ? z[k] / a2 : z[k]
        nw = zeta - Dz(E, zeta, p) / dDz(E, zeta, p)
        if isfinite(nw) && abs(Dz(E, nw, p)) < abs(Dz(E, zeta, p))
            zeta = nw
        end
        dD = dDz(E, zeta, p)
        abs(dD) > 0.0 || continue
        R = -Nz(E, zeta) / (im * dD)
        if far
            b2, b1 = 1.0 - w2 * conj(zeta), 1.0 - w1 * conj(zeta)
            ell += conj(R) * C64(log(abs(b2)) - log(abs(b1)), angle(b2) - angle(b1))
        else
            ell += R * dlog(zeta, lo, hi, w1, w2)
        end
    end
    # near-double pairs, each as one quadratic factor q = z² - σ z + π:
    # ∫ (A z + B)/q dz = (A/2) log q + C J(z - σ/2), the winding of log q
    # taken from the two roots' own continued logs
    for (j, l) in pairs
        z1, z2 = z[j], z[l]
        sig, pq = bairstow(E, p, z1 + z2, z1 * z2)
        A_sum, C_num = pair_terms(E, p, sig, pq)
        qh = w2 * w2 - sig * w2 + pq
        ql = w1 * w1 - sig * w1 + pq
        darg = angle(qh / ql)
        rough = imag(dlog(z1, lo, hi, w1, w2)) + imag(dlog(z2, lo, hi, w1, w2))
        wind = round((rough - darg) / TWO_PI)
        LQ = C64(log(abs(qh)) - log(abs(ql)), darg + TWO_PI * wind)
        Delta = 0.25 * sig * sig - pq
        cen = 0.5 * sig
        Jd = Jser(w2 - cen, Delta) - Jser(w1 - cen, Delta)
        ell += 0.5 * (im * A_sum) * LQ + (im * C_num) * Jd
    end
    return real(ell) / (2.0 * p)
end

# occulted flux of [μ⁰, μ², μ⁴, P_ε...] for the ellipse with centre (x0, y0)
# in its principal frame and semi-axes a > b
function ellipse_occ(x0, y0, a, b, L::Law)
    zero6 = ntuple(_ -> 0.0, 6)
    d = hypot(x0, y0)
    fast_in = d + a < 1.0
    general = !fast_in && (d - a < 1.0)
    !fast_in && !general && return zero6            # disjoint
    cnt, phi_a, phi_b = general ? crossings(x0, y0, a, b) : (0, 0.0, 0.0)
    partial = general && cnt == 2
    inside = fast_in || (general && cnt != 2 && d < 1.0)
    !partial && !inside && return zero6

    E = Ellipse(x0, y0, a, b)
    q = q_coeffs(E)
    if inside
        poles = ntuple(k -> k > L.K ? 0.0 : pole_inside(E, 1.0 + L.eps[k],
                                                        roots(E, 1.0 + L.eps[k])), 3)
        return (TWO_PI * real(q[1][1]), TWO_PI * real(q[2][1]), TWO_PI * real(q[3][1]), poles...)
    end
    lo, hi, dth = arcs(x0, y0, a, b, phi_a, phi_b)
    evens = ntuple(n -> even_arc(q[n], lo, hi) + dth / (2.0 * n), 3)
    poles = ntuple(3) do k
        k > L.K && return 0.0
        e = L.eps[k]
        p = 1.0 + e
        pole_arc(E, p, roots(E, p), lo, hi) + dth / (2.0 * p * e)
    end
    return (evens..., poles...)
end

# ===========================================================================
# public interface
# ===========================================================================

"""
    hybrid_flux(z, r, w, law=:hybrid4)

Normalised flux of a star with a hybrid law, occulted by a spherical planet
of radius `r` (stellar radii) at sky separation `z` (a number or an array).
`w` are the law's shape weights; `law` is `:hybrid2`, `:hybrid4` or
`:hybrid5` (or `HYBRID2`, ...). Exact for r < 1; for a companion larger than
the star the even terms lose about r⁵ × 1e-16.
"""
function hybrid_flux(z::Real, r::Real, w, law=HYBRID4)
    L = getlaw(law)
    checkweights(w, L)
    r >= 0 || throw(DomainError(r, "need r >= 0"))
    return flux_from_occ(sphere_occ(Float64(z), Float64(r), L), w, L)
end

hybrid_flux(z::AbstractArray, r::Real, w, law=HYBRID4) =
    map(zi -> hybrid_flux(zi, r, w, law), z)

"""
    hybrid_flux_oblate(X, Y, r, f, theta, w, law=:hybrid4)

Normalised flux of a star with a hybrid law, occulted by an oblate planet:
sky position `(X, Y)` (numbers or arrays of one shape), X along its
direction of motion; area-equivalent radius `r` = √(AB), with
A = r/√(1 - f) and B = A(1 - f); flattening `f` in [0, 1); long axis at
angle `theta` (radians) from X. Requires r ≤ (1 - f)^(3/2).
"""
function hybrid_flux_oblate(X::Real, Y::Real, r::Real, f::Real, theta::Real, w,
                            law=HYBRID4)
    L = getlaw(law)
    checkweights(w, L)
    0 <= f < 1 || throw(DomainError(f, "need 0 <= f < 1"))
    0 <= r <= (1 - f)^1.5 * (1 + 1e-12) ||
        throw(DomainError(r, "need 0 <= r <= (1 - f)^(3/2): the planet must " *
                             "cross the stellar limb at most twice"))
    s, c = sincos(Float64(theta))
    x0, y0 = X * c + Y * s, -X * s + Y * c        # principal frame
    occ = if f < F_SW
        sphere_occ(hypot(x0, y0), Float64(r), L)
    else
        a = r / sqrt(1.0 - f)
        ellipse_occ(x0, y0, a, a * (1.0 - f), L)
    end
    return flux_from_occ(occ, w, L)
end

function hybrid_flux_oblate(X::AbstractArray, Y::AbstractArray, r::Real, f::Real,
                            theta::Real, w, law=HYBRID4)
    size(X) == size(Y) || throw(DimensionMismatch("X and Y differ in size"))
    return map((x, y) -> hybrid_flux_oblate(x, y, r, f, theta, w, law), X, Y)
end

end # module
