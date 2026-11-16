# Check the Julia version against the shared test vectors:
#     julia --project=. -e 'using Pkg; Pkg.test()'      (or: julia --project=. test/runtests.jl)
using LimbDarkFast
using Test

const VECTORS = joinpath(@__DIR__, "..", "..", "tests", "vectors")
const TOL = 1e-12

function read_rows(path)
    rows = Vector{Float64}[]
    for line in eachline(path)
        (isempty(line) || startswith(line, "#")) && continue
        push!(rows, parse.(Float64, split(line)))
    end
    return rows
end

function max_error(kind, law)
    nw = nweights(law)
    rows = read_rows(joinpath(VECTORS, "$(kind)_$(law).txt"))
    worst = 0.0
    for d in rows
        F = kind == "sphere" ? hybrid_flux(d[1], d[2], d[3:2+nw], law) :
                               hybrid_flux_oblate(d[1], d[2], d[3], d[4], d[5], d[6:5+nw], law)
        err = abs(F - d[end])
        err <= worst || (worst = err)          # NaN counts as worst
    end
    return length(rows), worst
end

# ---------------------------------------------------------------------------
# brute-force integration, as in Agol, Luger & Foreman-Mackey (2020)
# (test/test_transit_poly.jl of Limbdark.jl): sum the intensity over NS
# annuli of the star, each weighted by the angle Phi(s) it spends inside
# the planet, with the midpoint rule
# ---------------------------------------------------------------------------

const NS = 5000                 # annuli
const M_ARC = 48                # samples in angle before bisecting for crossings
const NUM_TOL = 1e-6
const NUM_W = Dict(:hybrid2 => [0.12123702682335066, 0.22967167533748584],
                   :hybrid4 => [0.3709047608728925, 0.00711438176464726,
                                0.3162293789763985, 0.01888876666845794],
                   :hybrid5 => [0.136725258105043, 0.25114827117024674, 0.1088217352833381,
                                0.030866335871229982, 0.32786772243688317])
# their table of (r, b), every contact included, with r <= 2
const NUM_RB = [(0.01, 0.0), (0.01, 0.01), (0.01, 0.99), (0.01, 1.0), (0.01, 1.01),
                (0.1, 0.0), (0.1, 0.1), (0.1, 0.9), (0.1, 1.0), (0.1, 1.1),
                (1.0, 0.0), (1.0, 1.0), (1.0, 2.0), (2.0, 1.0), (2.0, 2.0), (2.0, 3.0)]
# the oblate table: the same contacts for r <= 0.5, at f = 0.2, the planet
# displaced along the direction 0.6 rad and tilted by 0.4 rad
const NUM_OBL = [filter(rb -> rb[1] < 1.0, NUM_RB);
                 [(0.5, 0.0), (0.5, 0.5), (0.5, 1.0), (0.5, 1.5)]]
const NUM_F, NUM_BETA, NUM_THETA = 0.2, 0.6, 0.4

num_total(w, law) =
    sum(2π * s / NS * intensity(sqrt(1 - s * s), w, law) for s in ((0:NS-1) .+ 0.5) ./ NS)

"Spherical planet, by brute force."
function numeric_flux(z, r, w, law)
    z >= 1 + r && return 1.0
    s1, s2 = max(0.0, z - r), min(1.0, z + r)
    ds = (s2 - s1) / NS
    occ = 0.0
    for j in 0:NS-1
        s = s1 + (j + 0.5) * ds
        phi = s < r - z ? 2π : 2 * acos(clamp((s * s + z * z - r * r) / (2 * z * s + 1e-300), -1, 1))
        occ += s * ds * phi * intensity(sqrt(1 - s * s), w, law)
    end
    return 1 - occ / num_total(w, law)
end

"""
The angle of the circle of radius s inside the ellipse
((x - x0)/A)² + ((y - y0)/B)² ≤ 1: sample M_ARC + 1 angles across the window
the ellipse can occupy, bisect each sign change, and add up the parts of the
intervals inside.
"""
function arc_inside(s, x0, y0, A, B)
    d = hypot(x0, y0)
    lo, hi = 0.0, 2π
    if d > A
        half, c = asin(A / d), atan(y0, x0)
        lo, hi = c - half, c + half
    end
    inside(p) = ((s * cos(p) - x0) / A)^2 + ((s * sin(p) - y0) / B)^2 < 1
    h = (hi - lo) / M_ARC
    total = 0.0
    in_a = inside(lo)
    for k in 0:M_ARC-1
        pa, pb = lo + h * k, lo + h * (k + 1)
        in_b = inside(pb)
        if in_a && in_b
            total += h
        elseif in_a != in_b
            a, b = pa, pb
            for _ in 1:50
                m = 0.5 * (a + b)
                inside(m) == in_a ? (a = m) : (b = m)
            end
            cross = 0.5 * (a + b)
            total += in_a ? cross - pa : pb - cross
        end
        in_a = in_b
    end
    return total
end

"Oblate planet, by brute force."
function numeric_flux_oblate(X, Y, r, f, theta, w, law)
    A = r / sqrt(1 - f)
    B = A * (1 - f)
    x0 = X * cos(theta) + Y * sin(theta)
    y0 = -X * sin(theta) + Y * cos(theta)
    d = hypot(x0, y0)
    s1, s2 = max(0.0, d - A), min(1.0, d + A)
    s1 >= s2 && return 1.0
    ds = (s2 - s1) / NS
    occ = 0.0
    for j in 0:NS-1
        s = s1 + (j + 0.5) * ds
        occ += s * ds * arc_inside(s, x0, y0, A, B) * intensity(sqrt(1 - s * s), w, law)
    end
    return 1 - occ / num_total(w, law)
end

"""
Largest difference from brute force: the table of contacts, then a light
curve (1000 points spherical, 201 oblate) for r = 0.1.
"""
function numeric_checks()
    worst = 0.0
    for law in (:hybrid2, :hybrid4, :hybrid5)
        w = NUM_W[law]
        for (r, b) in NUM_RB
            worst = max(worst, abs(hybrid_flux(b, r, w, law) - numeric_flux(b, r, w, law)))
        end
        for (r, b) in NUM_OBL
            X, Y = b * cos(NUM_BETA), b * sin(NUM_BETA)
            worst = max(worst, abs(hybrid_flux_oblate(X, Y, r, NUM_F, NUM_THETA, w, law) -
                                   numeric_flux_oblate(X, Y, r, NUM_F, NUM_THETA, w, law)))
        end
    end
    w, law = NUM_W[:hybrid4], :hybrid4
    for i in 1:1000
        b = abs(3.0 * (i / 1000.0 - 0.5))
        worst = max(worst, abs(hybrid_flux(b, 0.1, w, law) - numeric_flux(b, 0.1, w, law)))
    end
    for i in 0:200
        X = -1.5 + 3.0 * i / 200.0
        worst = max(worst, abs(hybrid_flux_oblate(X, 0.3, 0.1, NUM_F, NUM_THETA, w, law) -
                               numeric_flux_oblate(X, 0.3, 0.1, NUM_F, NUM_THETA, w, law)))
    end
    return worst
end

@testset "LimbDarkFast" begin
    @testset "vectors" begin
        for kind in ("sphere", "oblate"), law in (:hybrid2, :hybrid4, :hybrid5)
            n, err = max_error(kind, law)
            println(rpad(kind, 7), law, ": ", lpad(n, 5), " points, max |dF| = ",
                    round(err; sigdigits=3))
            @test err <= TOL
        end
    end
    @testset "law shapes" begin
        for law in (HYBRID2, HYBRID4, HYBRID5), j in 1:nweights(law)
            w = zeros(nweights(law)); w[j] = 1.0
            @test abs(intensity(0.0, w, law)) < 1e-12        # 1 at the limb
            @test abs(intensity(1.0, w, law) - 1.0) < 1e-12  # 0 at the centre
        end
    end
    @testset "brute-force integration" begin
        worst = numeric_checks()
        println("brute-force integration: max |dF| = ", round(worst; sigdigits=3))
        @test worst <= NUM_TOL
    end
    @testset "interface" begin
        w = [0.3, 0.2, 0.1, 0.05]
        z = range(0, 1.2, length=7)
        @test hybrid_flux(z, 0.1, w, :hybrid4) == [hybrid_flux(zi, 0.1, w, :hybrid4) for zi in z]
        @test hybrid_flux(0.3, 0.1, w, "hybrid4") == hybrid_flux(0.3, 0.1, w, HYBRID4)
        @test hybrid_flux_oblate(0.3, 0.2, 0.1, 0.0, 0.5, w) ≈ hybrid_flux(hypot(0.3, 0.2), 0.1, w)
        @test_throws ArgumentError hybrid_flux(0.3, 0.1, [0.1, 0.2], :hybrid4)
        @test_throws ArgumentError hybrid_flux(0.3, 0.1, w, :quadratic)
        @test_throws DomainError hybrid_flux_oblate(0.0, 0.0, 0.5, 0.5, 0.0, w)
    end
end
