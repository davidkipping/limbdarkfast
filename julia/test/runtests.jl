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
