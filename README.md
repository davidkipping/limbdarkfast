# limbdarkfast

Fast transit fluxes for the hybrid limb-darkening laws, in the languages
astronomers use, in the spirit of the Agol, Luger & Foreman-Mackey (2020)
code release. Each version computes one thing: the normalised flux of a star with
**hybrid2, hybrid4 or hybrid5** limb darkening, occulted by a **spherical**
or an **oblate** planet, in closed form and double precision. There are
no orbits, exposure integration or derivatives; the caller supplies the
planet's sky position.

If you use this code, please cite Kipping (2027); see [Citation](#citation).

## Interface (the same in every language)

```
F = hybrid_flux(z, r, w, law)                        spherical planet
F = hybrid_flux_oblate(X, Y, r, f, theta, w, law)    oblate planet
```

| argument | meaning |
|---|---|
| `z` | sky separation of planet and star centres, in stellar radii |
| `X, Y` | the planet's sky position, X along its direction of motion |
| `r` | planet radius; for an oblate planet the area-equivalent radius √(AB) |
| `f` | projected flattening 1 − B/A, in [0, 1) |
| `theta` | angle of the long axis from X, **radians** |
| `w` | the law's shape weights (2, 4 or 5) |
| `law` | `hybrid2`, `hybrid4` or `hybrid5` |

The semi-axes are A = r/√(1−f) and B = A(1−f). The oblate planet must
satisfy r ≤ (1−f)^{3/2}, so that its outline crosses the limb at most
twice (0.59 at f = 0.3, 0.35 at f = 0.5).

On a circular orbit with phase φ = 2π(t − t₀)/P, semi-major axis a and
impact parameter b, `X = a sin φ`, `Y = −b cos φ`, with the planet in front
of the star when cos φ > 0. Each version's example does this.

## The laws

I(μ)/I(1) = 1 − Σⱼ wⱼ Sⱼ(μ), each shape falling from 1 at the limb to 0 at
the centre, built from μ⁰, μ², μ⁴ and double poles (μ² + ε)⁻²:

| law | weights | shapes | poles ε |
|---|---:|---|---|
| hybrid2 | 2 | 1 − μ², Π_ε | 0.36 |
| hybrid4 | 4 | curvature basis S₁…S₄ | 0.0286, 0.3858 |
| hybrid5 | 5 | 1 − μ⁴, (1 − μ²)², Π_ε₀, Π_ε₁, Π_ε₂ | 0.0016, 0.0413, 0.3874 |

with Π_ε = [(μ²+ε)⁻² − (1+ε)⁻²] / [ε⁻² − (1+ε)⁻²]. Every version holds
the same constants, generated from the Python one by
`tools/law_constants.py`.

## Layout and status

| directory | contents | status |
|---|---|---|
| `python/` | `limbdarkfast.py` (NumPy), `test_limbdarkfast.py`, `example.py` | done: vectors at ≤ 2e-14 |
| `cpp/` | header-only `limbdarkfast.hpp` (C++11), test, example, Makefile | done: vectors at ≤ 2e-14 |
| `fortran/` | module `limbdarkfast.f90` (Fortran 2008), test, example, Makefile | done: vectors at ≤ 2e-14 |
| `julia/` | package `LimbDarkFast` (`src/LimbDarkFast.jl`), `test/runtests.jl`, `examples/` | done: vectors at ≤ 6e-14 |
| `idl/` | `limbdarkfast_hybrid_flux.pro`, `limbdarkfast_hybrid_flux_oblate.pro`, `limbdarkfast_intensity.pro`, `limbdarkfast_core.pro`, test, example | done in GDL: vectors at ≤ 2e-14 |
| `tests/` | `vectors/*.txt`, the shared reference fluxes | done |

## Using each version

**Python** (NumPy only; arguments broadcast as arrays):

```python
import limbdarkfast
F = limbdarkfast.hybrid_flux(z, 0.1, [0.3, 0.2, 0.1, 0.05], "hybrid4")
F = limbdarkfast.hybrid_flux_oblate(X, Y, 0.1, 0.2, 0.52, [0.3, 0.2, 0.1, 0.05], "hybrid4")
```

**C++** (one header, C++11, no dependencies; one point per call, or a
`std::vector` of positions; out-of-domain input throws `std::domain_error`):

```cpp
#include "limbdarkfast.hpp"
std::vector<double> w = {0.3, 0.2, 0.1, 0.05};
double F  = limbdarkfast::hybrid_flux(z, 0.1, w, limbdarkfast::Law::hybrid4);
double Fo = limbdarkfast::hybrid_flux_oblate(X, Y, 0.1, 0.2, 0.52, w, limbdarkfast::Law::hybrid4);
```

`make test` in `cpp/` builds and runs the checks. A spherical point costs
about 30 ns and an oblate one 8–11 µs (Apple M2 Max, `-O2`); adding
`-fcx-limited-range` saves 10–30% with clang and is safe here, since no
complex number in the calculation comes near overflow.

**Fortran** (one module, Fortran 2008, no dependencies; `z`, or `X` and
`Y`, scalars or rank-1 arrays; the law codes `HYBRID2`, `HYBRID4`,
`HYBRID5` equal their number of weights; out-of-domain input or a wrong
number of weights returns NaN, so the functions stay `pure`):

```fortran
use limbdarkfast
real(dp) :: w(4) = [0.3d0, 0.2d0, 0.1d0, 0.05d0]
F  = hybrid_flux(z, 0.1d0, w, HYBRID4)
Fo = hybrid_flux_oblate(X, Y, 0.1d0, 0.2d0, 0.52d0, w, HYBRID4)
```

`make test` in `fortran/` builds and runs the checks (gfortran; set `FC` for
another compiler). About 50 ns per spherical point and 8 µs per oblate one.

**Julia** (a package with no dependencies; `z`, or `X` and `Y`, numbers or
arrays; the law as `:hybrid4`, `"hybrid4"` or `HYBRID4`; bad input throws
`DomainError` or `ArgumentError`):

```julia
using Pkg; Pkg.develop(path="path/to/limbdarkfast/julia")
using LimbDarkFast
w = [0.3, 0.2, 0.1, 0.05]
F  = hybrid_flux(z, 0.1, w, :hybrid4)
Fo = hybrid_flux_oblate(X, Y, 0.1, 0.2, 0.52, w, :hybrid4)
```

`julia --project=. -e 'using Pkg; Pkg.test()'` in `julia/` runs the checks.
About 50 ns per spherical point and 13 µs per oblate one.

**IDL** (vectorised like EXOFAST's `occultquad`; `z`, or `X` and `Y`,
scalars or arrays; the law as a string; bad input stops with `MESSAGE`). Put
`idl/` on `!PATH`:

```idl
w  = [0.3d, 0.2d, 0.1d, 0.05d]
F  = limbdarkfast_hybrid_flux(z, 0.1d, w, 'hybrid4')
Fo = limbdarkfast_hybrid_flux_oblate(X, Y, 0.1d, 0.2d, 0.52d, w, 'hybrid4')
```

`test_limbdarkfast` (from `idl/`) runs the checks. Tested in GDL (the GNU Data
Language, weekly arm64 build of November 2026), not in IDL itself; about
2 µs per spherical point and 28 µs per oblate one there. That build has two
bugs the code works around, harmlessly for IDL: arithmetic between
one-element arrays is wrong (`r = a / b` leaves the quotient in `a`), so no
internal array ever has exactly one element; and
`MIN(..., DIMENSION=2)` returns subscript 0 for a row that is all `Infinity`,
so the code uses finite sentinels.

## Validation

`tests/vectors/` holds reference fluxes from an independent float64
implementation of the same closed forms: for each law, 2,226 spherical points (every regime, and the contacts, the line
z = r and each pole's switch line approached to 1e-13) and 5,049 oblate
points (chords at flattenings from 0 to 0.7, the edge of the domain, the
outline tangent to the limb and to every pole circle approached to 1e-13,
and centred, nearly round planets). The oblate references agree with
SquishierPlanet, a third implementation, to 2.4e-15. Every version must
reproduce them to 1e-12.

The files are plain text: `#` header lines, then rows of
`z r w… F` (spherical) or `X Y r f theta w… F` (oblate).

Accuracy: rounding (≲ 1e-15) for planets. For a companion larger than the
star the spherical even terms lose about r⁵ × 1e-16 (2e-14 at r = 3).

## Citation

If you use this code, please cite

> Kipping, D. 2027, "Limb Darkened Transits of Oblate Planets as Elementary
> Functions", AJ, submitted.

```bibtex
@article{kipping2027,
  author  = {Kipping, David},
  title   = {Limb Darkened Transits of Oblate Planets as Elementary Functions},
  journal = {AJ},
  year    = {2027},
  note    = {submitted}
}
```

## License

MIT; see [LICENSE](LICENSE).
