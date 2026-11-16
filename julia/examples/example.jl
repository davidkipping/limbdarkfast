# Example: a spherical and an oblate planet on a circular orbit.
#
# LimbDarkFast computes the flux at given sky positions; the orbit is up to the
# caller. On a circular orbit with orbital phase phi = 2 pi (t - t0) / P,
# semi-major axis a and impact parameter b (both in stellar radii), the
# planet's sky position, X along its direction of motion, is
#
#     X = a sin(phi),   Y = -b cos(phi),   in front of the star when cos(phi) > 0.
#
# Run from julia/:  julia --project=. examples/example.jl
using LimbDarkFast

P, t0, a, b = 3.0, 0.0, 10.0, 0.3          # days, days, R*, R*
r = 0.1                                    # radius (area-equivalent if oblate)
f, theta = 0.2, deg2rad(30.0)              # flattening; long axis from X
law, w = :hybrid4, [0.30, 0.20, 0.10, 0.05]

t = range(-0.1, 0.1, length=2001)
phi = @. 2pi * (t - t0) / P
X, Y = a .* sin.(phi), -b .* cos.(phi)
front = cos.(phi) .> 0

F_sphere = ifelse.(front, hybrid_flux(hypot.(X, Y), r, w, law), 1.0)
F_oblate = ifelse.(front, hybrid_flux_oblate(X, Y, r, f, theta, w, law), 1.0)

println("$law, w = $w")
println("depth at mid-transit: sphere ", 1 - F_sphere[1001], ", oblate ", 1 - F_oblate[1001])
println("largest oblate - sphere difference: ",
        round(maximum(abs.(F_oblate .- F_sphere)) * 1e6; digits=1), " ppm")
