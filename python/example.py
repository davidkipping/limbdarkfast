"""Example: a spherical and an oblate planet on a circular orbit.

limbdarkfast computes the flux at given sky positions; the orbit is up to the
caller. On a circular orbit with orbital phase phi = 2 pi (t - t0) / P,
semi-major axis a and impact parameter b (both in stellar radii), the
planet's sky position, X along its direction of motion, is

    X = a sin(phi),   Y = -b cos(phi),   in front of the star when cos(phi) > 0.
"""

import numpy as np

import limbdarkfast

P, t0, a, b = 3.0, 0.0, 10.0, 0.3          # days, days, R*, R*
r = 0.1                                    # radius (area-equivalent if oblate)
f, theta = 0.2, np.radians(30.0)           # flattening; long axis from X
law, w = "hybrid4", [0.30, 0.20, 0.10, 0.05]

t = np.linspace(-0.1, 0.1, 2001)
phi = 2.0 * np.pi * (t - t0) / P
X, Y = a * np.sin(phi), -b * np.cos(phi)
front = np.cos(phi) > 0.0

F_sphere = np.where(front, limbdarkfast.hybrid_flux(np.hypot(X, Y), r, w, law), 1.0)
F_oblate = np.where(front, limbdarkfast.hybrid_flux_oblate(X, Y, r, f, theta, w, law),
                    1.0)

print(f"{law}, w = {w}")
print(f"depth at mid-transit: sphere {1 - F_sphere[1000]:.6e}, "
      f"oblate {1 - F_oblate[1000]:.6e}")
print(f"largest oblate - sphere difference: "
      f"{np.max(np.abs(F_oblate - F_sphere)) * 1e6:.1f} ppm")

try:
    import matplotlib.pyplot as plt
except ImportError:
    pass
else:
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 5))
    ax1.plot(t * 24, F_sphere, label="spherical")
    ax1.plot(t * 24, F_oblate, "--", label=f"oblate, f = {f}")
    ax1.set_ylabel("flux")
    ax1.legend()
    ax2.plot(t * 24, (F_oblate - F_sphere) * 1e6)
    ax2.set_xlabel("time from mid-transit [h]")
    ax2.set_ylabel("oblate - sphere [ppm]")
    fig.tight_layout()
    fig.savefig("example.png", dpi=150)
    print("wrote example.png")
