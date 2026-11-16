!> Example: a spherical and an oblate planet on a circular orbit.
!>
!> limbdarkfast computes the flux at given sky positions; the orbit is up to the
!> caller. On a circular orbit with orbital phase phi = 2 pi (t - t0) / P,
!> semi-major axis a and impact parameter b (both in stellar radii), the
!> planet's sky position, X along its direction of motion, is
!>
!>     X = a sin(phi),   Y = -b cos(phi),   in front of the star when cos(phi) > 0.
!>
!> Writes t, F_sphere, F_oblate as three columns to standard output.
program example
  use limbdarkfast
  implicit none
  real(dp), parameter :: PI = 3.141592653589793d0
  real(dp), parameter :: P = 3.0d0, t0 = 0.0d0, a = 10.0d0, b = 0.3d0   ! days, days, R*, R*
  real(dp), parameter :: r = 0.1d0                     ! area-equivalent if oblate
  real(dp), parameter :: f = 0.2d0, theta = 30.0d0 * PI / 180.0d0
  real(dp), parameter :: w(4) = [0.30d0, 0.20d0, 0.10d0, 0.05d0]
  integer, parameter :: n = 2001
  real(dp) :: t, phi, X, Y, Fs, Fo
  integer :: i

  do i = 0, n - 1
    t = -0.1d0 + 0.2d0 * i / (n - 1)
    phi = 2.0d0 * PI * (t - t0) / P
    X = a * sin(phi)
    Y = -b * cos(phi)
    Fs = 1.0d0
    Fo = 1.0d0
    if (cos(phi) > 0.0d0) then
      Fs = hybrid_flux(hypot(X, Y), r, w, HYBRID4)
      Fo = hybrid_flux_oblate(X, Y, r, f, theta, w, HYBRID4)
    end if
    write (*, '(f13.10, 2(1x, es24.17))') t, Fs, Fo
  end do
end program example
