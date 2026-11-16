!> limbdarkfast: transit light curves for the hybrid limb-darkening laws.
!>
!> The Fortran version: one module, Fortran 2008. It gives the flux of a
!> star with hybrid2, hybrid4 or hybrid5 limb darkening, occulted by a
!> spherical or an oblate planet, in closed form and in double precision.
!> There are no orbits, exposures or derivatives: give it the planet's
!> position on the sky and it returns the normalised flux.
!>
!> If you use this code, please cite
!>
!>     Kipping, D. 2027, "Limb Darkened Transits of Oblate Planets as
!>     Elementary Functions", AJ, submitted.
!>
!> Copyright (c) 2026 David Kipping. Released under the MIT License: see
!> the LICENSE file, or https://opensource.org/licenses/MIT.
!>
!>     use limbdarkfast
!>     F = hybrid_flux(z, r, w, HYBRID4)                       ! spherical
!>     F = hybrid_flux_oblate(X, Y, r, f, theta, w, HYBRID4)   ! oblate
!>
!> Lengths are in stellar radii. z is the sky separation of the planet's
!> centre from the star's; (X, Y) is its sky position, X along the
!> direction of motion. r is the planet's radius, for an oblate planet its
!> area-equivalent radius sqrt(A B); f = 1 - B/A is the projected
!> flattening and theta (radians) the angle of the long axis from X. z (or
!> X and Y) may be scalars or rank-1 arrays. w are the law's shape weights:
!>
!>     I(mu) / I(1) = 1 - sum_j w_j S_j(mu),
!>
!> each shape falling from 1 at the limb to 0 at the centre. The law codes
!> HYBRID2, HYBRID4 and HYBRID5 equal their number of weights (2, 4, 5).
!> The shapes combine mu^0, mu^2, mu^4 with double poles (mu^2 + eps)^-2,
!> so that the occulted flux of every term is elementary for any
!> elliptical planet. The oblate planet must satisfy r <= (1 - f)^(3/2), so
!> that its outline crosses the stellar limb at most twice. Out-of-domain
!> input, or a wrong number of weights, gives NaN.
!>
!> Method (python/limbdarkfast.py has the same structure, with more comments).
!> Green's theorem with the radial potential h(s) = (1/2s) int_0^s I,
!> s = x^2 + y^2, turns the occulted flux into h(1) times the limb arc
!> inside the planet plus the integral of h(s) along the planet's outline
!> inside the star. For a circle that reduces to the lens geometry; for an
!> ellipse the even terms are finite Fourier sums and each pole a
!> partial-fraction sum over the roots of a quartic, near-double roots
!> being integrated as one quadratic factor.
module limbdarkfast
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
  implicit none
  private

  integer, parameter, public :: dp = kind(1.0d0)
  integer, parameter, public :: HYBRID2 = 2, HYBRID4 = 4, HYBRID5 = 5

  public :: hybrid_flux, hybrid_flux_oblate, intensity, n_weights

  !> flux of a spherical planet: z scalar or rank-1
  interface hybrid_flux
    module procedure hybrid_flux_0, hybrid_flux_1
  end interface hybrid_flux
  !> flux of an oblate planet: X, Y scalars or rank-1
  interface hybrid_flux_oblate
    module procedure hybrid_flux_oblate_0, hybrid_flux_oblate_1
  end interface hybrid_flux_oblate

  real(dp), parameter :: PI = 3.141592653589793d0, TWO_PI = 6.283185307179586d0
  real(dp), parameter :: EPS = 2.220446049250313d-16, TINY = 1d-280
  real(dp), parameter :: Y_SW = 1d-3        ! spherical pole: series half-width in y
  real(dp), parameter :: F_SW = 1d-10       ! below this flattening, the spherical forms
  real(dp), parameter :: Z_TOL = 1d-4       ! |z| ~ 1: a candidate limb crossing
  real(dp), parameter :: G_TOL = 1d-12      ! accepted crossing residual, relative
  real(dp), parameter :: PAIR_TOL = 1d-6    ! crossings closer than this are a tangency
  real(dp), parameter :: PAIR_RATIO = 0.1d0 ! near-double pole roots
  integer, parameter :: N_SERIES = 18, N_ABERTH = 8, N_BAIRSTOW = 2
  integer, parameter :: QR_SWEEPS4 = 6, QR_SWEEPS3 = 5, QR_BALANCE = 4
  integer, parameter :: PAIRS(2, 6) = reshape([1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4], [2, 6])

  !> A law: n_w weights and K poles eps; row j of G gives shape S_j as a
  !> combination of the generators [mu^0, mu^2, mu^4, (mu^2 + eps_k)^-2 ...];
  !> norms(j) is the unocculted flux of S_j. From tools/law_constants.py.
  type :: law_t
    integer :: n_w, K
    real(dp) :: eps(3)
    real(dp) :: G(5, 6)
    real(dp) :: norms(5)
  end type law_t

  ! hybrid2: G row by row (shape j over generators k)
  type(law_t), parameter :: law_hybrid2 = law_t(2, 1, &
      [ &
        0.35999999999999999d0, 0d0, 0d0], &
      reshape([ &
        1d0, -1d0, 0d0, &
        0d0, 0d0, 0d0, &
        -0.075348837209302341d0, 0d0, 0d0, &
        0.13936520930232557d0, 0d0, 0d0, &
        0d0, 0d0, 0d0, &
        0d0, 0d0, 0d0, &
        0d0, 0d0, 0d0, &
        0d0, 0d0, 0d0, &
        0d0, 0d0, 0d0, &
        0d0, 0d0, 0d0], [5, 6], order=[2, 1]), &
      [ &
        1.5707963267948966d0, 0.6575426484257707d0, 0d0, &
        0d0, 0d0])
  ! hybrid4: G row by row (shape j over generators k)
  type(law_t), parameter :: law_hybrid4 = law_t(4, 2, &
      [ &
        0.028601207249486824d0, 0.38582417584735224d0, 0d0], &
      reshape([ &
        1d0, -1d0, 0d0, &
        0d0, 0d0, 0d0, &
        1d0, -2d0, 1d0, &
        0d0, 0d0, 0d0, &
        0.0031720957997626959d0, -0.0028080275108228465d0, 0.0007779295064964302d0, &
        0.0008358238567779617d0, -0.003710396155608333d0, 0d0, &
        -0.45048769921539672d0, 0.51370715445091153d0, -0.17569337890652659d0, &
        -4.8355348357818505d-07, 0.21600802064223382d0, 0d0, &
        0d0, 0d0, 0d0, &
        0d0, 0d0, 0d0], [5, 6], order=[2, 1]), &
      [ &
        1.5707963267948966d0, 1.0471975511965974d0, 0.073823572890099914d0, &
        0.47682011132930302d0, 0d0])
  ! hybrid5: G row by row (shape j over generators k)
  type(law_t), parameter :: law_hybrid5 = law_t(5, 3, &
      [ &
        0.0015965451792942786d0, 0.041275044772743019d0, 0.38743471584956696d0], &
      reshape([ &
        1d0, 0d0, -1d0, &
        0d0, 0d0, 0d0, &
        1d0, -2d0, 1d0, &
        0d0, 0d0, 0d0, &
        -2.5408433670697262d-06, 0d0, 0d0, &
        2.5489629860270405d-06, 0d0, 0d0, &
        -0.0015737187012817282d0, 0d0, 0d0, &
        0d0, 0.0017063103543144378d0, 0d0, &
        -0.084572789617394117d0, 0d0, 0d0, &
        0d0, 0d0, 0.16280051336826454d0], [5, 6], order=[2, 1]), &
      [ &
        2.0943951023931957d0, 1.0471975511965974d0, 0.0049997300166832875d0, &
        0.11978141121311216d0, 0.68577554794744777d0])

  !> the ellipse x = x0 + a cos phi, y = y0 + b sin phi: with z = e^{i phi},
  !> s(phi) = x^2 + y^2 = S0 + 2 Re(S1 z) + 2 S2 cos 2 phi and
  !> (x dy - y dx)/dphi = K0 + 2 Re(K1 z)
  type :: ellipse_t
    real(dp) :: S0, S2, K0
    complex(dp) :: S1, K1
  end type ellipse_t

contains

  ! ==========================================================================
  ! public interface
  ! ==========================================================================

  !> number of shape weights a law takes (0 for an unknown law)
  pure integer function n_weights(law)
    integer, intent(in) :: law
    select case (law)
    case (HYBRID2, HYBRID4, HYBRID5)
      n_weights = law
    case default
      n_weights = 0
    end select
  end function n_weights

  pure function get_law(law) result(L)
    integer, intent(in) :: law
    type(law_t) :: L
    select case (law)
    case (HYBRID2)
      L = law_hybrid2
    case (HYBRID4)
      L = law_hybrid4
    case default
      L = law_hybrid5
    end select
  end function get_law

  pure real(dp) function nan()
    nan = ieee_value(1.0_dp, ieee_quiet_nan)
  end function nan

  !> I(mu) / I(1) = 1 - sum_j w_j S_j(mu)
  pure real(dp) function intensity(mu, w, law)
    real(dp), intent(in) :: mu, w(:)
    integer, intent(in) :: law
    type(law_t) :: L
    real(dp) :: m2, gen(6)
    integer :: j, k
    if (n_weights(law) == 0 .or. size(w) /= n_weights(law)) then
      intensity = nan()
      return
    end if
    L = get_law(law)
    m2 = mu * mu
    gen = 0.0_dp
    gen(1:3) = [1.0_dp, m2, m2 * m2]
    do k = 1, L%K
      gen(3 + k) = 1.0_dp / (m2 + L%eps(k))**2
    end do
    intensity = 1.0_dp
    do j = 1, L%n_w
      intensity = intensity - w(j) * sum(L%G(j, 1:3 + L%K) * gen(1:3 + L%K))
    end do
  end function intensity

  !> Normalised flux of a star with a hybrid law, occulted by a spherical
  !> planet of radius r (stellar radii) at sky separation z. Exact for
  !> r < 1; for a companion larger than the star the even terms lose about
  !> r^5 x 1e-16.
  pure real(dp) function hybrid_flux_0(z, r, w, law) result(flux)
    real(dp), intent(in) :: z, r, w(:)
    integer, intent(in) :: law
    type(law_t) :: L
    real(dp) :: occ(6)
    if (n_weights(law) == 0 .or. size(w) /= n_weights(law) .or. .not. r >= 0.0_dp) then
      flux = nan()
      return
    end if
    L = get_law(law)
    call sphere_occ(z, r, L, occ)
    flux = flux_from_occ(occ, w, L)
  end function hybrid_flux_0

  pure function hybrid_flux_1(z, r, w, law) result(flux)
    real(dp), intent(in) :: z(:), r, w(:)
    integer, intent(in) :: law
    real(dp) :: flux(size(z))
    integer :: i
    do i = 1, size(z)
      flux(i) = hybrid_flux_0(z(i), r, w, law)
    end do
  end function hybrid_flux_1

  !> Normalised flux of a star with a hybrid law, occulted by an oblate
  !> planet: sky position (X, Y), X along its direction of motion;
  !> area-equivalent radius r = sqrt(A B), with A = r / sqrt(1 - f) and
  !> B = A (1 - f); flattening f in [0, 1); long axis at angle theta
  !> (radians) from X. Requires r <= (1 - f)^(3/2).
  pure real(dp) function hybrid_flux_oblate_0(X, Y, r, f, theta, w, law) result(flux)
    real(dp), intent(in) :: X, Y, r, f, theta, w(:)
    integer, intent(in) :: law
    type(law_t) :: L
    real(dp) :: occ(6), x0, y0, a
    if (n_weights(law) == 0 .or. size(w) /= n_weights(law)) then
      flux = nan()
      return
    end if
    if (.not. (f >= 0.0_dp .and. f < 1.0_dp)) then
      flux = nan()
      return
    end if
    if (.not. (r >= 0.0_dp .and. r <= (1.0_dp - f)**1.5_dp * (1.0_dp + 1d-12))) then
      flux = nan()
      return
    end if
    L = get_law(law)
    x0 = X * cos(theta) + Y * sin(theta)       ! principal frame
    y0 = -X * sin(theta) + Y * cos(theta)
    if (f < F_SW) then
      call sphere_occ(hypot(x0, y0), r, L, occ)
    else
      a = r / sqrt(1.0_dp - f)
      call ellipse_occ(x0, y0, a, a * (1.0_dp - f), L, occ)
    end if
    flux = flux_from_occ(occ, w, L)
  end function hybrid_flux_oblate_0

  pure function hybrid_flux_oblate_1(X, Y, r, f, theta, w, law) result(flux)
    real(dp), intent(in) :: X(:), Y(:), r, f, theta, w(:)
    integer, intent(in) :: law
    real(dp) :: flux(size(X))
    integer :: i
    do i = 1, size(X)
      flux(i) = hybrid_flux_oblate_0(X(i), Y(i), r, f, theta, w, law)
    end do
  end function hybrid_flux_oblate_1

  !> F = 1 - (O_0 - sum_j w_j O_Sj) / (pi - sum_j w_j N_j), from the
  !> occulted flux occ(k) of each generator
  pure real(dp) function flux_from_occ(occ, w, L) result(F)
    real(dp), intent(in) :: occ(6), w(:)
    type(law_t), intent(in) :: L
    real(dp) :: blocked, norm
    integer :: j
    blocked = occ(1)
    norm = PI
    do j = 1, L%n_w
      blocked = blocked - w(j) * sum(L%G(j, 1:3 + L%K) * occ(1:3 + L%K))
      norm = norm - w(j) * L%norms(j)
    end do
    F = 1.0_dp - blocked / norm
  end function flux_from_occ

  ! ==========================================================================
  ! spherical planet
  ! ==========================================================================

  !> 16 x (area of the triangle with sides 1, z, r)^2, sorted (Kahan) form;
  !> negative when no such triangle exists
  pure real(dp) function kite_sqarea(z, r)
    real(dp), intent(in) :: z, r
    real(dp) :: a, b, c, bb, cc
    a = max(z, r)
    c = min(z, r)
    b = min(a, 1.0_dp)
    a = max(a, 1.0_dp)
    bb = max(b, c)
    cc = min(b, c)
    kite_sqarea = (a + (bb + cc)) * (cc - (a - bb)) * (cc + (a - bb)) * (a + (bb - cc))
  end function kite_sqarea

  !> occulted flux of [mu^0, mu^2, mu^4, P_eps...] for a circle of radius r
  !> at separation z
  pure subroutine sphere_occ(zin, r, L, occ)
    real(dp), intent(in) :: zin, r
    type(law_t), intent(in) :: L
    real(dp), intent(out) :: occ(6)
    real(dp) :: z, r2, z2, alpha, sqarea, s0, s2, s4, M0, M2, M4, A, Bp
    real(dp) :: kite, kap0, kap1, eta2, nBp, e, p, Kp, Q, sqQ, apb, U, ratio, y, y2, T, G, J
    logical :: tot
    integer :: k

    z = abs(zin)
    occ = 0.0_dp
    if (z >= 1.0_dp + r) return
    r2 = r * r
    z2 = z * z
    alpha = 1.0_dp - r2 - z2
    sqarea = kite_sqarea(z, r)
    A = max((r + 1.0_dp - z) * (1.0_dp - r + z), TINY)   ! 1 - (z - r)^2
    Bp = (1.0_dp - z - r) * (1.0_dp + z + r)              ! 1 - (z + r)^2

    if (z <= 1.0_dp - r) then                             ! wholly on the disc
      s0 = -PI * r2
      s2 = TWO_PI * r2 * (r2 + 2.0_dp * z2 - 1.0_dp)
      M0 = PI
      M2 = PI * alpha
      M4 = (6.0_dp * alpha * M2 + 2.0_dp * sqarea * M0) / 4.0_dp
      s4 = -(2.0_dp * r2 * M4 - (2.0_dp / 3.0_dp) * (alpha * M4 + sqarea * M2))
      occ(1) = -s0
      occ(2) = -(0.5_dp * s0 + 0.25_dp * s2)
      occ(3) = -(s0 / 3.0_dp + s2 / 6.0_dp + s4 / 6.0_dp)
      do k = 1, L%K
        e = L%eps(k)
        p = 1.0_dp + e
        Kp = p + r2 - z2
        sqQ = sqrt((e + A) * (e + Bp))
        occ(3 + k) = TWO_PI * r2 / (sqQ * (Kp + sqQ))
      end do
      return
    end if

    ! partial overlap (or, for r > 1, total occultation: kite = 0 there, and
    ! every term lands on its full-disc value)
    tot = z <= r - 1.0_dp
    if (tot) then
      kite = 0.0_dp
    else
      kite = sqrt(max(sqarea, (10.0_dp * EPS)**2))
    end if
    kap0 = atan2(kite, r2 + z2 - 1.0_dp)     ! half-angle at the planet
    kap1 = atan2(kite, 1.0_dp + z2 - r2)     ! half-angle at the star

    s0 = -(kap1 + r2 * kap0 - 0.5_dp * kite)
    eta2 = r2 * (r2 + 2.0_dp * z2)
    s2 = 2.0_dp * s0 + 2.0_dp * (kap1 + eta2 * kap0 - 0.25_dp * kite * (1.0_dp + 5.0_dp * r2 + z2))
    M0 = kap0
    M2 = kap0 * alpha + kite
    M4 = (6.0_dp * alpha * M2 + 2.0_dp * sqarea * M0) / 4.0_dp
    s4 = -(2.0_dp * r2 * M4 - (2.0_dp / 3.0_dp) * (alpha * M4 + sqarea * M2))
    occ(1) = -s0
    occ(2) = -(0.5_dp * s0 + 0.25_dp * s2)
    occ(3) = -(s0 / 3.0_dp + s2 / 6.0_dp + s4 / 6.0_dp)

    ! the smaller of A and -Bp is taken from kite^2 = A (-Bp), exactly
    if (.not. tot .and. A >= -Bp) then
      nBp = kite * kite / A
    else if (-Bp > 0.0_dp) then
      nBp = -Bp
    else
      nBp = kite * kite / A
    end if
    do k = 1, L%K
      e = L%eps(k)
      p = 1.0_dp + e
      Kp = p + r2 - z2
      Q = (e + A) * (e + Bp)
      ! J = int_{-kap0}^{kap0} dpsi / (p - z^2 - r^2 + 2 z r cos psi)
      apb = e + A
      U = apb * nBp
      ratio = (e + Bp) / nBp
      y = ratio * kite * kite / U
      if (y > Y_SW) then
        J = 4.0_dp / sqrt(Q) * atan2(kite * sqrt(ratio), sqrt(U))
      else if (y < -Y_SW) then
        J = 4.0_dp / sqrt(-Q) * log((sqrt(U) + kite * sqrt(-ratio)) &
            / (2.0_dp * sqrt(max(z * r * e, TINY))))
      else                                   ! across Q = 0: a series
        T = kite / nBp                       ! tan(kap0 / 2)
        y2 = y * y
        G = 1.0_dp - y / 3.0_dp + y2 / 5.0_dp - y * y2 / 7.0_dp + y2 * y2 / 9.0_dp &
            - y * y2 * y2 / 11.0_dp
        J = 4.0_dp * T * G / apb
      end if
      occ(3 + k) = kap1 / (p * e) + (Kp * J - 2.0_dp * kap0) / (4.0_dp * p)
    end do
  end subroutine sphere_occ

  ! ==========================================================================
  ! oblate planet: the points where s = p are the roots of the self-inversive
  ! quartic D(z) = S2 z^4 + S1 z^3 + (S0 - p) z^2 + conj(S1) z + S2
  ! ==========================================================================

  pure function make_ellipse(x0, y0, a, b) result(E)
    real(dp), intent(in) :: x0, y0, a, b
    type(ellipse_t) :: E
    E%S0 = x0 * x0 + y0 * y0 + 0.5_dp * (a * a + b * b)
    E%S2 = 0.25_dp * (a * a - b * b)
    E%K0 = a * b
    E%S1 = cmplx(a * x0, -b * y0, kind=dp)
    E%K1 = cmplx(0.5_dp * b * x0, -0.5_dp * a * y0, kind=dp)
  end function make_ellipse

  pure complex(dp) function Dz(E, z, p)
    type(ellipse_t), intent(in) :: E
    complex(dp), intent(in) :: z
    real(dp), intent(in) :: p
    Dz = (((E%S2 * z + E%S1) * z + (E%S0 - p)) * z + conjg(E%S1)) * z + E%S2
  end function Dz

  pure complex(dp) function dDz(E, z, p)
    type(ellipse_t), intent(in) :: E
    complex(dp), intent(in) :: z
    real(dp), intent(in) :: p
    dDz = ((4.0_dp * E%S2 * z + 3.0_dp * E%S1) * z + 2.0_dp * (E%S0 - p)) * z + conjg(E%S1)
  end function dDz

  pure complex(dp) function Nz(E, z)
    type(ellipse_t), intent(in) :: E
    complex(dp), intent(in) :: z
    Nz = conjg(E%K1) + E%K0 * z + E%K1 * z * z
  end function Nz

  pure real(dp) function carg(z)
    complex(dp), intent(in) :: z
    carg = atan2(aimag(z), real(z))
  end function carg

  pure logical function cfinite(z)
    complex(dp), intent(in) :: z
    cfinite = ieee_is_finite(real(z)) .and. ieee_is_finite(aimag(z))
  end function cfinite

  pure complex(dp) function expi(t)
    real(dp), intent(in) :: t
    expi = cmplx(cos(t), sin(t), kind=dp)
  end function expi

  ! ---- the quartic's roots: a balanced, shifted 4x4 QR, then Aberth --------

  !> (c, s), c real: [[c, s], [-conj(s), c]] [x, y] = [rho, 0]
  pure subroutine givens(x, y, c, s)
    complex(dp), intent(in) :: x, y
    real(dp), intent(out) :: c
    complex(dp), intent(out) :: s
    real(dp) :: xa, rho
    complex(dp) :: ph
    xa = abs(x)
    rho = sqrt(xa * xa + real(y)**2 + aimag(y)**2)
    if (xa > 0.0_dp) then
      ph = x / xa
    else
      ph = (1.0_dp, 0.0_dp)
    end if
    if (rho > 0.0_dp) then
      c = xa / rho
      s = ph * conjg(y) / rho
    else
      c = 1.0_dp
      s = (0.0_dp, 0.0_dp)
    end if
  end subroutine givens

  !> the eigenvalue of [[a, b], [c, d]] nearer d
  pure complex(dp) function wilkinson(a, b, c, d)
    complex(dp), intent(in) :: a, b, c, d
    complex(dp) :: h, disc, den
    h = 0.5_dp * (a - d)
    disc = sqrt(h * h + b * c)
    if (abs(h + disc) >= abs(h - disc)) then
      den = h + disc
    else
      den = h - disc
    end if
    if (abs(den) > 0.0_dp) then
      wilkinson = d - b * c / den
    else
      wilkinson = d
    end if
  end function wilkinson

  !> one shifted QR step on the leading m x m block of H
  pure subroutine qr_sweep(H, m)
    complex(dp), intent(inout) :: H(4, 4)
    integer, intent(in) :: m
    complex(dp) :: sig, x, y, sn(3)
    real(dp) :: cs(3)
    integer :: i, j, k
    sig = wilkinson(H(m - 1, m - 1), H(m - 1, m), H(m, m - 1), H(m, m))
    do i = 1, m
      H(i, i) = H(i, i) - sig
    end do
    do k = 1, m - 1
      call givens(H(k, k), H(k + 1, k), cs(k), sn(k))
      do j = k, m
        x = H(k, j)
        y = H(k + 1, j)
        H(k, j) = cs(k) * x + sn(k) * y
        H(k + 1, j) = cs(k) * y - conjg(sn(k)) * x
      end do
    end do
    do k = 1, m - 1
      do i = 1, min(k + 2, m)
        x = H(i, k)
        y = H(i, k + 1)
        H(i, k) = cs(k) * x + conjg(sn(k)) * y
        H(i, k + 1) = cs(k) * y - sn(k) * x
      end do
    end do
    do i = 1, m
      H(i, i) = H(i, i) + sig
    end do
  end subroutine qr_sweep

  !> eigenvalues of D's companion matrix
  pure subroutine qr_seeds(E, p, z)
    type(ellipse_t), intent(in) :: E
    real(dp), intent(in) :: p
    complex(dp), intent(out) :: z(4)
    complex(dp) :: H(4, 4), a, b, c, d, h2, disc, mid
    real(dp) :: inv, col, row, f
    integer :: it, i, j
    H = (0.0_dp, 0.0_dp)
    inv = 1.0_dp / E%S2
    H(1, 1) = -E%S1 * inv
    H(1, 2) = -(E%S0 - p) * inv
    H(1, 3) = -conjg(E%S1) * inv
    H(1, 4) = (-1.0_dp, 0.0_dp)
    H(2, 1) = (1.0_dp, 0.0_dp)
    H(3, 2) = (1.0_dp, 0.0_dp)
    H(4, 3) = (1.0_dp, 0.0_dp)
    do it = 1, QR_BALANCE                    ! diagonal similarity
      do i = 1, 4
        col = 0.0_dp
        row = 0.0_dp
        do j = 1, 4
          if (j == i) cycle
          col = col + abs(H(j, i))
          row = row + abs(H(i, j))
        end do
        f = 1.0_dp
        if (col > 0.0_dp .and. row > 0.0_dp) f = sqrt(row / col)
        H(:, i) = H(:, i) * f
        H(i, :) = H(i, :) / f
      end do
    end do
    do it = 1, QR_SWEEPS4
      call qr_sweep(H, 4)
    end do
    do it = 1, QR_SWEEPS3
      call qr_sweep(H, 3)
    end do
    a = H(1, 1)
    b = H(1, 2)
    c = H(2, 1)
    d = H(2, 2)
    h2 = 0.5_dp * (a - d)
    disc = sqrt(h2 * h2 + b * c)
    mid = 0.5_dp * (a + d)
    z = [mid + disc, mid - disc, H(3, 3), H(4, 4)]
  end subroutine qr_seeds

  !> the four roots of D at level p: eigenvalue seeds, then Aberth
  !> iterations (at most N_ABERTH, stopping once every correction is at
  !> rounding) and two Newton polishes; a correction that is not finite is
  !> skipped
  pure subroutine roots(E, p, z)
    type(ellipse_t), intent(in) :: E
    real(dp), intent(in) :: p
    complex(dp), intent(out) :: z(4)
    complex(dp) :: rat(4), nz(4), s, step, nw
    logical :: done
    integer :: it, i, j
    call qr_seeds(E, p, z)
    do it = 1, N_ABERTH
      do i = 1, 4
        rat(i) = Dz(E, z(i), p) / dDz(E, z(i), p)
      end do
      done = .true.
      do i = 1, 4
        s = (0.0_dp, 0.0_dp)
        do j = 1, 4
          if (j /= i .and. abs(z(i) - z(j)) > 0.0_dp) s = s + 1.0_dp / (z(i) - z(j))
        end do
        step = rat(i) / (1.0_dp - rat(i) * s)
        if (cfinite(step)) then
          nz(i) = z(i) - step
          if (abs(step) > 4.0_dp * EPS * abs(z(i))) done = .false.
        else
          nz(i) = z(i)
        end if
      end do
      z = nz
      if (done) exit
    end do
    do it = 1, 2
      do i = 1, 4
        nw = z(i) - Dz(E, z(i), p) / dDz(E, z(i), p)
        if (cfinite(nw)) z(i) = nw
      end do
    end do
  end subroutine roots

  ! ---- limb crossings and the arcs ------------------------------------------

  pure real(dp) function mod2pi(t)
    real(dp), intent(in) :: t
    mod2pi = modulo(t, TWO_PI)           ! real-valued, in [0, 2 pi)
  end function mod2pi

  pure subroutine g_dg(phi, x0, y0, a, b, g, dg)
    real(dp), intent(in) :: phi, x0, y0, a, b
    real(dp), intent(out) :: g, dg
    real(dp) :: s, c, x, y
    s = sin(phi)
    c = cos(phi)
    x = x0 + a * c
    y = y0 + b * s
    g = x * x + y * y - 1.0_dp
    dg = 2.0_dp * (-x * a * s + y * b * c)
  end subroutine g_dg

  !> number of limb crossings, and the first two (phi_a <= phi_b) on the
  !> ellipse parameter; crossings closer than PAIR_TOL are a tangency and
  !> are dropped together
  pure subroutine crossings(x0, y0, a, b, n, phi_a, phi_b)
    real(dp), intent(in) :: x0, y0, a, b
    integer, intent(out) :: n
    real(dp), intent(out) :: phi_a, phi_b
    type(ellipse_t) :: E
    complex(dp) :: z(4)
    real(dp) :: sc, phi(4), keep(4), ph, g, dg, trial, gt, dgt, d
    logical :: valid(4), drop(4), cand
    integer :: i, j, it
    E = make_ellipse(x0, y0, a, b)
    call roots(E, 1.0_dp, z)
    sc = 1.0_dp + (abs(x0) + a)**2 + (abs(y0) + b)**2
    do i = 1, 4
      cand = abs(abs(z(i)) - 1.0_dp) < Z_TOL
      ph = 0.0_dp
      if (cand) ph = carg(z(i))
      call g_dg(ph, x0, y0, a, b, g, dg)
      do it = 1, 3
        trial = ph
        if (abs(dg) > 0.0_dp) trial = ph - g / dg
        call g_dg(trial, x0, y0, a, b, gt, dgt)
        if (abs(gt) < abs(g)) then
          ph = trial
          g = gt
          dg = dgt
        end if
      end do
      valid(i) = cand .and. abs(g) <= G_TOL * sc
      phi(i) = mod2pi(ph)
    end do
    drop = .false.
    do i = 1, 4
      do j = i + 1, 4
        d = abs(mod2pi(phi(i) - phi(j) + PI) - PI)
        if (valid(i) .and. valid(j) .and. d < PAIR_TOL) then
          drop(i) = .true.
          drop(j) = .true.
        end if
      end do
    end do
    n = 0
    do i = 1, 4                              ! insertion sort, ascending
      if (.not. valid(i) .or. drop(i)) cycle
      n = n + 1
      j = n
      do while (j > 1)
        if (keep(j - 1) <= phi(i)) exit
        keep(j) = keep(j - 1)
        j = j - 1
      end do
      keep(j) = phi(i)
    end do
    phi_a = 0.0_dp
    phi_b = 0.0_dp
    if (n > 0) phi_a = keep(1)
    if (n > 1) phi_b = keep(2)
  end subroutine crossings

  !> the ellipse arc [lo, hi] inside the star and the length dth of the limb
  !> arc inside the ellipse, for two crossings. The ellipse arc is decided at
  !> the midpoint of the longer of its two arcs; the limb arc inside the
  !> ellipse is always the shorter one (the planet lies within radius A < 1).
  pure subroutine arcs(x0, y0, a, b, phi_a, phi_b, lo, hi, dth)
    real(dp), intent(in) :: x0, y0, a, b, phi_a, phi_b
    real(dp), intent(out) :: lo, hi, dth
    logical :: long_first, mid_in, first
    real(dp) :: mid, xm, ym, ta, tb, dt
    long_first = (phi_b - phi_a) >= PI
    if (long_first) then
      mid = 0.5_dp * (phi_a + phi_b)
    else
      mid = 0.5_dp * (phi_a + phi_b) + PI
    end if
    xm = x0 + a * cos(mid)
    ym = y0 + b * sin(mid)
    mid_in = xm * xm + ym * ym < 1.0_dp
    if (long_first) then                     ! is [phi_a, phi_b] inside?
      first = mid_in
    else
      first = .not. mid_in
    end if
    if (first) then
      lo = phi_a
      hi = phi_b
    else
      lo = phi_b
      hi = phi_a + TWO_PI
    end if
    ta = atan2(y0 + b * sin(phi_a), x0 + a * cos(phi_a))
    tb = atan2(y0 + b * sin(phi_b), x0 + a * cos(phi_b))
    dt = abs(mod2pi(ta) - mod2pi(tb))
    dth = min(dt, TWO_PI - dt)
  end subroutine arcs

  ! ---- even terms: Fourier coefficients of h_n(s) K ---------------------------

  !> q(n, m), n = 0..2, m = 0..5: h_n(s) K = Re sum_m q_m e^{i m phi}, with
  !> h_n(s) = sum_p alpha(n, p) s^p, alpha(n, p) = (-1)^p C(n, p) / (2 (p + 1))
  pure subroutine q_coeffs(E, q)
    type(ellipse_t), intent(in) :: E
    complex(dp), intent(out) :: q(0:2, 0:5)
    real(dp), parameter :: alpha(0:2, 0:2) = reshape( &
        [0.5_dp, 0.5_dp, 0.5_dp, 0.0_dp, -0.25_dp, -0.5_dp, 0.0_dp, 0.0_dp, 1.0_dp / 6.0_dp], &
        [3, 3])
    complex(dp) :: T(-5:5), Tn(-5:5), Sk(-2:2)
    integer :: n, m, p, i, j
    T = (0.0_dp, 0.0_dp)
    T(-1) = conjg(E%K1)
    T(0) = E%K0
    T(1) = E%K1
    Sk = [cmplx(E%S2, 0.0_dp, kind=dp), conjg(E%S1), cmplx(E%S0, 0.0_dp, kind=dp), E%S1, &
          cmplx(E%S2, 0.0_dp, kind=dp)]
    q = (0.0_dp, 0.0_dp)
    do p = 0, 2
      if (p > 0) then                        ! T <- S * T
        do i = -5, 5
          Tn(i) = (0.0_dp, 0.0_dp)
          do j = -2, 2
            if (i - j >= -5 .and. i - j <= 5) Tn(i) = Tn(i) + Sk(j) * T(i - j)
          end do
        end do
        T = Tn
      end if
      do n = p, 2
        do m = 0, 5
          q(n, m) = q(n, m) + alpha(n, p) * T(m)
        end do
      end do
    end do
  end subroutine q_coeffs

  !> Psi_n(hi) - Psi_n(lo), Psi_n the exact antiderivative, in difference form
  pure real(dp) function even_arc(qn, lo, hi)
    complex(dp), intent(in) :: qn(0:5)
    real(dp), intent(in) :: lo, hi
    real(dp) :: mid, half
    integer :: m
    mid = 0.5_dp * (lo + hi)
    half = 0.5_dp * (hi - lo)
    even_arc = real(qn(0)) * (hi - lo)
    do m = 1, 5
      even_arc = even_arc + (4.0_dp / m) * real(qn(m) * expi(m * mid)) * sin(m * half)
    end do
  end function even_arc

  ! ---- poles: partial fractions over the roots of D ---------------------------

  !> Newton steps on (sigma, pi) for the quadratic factor z^2 - sigma z + pi
  pure subroutine bairstow(E, p, sig, pq)
    type(ellipse_t), intent(in) :: E
    real(dp), intent(in) :: p
    complex(dp), intent(inout) :: sig, pq
    complex(dp) :: dc(0:4), al, be, das, dbs, dap, dbp, r1, r0, J11, J12, J21, J22
    complex(dp) :: nal, nbe, ndas, ndbs, ndap, ndbp, det
    integer :: it, k
    dc = [cmplx(E%S2, 0.0_dp, kind=dp), conjg(E%S1), cmplx(E%S0 - p, 0.0_dp, kind=dp), E%S1, &
          cmplx(E%S2, 0.0_dp, kind=dp)]              ! low -> high
    do it = 1, N_BAIRSTOW
      al = (0.0_dp, 0.0_dp)                  ! z^k mod q = al z + be
      be = (1.0_dp, 0.0_dp)
      das = (0.0_dp, 0.0_dp)
      dbs = das
      dap = das
      dbp = das
      r1 = dc(0) * al
      r0 = dc(0) * be
      J11 = das
      J12 = das
      J21 = das
      J22 = das
      do k = 1, 4
        nal = al * sig + be
        nbe = -al * pq
        ndas = das * sig + al + dbs
        ndbs = -das * pq
        ndap = dap * sig + dbp
        ndbp = -dap * pq - al
        al = nal
        be = nbe
        das = ndas
        dbs = ndbs
        dap = ndap
        dbp = ndbp
        r1 = r1 + dc(k) * al
        r0 = r0 + dc(k) * be
        J11 = J11 + dc(k) * das
        J12 = J12 + dc(k) * dap
        J21 = J21 + dc(k) * dbs
        J22 = J22 + dc(k) * dbp
      end do
      det = J11 * J22 - J12 * J21
      if (.not. abs(det) > 0.0_dp) exit
      sig = sig - (J22 * r1 - J12 * r0) / det
      pq = pq - (-J21 * r1 + J11 * r0) / det
    end do
  end subroutine bairstow

  !> for the quadratic factor q of D, with N = aN z + bN and Q = D/q = aQ z + bQ
  !> (mod q): the sum of N/D' over its two roots (A_sum) and the C of
  !> int (A z + B)/q dz (C_num, before the factor i); symmetric in the roots
  pure subroutine pair_terms(E, p, sig, pq, A_sum, C_num)
    type(ellipse_t), intent(in) :: E
    real(dp), intent(in) :: p
    complex(dp), intent(in) :: sig, pq
    complex(dp), intent(out) :: A_sum, C_num
    complex(dp) :: Q2, Q1, Q0, aN, bN, aQ, bQ, den
    Q2 = E%S2
    Q1 = E%S1 + sig * Q2
    Q0 = sig * Q1 + (E%S0 - p) - pq * Q2
    aN = E%K1 * sig + E%K0
    bN = conjg(E%K1) - E%K1 * pq
    aQ = Q2 * sig + Q1
    bQ = Q0 - Q2 * pq
    den = aQ * aQ * pq + aQ * bQ * sig + bQ * bQ
    if (.not. abs(den) > 0.0_dp) den = (1.0_dp, 0.0_dp)
    A_sum = (aN * bQ - bN * aQ) / den
    C_num = 0.5_dp * (2.0_dp * aN * aQ * pq + (aN * bQ + bN * aQ) * sig + 2.0_dp * bN * bQ) / den
  end subroutine pair_terms

  !> a whole ellipse on the disc: -(pi/p) Re sum_{|z_r| < 1} N/D'(z_r), the two
  !> inner roots summed as one quadratic factor
  pure real(dp) function pole_inside(E, p, z)
    type(ellipse_t), intent(in) :: E
    real(dp), intent(in) :: p
    complex(dp), intent(in) :: z(4)
    complex(dp) :: sig, pq, A_sum, C_num
    integer :: i
    sig = (0.0_dp, 0.0_dp)
    pq = (1.0_dp, 0.0_dp)
    do i = 1, 4
      if (abs(z(i)) < 1.0_dp) then
        sig = sig + z(i)
        pq = pq * z(i)
      end if
    end do
    call bairstow(E, p, sig, pq)
    call pair_terms(E, p, sig, pq, A_sum, C_num)
    pole_inside = -PI * real(A_sum) / p
  end function pole_inside

  !> change of log(e^{i phi} - z) for phi from lo to hi, continuous on the arc
  pure complex(dp) function dlog_arc(z, lo, hi, w1, w2)
    complex(dp), intent(in) :: z, w1, w2
    real(dp), intent(in) :: lo, hi
    real(dp) :: dm, da
    dm = log(abs(w2 - z)) - log(abs(w1 - z))
    if (abs(z) < 1.0_dp) then
      da = (hi - lo) + carg(1.0_dp - z / w2) - carg(1.0_dp - z / w1)
    else
      da = carg(1.0_dp - w2 / z) - carg(1.0_dp - w1 / z)
    end if
    dlog_arc = cmplx(dm, da, kind=dp)
  end function dlog_arc

  !> int du / (u^2 - Delta) = -(1/u) sum_k (Delta/u^2)^k / (2k + 1)
  pure complex(dp) function Jser(u, Delta)
    complex(dp), intent(in) :: u, Delta
    complex(dp) :: x, acc
    integer :: k
    x = Delta / (u * u)
    acc = (0.0_dp, 0.0_dp)
    do k = N_SERIES - 1, 0, -1
      acc = acc * x + 1.0_dp / (2.0_dp * k + 1.0_dp)
    end do
    Jser = -acc / u
  end function Jser

  !> (1/2p) int_lo^hi K / (p - s) dphi along one ellipse arc inside the star:
  !> (1/2p) Re sum_r R_r [log(e^{i phi} - z_r)], R_r = -N(z_r) / (i D'(z_r))
  pure real(dp) function pole_arc(E, p, z, lo, hi)
    type(ellipse_t), intent(in) :: E
    real(dp), intent(in) :: p, lo, hi
    complex(dp), intent(in) :: z(4)
    complex(dp), parameter :: I1 = (0.0_dp, 1.0_dp)
    complex(dp) :: w1, w2, Dl, c, ell, zeta, nw, dDv, R, b1, b2, z1, z2, sig, pq
    complex(dp) :: A_sum, C_num, qh, ql, LQ, Delta, cen, Jd
    real(dp) :: best, score, ang, d, a2, darg, rough, wind
    logical :: used(4), on, far, pgood(2)
    integer :: rnd, m, j, l, bi, k, pj(2), pl(2)

    w1 = expi(lo)
    w2 = expi(hi)
    ! near-double pairs: up to two disjoint, scored by |Delta| / d^2, d the
    ! distance from the pair's centre to the arc
    used = .false.
    do rnd = 1, 2
      best = huge(1.0_dp)
      bi = 1
      do m = 1, 6
        j = PAIRS(1, m)
        l = PAIRS(2, m)
        if (used(j) .or. used(l)) cycle
        Dl = 0.25_dp * (z(j) - z(l)) * (z(j) - z(l))
        c = 0.5_dp * (z(j) + z(l))
        ang = mod2pi(carg(c))
        on = (ang >= lo .and. ang <= hi) .or. (ang + TWO_PI >= lo .and. ang + TWO_PI <= hi)
        if (on) then
          d = abs(abs(c) - 1.0_dp)
        else
          d = min(abs(w1 - c), abs(w2 - c))
        end if
        score = abs(Dl) / max(d * d, 1d-300)
        if (score < best) then
          best = score
          bi = m
        end if
      end do
      pj(rnd) = PAIRS(1, bi)
      pl(rnd) = PAIRS(2, bi)
      pgood(rnd) = best <= PAIR_RATIO
      if (pgood(rnd)) then
        used(pj(rnd)) = .true.
        used(pl(rnd)) = .true.
      end if
    end do

    ell = (0.0_dp, 0.0_dp)
    ! isolated roots; one beyond |z| = 2 is taken through its partner
    ! zeta = 1/conj(z) inside the disc, where R(z) = conj(R(zeta)) and the
    ! log's change is that of log(1 - e^{i phi} conj(zeta))
    do k = 1, 4
      if (used(k)) cycle
      a2 = real(z(k))**2 + aimag(z(k))**2
      far = a2 > 4.0_dp
      if (far) then
        zeta = z(k) / a2
      else
        zeta = z(k)
      end if
      nw = zeta - Dz(E, zeta, p) / dDz(E, zeta, p)
      if (cfinite(nw)) then
        if (abs(Dz(E, nw, p)) < abs(Dz(E, zeta, p))) zeta = nw
      end if
      dDv = dDz(E, zeta, p)
      if (.not. abs(dDv) > 0.0_dp) cycle
      R = -Nz(E, zeta) / (I1 * dDv)
      if (far) then
        b2 = 1.0_dp - w2 * conjg(zeta)
        b1 = 1.0_dp - w1 * conjg(zeta)
        ell = ell + conjg(R) * cmplx(log(abs(b2)) - log(abs(b1)), carg(b2) - carg(b1), kind=dp)
      else
        ell = ell + R * dlog_arc(zeta, lo, hi, w1, w2)
      end if
    end do
    ! near-double pairs, each as one quadratic factor q = z^2 - sigma z + pi:
    ! int (A z + B)/q dz = (A/2) log q + C J(z - sigma/2), the winding of
    ! log q taken from the two roots' own continued logs
    do rnd = 1, 2
      if (.not. pgood(rnd)) cycle
      z1 = z(pj(rnd))
      z2 = z(pl(rnd))
      sig = z1 + z2
      pq = z1 * z2
      call bairstow(E, p, sig, pq)
      call pair_terms(E, p, sig, pq, A_sum, C_num)
      qh = w2 * w2 - sig * w2 + pq
      ql = w1 * w1 - sig * w1 + pq
      darg = carg(qh / ql)
      rough = aimag(dlog_arc(z1, lo, hi, w1, w2)) + aimag(dlog_arc(z2, lo, hi, w1, w2))
      wind = anint((rough - darg) / TWO_PI)
      LQ = cmplx(log(abs(qh)) - log(abs(ql)), darg + TWO_PI * wind, kind=dp)
      Delta = 0.25_dp * sig * sig - pq
      cen = 0.5_dp * sig
      Jd = Jser(w2 - cen, Delta) - Jser(w1 - cen, Delta)
      ell = ell + 0.5_dp * (I1 * A_sum) * LQ + (I1 * C_num) * Jd
    end do
    pole_arc = real(ell) / (2.0_dp * p)
  end function pole_arc

  !> occulted flux of [mu^0, mu^2, mu^4, P_eps...] for the ellipse with centre
  !> (x0, y0) in its principal frame and semi-axes a > b
  pure subroutine ellipse_occ(x0, y0, a, b, L, occ)
    real(dp), intent(in) :: x0, y0, a, b
    type(law_t), intent(in) :: L
    real(dp), intent(out) :: occ(6)
    type(ellipse_t) :: E
    complex(dp) :: q(0:2, 0:5), z(4)
    real(dp) :: d, phi_a, phi_b, lo, hi, dth, ek, p
    logical :: fast_in, general, partial, inside
    integer :: cnt, n, k
    occ = 0.0_dp
    d = hypot(x0, y0)
    fast_in = d + a < 1.0_dp
    general = .not. fast_in .and. (d - a < 1.0_dp)
    if (.not. fast_in .and. .not. general) return     ! disjoint
    cnt = 0
    phi_a = 0.0_dp
    phi_b = 0.0_dp
    if (general) call crossings(x0, y0, a, b, cnt, phi_a, phi_b)
    partial = general .and. cnt == 2
    inside = fast_in .or. (general .and. cnt /= 2 .and. d < 1.0_dp)
    if (.not. partial .and. .not. inside) return

    E = make_ellipse(x0, y0, a, b)
    call q_coeffs(E, q)
    if (inside) then
      do n = 0, 2
        occ(n + 1) = TWO_PI * real(q(n, 0))
      end do
      do k = 1, L%K
        p = 1.0_dp + L%eps(k)
        call roots(E, p, z)
        occ(3 + k) = pole_inside(E, p, z)
      end do
      return
    end if
    call arcs(x0, y0, a, b, phi_a, phi_b, lo, hi, dth)
    do n = 0, 2
      occ(n + 1) = even_arc(q(n, :), lo, hi) + dth / (2.0_dp * (n + 1))
    end do
    do k = 1, L%K
      ek = L%eps(k)
      p = 1.0_dp + ek
      call roots(E, p, z)
      occ(3 + k) = pole_arc(E, p, z, lo, hi) + dth / (2.0_dp * p * ek)
    end do
  end subroutine ellipse_occ

end module limbdarkfast
