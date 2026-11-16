!> Check the Fortran version against the shared test vectors.
!>
!>     ./test_limbdarkfast [path/to/tests/vectors]
!>
!> Stops with a non-zero code if any point differs from the reference by
!> more than TOL.
program test_limbdarkfast
  use limbdarkfast
  implicit none
  real(dp), parameter :: TOL = 1d-12
  ! brute-force integration, as in Agol, Luger & Foreman-Mackey (2020)
  ! (test/test_transit_poly.jl of Limbdark.jl): sum the intensity over NS
  ! annuli of the star, each weighted by the angle Phi(s) it spends inside
  ! the planet, with the midpoint rule
  integer, parameter :: NS = 5000                 ! annuli
  integer, parameter :: M_ARC = 48                ! samples in angle before bisecting
  real(dp), parameter :: NUM_TOL = 1d-6, PI = 3.141592653589793d0
  real(dp), parameter :: NUM_F = 0.2d0, NUM_BETA = 0.6d0, NUM_THETA = 0.4d0
  character(len=*), parameter :: laws(3) = ["hybrid2", "hybrid4", "hybrid5"]
  integer, parameter :: codes(3) = [HYBRID2, HYBRID4, HYBRID5]
  character(len=6), parameter :: kinds(2) = ["sphere", "oblate"]
  character(len=512) :: dir
  integer :: ik, il, j, failures, npts
  real(dp) :: worst, w(5)

  dir = "../tests/vectors"
  if (command_argument_count() > 0) call get_command_argument(1, dir)
  failures = 0
  do ik = 1, 2
    do il = 1, 3
      call check(trim(dir) // "/" // trim(kinds(ik)) // "_" // laws(il) // ".txt", &
                 ik == 1, codes(il), npts, worst)
      if (worst <= TOL) then
        write (*, '(a6, 1x, a, ": ", i5, " points, max |dF| = ", es8.2, "  ok")') &
            kinds(ik), laws(il), npts, worst
      else
        write (*, '(a6, 1x, a, ": ", i5, " points, max |dF| = ", es8.2, "  FAIL")') &
            kinds(ik), laws(il), npts, worst
        failures = failures + 1
      end if
    end do
  end do
  ! every shape is 1 at the limb and 0 at the centre
  do il = 1, 3
    do j = 1, n_weights(codes(il))
      w = 0.0_dp
      w(j) = 1.0_dp
      if (abs(intensity(0.0_dp, w(1:n_weights(codes(il))), codes(il))) > 1d-12 .or. &
          abs(intensity(1.0_dp, w(1:n_weights(codes(il))), codes(il)) - 1.0_dp) > 1d-12) then
        write (*, '(a, i0, a, a)') "shape ", j, " is not 1 at the limb and 0 at the centre: ", laws(il)
        failures = failures + 1
      end if
    end do
  end do
  if (failures == 0) write (*, '(a)') "law shapes: ok"
  worst = numeric_checks()
  if (worst <= NUM_TOL) then
    write (*, '("brute-force integration: max |dF| = ", es8.2, "  ok")') worst
  else
    write (*, '("brute-force integration: max |dF| = ", es8.2, "  FAIL")') worst
    failures = failures + 1
  end if
  if (failures > 0) error stop 1

contains

  subroutine check(path, sphere, law, npts, worst)
    character(len=*), intent(in) :: path
    logical, intent(in) :: sphere
    integer, intent(in) :: law
    integer, intent(out) :: npts
    real(dp), intent(out) :: worst
    character(len=1024) :: line
    real(dp) :: row(11), F, err
    integer :: unit, ios, ngeo, nw
    nw = n_weights(law)
    ngeo = merge(2, 5, sphere)
    open (newunit=unit, file=path, status="old", action="read", iostat=ios)
    if (ios /= 0) then
      write (*, '(2a)') "cannot open ", path
      error stop 2
    end if
    npts = 0
    worst = 0.0_dp
    do
      read (unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0 .or. line(1:1) == "#") cycle
      read (line, *) row(1:ngeo + nw + 1)
      if (sphere) then
        F = hybrid_flux(row(1), row(2), row(3:2 + nw), law)
      else
        F = hybrid_flux_oblate(row(1), row(2), row(3), row(4), row(5), row(6:5 + nw), law)
      end if
      err = abs(F - row(ngeo + nw + 1))
      if (.not. err <= worst) worst = err        ! NaN counts as worst
      npts = npts + 1
    end do
    close (unit)
  end subroutine check

  real(dp) function num_total(w, law)
    real(dp), intent(in) :: w(:)
    integer, intent(in) :: law
    real(dp) :: s
    integer :: j
    num_total = 0.0_dp
    do j = 0, NS - 1
      s = (j + 0.5_dp) / NS
      num_total = num_total + 2.0_dp * PI * s / NS * intensity(sqrt(1.0_dp - s * s), w, law)
    end do
  end function num_total

  !> spherical planet, by brute force
  real(dp) function numeric_flux(z, r, w, law)
    real(dp), intent(in) :: z, r, w(:)
    integer, intent(in) :: law
    real(dp) :: s1, s2, ds, s, phi, c, occ
    integer :: j
    numeric_flux = 1.0_dp
    if (z >= 1.0_dp + r) return
    s1 = max(0.0_dp, z - r)
    s2 = min(1.0_dp, z + r)
    ds = (s2 - s1) / NS
    occ = 0.0_dp
    do j = 0, NS - 1
      s = s1 + (j + 0.5_dp) * ds
      phi = 2.0_dp * PI
      if (.not. s < r - z) then
        c = (s * s + z * z - r * r) / (2.0_dp * z * s + 1d-300)
        phi = 2.0_dp * acos(max(-1.0_dp, min(1.0_dp, c)))
      end if
      occ = occ + s * ds * phi * intensity(sqrt(1.0_dp - s * s), w, law)
    end do
    numeric_flux = 1.0_dp - occ / num_total(w, law)
  end function numeric_flux

  logical function inside(p, s, x0, y0, A, B)
    real(dp), intent(in) :: p, s, x0, y0, A, B
    inside = ((s * cos(p) - x0) / A)**2 + ((s * sin(p) - y0) / B)**2 < 1.0_dp
  end function inside

  !> The angle of the circle of radius s inside the ellipse
  !> ((x - x0)/A)^2 + ((y - y0)/B)^2 <= 1: sample M_ARC + 1 angles across
  !> the window the ellipse can occupy, bisect each sign change, and add up
  !> the parts of the intervals inside.
  real(dp) function arc_inside(s, x0, y0, A, B)
    real(dp), intent(in) :: s, x0, y0, A, B
    real(dp) :: d, lo, hi, h, half, pa, pb, a_, b_, m, cross
    logical :: in_a, in_b
    integer :: k, it
    d = hypot(x0, y0)
    lo = 0.0_dp
    hi = 2.0_dp * PI
    if (d > A) then
      half = asin(A / d)
      lo = atan2(y0, x0) - half
      hi = atan2(y0, x0) + half
    end if
    h = (hi - lo) / M_ARC
    arc_inside = 0.0_dp
    in_a = inside(lo, s, x0, y0, A, B)
    do k = 0, M_ARC - 1
      pa = lo + h * k
      pb = lo + h * (k + 1)
      in_b = inside(pb, s, x0, y0, A, B)
      if (in_a .and. in_b) then
        arc_inside = arc_inside + h
      else if (in_a .neqv. in_b) then
        a_ = pa
        b_ = pb
        do it = 1, 50
          m = 0.5_dp * (a_ + b_)
          if (inside(m, s, x0, y0, A, B) .eqv. in_a) then
            a_ = m
          else
            b_ = m
          end if
        end do
        cross = 0.5_dp * (a_ + b_)
        arc_inside = arc_inside + merge(cross - pa, pb - cross, in_a)
      end if
      in_a = in_b
    end do
  end function arc_inside

  !> oblate planet, by brute force
  real(dp) function numeric_flux_oblate(X, Y, r, f, theta, w, law)
    real(dp), intent(in) :: X, Y, r, f, theta, w(:)
    integer, intent(in) :: law
    real(dp) :: A, B, x0, y0, d, s1, s2, ds, s, occ
    integer :: j
    A = r / sqrt(1.0_dp - f)
    B = A * (1.0_dp - f)
    x0 = X * cos(theta) + Y * sin(theta)
    y0 = -X * sin(theta) + Y * cos(theta)
    d = hypot(x0, y0)
    s1 = max(0.0_dp, d - A)
    s2 = min(1.0_dp, d + A)
    numeric_flux_oblate = 1.0_dp
    if (s1 >= s2) return
    ds = (s2 - s1) / NS
    occ = 0.0_dp
    do j = 0, NS - 1
      s = s1 + (j + 0.5_dp) * ds
      occ = occ + s * ds * arc_inside(s, x0, y0, A, B) * intensity(sqrt(1.0_dp - s * s), w, law)
    end do
    numeric_flux_oblate = 1.0_dp - occ / num_total(w, law)
  end function numeric_flux_oblate

  !> Largest difference from brute force: their table of (r, b) with every
  !> contact (r <= 2), its oblate twin, then a light curve for r = 0.1
  !> (1000 points spherical, 201 oblate).
  real(dp) function numeric_checks() result(worst)
    real(dp), parameter :: W2(2) = [0.12123702682335066d0, 0.22967167533748584d0]
    real(dp), parameter :: W4(4) = [0.3709047608728925d0, 0.00711438176464726d0, &
        0.3162293789763985d0, 0.01888876666845794d0]
    real(dp), parameter :: W5(5) = [0.136725258105043d0, 0.25114827117024674d0, &
        0.1088217352833381d0, 0.030866335871229982d0, 0.32786772243688317d0]
    worst = max(num_table(W2, HYBRID2), num_table(W4, HYBRID4), num_table(W5, HYBRID5), &
                num_curve(W4, HYBRID4))
  end function numeric_checks

  real(dp) function num_table(w, law) result(worst)
    real(dp), intent(in) :: w(:)
    integer, intent(in) :: law
    real(dp), parameter :: RB(2, 20) = reshape([ &
        0.01d0, 0.0d0, 0.01d0, 0.01d0, 0.01d0, 0.99d0, 0.01d0, 1.0d0, 0.01d0, 1.01d0, &
        0.1d0, 0.0d0, 0.1d0, 0.1d0, 0.1d0, 0.9d0, 0.1d0, 1.0d0, 0.1d0, 1.1d0, &
        1.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0, 2.0d0, 2.0d0, 1.0d0, 2.0d0, 2.0d0, 2.0d0, 3.0d0, &
        0.5d0, 0.0d0, 0.5d0, 0.5d0, 0.5d0, 1.0d0, 0.5d0, 1.5d0], [2, 20])
    real(dp) :: X, Y
    integer :: i
    worst = 0.0_dp
    do i = 1, 16                    ! spherical: the first 16
      worst = max(worst, abs(hybrid_flux(RB(2, i), RB(1, i), w, law) &
                             - numeric_flux(RB(2, i), RB(1, i), w, law)))
    end do
    do i = 1, 20                    ! oblate: r < 1
      if (RB(1, i) >= 1.0_dp) cycle
      X = RB(2, i) * cos(NUM_BETA)
      Y = RB(2, i) * sin(NUM_BETA)
      worst = max(worst, abs(hybrid_flux_oblate(X, Y, RB(1, i), NUM_F, NUM_THETA, w, law) &
                             - numeric_flux_oblate(X, Y, RB(1, i), NUM_F, NUM_THETA, w, law)))
    end do
  end function num_table

  real(dp) function num_curve(w, law) result(worst)
    real(dp), intent(in) :: w(:)
    integer, intent(in) :: law
    real(dp) :: b, X
    integer :: i
    worst = 0.0_dp
    do i = 1, 1000
      b = abs(3.0_dp * (i / 1000.0_dp - 0.5_dp))
      worst = max(worst, abs(hybrid_flux(b, 0.1_dp, w, law) - numeric_flux(b, 0.1_dp, w, law)))
    end do
    do i = 0, 200
      X = -1.5_dp + 3.0_dp * i / 200.0_dp
      worst = max(worst, abs(hybrid_flux_oblate(X, 0.3_dp, 0.1_dp, NUM_F, NUM_THETA, w, law) &
                             - numeric_flux_oblate(X, 0.3_dp, 0.1_dp, NUM_F, NUM_THETA, w, law)))
    end do
  end function num_curve

end program test_limbdarkfast
