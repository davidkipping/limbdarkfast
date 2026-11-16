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
  if (failures == 0) then
    write (*, '(a)') "law shapes: ok"
  else
    error stop 1
  end if

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

end program test_limbdarkfast
