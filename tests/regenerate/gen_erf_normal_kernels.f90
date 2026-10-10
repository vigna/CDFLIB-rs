! Reference tables for the erf and standard-normal kernels:
! error_f, error_fc, error_fc_scaled, cumnor, dinvnr. The tables end
! with rows at NaN arguments and, except for dinvnr, at infinite and
! huge ones, where the documentation states the F90 result.
!
! Build and run via tests/regenerate/regenerate.sh.

program gen_erf_normal_kernels
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer, parameter :: rk = kind(1.0d0)

  real(kind=rk), external :: error_f, error_fc, dinvnr
  external :: cumnor

  call gen_error_f()
  call gen_error_fc()
  call gen_error_fc_scaled()
  call gen_cumnor()
  call gen_dinvnr()

  write(0, '(a)') 'wrote 5 tables under tests/data/'

contains

  ! Write a single double-precision value to unit followed by a comma,
  ! or by a newline if last is true. Uses 18 significant digits so each
  ! distinct f64 round-trips through decimal.
  subroutine putval(unit, v, last)
    integer, intent(in) :: unit
    real(kind=rk), intent(in) :: v
    logical, intent(in) :: last
    character(len=32) :: buf
    ! es25.17e3 = sign + digit + '.' + 17 frac + 'E' + signed-3 = 25 chars
    write(buf, '(es25.17e3)') v
    if (last) then
      write(unit, '(a)') trim(adjustl(buf))
    else
      write(unit, '(a, a)', advance='no') trim(adjustl(buf)), ','
    end if
  end subroutine putval

  ! NaN, plus and minus infinity, and plus and minus huge.
  subroutine specials(v)
    real(kind=rk), intent(out) :: v(5)
    v(1) = ieee_value(v(1), ieee_quiet_nan)
    v(2) = ieee_value(v(2), ieee_positive_inf)
    v(3) = ieee_value(v(3), ieee_negative_inf)
    v(4) = huge(1.0_rk)
    v(5) = -huge(1.0_rk)
  end subroutine specials

  subroutine gen_error_f()
    integer :: unit, i
    real(kind=rk) :: x, v(5)
    open(newunit=unit, file='tests/data/error_f.csv', status='replace', action='write')
    write(unit, '(a)') '# x, erf(x)'
    x = -6.0_rk
    do while (x <= 6.0_rk + 1.0e-12_rk)
      call putval(unit, x, .false.)
      call putval(unit, error_f(x), .true.)
      x = x + 0.0625_rk
    end do
    ! A NaN x falls through every branch to 1 (cdflib.f90:9440) and is not
    ! negated.
    call specials(v)
    do i = 1, size(v)
      call putval(unit, v(i), .false.)
      call putval(unit, error_f(v(i)), .true.)
    end do
    close(unit)
  end subroutine gen_error_f

  subroutine gen_error_fc()
    integer :: unit, ind0, i
    real(kind=rk) :: x, v(5)
    ind0 = 0
    open(newunit=unit, file='tests/data/error_fc.csv', status='replace', action='write')
    write(unit, '(a)') '# x, erfc(x)'
    x = -6.0_rk
    do while (x <= 30.0_rk + 1.0e-12_rk)
      call putval(unit, x, .false.)
      call putval(unit, error_fc(ind0, x), .true.)
      x = x + 0.0625_rk
    end do
    call specials(v)
    do i = 1, size(v)
      call putval(unit, v(i), .false.)
      call putval(unit, error_fc(ind0, v(i)), .true.)
    end do
    close(unit)
  end subroutine gen_error_fc

  subroutine gen_error_fc_scaled()
    integer :: unit, ind1, i
    real(kind=rk) :: x, v(5)
    ind1 = 1
    open(newunit=unit, file='tests/data/error_fc_scaled.csv', status='replace', action='write')
    write(unit, '(a)') '# x, erfc(x)*exp(x^2)'
    x = -6.0_rk
    do while (x <= 60.0_rk + 1.0e-12_rk)
      call putval(unit, x, .false.)
      call putval(unit, error_fc(ind1, x), .true.)
      x = x + 0.0625_rk
    end do
    call specials(v)
    do i = 1, size(v)
      call putval(unit, v(i), .false.)
      call putval(unit, error_fc(ind1, v(i)), .true.)
    end do
    close(unit)
  end subroutine gen_error_fc_scaled

  subroutine gen_cumnor()
    integer :: unit, i
    real(kind=rk) :: x, cum, ccum, v(11), m16
    open(newunit=unit, file='tests/data/cumnor.csv', status='replace', action='write')
    write(unit, '(a)') '# x, cum, ccum'
    x = -38.0_rk
    do while (x <= 38.0_rk + 1.0e-12_rk)
      call cumnor(x, cum, ccum)
      call putval(unit, x, .false.)
      call putval(unit, cum, .false.)
      call putval(unit, ccum, .true.)
      x = x + 0.0625_rk
    end do
    ! The result is (NaN, NaN) for a NaN x and where 16 |x| overflows in
    ! aint ( x * sixten ) (cdflib.f90:7774): around huge/16 and beyond.
    call specials(v(1:5))
    m16 = huge(1.0_rk) / 16.0_rk
    v(6) = m16
    v(7) = nearest(m16, 1.0_rk)
    v(8) = nearest(m16, -1.0_rk)
    v(9) = -m16
    v(10) = -nearest(m16, 1.0_rk)
    v(11) = -nearest(m16, -1.0_rk)
    do i = 1, size(v)
      call cumnor(v(i), cum, ccum)
      call putval(unit, v(i), .false.)
      call putval(unit, cum, .false.)
      call putval(unit, ccum, .true.)
    end do
    close(unit)
  end subroutine gen_cumnor

  subroutine gen_dinvnr()
    ! Parametrize by x rather than p; store (cum, ccum, back). Skips rows
    ! where cum/ccum has saturated to 0 or 1.
    integer :: unit
    real(kind=rk) :: x, cum, ccum, back, nan
    open(newunit=unit, file='tests/data/dinvnr.csv', status='replace', action='write')
    write(unit, '(a)') '# p, q, x'
    x = -7.0_rk
    do while (x <= 7.0_rk + 1.0e-12_rk)
      call cumnor(x, cum, ccum)
      if (cum > 0.0_rk .and. cum < 1.0_rk .and. ccum > 0.0_rk .and. ccum < 1.0_rk) then
        back = dinvnr(cum, ccum)
        call putval(unit, cum, .false.)
        call putval(unit, ccum, .false.)
        call putval(unit, back, .true.)
      end if
      x = x + 0.0625_rk
    end do
    ! The endpoints, where the Newton steps never converge and the result
    ! is the NaN start value (cdflib.f90:8128, :8130).
    call putval(unit, 0.0_rk, .false.)
    call putval(unit, 1.0_rk, .false.)
    call putval(unit, dinvnr(0.0_rk, 1.0_rk), .true.)
    call putval(unit, 1.0_rk, .false.)
    call putval(unit, 0.0_rk, .false.)
    call putval(unit, dinvnr(1.0_rk, 0.0_rk), .true.)
    ! NaN arguments. With both NaN the start value and every Newton step
    ! are NaN. With one NaN, gfortran's MIN returns the other argument
    ! (the Fortran standard leaves it to the processor), and the final
    ! test p <= q is false.
    nan = ieee_value(nan, ieee_quiet_nan)
    call putval(unit, nan, .false.)
    call putval(unit, nan, .false.)
    call putval(unit, dinvnr(nan, nan), .true.)
    call putval(unit, nan, .false.)
    call putval(unit, 0.3_rk, .false.)
    call putval(unit, dinvnr(nan, 0.3_rk), .true.)
    call putval(unit, 0.3_rk, .false.)
    call putval(unit, nan, .false.)
    call putval(unit, dinvnr(0.3_rk, nan), .true.)
    close(unit)
  end subroutine gen_dinvnr

end program gen_erf_normal_kernels
