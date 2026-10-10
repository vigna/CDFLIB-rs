! Reference tables for the remaining special-function kernels:
! psi, dlanor, dt1, stvaln.
!
! Build and run via tests/regenerate/regenerate.sh.

program gen_psi_dt1_kernels
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer, parameter :: rk = kind(1.0d0)

  real(kind=rk), external :: psi, dlanor, dt1, stvaln

  call gen_psi()
  call gen_dlanor()
  call gen_dt1()
  call gen_stvaln()

  write(0, '(a)') 'wrote 4 tables under tests/data/'

contains

  subroutine putval(unit, v, last)
    integer, intent(in) :: unit
    real(kind=rk), intent(in) :: v
    logical, intent(in) :: last
    character(len=32) :: buf
    write(buf, '(es25.17e3)') v
    if (last) then
      write(unit, '(a)') trim(adjustl(buf))
    else
      write(unit, '(a, a)', advance='no') trim(adjustl(buf)), ','
    end if
  end subroutine putval

  subroutine gen_psi()
    integer :: unit
    real(kind=rk) :: x
    open(newunit=unit, file='tests/data/psi.csv', status='replace', action='write')
    write(unit, '(a)') '# x, psi(x)'

    ! Positive arguments: cover rational-approx regime (x < 3) and
    ! reduction-to-asymptotic regime (x >= 3).
    x = 0.0625_rk
    do while (x <= 20.0_rk + 1.0e-12_rk)
      call putval(unit, x, .false.)
      call putval(unit, psi(x), .true.)
      x = x + 0.0625_rk
    end do
    ! Large x: asymptotic regime.
    x = 50.0_rk
    do while (x <= 1000.0_rk + 1.0e-12_rk)
      call putval(unit, x, .false.)
      call putval(unit, psi(x), .true.)
      x = x + 50.0_rk
    end do
    ! Negative non-integer: reflection formula. Skip integers (poles).
    x = -0.9375_rk
    do while (x >= -10.0_rk - 1.0e-12_rk)
      call putval(unit, x, .false.)
      call putval(unit, psi(x), .true.)
      x = x - 0.0625_rk
      ! Skip values too close to non-positive integers (poles).
      if (abs(x - nint(x)) < 0.001_rk) x = x - 0.0625_rk
    end do
    ! Arguments below 0.5 that are not multiples of 1/16, for which the
    ! reflection formula takes the sine and cosine of arbitrary arguments.
    x = 0.4999_rk
    do while (x >= -10.0_rk)
      if (abs(x - nint(x)) >= 0.001_rk) then
        call putval(unit, x, .false.)
        call putval(unit, psi(x), .true.)
      end if
      x = x - 0.0137_rk
    end do

    close(unit)
  end subroutine gen_psi

  subroutine gen_dlanor()
    integer :: unit
    real(kind=rk) :: x
    open(newunit=unit, file='tests/data/dlanor.csv', status='replace', action='write')
    write(unit, '(a)') '# x, dlanor(x)'
    ! dlanor requires |x| >= 5.
    x = 5.0_rk
    do while (x <= 38.0_rk + 1.0e-12_rk)
      call putval(unit, x, .false.)
      call putval(unit, dlanor(x), .true.)
      x = x + 0.25_rk
    end do
    ! Also test negative x.
    x = -5.0_rk
    do while (x >= -38.0_rk - 1.0e-12_rk)
      call putval(unit, x, .false.)
      call putval(unit, dlanor(x), .true.)
      x = x - 0.25_rk
    end do
    close(unit)
  end subroutine gen_dlanor

  subroutine gen_dt1()
    integer :: unit, ip, idf
    real(kind=rk) :: p, q, df, nan, inf
    real(kind=rk), parameter :: ps(*) = (/ &
      0.001_rk, 0.01_rk, 0.025_rk, 0.05_rk, 0.1_rk, 0.2_rk, &
      0.3_rk, 0.4_rk, 0.5_rk, 0.6_rk, 0.7_rk, 0.8_rk, 0.9_rk, &
      0.95_rk, 0.975_rk, 0.99_rk, 0.999_rk /)
    real(kind=rk), parameter :: dfs(*) = (/ &
      1.0_rk, 2.0_rk, 3.0_rk, 5.0_rk, 10.0_rk, 20.0_rk, &
      50.0_rk, 100.0_rk, 1000.0_rk, 1.0e6_rk /)

    open(newunit=unit, file='tests/data/dt1.csv', status='replace', action='write')
    write(unit, '(a)') '# p, q, df, dt1(p, q, df)'
    do idf = 1, size(dfs)
      df = dfs(idf)
      do ip = 1, size(ps)
        p = ps(ip)
        q = 1.0_rk - p
        call putval(unit, p, .false.)
        call putval(unit, q, .false.)
        call putval(unit, df, .false.)
        call putval(unit, dt1(p, q, df), .true.)
      end do
    end do
    ! NaN at the endpoints, where dinvnr is NaN, and for NaN arguments;
    ! the normal quantile for an infinite df; the small-df examples of the
    ! documentation (a wrong sign at df = 0.1, and NaN from infinite terms
    ! of opposite sign at df = 1e-300).
    nan = ieee_value(nan, ieee_quiet_nan)
    inf = ieee_value(inf, ieee_positive_inf)
    call dt1_row(unit, 0.0_rk, 1.0_rk, 10.0_rk)
    call dt1_row(unit, 1.0_rk, 0.0_rk, 10.0_rk)
    call dt1_row(unit, nan, nan, 10.0_rk)
    call dt1_row(unit, 0.3_rk, 0.7_rk, nan)
    call dt1_row(unit, 0.3_rk, 0.7_rk, inf)
    call dt1_row(unit, 0.6_rk, 0.4_rk, 0.1_rk)
    call dt1_row(unit, 0.3_rk, 0.7_rk, 1.0e-300_rk)
    close(unit)
  end subroutine gen_dt1

  subroutine dt1_row(unit, p, q, df)
    integer, intent(in) :: unit
    real(kind=rk), intent(in) :: p, q, df
    call putval(unit, p, .false.)
    call putval(unit, q, .false.)
    call putval(unit, df, .false.)
    call putval(unit, dt1(p, q, df), .true.)
  end subroutine dt1_row

  subroutine gen_stvaln()
    integer :: unit
    real(kind=rk) :: p, nan
    open(newunit=unit, file='tests/data/stvaln.csv', status='replace', action='write')
    write(unit, '(a)') '# p, stvaln(p)'
    p = 0.001_rk
    do while (p <= 0.999_rk + 1.0e-12_rk)
      call putval(unit, p, .false.)
      call putval(unit, stvaln(p), .true.)
      p = p + 0.001_rk
    end do
    ! NaN at 0 and 1, where log(0) makes the rational function inf/inf,
    ! and at NaN.
    nan = ieee_value(nan, ieee_quiet_nan)
    call putval(unit, 0.0_rk, .false.)
    call putval(unit, stvaln(0.0_rk), .true.)
    call putval(unit, 1.0_rk, .false.)
    call putval(unit, stvaln(1.0_rk), .true.)
    call putval(unit, nan, .false.)
    call putval(unit, stvaln(nan), .true.)
    close(unit)
  end subroutine gen_stvaln

end program gen_psi_dt1_kernels
