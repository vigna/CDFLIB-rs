! Reference tables for the gamma kernels: gamma_log, gamma, gamma_inc,
! gamma_inc_inv, dstrem.

program gen_gamma_kernels
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), external :: gamma_log, gamma_user, dstrem
  external :: gamma_inc, gamma_inc_inv

  real(kind=rk), parameter :: as(*) = (/ &
    0.25_rk, 0.5_rk, 0.9_rk, 1.0_rk, 1.5_rk, 2.5_rk, 5.0_rk, &
    9.5_rk, 12.0_rk, 17.0_rk, 25.0_rk, 50.0_rk, 100.0_rk, 500.0_rk /)
  real(kind=rk), parameter :: log_special(*) = (/ &
    0.149999999999_rk, 0.15_rk, 0.150000000001_rk, &
    0.374999999999_rk, 0.375_rk, 0.375000000001_rk, &
    0.799999999999_rk, 0.8_rk, 0.800000000001_rk, &
    2.249999999999_rk, 2.25_rk, 2.250000000001_rk, &
    9.999999999999_rk, 10.0_rk, 10.000000000001_rk /)
  real(kind=rk), parameter :: gamma_special(*) = (/ &
    14.999999999999_rk, 15.0_rk, 15.000000000001_rk /)
  ! Negative non-integer arguments: the |a| < 15 rows drive the
  ! fractional-part peel loop and the (x + 0.5) + 0.5 rounding (the
  ! -1.0e-8 row distinguishes it from a single x + 1.0 rounding); the
  ! |a| >= 15 rows drive the reflection branch.
  real(kind=rk), parameter :: gamma_negative(*) = (/ &
    -1.0e-8_rk, -0.3_rk, -0.5_rk, -0.9_rk, -0.999999_rk, &
    -1.5_rk, -2.000001_rk, -2.5_rk, -5.5_rk, -9.5_rk, &
    -13.5_rk, -14.5_rk, -15.5_rk, -20.5_rk, -42.5_rk /)
  real(kind=rk), parameter :: inc_special_a(*) = (/ &
    0.149999999999_rk, 0.15_rk, 0.150000000001_rk, &
    0.374999999999_rk, 0.375_rk, 0.375000000001_rk, &
    0.799999999999_rk, 0.8_rk, 0.800000000001_rk, &
    2.249999999999_rk, 2.25_rk, 2.250000000001_rk, &
    9.999999999999_rk, 10.0_rk, 10.000000000001_rk, &
    19.999999999999_rk, 20.0_rk, 20.000000000001_rk /)
  real(kind=rk), parameter :: inc_special_x(*) = (/ &
    0.149999999999_rk, 0.15_rk, 0.150000000001_rk, &
    0.374999999999_rk, 0.375_rk, 0.375000000001_rk, &
    0.799999999999_rk, 0.8_rk, 0.800000000001_rk, &
    1.999999999999_rk, 2.0_rk, 2.000000000001_rk, &
    9.999999999999_rk, 10.0_rk, 10.000000000001_rk, &
    19.999999999999_rk, 20.0_rk, 20.000000000001_rk /)
  integer :: i, j, unit, ind
  real(kind=rk) :: a, x, p, q, max_x, step

  ! gamma_log
  open(newunit=unit, file='tests/data/gamma_log.csv', status='replace', action='write')
  write(unit, '(a)') '# a, ln_gamma(a)'
  a = 0.05_rk
  do while (a <= 20.0_rk + 1.0e-12_rk)
    call putval(unit, a, .false.)
    call putval(unit, gamma_log(a), .true.)
    a = a + 0.05_rk
  end do
  do i = 1, size(log_special)
    a = log_special(i)
    call putval(unit, a, .false.)
    call putval(unit, gamma_log(a), .true.)
  end do
  close(unit)

  ! gamma
  open(newunit=unit, file='tests/data/gamma.csv', status='replace', action='write')
  write(unit, '(a)') '# a, gamma(a)'
  a = 0.1_rk
  do while (a <= 14.0_rk + 1.0e-12_rk)
    call putval(unit, a, .false.)
    call putval(unit, gamma_user(a), .true.)
    a = a + 0.1_rk
  end do
  do i = 1, size(gamma_special)
    a = gamma_special(i)
    call putval(unit, a, .false.)
    call putval(unit, gamma_user(a), .true.)
  end do
  do i = 1, size(gamma_negative)
    a = gamma_negative(i)
    call putval(unit, a, .false.)
    call putval(unit, gamma_user(a), .true.)
  end do
  close(unit)

  ! gamma_inc at all three accuracy regimes (ind = 0, 1, 2 in CDFLIB);
  ! each ind writes to its own CSV so Rust can read the right fixture per
  ! GammaIncAcc variant.
  call write_gamma_inc('tests/data/gamma_inc.csv', 0)
  call write_gamma_inc('tests/data/gamma_inc_d6.csv', 1)
  call write_gamma_inc('tests/data/gamma_inc_d3.csv', 2)

  call write_gamma_inc_inv()
  call write_dstrem()

  write(0, '(a)') 'wrote 7 tables under tests/data/'

contains
  subroutine write_gamma_inc(path, ind_arg)
    character(len=*), intent(in) :: path
    integer, intent(in) :: ind_arg
    integer :: ii, jj, u
    real(kind=rk) :: aa, xx, pp, qq, mx, st
    open(newunit=u, file=path, status='replace', action='write')
    write(u, '(a)') '# a, x, P(a, x), Q(a, x)'
    do ii = 1, size(as)
      aa = as(ii)
      if (aa > 1.0_rk) then
        mx = 3.0_rk * aa
      else
        mx = 8.0_rk
      end if
      st = mx / 50.0_rk
      xx = st / 2.0_rk
      do while (xx <= mx + 1.0e-12_rk)
        pp = 0.0_rk
        qq = 0.0_rk
        call gamma_inc(aa, xx, pp, qq, ind_arg)
        if (pp /= 2.0_rk) then  ! 2.0 is CDFLIB's error sentinel
          call putval(u, aa, .false.)
          call putval(u, xx, .false.)
          call putval(u, pp, .false.)
          call putval(u, qq, .true.)
        end if
        xx = xx + st
      end do
    end do
    do ii = 1, size(inc_special_a)
      aa = inc_special_a(ii)
      do jj = 1, size(inc_special_x)
        xx = inc_special_x(jj)
        pp = 0.0_rk
        qq = 0.0_rk
        call gamma_inc(aa, xx, pp, qq, ind_arg)
        if (pp /= 2.0_rk) then
          call putval(u, aa, .false.)
          call putval(u, xx, .false.)
          call putval(u, pp, .false.)
          call putval(u, qq, .true.)
        end if
      end do
    end do
    close(u)
  end subroutine write_gamma_inc

  subroutine write_gamma_inc_inv()
    integer :: u, ii, ip, ierr
    real(kind=rk) :: a, p, qq, xx
    ! a values spanning the initial-approximation branches:
    ! small a (< 1), moderate a (1..20), large a (> 20).
    real(kind=rk), parameter :: inv_as(*) = (/ &
      0.1_rk, 0.25_rk, 0.5_rk, 0.75_rk, 0.9_rk, &
      1.0_rk, 1.5_rk, 2.0_rk, 3.0_rk, 5.0_rk, &
      10.0_rk, 15.0_rk, 20.0_rk, 25.0_rk, 50.0_rk, &
      100.0_rk, 500.0_rk /)
    ! p values covering both tails and the body.
    real(kind=rk), parameter :: ps(*) = (/ &
      1.0e-12_rk, 1.0e-8_rk, 1.0e-4_rk, 0.001_rk, 0.01_rk, &
      0.05_rk, 0.1_rk, 0.2_rk, 0.3_rk, 0.4_rk, 0.5_rk, &
      0.6_rk, 0.7_rk, 0.8_rk, 0.9_rk, 0.95_rk, 0.99_rk, &
      0.999_rk, 0.9999_rk, 0.99999999_rk, 0.999999999999_rk /)

    open(newunit=u, file='tests/data/gamma_inc_inv.csv', &
         status='replace', action='write')
    write(u, '(a)') '# a, p, q, x, ierr'
    do ii = 1, size(inv_as)
      a = inv_as(ii)
      do ip = 1, size(ps)
        p = ps(ip)
        qq = 1.0_rk - p
        xx = 0.0_rk
        ierr = 0
        call gamma_inc_inv(a, xx, -1.0_rk, p, qq, ierr)
        if (ierr >= 0) then
          call putval(u, a, .false.)
          call putval(u, p, .false.)
          call putval(u, qq, .false.)
          call putval(u, xx, .false.)
          call putint(u, ierr, .true.)
        end if
      end do
    end do
    close(u)
  end subroutine write_gamma_inc_inv

  subroutine write_dstrem()
    integer :: u
    real(kind=rk) :: z
    open(newunit=u, file='tests/data/dstrem.csv', &
         status='replace', action='write')
    write(u, '(a)') '# z, dstrem(z)'
    ! Fine grid across the branch point at z = 6.
    z = 0.1_rk
    do while (z <= 20.0_rk + 1.0e-12_rk)
      call putval(u, z, .false.)
      call putval(u, dstrem(z), .true.)
      z = z + 0.1_rk
    end do
    ! Large z: asymptotic regime.
    z = 50.0_rk
    do while (z <= 500.0_rk + 1.0e-12_rk)
      call putval(u, z, .false.)
      call putval(u, dstrem(z), .true.)
      z = z + 50.0_rk
    end do
    close(u)
  end subroutine write_dstrem

  subroutine putint(unit, v, last)
    integer, intent(in) :: unit, v
    logical, intent(in) :: last
    character(len=16) :: buf
    write(buf, '(i0)') v
    if (last) then
      write(unit, '(a)') trim(adjustl(buf))
    else
      write(unit, '(a, a)', advance='no') trim(adjustl(buf)), ','
    end if
  end subroutine putint

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
end program gen_gamma_kernels
