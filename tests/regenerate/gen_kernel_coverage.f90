! Reference tables that drive every branch of the low-level kernels of
! cdflib.f90, plus the error exits of gamma_user, psi, gamma_inc,
! gamma_inc_inv and beta_inc.
!
! Each grid is chosen from the branch structure of the routine it tests;
! the comment above each block names the branches it reaches.
! Companion Rust test: tests/kernel_coverage.rs.

program gen_kernel_coverage
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), external :: algdiv, alnrel, apser, bcorr, beta, beta_asym, &
    beta_frac, beta_pser, beta_rcomp, beta_rcomp1, beta_up, dexpm1, esum, &
    exparg, fpser, gam1, gamma_ln1, gamma_user, gsumln, psi, rcomp, rexp, &
    rlog, rlog1
  external :: beta_grat, beta_inc, gamma_inc, gamma_inc_inv, gamma_rat1

  ! (a, b) pairs and x values shared by beta_rcomp and beta_rcomp1.
  integer, parameter :: nab = 15
  real(kind=rk), parameter :: rab(2, nab) = reshape((/ &
    0.3_rk, 0.5_rk,   0.6_rk, 0.8_rk,   0.5_rk, 3.0_rk,   0.2_rk, 1.5_rk, &
    0.3_rk, 20.0_rk,  0.9_rk, 7.9_rk,   2.0_rk, 5.0_rk,   5.0_rk, 2.0_rk, &
    10.0_rk, 20.0_rk, 20.0_rk, 10.0_rk, 8.0_rk, 8.0_rk,   100.0_rk, 1000.0_rk, &
    1000.0_rk, 100.0_rk, 1.0e-5_rk, 0.5_rk, 0.8_rk, 0.9_rk /), (/ 2, nab /))
  integer, parameter :: nx = 9
  real(kind=rk), parameter :: rx(nx) = (/ 0.0_rk, 0.001_rk, 0.1_rk, 0.375_rk, &
    0.5_rk, 0.625_rk, 0.7_rk, 0.99_rk, 1.0_rk /)
  ! Arguments shared by dexpm1 and rexp.
  real(kind=rk), parameter :: ev(13) = (/ -50.0_rk, -1.0_rk, -0.15_rk, -0.1_rk, &
    -1.0e-10_rk, 0.0_rk, 1.0e-10_rk, 0.1_rk, 0.15_rk, 0.16_rk, 1.0_rk, 50.0_rk, &
    700.0_rk /)

  call gen_algdiv()
  call gen_alnrel()
  call gen_apser()
  call gen_bcorr()
  call gen_beta()
  call gen_beta_asym()
  call gen_beta_frac()
  call gen_beta_grat()
  call gen_beta_pser()
  call gen_beta_rcomp()
  call gen_beta_rcomp1()
  call gen_beta_up()
  call gen_dexpm1()
  call gen_esum()
  call gen_exparg()
  call gen_fpser()
  call gen_gam1()
  call gen_gamma_ln1()
  call gen_gamma_rat1()
  call gen_gsumln()
  call gen_rcomp()
  call gen_rexp()
  call gen_rlog()
  call gen_rlog1()
  call gen_gamma_edge()
  call gen_psi_edge()
  call gen_gamma_inc_edge()
  call gen_gamma_inc_inv_edge()
  call gen_beta_inc_regimes()

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

  subroutine openf(unit, name, header)
    integer, intent(out) :: unit
    character(len=*), intent(in) :: name, header
    open(newunit=unit, file='tests/data/'//name, status='replace', action='write')
    write(unit, '(a)') header
  end subroutine openf

  ! algdiv requires 8 <= b; covers b < a and a <= b, and v < u.
  subroutine gen_algdiv()
    real(kind=rk), parameter :: av(10) = (/ 0.5_rk, 1.0_rk, 2.0_rk, 5.0_rk, &
      7.9_rk, 8.0_rk, 10.0_rk, 30.0_rk, 100.0_rk, 1.0e5_rk /)
    real(kind=rk), parameter :: bv(7) = (/ 8.0_rk, 9.5_rk, 10.0_rk, 20.0_rk, &
      50.0_rk, 1.0e3_rk, 1.0e8_rk /)
    integer :: unit, i, j
    call openf(unit, 'algdiv.csv', '# a, b, algdiv(a, b)')
    do i = 1, size(av)
      do j = 1, size(bv)
        call putval(unit, av(i), .false.)
        call putval(unit, bv(j), .false.)
        call putval(unit, algdiv(av(i), bv(j)), .true.)
      end do
    end do
    close(unit)
  end subroutine gen_algdiv

  ! alnrel: |a| <= 0.375 and its complement.
  subroutine gen_alnrel()
    real(kind=rk), parameter :: av(14) = (/ -0.999_rk, -0.5_rk, -0.375_rk, &
      -0.2_rk, -1.0e-10_rk, 0.0_rk, 1.0e-10_rk, 0.1_rk, 0.375_rk, 0.376_rk, &
      1.0_rk, 10.0_rk, 1.0e10_rk, 1.0e300_rk /)
    integer :: unit, i
    call openf(unit, 'alnrel.csv', '# a, alnrel(a)')
    do i = 1, size(av)
      call putval(unit, av(i), .false.)
      call putval(unit, alnrel(av(i)), .true.)
    end do
    close(unit)
  end subroutine gen_alnrel

  ! apser, used by beta_inc when a <= min(eps, eps*b) and b*x <= 1:
  ! covers b*eps <= 0.02 and 0.02 < b*eps.
  subroutine gen_apser()
    real(kind=rk), parameter :: av(3) = (/ 1.0e-20_rk, 1.0e-17_rk, 1.0e-16_rk /)
    real(kind=rk), parameter :: bv(8) = (/ 0.5_rk, 1.0_rk, 2.0_rk, 10.0_rk, &
      100.0_rk, 1.0e5_rk, 1.0e14_rk, 1.0e16_rk /)
    real(kind=rk), parameter :: xv(7) = (/ 1.0e-17_rk, 1.0e-12_rk, 1.0e-6_rk, &
      1.0e-3_rk, 0.01_rk, 0.1_rk, 0.5_rk /)
    real(kind=rk), parameter :: eps = 1.0e-15_rk
    integer :: unit, i, j, k
    call openf(unit, 'apser.csv', '# a, b, x, eps, apser(a, b, x, eps)')
    do i = 1, size(av)
      do j = 1, size(bv)
        do k = 1, size(xv)
          if (1.0_rk < bv(j) * xv(k)) cycle
          call putval(unit, av(i), .false.)
          call putval(unit, bv(j), .false.)
          call putval(unit, xv(k), .false.)
          call putval(unit, eps, .false.)
          call putval(unit, apser(av(i), bv(j), xv(k), eps), .true.)
        end do
      end do
    end do
    close(unit)
  end subroutine gen_apser

  ! bcorr requires 8 <= a0 and 8 <= b0.
  subroutine gen_bcorr()
    real(kind=rk), parameter :: v(8) = (/ 8.0_rk, 9.0_rk, 10.0_rk, 15.0_rk, &
      50.0_rk, 100.0_rk, 1.0e5_rk, 1.0e10_rk /)
    integer :: unit, i, j
    call openf(unit, 'bcorr.csv', '# a0, b0, bcorr(a0, b0)')
    do i = 1, size(v)
      do j = 1, size(v)
        call putval(unit, v(i), .false.)
        call putval(unit, v(j), .false.)
        call putval(unit, bcorr(v(i), v(j)), .true.)
      end do
    end do
    close(unit)
  end subroutine gen_bcorr

  ! The complete beta function.
  subroutine gen_beta()
    real(kind=rk), parameter :: v(8) = (/ 0.1_rk, 0.5_rk, 1.0_rk, 2.5_rk, &
      10.0_rk, 50.0_rk, 100.0_rk, 200.0_rk /)
    integer :: unit, i, j
    call openf(unit, 'beta.csv', '# a, b, beta(a, b)')
    do i = 1, size(v)
      do j = 1, size(v)
        call putval(unit, v(i), .false.)
        call putval(unit, v(j), .false.)
        call putval(unit, beta(v(i), v(j)), .true.)
      end do
    end do
    close(unit)
  end subroutine gen_beta

  ! beta_asym requires 15 <= a, b: covers a < b and b <= a, the t == 0
  ! early exit, convergence inside the loop and after it.
  subroutine gen_beta_asym()
    real(kind=rk), parameter :: ab(2, 6) = reshape((/ 15.0_rk, 15.0_rk, &
      20.0_rk, 100.0_rk, 100.0_rk, 20.0_rk, 1000.0_rk, 1000.0_rk, &
      1.0e5_rk, 2.0e5_rk, 1.0e8_rk, 1.0e8_rk /), (/ 2, 6 /))
    real(kind=rk), parameter :: fr(6) = (/ 0.0_rk, 0.001_rk, 0.01_rk, 0.1_rk, &
      0.5_rk, 0.9_rk /)
    real(kind=rk), parameter :: eps = 100.0e-15_rk
    real(kind=rk) :: a, b, lambda
    integer :: unit, i, k
    call openf(unit, 'beta_asym.csv', '# a, b, lambda, eps, beta_asym(a, b, lambda, eps)')
    do i = 1, size(ab, 2)
      a = ab(1, i)
      b = ab(2, i)
      do k = 1, size(fr)
        lambda = fr(k) * min(a, b)
        call putval(unit, a, .false.)
        call putval(unit, b, .false.)
        call putval(unit, lambda, .false.)
        call putval(unit, eps, .false.)
        call putval(unit, beta_asym(a, b, lambda, eps), .true.)
      end do
    end do
    close(unit)
  end subroutine gen_beta_asym

  ! beta_frac with the arguments beta_inc passes at label 130: 1 < a0, b0,
  ! lambda >= 0 after the swap at label 40, and 40 <= b0. Covers the
  ! beta_rcomp == 0 early exit.
  subroutine gen_beta_frac()
    real(kind=rk), parameter :: v(5) = (/ 1.5_rk, 5.0_rk, 40.0_rk, 100.0_rk, 500.0_rk /)
    real(kind=rk), parameter :: xv(5) = (/ 0.001_rk, 0.1_rk, 0.5_rk, 0.9_rk, 0.999_rk /)
    real(kind=rk), parameter :: eps = 15.0e-15_rk
    real(kind=rk) :: a, b, x, y, lambda, t
    integer :: unit, i, j, k
    call openf(unit, 'beta_frac.csv', '# a, b, x, y, lambda, eps, beta_frac(a, b, x, y, lambda, eps)')
    do i = 1, size(v)
      do j = 1, size(v)
        do k = 1, size(xv)
          a = v(i)
          b = v(j)
          x = xv(k)
          y = 1.0_rk - x
          if (a <= b) then
            lambda = a - (a + b) * x
          else
            lambda = (a + b) * y - b
          end if
          if (lambda < 0.0_rk) then
            t = a
            a = b
            b = t
            t = x
            x = y
            y = t
            lambda = abs(lambda)
          end if
          if (b < 40.0_rk) cycle
          call putval(unit, a, .false.)
          call putval(unit, b, .false.)
          call putval(unit, x, .false.)
          call putval(unit, y, .false.)
          call putval(unit, lambda, .false.)
          call putval(unit, eps, .false.)
          call putval(unit, beta_frac(a, b, x, y, lambda, eps), .true.)
        end do
      end do
    end do
    close(unit)
  end subroutine gen_beta_frac

  ! beta_grat requires 15 <= a and b <= 1: covers y <= 0.375 and its
  ! complement, convergence inside and after the loop, and the three
  ! ierr = 1 exits (b*z == 0 at y == 0, u == 0 for huge a, and a
  ! nonpositive sum).
  subroutine gen_beta_grat()
    real(kind=rk), parameter :: av(6) = (/ 15.0_rk, 20.0_rk, 50.0_rk, 100.0_rk, &
      1000.0_rk, 1.0e5_rk /)
    real(kind=rk), parameter :: bv(5) = (/ 0.01_rk, 0.1_rk, 0.5_rk, 0.9_rk, 1.0_rk /)
    real(kind=rk), parameter :: xv(9) = (/ 0.01_rk, 0.1_rk, 0.3_rk, 0.5_rk, &
      0.625_rk, 0.8_rk, 0.95_rk, 0.99_rk, 0.999_rk /)
    real(kind=rk), parameter :: wv(2) = (/ 0.0_rk, 0.3_rk /)
    real(kind=rk), parameter :: eps = 15.0e-15_rk
    integer :: unit, i, j, k, l
    call openf(unit, 'beta_grat.csv', '# a, b, x, y, w, eps, w_out, ierr')
    do i = 1, size(av)
      do j = 1, size(bv)
        do k = 1, size(xv)
          do l = 1, size(wv)
            call grat_row(unit, av(i), bv(j), xv(k), 1.0_rk - xv(k), wv(l), eps)
          end do
        end do
      end do
    end do
    call grat_row(unit, 20.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.0_rk, eps)
    call grat_row(unit, 1.0e5_rk, 1.0_rk, 0.01_rk, 0.99_rk, 0.0_rk, eps)
    call grat_row(unit, 1.0e5_rk, 0.5_rk, 0.001_rk, 0.999_rk, 0.25_rk, eps)
    call grat_row(unit, 15.0_rk, 0.999_rk, 1.0e-300_rk, 1.0_rk, 0.0_rk, eps)
    close(unit)
  end subroutine gen_beta_grat

  subroutine grat_row(unit, a, b, x, y, w, eps)
    integer, intent(in) :: unit
    real(kind=rk), intent(in) :: a, b, x, y, w, eps
    real(kind=rk) :: wo
    integer :: ierr
    wo = w
    ierr = 0
    call beta_grat(a, b, x, y, wo, eps, ierr)
    call putval(unit, a, .false.)
    call putval(unit, b, .false.)
    call putval(unit, x, .false.)
    call putval(unit, y, .false.)
    call putval(unit, w, .false.)
    call putval(unit, eps, .false.)
    call putval(unit, wo, .false.)
    call putval(unit, real(ierr, kind=rk), .true.)
  end subroutine grat_row

  ! beta_pser: x == 0, 1 <= a0, and for a0 < 1 the three b0 regimes with
  ! both apb branches, a zero factor and the a <= 0.1*eps exit.
  subroutine gen_beta_pser()
    real(kind=rk), parameter :: ab(2, 16) = reshape((/ &
      0.5_rk, 0.5_rk,   0.3_rk, 0.9_rk,   0.2_rk, 0.6_rk,   0.7_rk, 0.8_rk, &
      0.5_rk, 3.0_rk,   0.8_rk, 5.0_rk,   0.3_rk, 1.5_rk,   0.5_rk, 10.0_rk, &
      0.5_rk, 100.0_rk, 2.0_rk, 3.0_rk,   1.0_rk, 1.0_rk,   5.0_rk, 50.0_rk, &
      1.0e-17_rk, 0.5_rk, 500.0_rk, 2.0_rk, 3.0_rk, 0.5_rk, 0.9_rk, 7.9_rk /), &
      (/ 2, 16 /))
    real(kind=rk), parameter :: xv(7) = (/ 0.0_rk, 1.0e-300_rk, 1.0e-10_rk, &
      0.01_rk, 0.1_rk, 0.3_rk, 0.7_rk /)
    real(kind=rk), parameter :: eps = 1.0e-15_rk
    integer :: unit, i, k
    call openf(unit, 'beta_pser.csv', '# a, b, x, eps, beta_pser(a, b, x, eps)')
    do i = 1, size(ab, 2)
      do k = 1, size(xv)
        call putval(unit, ab(1, i), .false.)
        call putval(unit, ab(2, i), .false.)
        call putval(unit, xv(k), .false.)
        call putval(unit, eps, .false.)
        call putval(unit, beta_pser(ab(1, i), ab(2, i), xv(k), eps), .true.)
      end do
    end do
    close(unit)
  end subroutine gen_beta_pser

  ! beta_rcomp and beta_rcomp1: x or y zero; a0 < 8 with the three lnx/lny
  ! choices and every a0/b0 regime; 8 <= a0 with a <= b and b < a and
  ! both sides of |e| <= 0.6.
  subroutine gen_beta_rcomp()
    integer :: unit, i, k
    call openf(unit, 'beta_rcomp.csv', '# a, b, x, y, beta_rcomp(a, b, x, y)')
    do i = 1, nab
      do k = 1, nx
        call putval(unit, rab(1, i), .false.)
        call putval(unit, rab(2, i), .false.)
        call putval(unit, rx(k), .false.)
        call putval(unit, 1.0_rk - rx(k), .false.)
        call putval(unit, beta_rcomp(rab(1, i), rab(2, i), rx(k), 1.0_rk - rx(k)), .true.)
      end do
    end do
    close(unit)
  end subroutine gen_beta_rcomp

  subroutine gen_beta_rcomp1()
    integer, parameter :: muv(4) = (/ -10, 0, 5, 700 /)
    integer :: unit, i, k, m
    call openf(unit, 'beta_rcomp1.csv', '# mu, a, b, x, y, beta_rcomp1(mu, a, b, x, y)')
    do m = 1, size(muv)
      do i = 1, nab
        do k = 1, nx
          call putval(unit, real(muv(m), kind=rk), .false.)
          call putval(unit, rab(1, i), .false.)
          call putval(unit, rab(2, i), .false.)
          call putval(unit, rx(k), .false.)
          call putval(unit, 1.0_rk - rx(k), .false.)
          call putval(unit, beta_rcomp1(muv(m), rab(1, i), rab(2, i), rx(k), &
            1.0_rk - rx(k)), .true.)
        end do
      end do
    end do
    close(unit)
  end subroutine gen_beta_rcomp1

  ! beta_up: n == 1 and 1 < n, a < 1 and 1 <= a, both apb branches, b <= 1
  ! and 1 < b with y <= 0.0001, r < 1, 1 <= r < n - 1 and n - 1 <= r, and
  ! convergence of the remaining terms.
  subroutine gen_beta_up()
    real(kind=rk), parameter :: av(5) = (/ 0.5_rk, 1.0_rk, 5.0_rk, 20.0_rk, 600.0_rk /)
    real(kind=rk), parameter :: bv(5) = (/ 0.5_rk, 1.5_rk, 10.0_rk, 100.0_rk, 1000.0_rk /)
    real(kind=rk), parameter :: xv(5) = (/ 0.001_rk, 0.3_rk, 0.6_rk, 0.9_rk, 0.99995_rk /)
    integer, parameter :: nv(3) = (/ 1, 5, 20 /)
    real(kind=rk), parameter :: eps = 1.0e-15_rk
    integer :: unit, i, j, k, l
    call openf(unit, 'beta_up.csv', '# a, b, x, y, n, eps, beta_up(a, b, x, y, n, eps)')
    do i = 1, size(av)
      do j = 1, size(bv)
        do k = 1, size(xv)
          do l = 1, size(nv)
            call putval(unit, av(i), .false.)
            call putval(unit, bv(j), .false.)
            call putval(unit, xv(k), .false.)
            call putval(unit, 1.0_rk - xv(k), .false.)
            call putval(unit, real(nv(l), kind=rk), .false.)
            call putval(unit, eps, .false.)
            call putval(unit, beta_up(av(i), bv(j), xv(k), 1.0_rk - xv(k), nv(l), &
              eps), .true.)
          end do
        end do
      end do
    end do
    close(unit)
  end subroutine gen_beta_up

  subroutine gen_dexpm1()
    integer :: unit, i
    call openf(unit, 'dexpm1.csv', '# x, dexpm1(x)')
    do i = 1, size(ev)
      call putval(unit, ev(i), .false.)
      call putval(unit, dexpm1(ev(i)), .true.)
    end do
    close(unit)
  end subroutine gen_dexpm1

  ! esum: every sign combination of mu and x and both signs of w.
  subroutine gen_esum()
    integer, parameter :: muv(7) = (/ -800, -10, -1, 0, 1, 10, 800 /)
    real(kind=rk), parameter :: xv(7) = (/ -800.0_rk, -10.0_rk, -0.5_rk, 0.0_rk, &
      0.5_rk, 10.0_rk, 800.0_rk /)
    integer :: unit, i, j
    call openf(unit, 'esum.csv', '# mu, x, esum(mu, x)')
    do i = 1, size(muv)
      do j = 1, size(xv)
        call putval(unit, real(muv(i), kind=rk), .false.)
        call putval(unit, xv(j), .false.)
        call putval(unit, esum(muv(i), xv(j)), .true.)
      end do
    end do
    close(unit)
  end subroutine gen_esum

  subroutine gen_exparg()
    integer :: unit, l
    call openf(unit, 'exparg.csv', '# l, exparg(l)')
    do l = -1, 1
      call putval(unit, real(l, kind=rk), .false.)
      call putval(unit, exparg(l), .true.)
    end do
    close(unit)
  end subroutine gen_exparg

  ! fpser, used by beta_inc when b <= min(eps, eps*a): covers
  ! a <= 0.001*eps, 0.001*eps < a, and the t < exparg(1) exit.
  subroutine gen_fpser()
    real(kind=rk), parameter :: av(7) = (/ 1.0e-20_rk, 1.0e-17_rk, 1.0e-10_rk, &
      0.5_rk, 1.0_rk, 5.0_rk, 100.0_rk /)
    real(kind=rk), parameter :: bv(2) = (/ 1.0e-20_rk, 1.0e-17_rk /)
    real(kind=rk), parameter :: xv(5) = (/ 1.0e-300_rk, 1.0e-10_rk, 0.01_rk, &
      0.3_rk, 0.5_rk /)
    real(kind=rk), parameter :: eps = 1.0e-15_rk
    integer :: unit, i, j, k
    call openf(unit, 'fpser.csv', '# a, b, x, eps, fpser(a, b, x, eps)')
    do i = 1, size(av)
      do j = 1, size(bv)
        do k = 1, size(xv)
          call putval(unit, av(i), .false.)
          call putval(unit, bv(j), .false.)
          call putval(unit, xv(k), .false.)
          call putval(unit, eps, .false.)
          call putval(unit, fpser(av(i), bv(j), xv(k), eps), .true.)
        end do
      end do
    end do
    close(unit)
  end subroutine gen_fpser

  ! gam1 on [-0.5..1.5]: t == 0, 0 < t, t < 0, both signs of d.
  subroutine gen_gam1()
    integer :: unit, i
    real(kind=rk) :: a
    call openf(unit, 'gam1.csv', '# a, gam1(a)')
    do i = -20, 60
      a = real(i, kind=rk) * 0.025_rk
      call putval(unit, a, .false.)
      call putval(unit, gam1(a), .true.)
    end do
    close(unit)
  end subroutine gen_gam1

  ! gamma_ln1 on [-0.2..1.25]: a < 0.6 and 0.6 <= a.
  subroutine gen_gamma_ln1()
    integer :: unit, i
    real(kind=rk) :: a
    call openf(unit, 'gamma_ln1.csv', '# a, gamma_ln1(a)')
    do i = -8, 50
      a = real(i, kind=rk) * 0.025_rk
      call putval(unit, a, .false.)
      call putval(unit, gamma_ln1(a), .true.)
    end do
    close(unit)
  end subroutine gen_gamma_ln1

  ! gamma_rat1 for a <= 1: a*x == 0 both ways, a == 0.5 both ways, the
  ! Taylor series with labels 30, 40, 50 and the q < 0 exit, and the
  ! continued fraction. r is the value exp(-x)*x^a/Gamma(a) passed by
  ! beta_grat.
  subroutine gen_gamma_rat1()
    real(kind=rk), parameter :: av(13) = (/ 0.0_rk, 1.0e-300_rk, 1.0e-100_rk, &
      1.0e-20_rk, 1.0e-16_rk, 1.0e-10_rk, 0.01_rk, 0.1_rk, 0.5_rk, 0.7_rk, &
      0.9_rk, 1.0_rk, 1.0e-15_rk /)
    real(kind=rk), parameter :: xv(11) = (/ 0.0_rk, 1.0e-5_rk, 0.1_rk, 0.2_rk, &
      0.24_rk, 0.3_rk, 0.5_rk, 0.9_rk, 1.0_rk, 1.5_rk, 20.0_rk /)
    real(kind=rk), parameter :: eps = 15.0e-15_rk
    real(kind=rk) :: a, x, r, p, q
    integer :: unit, i, k
    call openf(unit, 'gamma_rat1.csv', '# a, x, r, eps, p, q')
    do i = 1, size(av)
      do k = 1, size(xv)
        a = av(i)
        x = xv(k)
        if (0.0_rk < a .and. 0.0_rk < x) then
          r = exp(-x) * x**a / gamma_user(a)
        else
          r = 0.0_rk
        end if
        call gamma_rat1(a, x, r, p, q, eps)
        call putval(unit, a, .false.)
        call putval(unit, x, .false.)
        call putval(unit, r, .false.)
        call putval(unit, eps, .false.)
        call putval(unit, p, .false.)
        call putval(unit, q, .true.)
      end do
    end do
    close(unit)
  end subroutine gen_gamma_rat1

  ! gsumln for 1 <= a, b <= 2: x <= 0.25, x <= 1.25, 1.25 < x.
  subroutine gen_gsumln()
    real(kind=rk), parameter :: v(7) = (/ 1.0_rk, 1.1_rk, 1.25_rk, 1.5_rk, &
      1.625_rk, 1.75_rk, 2.0_rk /)
    integer :: unit, i, j
    call openf(unit, 'gsumln.csv', '# a, b, gsumln(a, b)')
    do i = 1, size(v)
      do j = 1, size(v)
        call putval(unit, v(i), .false.)
        call putval(unit, v(j), .false.)
        call putval(unit, gsumln(v(i), v(j)), .true.)
      end do
    end do
    close(unit)
  end subroutine gen_gsumln

  ! rcomp: a < 1, 1 <= a < 20, 20 <= a with u == 0 and u /= 0, and an
  ! underflowing result.
  subroutine gen_rcomp()
    real(kind=rk), parameter :: av(8) = (/ 0.1_rk, 0.5_rk, 1.0_rk, 5.0_rk, &
      19.9_rk, 20.0_rk, 50.0_rk, 1000.0_rk /)
    real(kind=rk), parameter :: xv(7) = (/ 0.0_rk, 1.0e-300_rk, 0.1_rk, 1.0_rk, &
      10.0_rk, 100.0_rk, 1000.0_rk /)
    integer :: unit, i, j
    call openf(unit, 'rcomp.csv', '# a, x, rcomp(a, x)')
    do i = 1, size(av)
      do j = 1, size(xv)
        call putval(unit, av(i), .false.)
        call putval(unit, xv(j), .false.)
        call putval(unit, rcomp(av(i), xv(j)), .true.)
      end do
    end do
    close(unit)
  end subroutine gen_rcomp

  subroutine gen_rexp()
    integer :: unit, i
    call openf(unit, 'rexp.csv', '# x, rexp(x)')
    do i = 1, size(ev)
      call putval(unit, ev(i), .false.)
      call putval(unit, rexp(ev(i)), .true.)
    end do
    close(unit)
  end subroutine gen_rexp

  ! rlog: every branch boundary.
  subroutine gen_rlog()
    real(kind=rk), parameter :: xv(15) = (/ 1.0e-300_rk, 0.3_rk, 0.6_rk, &
      0.61_rk, 0.7_rk, 0.82_rk, 0.9_rk, 1.0_rk, 1.1_rk, 1.18_rk, 1.3_rk, &
      1.57_rk, 2.0_rk, 10.0_rk, 1.0e10_rk /)
    integer :: unit, i
    call openf(unit, 'rlog.csv', '# x, rlog(x)')
    do i = 1, size(xv)
      call putval(unit, xv(i), .false.)
      call putval(unit, rlog(xv(i)), .true.)
    end do
    close(unit)
  end subroutine gen_rlog

  ! rlog1: every branch boundary.
  subroutine gen_rlog1()
    real(kind=rk), parameter :: xv(15) = (/ -0.999_rk, -0.9_rk, -0.39_rk, &
      -0.3_rk, -0.18_rk, -0.1_rk, 0.0_rk, 0.1_rk, 0.18_rk, 0.3_rk, 0.57_rk, &
      0.6_rk, 2.0_rk, 1.0e10_rk, 1.0e300_rk /)
    integer :: unit, i
    call openf(unit, 'rlog1.csv', '# x, rlog1(x)')
    do i = 1, size(xv)
      call putval(unit, xv(i), .false.)
      call putval(unit, rlog1(xv(i)), .true.)
    end do
    close(unit)
  end subroutine gen_rlog1

  ! gamma_user edge cases: poles at 0 and negative integers, the overflow
  ! of 1/t, 1000 <= |a|, negative a with 15 <= |a|, s == 0, and overflow
  ! of exp(w). A zero result is the F90 error value.
  subroutine gen_gamma_edge()
    real(kind=rk), parameter :: av(24) = (/ -15.95_rk, -20.95_rk, 0.0_rk, -1.0_rk, -2.0_rk, -14.0_rk, &
      1.0e-308_rk, -1.0e-308_rk, 1.0e-31_rk, 4.0e-324_rk, -4.0e-324_rk, &
      -15.0_rk, -20.0_rk, -20.5_rk, -100.25_rk, -170.5_rk, -999.5_rk, &
      1000.0_rk, -1000.0_rk, 1.0e10_rk, 171.6_rk, 172.0_rk, -0.95_rk, &
      -15.05_rk /)
    integer :: unit, i
    call openf(unit, 'gamma_edge.csv', '# a, gamma_user(a)')
    do i = 1, size(av)
      call putval(unit, av(i), .false.)
      call putval(unit, gamma_user(av(i)), .true.)
    end do
    close(unit)
  end subroutine gen_gamma_edge

  ! psi edge cases: x == 0, |x| <= xsmall, x <= -xmax1, negative
  ! integers (z == 0), every quadrant of the cotangent reduction, and
  ! xmax1 <= x. A zero result is the F90 error value.
  subroutine gen_psi_edge()
    real(kind=rk), parameter :: xv(22) = (/ 0.0_rk, 1.0e-10_rk, -1.0e-10_rk, &
      -1.0_rk, -2.0_rk, -3.0_rk, -5.0e15_rk, -1.0e16_rk, -0.125_rk, -0.375_rk, &
      -0.625_rk, -0.875_rk, -1.125_rk, -1.375_rk, -1.625_rk, -1.875_rk, &
      0.25_rk, 0.49_rk, 3.0_rk, 3.5_rk, 5.0e15_rk, 1.0e300_rk /)
    integer :: unit, i
    call openf(unit, 'psi_edge.csv', '# x, psi(x)')
    do i = 1, size(xv)
      call putval(unit, xv(i), .false.)
      call putval(unit, psi(xv(i)), .true.)
    end do
    close(unit)
  end subroutine gen_psi_edge

  ! gamma_inc error and edge exits: a < 0, x < 0, a = x = 0 (ans = 2),
  ! a*x == 0, u == 0, l == 0, r == 0 with x <= a, the Taylor series exits
  ! at label 180, the indeterminate results at labels 270, 330 and 410,
  ! and label 410 without the error. The last column is 1 when ans = 2.
  subroutine gen_gamma_inc_edge()
    real(kind=rk), parameter :: ax(2, 24) = reshape((/ &
      -1.0_rk, 1.0_rk,        1.0_rk, -1.0_rk,        0.0_rk, 0.0_rk, &
      1.0_rk, 0.0_rk,         0.0_rk, 1.0_rk,         0.3_rk, 1000.0_rk, &
      100.0_rk, 4.9e-324_rk,  10.0_rk, 1.0e-300_rk,   10.0_rk, 1.0e4_rk, &
      0.01_rk, 0.1_rk,        1.0e-10_rk, 1.0_rk,     1.0e-10_rk, 0.9_rk, &
      1.0e-6_rk, 1.0_rk,      1.0e30_rk, 1.0e30_rk,   1.0e40_rk, 1.0e40_rk, &
      1.0e6_rk, 1.0e3_rk,     1.0e6_rk, 1.0e9_rk,     2.5_rk, 1.0e-300_rk, &
      0.5_rk, 0.1_rk,         0.5_rk, 2.0_rk,         25.5_rk, 26.0_rk, &
      30.0_rk, 30.0_rk,       1.0e3_rk, 1.0e3_rk,     1.0e300_rk, 1.0e-300_rk /), &
      (/ 2, 24 /))
    real(kind=rk), parameter :: tiny_a(6) = (/ 1.0e-300_rk, 1.0e-100_rk, 1.0e-20_rk, &
      1.0e-17_rk, 1.0e-16_rk, 1.0e-15_rk /)
    real(kind=rk), parameter :: tiny_x(7) = (/ 0.01_rk, 0.1_rk, 0.24_rk, 0.3_rk, &
      0.5_rk, 0.9_rk, 1.0_rk /)
    integer :: j
    real(kind=rk) :: a, x
    integer :: unit, i, ind
    call openf(unit, 'gamma_inc_edge.csv', '# a, x, ind, ans, qans, error')
    do ind = 0, 2
      do i = 1, size(ax, 2)
        call ginc_row(unit, ax(1, i), ax(2, i), ind)
      end do
      ! a just above and below x at the scale where gamma_inc fails.
      a = 1.0e40_rk
      x = a * (1.0_rk - epsilon(1.0_rk))
      call ginc_row(unit, a, x, ind)
      a = 1.0e30_rk
      x = a * (1.0_rk + epsilon(1.0_rk))
      call ginc_row(unit, a, x, ind)
      a = 1.0e36_rk
      x = a * (1.0_rk + 2.0_rk * epsilon(1.0_rk))
      call ginc_row(unit, a, x, ind)
      ! Tiny a with x < 1.1, where the Taylor series for P(A,X)/X^A can
      ! yield a negative qans.
      do i = 1, size(tiny_a)
        do j = 1, size(tiny_x)
          call ginc_row(unit, tiny_a(i), tiny_x(j), ind)
        end do
      end do
    end do
    close(unit)
  end subroutine gen_gamma_inc_edge

  subroutine ginc_row(unit, a, x, ind)
    integer, intent(in) :: unit
    real(kind=rk), intent(in) :: a, x
    integer, intent(in) :: ind
    real(kind=rk) :: ans, qans
    ans = 0.0_rk
    qans = 0.0_rk
    call gamma_inc(a, x, ans, qans, ind)
    call putval(unit, a, .false.)
    call putval(unit, x, .false.)
    call putval(unit, real(ind, kind=rk), .false.)
    if (ans == 2.0_rk) then
      call putval(unit, 0.0_rk, .false.)
      call putval(unit, 0.0_rk, .false.)
      call putval(unit, 1.0_rk, .true.)
    else
      call putval(unit, ans, .false.)
      call putval(unit, qans, .false.)
      call putval(unit, 0.0_rk, .true.)
    end if
  end subroutine ginc_row

  ! gamma_inc_inv with an initial approximation x0 and on its error exits:
  ! ierr = -2, -3, -4, -6, -7, -8, q == 0 (x = huge), and the label 30
  ! early return. When ierr = -2, -4, -3, -7 the F90 x is not meaningful;
  ! the CSV writes it anyway.
  subroutine gen_gamma_inc_inv_edge()
    real(kind=rk), parameter :: rows(4, 30) = reshape((/ &
      ! a, x0, p, q
      0.0_rk, -1.0_rk, 0.5_rk, 0.5_rk, &
      -1.0_rk, -1.0_rk, 0.5_rk, 0.5_rk, &
      2.0_rk, -1.0_rk, 0.5_rk, 0.6_rk, &
      2.0_rk, -1.0_rk, 1.0_rk, 0.0_rk, &
      2.0_rk, -1.0_rk, 0.0_rk, 1.0_rk, &
      1.0_rk, -1.0_rk, 0.05_rk, 0.95_rk, &
      1.0_rk, -1.0_rk, 0.5_rk, 0.5_rk, &
      0.5_rk, -1.0_rk, 1.0_rk, 4.9e-324_rk, &
      1.0e-300_rk, -1.0_rk, 1.0e-10_rk, 0.9999999999_rk, &
      1.0e-6_rk, -1.0_rk, 0.9_rk, 0.1_rk, &
      0.01_rk, -1.0_rk, 0.99_rk, 0.01_rk, &
      0.2_rk, -1.0_rk, 0.1_rk, 0.9_rk, &
      0.2_rk, -1.0_rk, 1.0e-300_rk, 1.0_rk, &
      0.9_rk, -1.0_rk, 1.0e-30_rk, 1.0_rk, &
      600.0_rk, -1.0_rk, 0.5_rk, 0.5_rk, &
      600.0_rk, -1.0_rk, 0.4_rk, 0.6_rk, &
      1.0e10_rk, -1.0_rk, 0.5_rk, 0.5_rk, &
      5.0_rk, -1.0_rk, 1.0e-300_rk, 1.0_rk, &
      5.0_rk, -1.0_rk, 1.0_rk, 1.0e-300_rk, &
      1.5_rk, -1.0_rk, 0.999999_rk, 1.0e-6_rk, &
      2.0_rk, 1.0_rk, 0.3_rk, 0.7_rk, &
      2.0_rk, 3.0_rk, 0.7_rk, 0.3_rk, &
      2.0_rk, 1.0e-300_rk, 0.3_rk, 0.7_rk, &
      2.0_rk, 1.0e300_rk, 0.7_rk, 0.3_rk, &
      2.0_rk, 1.0e-300_rk, 0.9_rk, 0.1_rk, &
      1.0e20_rk, 1.0e20_rk, 0.5_rk, 0.5_rk, &
      1.0e20_rk, 1.0e20_rk, 0.3_rk, 0.7_rk, &
      0.5_rk, 1.0e-200_rk, 1.0e-300_rk, 1.0_rk, &
      3.0_rk, 50.0_rk, 0.3_rk, 0.7_rk, &
      3.0_rk, 1.0e-5_rk, 0.7_rk, 0.3_rk /), (/ 4, 30 /))
    real(kind=rk), parameter :: sa(7) = (/ 0.5_rk, 2.0_rk, 5.0_rk, 50.0_rk, &
      1000.0_rk, 1.0e21_rk, 0.05_rk /)
    real(kind=rk), parameter :: sx(11) = (/ 1.0e-300_rk, 1.0e-10_rk, 1.0e-3_rk, &
      0.1_rk, 0.5_rk, 1.0_rk, 2.0_rk, 10.0_rk, 100.0_rk, 1.0e10_rk, 1.0e300_rk /)
    real(kind=rk), parameter :: sp(7) = (/ 1.0e-10_rk, 0.1_rk, 0.3_rk, 0.5_rk, &
      0.7_rk, 0.9_rk, 0.9999999999_rk /)
    integer :: unit, i, j, k
    call openf(unit, 'gamma_inc_inv_edge.csv', '# a, x0, p, q, x, ierr')
    do i = 1, size(rows, 2)
      call ginv_row(unit, rows(1, i), rows(2, i), rows(3, i), rows(4, i))
    end do
    ! A q of 2^-1074 (the smallest subnormal) cannot be written as a
    ! literal, which gfortran flushes to zero.
    call ginv_row(unit, 0.5_rk, -1.0_rk, 1.0_rk, tiny(1.0_rk) * epsilon(1.0_rk))
    call ginv_row(unit, 5.0_rk, -1.0_rk, 1.0_rk, tiny(1.0_rk) * epsilon(1.0_rk))
    ! Label 30 returning without iteration, and label 40 with b*q <= 1e-8.
    call ginv_row(unit, 0.5_rk, -1.0_rk, 1.0_rk, 1.0e-30_rk)
    call ginv_row(unit, 1.0e-10_rk, -1.0_rk, 1.0_rk - 1.0e-9_rk, 1.0e-9_rk)
    ! Schroder iterations from a caller-supplied x0 in both branches.
    do i = 1, size(sa)
      do j = 1, size(sx)
        do k = 1, size(sp)
          call ginv_row(unit, sa(i), min(sx(j) * max(1.0_rk, sa(i)), 1.0e300_rk), sp(k), &
            1.0_rk - sp(k))
        end do
      end do
    end do
    close(unit)
  end subroutine gen_gamma_inc_inv_edge

  subroutine ginv_row(unit, a, x0, p, q)
    integer, intent(in) :: unit
    real(kind=rk), intent(in) :: a, x0, p, q
    real(kind=rk) :: x
    integer :: ierr
    x = 0.0_rk
    ierr = 0
    call gamma_inc_inv(a, x, x0, p, q, ierr)
    call putval(unit, a, .false.)
    call putval(unit, x0, .false.)
    call putval(unit, p, .false.)
    call putval(unit, q, .false.)
    call putval(unit, x, .false.)
    call putval(unit, real(ierr, kind=rk), .true.)
  end subroutine ginv_row

  ! beta_inc: every evaluation label (90 fpser, 100 apser, 110/120
  ! beta_pser, 130 beta_frac, 140/150/160 beta_up and beta_grat, 200
  ! beta_asym, 260 tiny a and b), the special values of x, y, a, b, and
  ! every ierr. The CSV writes w = w1 = 0 when ierr /= 0.
  subroutine gen_beta_inc_regimes()
    real(kind=rk), parameter :: rows(4, 32) = reshape((/ &
      ! a, b, x, y
      1.0e-20_rk, 1.0e-20_rk, 0.3_rk, 0.7_rk, &
      0.5_rk, 1.0e-20_rk, 0.3_rk, 0.7_rk, &
      1.0e-20_rk, 0.5_rk, 0.3_rk, 0.7_rk, &
      1.0e-20_rk, 0.5_rk, 0.7_rk, 0.3_rk, &
      0.01_rk, 0.5_rk, 0.1_rk, 0.9_rk, &
      0.01_rk, 0.5_rk, 0.4_rk, 0.6_rk, &
      0.5_rk, 0.01_rk, 0.9_rk, 0.1_rk, &
      0.3_rk, 0.5_rk, 0.1_rk, 0.9_rk, &
      0.5_rk, 3.0_rk, 0.2_rk, 0.8_rk, &
      0.5_rk, 3.0_rk, 0.01_rk, 0.99_rk, &
      0.5_rk, 20.0_rk, 0.2_rk, 0.8_rk, &
      0.5_rk, 20.0_rk, 0.05_rk, 0.95_rk, &
      200.0_rk, 300.0_rk, 0.4_rk, 0.6_rk, &
      200.0_rk, 300.0_rk, 0.38_rk, 0.62_rk, &
      300.0_rk, 200.0_rk, 0.6_rk, 0.4_rk, &
      300.0_rk, 200.0_rk, 0.58_rk, 0.42_rk, &
      50.0_rk, 60.0_rk, 0.45_rk, 0.55_rk, &
      30.0_rk, 5.0_rk, 0.8_rk, 0.2_rk, &
      10.0_rk, 3.0_rk, 0.75_rk, 0.25_rk, &
      2.0_rk, 50.0_rk, 0.01_rk, 0.99_rk, &
      2.0_rk, 0.0_rk, 0.5_rk, 0.5_rk, &
      0.0_rk, 2.0_rk, 0.5_rk, 0.5_rk, &
      2.0_rk, 3.0_rk, 0.0_rk, 1.0_rk, &
      2.0_rk, 3.0_rk, 1.0_rk, 0.0_rk, &
      -1.0_rk, 2.0_rk, 0.5_rk, 0.5_rk, &
      0.0_rk, 0.0_rk, 0.5_rk, 0.5_rk, &
      2.0_rk, 3.0_rk, -0.1_rk, 1.1_rk, &
      2.0_rk, 3.0_rk, 0.5_rk, 1.5_rk, &
      2.0_rk, 3.0_rk, 0.5_rk, 0.4_rk, &
      0.0_rk, 3.0_rk, 0.0_rk, 1.0_rk, &
      2.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, &
      1.0e5_rk, 1.0e5_rk, 0.5_rk, 0.5_rk /), (/ 4, 32 /))
    real(kind=rk) :: w, w1
    integer :: unit, i, ierr
    call openf(unit, 'beta_inc_regimes.csv', '# a, b, x, y, w, w1, ierr')
    do i = 1, size(rows, 2)
      w = 0.0_rk
      w1 = 0.0_rk
      ierr = 0
      call beta_inc(rows(1, i), rows(2, i), rows(3, i), rows(4, i), w, w1, ierr)
      if (ierr /= 0) then
        w = 0.0_rk
        w1 = 0.0_rk
      end if
      call putval(unit, rows(1, i), .false.)
      call putval(unit, rows(2, i), .false.)
      call putval(unit, rows(3, i), .false.)
      call putval(unit, rows(4, i), .false.)
      call putval(unit, w, .false.)
      call putval(unit, w1, .false.)
      call putval(unit, real(ierr, kind=rk), .true.)
    end do
    close(unit)
  end subroutine gen_beta_inc_regimes

end program gen_kernel_coverage

