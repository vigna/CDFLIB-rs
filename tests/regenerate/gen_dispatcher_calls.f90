! Call logs of every cdf* dispatcher in cdflib.f90.
!
! Each row of tests/data/<cdf>_calls.csv records one call: which, the
! arguments before the call, status and bound, and the arguments after
! the call. The rows cover every value of which, both sides of p <= q,
! invalid arguments for each negative status, P + Q /= 1, and searches
! that fail at either bound or meet the error value of gamma_inc (status
! 10).
!
! Conventions that let the Rust API, which derives complements, receive
! exactly the arguments the F90 received: a q of -99 means q = 1 - p and
! then p = 1 - q (so that each is the exact complement of the other); a
! y or ompr of -99 means 1 - x or 1 - pr.
!
! Companion Rust test: tests/dispatcher_calls.rs.

program gen_dispatcher_calls
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), parameter :: auto = -99.0_rk
  external :: cdfbet, cdfbin, cdfchi, cdfchn, cdff, cdffnc, cdfgam, cdfnbn, &
    cdfnor, cdfpoi, cdft

  call log_cdfbet()
  call log_cdfbin()
  call log_cdfchi()
  call log_cdfchn()
  call log_cdff()
  call log_cdffnc()
  call log_cdfgam()
  call log_cdfnbn()
  call log_cdfnor()
  call log_cdfpoi()
  call log_cdft()

  write(0, '(a)') 'wrote 11 tables under tests/data/'

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

  ! Writes which, the arguments before the call, status, bound, and the
  ! arguments after the call.
  subroutine emit(unit, which, before, status, bound, after)
    integer, intent(in) :: unit, which, status
    real(kind=rk), intent(in) :: before(:), bound, after(:)
    integer :: i
    call putval(unit, real(which, kind=rk), .false.)
    do i = 1, size(before)
      call putval(unit, before(i), .false.)
    end do
    call putval(unit, real(status, kind=rk), .false.)
    call putval(unit, bound, .false.)
    do i = 1, size(after) - 1
      call putval(unit, after(i), .false.)
    end do
    call putval(unit, after(size(after)), .true.)
  end subroutine emit

  ! Resolves the q = -99 convention.
  subroutine resolve_pq(p, q)
    real(kind=rk), intent(inout) :: p, q
    if (q == auto) then
      q = 1.0_rk - p
      p = 1.0_rk - q
    end if
  end subroutine resolve_pq

  subroutine log_cdfbet()
    ! which, p, q, x, y, a, b
    real(kind=rk), parameter :: r(7, 38) = reshape((/ &
      1.0_rk, 0.0_rk, 0.0_rk, 0.3_rk, auto, 2.0_rk, 5.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 0.999_rk, auto, 0.5_rk, 0.5_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, auto, 2.0_rk, 5.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, auto, 2.0_rk, 5.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 0.3_rk, auto, 0.0_rk, 5.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 0.3_rk, auto, 2.0_rk, -1.0_rk, &
      2.0_rk, 0.1_rk, auto, 0.0_rk, 0.0_rk, 2.0_rk, 5.0_rk, &
      2.0_rk, 0.9_rk, auto, 0.0_rk, 0.0_rk, 2.0_rk, 5.0_rk, &
      2.0_rk, 1.0e-10_rk, auto, 0.0_rk, 0.0_rk, 0.5_rk, 0.5_rk, &
      2.0_rk, 0.999999_rk, auto, 0.0_rk, 0.0_rk, 50.0_rk, 25.0_rk, &
      2.0_rk, -0.1_rk, 1.1_rk, 0.0_rk, 0.0_rk, 2.0_rk, 5.0_rk, &
      2.0_rk, 1.1_rk, -0.1_rk, 0.0_rk, 0.0_rk, 2.0_rk, 5.0_rk, &
      2.0_rk, 0.5_rk, -0.1_rk, 0.0_rk, 0.0_rk, 2.0_rk, 5.0_rk, &
      2.0_rk, 0.5_rk, 1.5_rk, 0.0_rk, 0.0_rk, 2.0_rk, 5.0_rk, &
      2.0_rk, 0.5_rk, 0.6_rk, 0.0_rk, 0.0_rk, 2.0_rk, 5.0_rk, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, 0.0_rk, 0.0_rk, 5.0_rk, &
      3.0_rk, 0.1_rk, auto, 0.3_rk, auto, 0.0_rk, 5.0_rk, &
      3.0_rk, 0.9_rk, auto, 0.3_rk, auto, 0.0_rk, 5.0_rk, &
      3.0_rk, 0.5_rk, auto, 0.5_rk, auto, 0.0_rk, 2.0_rk, &
      3.0_rk, 1.0e-300_rk, auto, 0.9_rk, auto, 0.0_rk, 1.0_rk, &
      3.0_rk, 0.999_rk, auto, 1.0e-300_rk, auto, 0.0_rk, 1.0_rk, &
      3.0_rk, 0.5_rk, auto, -0.1_rk, auto, 0.0_rk, 2.0_rk, &
      3.0_rk, 0.5_rk, auto, 1.1_rk, auto, 0.0_rk, 2.0_rk, &
      3.0_rk, 0.5_rk, auto, 0.5_rk, auto, 0.0_rk, 0.0_rk, &
      3.0_rk, 0.5_rk, 0.6_rk, 0.5_rk, auto, 0.0_rk, 2.0_rk, &
      4.0_rk, 0.1_rk, auto, 0.3_rk, auto, 2.0_rk, 0.0_rk, &
      4.0_rk, 0.9_rk, auto, 0.3_rk, auto, 2.0_rk, 0.0_rk, &
      4.0_rk, 0.5_rk, auto, 0.5_rk, auto, 3.0_rk, 0.0_rk, &
      4.0_rk, 1.0e-300_rk, auto, 0.1_rk, auto, 1.0_rk, 0.0_rk, &
      4.0_rk, 0.999_rk, auto, 0.999_rk, auto, 1.0_rk, 0.0_rk, &
      4.0_rk, 0.5_rk, auto, 0.5_rk, auto, -2.0_rk, 0.0_rk, &
      4.0_rk, -0.5_rk, auto, 0.5_rk, auto, 2.0_rk, 0.0_rk, &
      4.0_rk, 0.5_rk, 0.6_rk, 0.5_rk, auto, 2.0_rk, 0.0_rk, &
      4.0_rk, 0.5_rk, auto, 1.5_rk, auto, 2.0_rk, 0.0_rk, &
      3.0_rk, 0.5_rk, auto, 0.0_rk, auto, 0.0_rk, 2.0_rk, &
      3.0_rk, 0.5_rk, auto, 1.0_rk, auto, 0.0_rk, 2.0_rk, &
      4.0_rk, 0.5_rk, auto, 0.0_rk, auto, 2.0_rk, 0.0_rk, &
      4.0_rk, 0.5_rk, auto, 1.0_rk, auto, 2.0_rk, 0.0_rk /), (/ 7, 38 /))
    real(kind=rk) :: v(6), a(6), bound
    integer :: unit, i, which, status
    open(newunit=unit, file='tests/data/cdfbet_calls.csv', status='replace', action='write')
    write(unit, '(a)') '# which, p, q, x, y, a, b, status, bound, p, q, x, y, a, b (after)'
    do i = 1, size(r, 2)
      which = int(r(1, i))
      v = r(2:7, i)
      call resolve_pq(v(1), v(2))
      if (v(4) == auto) v(4) = 1.0_rk - v(3)
      a = v
      status = 0
      bound = 0.0_rk
      call cdfbet(which, a(1), a(2), a(3), a(4), a(5), a(6), status, bound)
      call emit(unit, which, v, status, bound, a)
    end do
    close(unit)
  end subroutine log_cdfbet

  subroutine log_cdfbin()
    ! which, p, q, s, xn, pr, ompr
    real(kind=rk), parameter :: r(7, 41) = reshape((/ &
      1.0_rk, 0.0_rk, 0.0_rk, 3.0_rk, 10.0_rk, 0.3_rk, auto, &
      1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 10.0_rk, 0.3_rk, auto, &
      1.0_rk, 0.0_rk, 0.0_rk, 10.0_rk, 10.0_rk, 0.3_rk, auto, &
      1.0_rk, 0.0_rk, 0.0_rk, 3.0_rk, 10.0_rk, 0.0_rk, auto, &
      1.0_rk, 0.0_rk, 0.0_rk, 3.0_rk, 10.0_rk, 1.0_rk, auto, &
      1.0_rk, 0.0_rk, 0.0_rk, 3.0_rk, 0.0_rk, 0.3_rk, auto, &
      1.0_rk, 0.0_rk, 0.0_rk, 3.0_rk, 10.0_rk, 1.5_rk, auto, &
      2.0_rk, 0.2_rk, auto, 0.0_rk, 10.0_rk, 0.3_rk, auto, &
      2.0_rk, 0.8_rk, auto, 0.0_rk, 10.0_rk, 0.3_rk, auto, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, 100.0_rk, 0.25_rk, auto, &
      2.0_rk, 1.0e-30_rk, auto, 0.0_rk, 10.0_rk, 0.5_rk, auto, &
      2.0_rk, 0.9999999_rk, auto, 0.0_rk, 10.0_rk, 0.5_rk, auto, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, 10.0_rk, -0.1_rk, auto, &
      2.0_rk, 0.5_rk, 0.7_rk, 0.0_rk, 10.0_rk, 0.3_rk, auto, &
      2.0_rk, -0.5_rk, auto, 0.0_rk, 10.0_rk, 0.3_rk, auto, &
      3.0_rk, 0.2_rk, auto, 3.0_rk, 0.0_rk, 0.3_rk, auto, &
      3.0_rk, 0.8_rk, auto, 3.0_rk, 0.0_rk, 0.3_rk, auto, &
      3.0_rk, 0.5_rk, auto, 25.0_rk, 0.0_rk, 0.5_rk, auto, &
      3.0_rk, 1.0e-300_rk, auto, 3.0_rk, 0.0_rk, 0.3_rk, auto, &
      3.0_rk, 0.9999_rk, auto, 3.0_rk, 0.0_rk, 0.3_rk, auto, &
      3.0_rk, 0.5_rk, auto, 3.0_rk, 0.0_rk, 2.0_rk, auto, &
      3.0_rk, 0.5_rk, 0.25_rk, 3.0_rk, 0.0_rk, 0.3_rk, auto, &
      4.0_rk, 0.2_rk, auto, 3.0_rk, 10.0_rk, 0.0_rk, auto, &
      4.0_rk, 0.8_rk, auto, 3.0_rk, 10.0_rk, 0.0_rk, auto, &
      4.0_rk, 0.5_rk, auto, 25.0_rk, 50.0_rk, 0.0_rk, auto, &
      4.0_rk, 0.5_rk, auto, 10.0_rk, 10.0_rk, 0.0_rk, auto, &
      4.0_rk, 1.0e-300_rk, auto, 0.0_rk, 10.0_rk, 0.0_rk, auto, &
      4.0_rk, 0.5_rk, auto, 11.0_rk, 10.0_rk, 0.0_rk, auto, &
      4.0_rk, 0.5_rk, auto, 3.0_rk, 0.0_rk, 0.0_rk, auto, &
      4.0_rk, 1.5_rk, auto, 3.0_rk, 10.0_rk, 0.0_rk, auto, &
      4.0_rk, 0.5_rk, 0.6_rk, 3.0_rk, 10.0_rk, 0.0_rk, auto , &
      2.0_rk, 0.5_rk, -0.5_rk, 0.0_rk, 10.0_rk, 0.3_rk, auto, &
      2.0_rk, 0.5_rk, 1.5_rk, 0.0_rk, 10.0_rk, 0.3_rk, auto, &
      3.0_rk, 0.5_rk, auto, -1.0_rk, 0.0_rk, 0.3_rk, auto, &
      4.0_rk, 0.8_rk, auto, 10.0_rk, 10.0_rk, 0.0_rk, auto, &
      3.0_rk, 0.5_rk, auto, 3.0_rk, 0.0_rk, 0.0_rk, auto, &
      3.0_rk, 0.01_rk, auto, 3.0_rk, 0.0_rk, 1.0e-300_rk, auto, &
      2.0_rk, 1.0_rk, auto, 0.0_rk, 10.0_rk, 0.3_rk, auto, &
      2.0_rk, 0.0_rk, auto, 0.0_rk, 10.0_rk, 0.3_rk, auto, &
      2.0_rk, 0.7_rk, auto, 0.0_rk, 10.0_rk, 0.0_rk, auto, &
      3.0_rk, 1.0_rk, auto, 0.0_rk, 0.0_rk, 0.5_rk, auto /), (/ 7, 41 /))
    real(kind=rk) :: v(6), a(6), bound
    integer :: unit, i, which, status
    open(newunit=unit, file='tests/data/cdfbin_calls.csv', status='replace', action='write')
    write(unit, '(a)') '# which, p, q, s, xn, pr, ompr, status, bound, p, q, s, xn, pr, ompr (after)'
    do i = 1, size(r, 2)
      which = int(r(1, i))
      v = r(2:7, i)
      call resolve_pq(v(1), v(2))
      if (v(6) == auto) v(6) = 1.0_rk - v(5)
      a = v
      status = 0
      bound = 0.0_rk
      call cdfbin(which, a(1), a(2), a(3), a(4), a(5), a(6), status, bound)
      call emit(unit, which, v, status, bound, a)
    end do
    close(unit)
  end subroutine log_cdfbin

  subroutine log_cdfchi()
    ! which, p, q, x, df
    real(kind=rk), parameter :: r(5, 29) = reshape((/ &
      1.0_rk, 0.0_rk, 0.0_rk, 3.84_rk, 1.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 5.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 1000.0_rk, 5.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 3.0_rk, 0.0_rk, &
      2.0_rk, 0.05_rk, auto, 0.0_rk, 5.0_rk, &
      2.0_rk, 0.95_rk, auto, 0.0_rk, 5.0_rk, &
      2.0_rk, 1.0e-12_rk, auto, 0.0_rk, 1000.0_rk, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, 1.0e-3_rk, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, -1.0_rk, &
      2.0_rk, -0.5_rk, auto, 0.0_rk, 5.0_rk, &
      2.0_rk, 0.5_rk, 0.6_rk, 0.0_rk, 5.0_rk, &
      3.0_rk, 0.05_rk, auto, 3.84_rk, 0.0_rk, &
      3.0_rk, 0.95_rk, auto, 3.84_rk, 0.0_rk, &
      3.0_rk, 0.5_rk, auto, 100.0_rk, 0.0_rk, &
      3.0_rk, 1.0e-300_rk, auto, 3.0_rk, 0.0_rk, &
      3.0_rk, 0.9999999_rk, auto, 1.0e-3_rk, 0.0_rk, &
      3.0_rk, 0.5_rk, auto, -1.0_rk, 0.0_rk, &
      3.0_rk, 0.5_rk, -0.5_rk, 3.0_rk, 0.0_rk, &
      3.0_rk, 0.5_rk, 0.4_rk, 3.0_rk, 0.0_rk, &
      3.0_rk, 1.5_rk, auto, 3.0_rk, 0.0_rk, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, 1.0e10_rk , &
      3.0_rk, 0.5_rk, 1.5_rk, 3.0_rk, 0.0_rk, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, 1.0e300_rk, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, 1.0e305_rk, &
      3.0_rk, 0.5_rk, auto, 1.0e300_rk, 0.0_rk, &
      3.0_rk, 0.5_rk, auto, 1.42108547152020022e29_rk, 0.0_rk, &
      2.0_rk, 0.7_rk, auto, 0.0_rk, 1.0e305_rk, &
      3.0_rk, 0.5_rk, auto, 1.5e300_rk, 0.0_rk, &
      3.0_rk, 0.7_rk, auto, 1.5e300_rk, 0.0_rk /), (/ 5, 29 /))
    real(kind=rk) :: v(4), a(4), bound
    integer :: unit, i, which, status
    open(newunit=unit, file='tests/data/cdfchi_calls.csv', status='replace', action='write')
    write(unit, '(a)') '# which, p, q, x, df, status, bound, p, q, x, df (after)'
    do i = 1, size(r, 2)
      which = int(r(1, i))
      v = r(2:5, i)
      call resolve_pq(v(1), v(2))
      a = v
      status = 0
      bound = 0.0_rk
      call cdfchi(which, a(1), a(2), a(3), a(4), status, bound)
      call emit(unit, which, v, status, bound, a)
    end do
    close(unit)
  end subroutine log_cdfchi

  subroutine log_cdfchn()
    ! which, p, q, x, df, pnonc
    real(kind=rk), parameter :: r(6, 28) = reshape((/ &
      1.0_rk, 0.0_rk, 0.0_rk, 15.0_rk, 5.0_rk, 10.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 3.0_rk, 4.0_rk, 0.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 4.0_rk, 2.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 3.0_rk, 0.0_rk, 2.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 3.0_rk, 4.0_rk, -1.0_rk, &
      2.0_rk, 0.25_rk, 0.0_rk, 0.0_rk, 5.0_rk, 2.0_rk, &
      2.0_rk, 0.9_rk, 0.0_rk, 0.0_rk, 5.0_rk, 2.0_rk, &
      2.0_rk, 0.5_rk, 0.0_rk, 0.0_rk, 2.0_rk, 0.5_rk, &
      2.0_rk, 0.5_rk, 0.0_rk, 0.0_rk, 2.0_rk, 0.0_rk, &
      2.0_rk, -0.1_rk, 0.0_rk, 0.0_rk, 2.0_rk, 1.0_rk, &
      2.0_rk, 0.5_rk, 0.0_rk, 0.0_rk, -2.0_rk, 1.0_rk, &
      3.0_rk, 0.25_rk, 0.0_rk, 15.0_rk, 0.0_rk, 10.0_rk, &
      3.0_rk, 0.75_rk, 0.0_rk, 15.0_rk, 0.0_rk, 10.0_rk, &
      3.0_rk, 1.0e-300_rk, 0.0_rk, 15.0_rk, 0.0_rk, 10.0_rk, &
      3.0_rk, 0.999999_rk, 0.0_rk, 0.01_rk, 0.0_rk, 10.0_rk, &
      3.0_rk, 0.5_rk, 0.0_rk, -1.0_rk, 0.0_rk, 10.0_rk, &
      3.0_rk, 0.5_rk, 0.0_rk, 15.0_rk, 0.0_rk, -1.0_rk, &
      4.0_rk, 0.25_rk, 0.0_rk, 15.0_rk, 5.0_rk, 0.0_rk, &
      4.0_rk, 0.75_rk, 0.0_rk, 15.0_rk, 5.0_rk, 0.0_rk, &
      4.0_rk, 1.0e-300_rk, 0.0_rk, 15.0_rk, 5.0_rk, 0.0_rk, &
      4.0_rk, 0.999_rk, 0.0_rk, 15.0_rk, 5.0_rk, 0.0_rk, &
      4.0_rk, 0.5_rk, 0.0_rk, 15.0_rk, -5.0_rk, 0.0_rk, &
      4.0_rk, 1.5_rk, 0.0_rk, 15.0_rk, 5.0_rk, 0.0_rk , &
      4.0_rk, 1.0e-30_rk, 0.0_rk, 15.0_rk, 5.0_rk, 0.0_rk, &
      3.0_rk, 1.0e-30_rk, 0.0_rk, 15.0_rk, 0.0_rk, 10.0_rk , &
      4.0_rk, 1.0e-10_rk, 0.0_rk, 9000.0_rk, 5.0_rk, 0.0_rk, &
      2.0_rk, 0.9999999_rk, 0.0_rk, 0.0_rk, 5.0_rk, 2.0_rk, &
      3.0_rk, 0.1_rk, 0.0_rk, 0.5_rk, 0.0_rk, 100.0_rk /), (/ 6, 28 /))
    real(kind=rk) :: v(5), a(5), bound
    integer :: unit, i, which, status
    open(newunit=unit, file='tests/data/cdfchn_calls.csv', status='replace', action='write')
    write(unit, '(a)') '# which, p, q, x, df, pnonc, status, bound, p, q, x, df, pnonc (after)'
    do i = 1, size(r, 2)
      which = int(r(1, i))
      v = r(2:6, i)
      a = v
      status = 0
      bound = 0.0_rk
      call cdfchn(which, a(1), a(2), a(3), a(4), a(5), status, bound)
      call emit(unit, which, v, status, bound, a)
    end do
    close(unit)
  end subroutine log_cdfchn

  subroutine log_cdff()
    ! which, p, q, f, dfn, dfd
    real(kind=rk), parameter :: r(6, 33) = reshape((/ &
      1.0_rk, 0.0_rk, 0.0_rk, 3.33_rk, 5.0_rk, 10.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 0.01_rk, 5.0_rk, 10.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 5.0_rk, 10.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 3.0_rk, 0.0_rk, 10.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 3.0_rk, 5.0_rk, 0.0_rk, &
      2.0_rk, 0.05_rk, auto, 0.0_rk, 5.0_rk, 10.0_rk, &
      2.0_rk, 0.95_rk, auto, 0.0_rk, 5.0_rk, 10.0_rk, &
      2.0_rk, 1.0e-20_rk, auto, 0.0_rk, 2.0_rk, 2.0_rk, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, -5.0_rk, 10.0_rk, &
      2.0_rk, 0.5_rk, 0.6_rk, 0.0_rk, 5.0_rk, 10.0_rk, &
      3.0_rk, 0.05_rk, auto, 3.33_rk, 0.0_rk, 10.0_rk, &
      3.0_rk, 0.95_rk, auto, 3.33_rk, 0.0_rk, 10.0_rk, &
      3.0_rk, 0.5_rk, auto, 1.0_rk, 0.0_rk, 10.0_rk, &
      3.0_rk, 1.0e-300_rk, auto, 3.33_rk, 0.0_rk, 10.0_rk, &
      3.0_rk, 0.9999999_rk, auto, 0.01_rk, 0.0_rk, 10.0_rk, &
      3.0_rk, 0.5_rk, auto, -1.0_rk, 0.0_rk, 10.0_rk, &
      3.0_rk, 0.5_rk, auto, 3.0_rk, 0.0_rk, 0.0_rk, &
      3.0_rk, 0.5_rk, -0.1_rk, 3.0_rk, 0.0_rk, 10.0_rk, &
      4.0_rk, 0.05_rk, auto, 3.33_rk, 5.0_rk, 0.0_rk, &
      4.0_rk, 0.95_rk, auto, 3.33_rk, 5.0_rk, 0.0_rk, &
      4.0_rk, 0.5_rk, auto, 1.0_rk, 5.0_rk, 0.0_rk, &
      4.0_rk, 1.0e-300_rk, auto, 3.33_rk, 5.0_rk, 0.0_rk, &
      4.0_rk, 0.9999999_rk, auto, 0.01_rk, 5.0_rk, 0.0_rk, &
      4.0_rk, 0.5_rk, auto, 3.0_rk, 0.0_rk, 0.0_rk, &
      4.0_rk, 0.5_rk, 0.6_rk, 3.0_rk, 5.0_rk, 0.0_rk , &
      2.0_rk, -0.5_rk, auto, 0.0_rk, 5.0_rk, 10.0_rk, &
      2.0_rk, 1.5_rk, auto, 0.0_rk, 5.0_rk, 10.0_rk, &
      2.0_rk, 0.5_rk, 1.5_rk, 0.0_rk, 5.0_rk, 10.0_rk, &
      3.0_rk, 1.0e-10_rk, auto, 3.33_rk, 0.0_rk, 10.0_rk, &
      3.0_rk, 0.9999999999_rk, auto, 3.33_rk, 0.0_rk, 10.0_rk, &
      4.0_rk, 1.0e-10_rk, auto, 3.33_rk, 5.0_rk, 0.0_rk, &
      4.0_rk, 0.9999999999_rk, auto, 3.33_rk, 5.0_rk, 0.0_rk, &
      2.0_rk, 0.75_rk, auto, 0.0_rk, 1.0e-10_rk, 1.0e-10_rk /), (/ 6, 33 /))
    real(kind=rk) :: v(5), a(5), bound
    integer :: unit, i, which, status
    open(newunit=unit, file='tests/data/cdff_calls.csv', status='replace', action='write')
    write(unit, '(a)') '# which, p, q, f, dfn, dfd, status, bound, p, q, f, dfn, dfd (after)'
    do i = 1, size(r, 2)
      which = int(r(1, i))
      v = r(2:6, i)
      call resolve_pq(v(1), v(2))
      a = v
      status = 0
      bound = 0.0_rk
      call cdff(which, a(1), a(2), a(3), a(4), a(5), status, bound)
      call emit(unit, which, v, status, bound, a)
    end do
    close(unit)
  end subroutine log_cdff

  subroutine log_cdffnc()
    ! which, p, q, f, dfn, dfd, pnonc
    real(kind=rk), parameter :: r(7, 32) = reshape((/ &
      1.0_rk, 0.0_rk, 0.0_rk, 4.0_rk, 5.0_rk, 10.0_rk, 2.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 4.0_rk, 5.0_rk, 10.0_rk, 0.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 5.0_rk, 10.0_rk, 2.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 4.0_rk, 0.0_rk, 10.0_rk, 2.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 4.0_rk, 5.0_rk, 0.0_rk, 2.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 4.0_rk, 5.0_rk, 10.0_rk, -1.0_rk, &
      2.0_rk, 0.25_rk, 0.0_rk, 0.0_rk, 5.0_rk, 10.0_rk, 2.0_rk, &
      2.0_rk, 0.9_rk, 0.0_rk, 0.0_rk, 5.0_rk, 10.0_rk, 2.0_rk, &
      2.0_rk, -0.1_rk, 0.0_rk, 0.0_rk, 5.0_rk, 10.0_rk, 2.0_rk, &
      2.0_rk, 0.5_rk, 0.0_rk, 0.0_rk, 5.0_rk, 10.0_rk, -2.0_rk, &
      3.0_rk, 0.25_rk, 0.0_rk, 4.0_rk, 0.0_rk, 10.0_rk, 2.0_rk, &
      3.0_rk, 0.75_rk, 0.0_rk, 4.0_rk, 0.0_rk, 10.0_rk, 2.0_rk, &
      3.0_rk, 1.0e-300_rk, 0.0_rk, 4.0_rk, 0.0_rk, 10.0_rk, 2.0_rk, &
      3.0_rk, 0.99_rk, 0.0_rk, 4.0_rk, 0.0_rk, 10.0_rk, 2.0_rk, &
      3.0_rk, 0.5_rk, 0.0_rk, -4.0_rk, 0.0_rk, 10.0_rk, 2.0_rk, &
      4.0_rk, 0.25_rk, 0.0_rk, 4.0_rk, 5.0_rk, 0.0_rk, 2.0_rk, &
      4.0_rk, 0.75_rk, 0.0_rk, 4.0_rk, 5.0_rk, 0.0_rk, 2.0_rk, &
      4.0_rk, 1.0e-300_rk, 0.0_rk, 4.0_rk, 5.0_rk, 0.0_rk, 2.0_rk, &
      4.0_rk, 0.999_rk, 0.0_rk, 4.0_rk, 5.0_rk, 0.0_rk, 2.0_rk, &
      4.0_rk, 0.5_rk, 0.0_rk, 4.0_rk, 5.0_rk, 0.0_rk, -2.0_rk, &
      5.0_rk, 0.25_rk, 0.0_rk, 4.0_rk, 5.0_rk, 10.0_rk, 0.0_rk, &
      5.0_rk, 0.5_rk, 0.0_rk, 4.0_rk, 5.0_rk, 10.0_rk, 0.0_rk, &
      5.0_rk, 1.0e-300_rk, 0.0_rk, 4.0_rk, 5.0_rk, 10.0_rk, 0.0_rk, &
      5.0_rk, 0.999_rk, 0.0_rk, 4.0_rk, 5.0_rk, 10.0_rk, 0.0_rk, &
      5.0_rk, 0.5_rk, 0.0_rk, 4.0_rk, 5.0_rk, -10.0_rk, 0.0_rk, &
      5.0_rk, 1.5_rk, 0.0_rk, 4.0_rk, 5.0_rk, 10.0_rk, 0.0_rk , &
      3.0_rk, 1.0e-30_rk, 0.0_rk, 4.0_rk, 0.0_rk, 10.0_rk, 2.0_rk, &
      5.0_rk, 1.0e-30_rk, 0.0_rk, 4.0_rk, 5.0_rk, 10.0_rk, 0.0_rk, &
      4.0_rk, 1.0e-30_rk, 0.0_rk, 4.0_rk, 5.0_rk, 0.0_rk, 2.0_rk , &
      5.0_rk, 1.0e-10_rk, 0.0_rk, 1000.0_rk, 5.0_rk, 10.0_rk, 0.0_rk, &
      3.0_rk, 0.9999999_rk, 0.0_rk, 4.0_rk, 0.0_rk, 10.0_rk, 2.0_rk, &
      2.0_rk, 0.999999_rk, 0.0_rk, 0.0_rk, 5.0_rk, 10.0_rk, 2.0_rk /), (/ 7, 32 /))
    real(kind=rk) :: v(6), a(6), bound
    integer :: unit, i, which, status
    open(newunit=unit, file='tests/data/cdffnc_calls.csv', status='replace', action='write')
    write(unit, '(a)') '# which, p, q, f, dfn, dfd, pnonc, status, bound, p, q, f, dfn, dfd, pnonc (after)'
    do i = 1, size(r, 2)
      which = int(r(1, i))
      v = r(2:7, i)
      a = v
      status = 0
      bound = 0.0_rk
      call cdffnc(which, a(1), a(2), a(3), a(4), a(5), a(6), status, bound)
      call emit(unit, which, v, status, bound, a)
    end do
    close(unit)
  end subroutine log_cdffnc

  subroutine log_cdfgam()
    ! which, p, q, x, shape, scale
    real(kind=rk), parameter :: r(6, 32) = reshape((/ &
      1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 2.0_rk, 1.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 1.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 0.0_rk, 1.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 2.0_rk, -1.0_rk, &
      2.0_rk, 0.1_rk, auto, 0.0_rk, 2.0_rk, 3.0_rk, &
      2.0_rk, 0.9_rk, auto, 0.0_rk, 2.0_rk, 3.0_rk, &
      2.0_rk, 1.0e-200_rk, auto, 0.0_rk, 0.1_rk, 1.0_rk, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, 1.0e6_rk, 1.0_rk, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, -2.0_rk, 1.0_rk, &
      2.0_rk, 0.5_rk, 0.6_rk, 0.0_rk, 2.0_rk, 1.0_rk, &
      3.0_rk, 0.1_rk, auto, 5.0_rk, 0.0_rk, 2.0_rk, &
      3.0_rk, 0.9_rk, auto, 5.0_rk, 0.0_rk, 2.0_rk, &
      3.0_rk, 0.5_rk, auto, 100.0_rk, 0.0_rk, 1.0_rk, &
      3.0_rk, 1.0e-300_rk, auto, 5.0_rk, 0.0_rk, 2.0_rk, &
      3.0_rk, 0.9999999_rk, auto, 1.0e-3_rk, 0.0_rk, 1.0_rk, &
      3.0_rk, 0.5_rk, auto, -5.0_rk, 0.0_rk, 2.0_rk, &
      3.0_rk, 0.5_rk, auto, 5.0_rk, 0.0_rk, 0.0_rk, &
      3.0_rk, -0.5_rk, auto, 5.0_rk, 0.0_rk, 2.0_rk, &
      3.0_rk, 0.5_rk, 0.6_rk, 5.0_rk, 0.0_rk, 2.0_rk, &
      4.0_rk, 0.1_rk, auto, 5.0_rk, 2.0_rk, 0.0_rk, &
      4.0_rk, 0.9_rk, auto, 5.0_rk, 2.0_rk, 0.0_rk, &
      4.0_rk, 0.5_rk, auto, 5.0_rk, 1.0e10_rk, 0.0_rk, &
      4.0_rk, 0.5_rk, auto, -5.0_rk, 2.0_rk, 0.0_rk, &
      4.0_rk, 0.5_rk, auto, 5.0_rk, 0.0_rk, 0.0_rk, &
      2.0_rk, 1.0e-10_rk, auto, 0.0_rk, 1.0e-3_rk, 1.0_rk, &
      4.0_rk, 1.0e-10_rk, auto, 5.0_rk, 1.0e-3_rk, 0.0_rk , &
      2.0_rk, 1.5_rk, auto, 0.0_rk, 2.0_rk, 3.0_rk, &
      2.0_rk, 0.5_rk, -0.5_rk, 0.0_rk, 2.0_rk, 3.0_rk, &
      2.0_rk, 0.5_rk, 1.5_rk, 0.0_rk, 2.0_rk, 3.0_rk, &
      3.0_rk, 0.5_rk, auto, 1.0e300_rk, 0.0_rk, 1.0_rk, &
      3.0_rk, 0.5_rk, auto, 1.5e300_rk, 0.0_rk, 1.0_rk, &
      3.0_rk, 0.7_rk, auto, 1.5e300_rk, 0.0_rk, 1.0_rk /), (/ 6, 32 /))
    real(kind=rk) :: v(5), a(5), bound
    integer :: unit, i, which, status
    open(newunit=unit, file='tests/data/cdfgam_calls.csv', status='replace', action='write')
    write(unit, '(a)') '# which, p, q, x, shape, scale, status, bound, p, q, x, shape, scale (after)'
    do i = 1, size(r, 2)
      which = int(r(1, i))
      v = r(2:6, i)
      call resolve_pq(v(1), v(2))
      a = v
      status = 0
      bound = 0.0_rk
      call cdfgam(which, a(1), a(2), a(3), a(4), a(5), status, bound)
      call emit(unit, which, v, status, bound, a)
    end do
    close(unit)
  end subroutine log_cdfgam

  subroutine log_cdfnbn()
    ! which, p, q, f, s, pr, ompr
    real(kind=rk), parameter :: r(7, 35) = reshape((/ &
      1.0_rk, 0.0_rk, 0.0_rk, 3.0_rk, 5.0_rk, 0.5_rk, auto, &
      1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 5.0_rk, 0.5_rk, auto, &
      1.0_rk, 0.0_rk, 0.0_rk, 3.0_rk, 5.0_rk, 1.0_rk, auto, &
      1.0_rk, 0.0_rk, 0.0_rk, 3.0_rk, 5.0_rk, 1.5_rk, auto, &
      2.0_rk, 0.2_rk, auto, 0.0_rk, 5.0_rk, 0.5_rk, auto, &
      2.0_rk, 0.8_rk, auto, 0.0_rk, 5.0_rk, 0.5_rk, auto, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, 50.0_rk, 0.2_rk, auto, &
      2.0_rk, 1.0e-30_rk, auto, 0.0_rk, 5.0_rk, 0.5_rk, auto, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, 5.0_rk, -0.5_rk, auto, &
      2.0_rk, 0.5_rk, 0.6_rk, 0.0_rk, 5.0_rk, 0.5_rk, auto, &
      3.0_rk, 0.2_rk, auto, 3.0_rk, 0.0_rk, 0.5_rk, auto, &
      3.0_rk, 0.8_rk, auto, 3.0_rk, 0.0_rk, 0.5_rk, auto, &
      3.0_rk, 0.5_rk, auto, 30.0_rk, 0.0_rk, 0.3_rk, auto, &
      3.0_rk, 1.0e-300_rk, auto, 3.0_rk, 0.0_rk, 0.5_rk, auto, &
      3.0_rk, 0.9999999_rk, auto, 3.0_rk, 0.0_rk, 0.5_rk, auto, &
      3.0_rk, 0.5_rk, auto, 3.0_rk, 0.0_rk, 1.5_rk, auto, &
      3.0_rk, -0.5_rk, auto, 3.0_rk, 0.0_rk, 0.5_rk, auto, &
      4.0_rk, 0.2_rk, auto, 3.0_rk, 5.0_rk, 0.0_rk, auto, &
      4.0_rk, 0.8_rk, auto, 3.0_rk, 5.0_rk, 0.0_rk, auto, &
      4.0_rk, 0.5_rk, auto, 30.0_rk, 10.0_rk, 0.0_rk, auto, &
      4.0_rk, 1.0e-300_rk, auto, 3.0_rk, 5.0_rk, 0.0_rk, auto, &
      4.0_rk, 0.9999999_rk, auto, 3.0_rk, 5.0_rk, 0.0_rk, auto, &
      4.0_rk, 0.5_rk, auto, 0.0_rk, 0.0_rk, 0.0_rk, auto, &
      4.0_rk, 0.5_rk, 0.6_rk, 3.0_rk, 5.0_rk, 0.0_rk, auto, &
      4.0_rk, 1.5_rk, auto, 3.0_rk, 5.0_rk, 0.0_rk, auto, &
      4.0_rk, 0.5_rk, auto, 3.0_rk, 1.0_rk, 0.0_rk, auto, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, 1.0_rk, 0.9_rk, auto , &
      2.0_rk, 0.5_rk, -0.5_rk, 0.0_rk, 5.0_rk, 0.5_rk, auto, &
      2.0_rk, 0.5_rk, 1.5_rk, 0.0_rk, 5.0_rk, 0.5_rk, auto, &
      3.0_rk, 0.5_rk, auto, -1.0_rk, 0.0_rk, 0.5_rk, auto, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, -1.0_rk, 0.5_rk, auto, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, 5.0_rk, 1.0e-300_rk, auto, &
      3.0_rk, 0.5_rk, auto, 3.0_rk, 0.0_rk, 1.0_rk, auto, &
      3.0_rk, 0.7_rk, auto, 5.0_rk, 0.0_rk, 1.0_rk, auto, &
      3.0_rk, 1.0_rk, auto, 3.0_rk, 0.0_rk, 0.5_rk, auto /), (/ 7, 35 /))
    real(kind=rk) :: v(6), a(6), bound
    integer :: unit, i, which, status
    open(newunit=unit, file='tests/data/cdfnbn_calls.csv', status='replace', action='write')
    write(unit, '(a)') '# which, p, q, f, s, pr, ompr, status, bound, p, q, f, s, pr, ompr (after)'
    do i = 1, size(r, 2)
      which = int(r(1, i))
      v = r(2:7, i)
      call resolve_pq(v(1), v(2))
      if (v(6) == auto) v(6) = 1.0_rk - v(5)
      a = v
      status = 0
      bound = 0.0_rk
      call cdfnbn(which, a(1), a(2), a(3), a(4), a(5), a(6), status, bound)
      call emit(unit, which, v, status, bound, a)
    end do
    close(unit)
  end subroutine log_cdfnbn

  subroutine log_cdfnor()
    ! which, p, q, x, mean, sd
    real(kind=rk), parameter :: r(6, 22) = reshape((/ &
      1.0_rk, 0.0_rk, 0.0_rk, 1.96_rk, 0.0_rk, 1.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, -40.0_rk, 0.0_rk, 1.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, &
      2.0_rk, 0.025_rk, auto, 0.0_rk, 1.0_rk, 2.0_rk, &
      2.0_rk, 0.975_rk, auto, 0.0_rk, 1.0_rk, 2.0_rk, &
      2.0_rk, 1.0e-300_rk, auto, 0.0_rk, 0.0_rk, 1.0_rk, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, 3.0_rk, -1.0_rk, &
      2.0_rk, -0.5_rk, auto, 0.0_rk, 3.0_rk, 1.0_rk, &
      2.0_rk, 0.5_rk, 1.5_rk, 0.0_rk, 3.0_rk, 1.0_rk, &
      2.0_rk, 0.5_rk, 0.6_rk, 0.0_rk, 3.0_rk, 1.0_rk, &
      3.0_rk, 0.025_rk, auto, 1.0_rk, 0.0_rk, 2.0_rk, &
      3.0_rk, 0.975_rk, auto, 1.0_rk, 0.0_rk, 2.0_rk, &
      3.0_rk, 0.5_rk, auto, 1.0_rk, 0.0_rk, 0.0_rk, &
      3.0_rk, 0.5_rk, 0.6_rk, 1.0_rk, 0.0_rk, 2.0_rk, &
      4.0_rk, 0.025_rk, auto, -3.0_rk, 1.0_rk, 0.0_rk, &
      4.0_rk, 0.975_rk, auto, 5.0_rk, 1.0_rk, 0.0_rk, &
      4.0_rk, 0.5_rk, auto, 1.0_rk, 1.0_rk, 0.0_rk, &
      4.0_rk, 0.975_rk, auto, -3.0_rk, 1.0_rk, 0.0_rk, &
      4.0_rk, 0.975_rk, auto, 1.0_rk, 1.0_rk, 0.0_rk, &
      4.0_rk, 1.5_rk, auto, 1.0_rk, 1.0_rk, 0.0_rk, &
      4.0_rk, 0.5_rk, 0.6_rk, 1.0_rk, 1.0_rk, 0.0_rk, &
      4.0_rk, 0.5_rk, -0.5_rk, 1.0_rk, 1.0_rk, 0.0_rk /), (/ 6, 22 /))
    real(kind=rk) :: v(5), a(5), bound
    integer :: unit, i, which, status
    open(newunit=unit, file='tests/data/cdfnor_calls.csv', status='replace', action='write')
    write(unit, '(a)') '# which, p, q, x, mean, sd, status, bound, p, q, x, mean, sd (after)'
    do i = 1, size(r, 2)
      which = int(r(1, i))
      v = r(2:6, i)
      call resolve_pq(v(1), v(2))
      a = v
      status = 0
      bound = 0.0_rk
      call cdfnor(which, a(1), a(2), a(3), a(4), a(5), status, bound)
      call emit(unit, which, v, status, bound, a)
    end do
    close(unit)
  end subroutine log_cdfnor

  subroutine log_cdfpoi()
    ! which, p, q, s, xlam
    real(kind=rk), parameter :: r(5, 23) = reshape((/ &
      1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, 3.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 285.0_rk, 200.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 2.0_rk, -3.0_rk, &
      2.0_rk, 0.2_rk, auto, 0.0_rk, 3.0_rk, &
      2.0_rk, 0.8_rk, auto, 0.0_rk, 3.0_rk, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, 1.0e6_rk, &
      2.0_rk, 1.0e-300_rk, auto, 0.0_rk, 3.0_rk, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, -3.0_rk, &
      2.0_rk, 0.5_rk, 0.6_rk, 0.0_rk, 3.0_rk, &
      3.0_rk, 0.2_rk, auto, 3.0_rk, 0.0_rk, &
      3.0_rk, 0.8_rk, auto, 3.0_rk, 0.0_rk, &
      3.0_rk, 0.5_rk, auto, 1000.0_rk, 0.0_rk, &
      3.0_rk, 1.0e-300_rk, auto, 3.0_rk, 0.0_rk, &
      3.0_rk, 0.9999999_rk, auto, 0.0_rk, 0.0_rk, &
      3.0_rk, 1.0_rk, auto, 0.0_rk, 0.0_rk, &
      3.0_rk, -0.5_rk, auto, 3.0_rk, 0.0_rk, &
      3.0_rk, 0.5_rk, -0.5_rk, 3.0_rk, 0.0_rk, &
      3.0_rk, 0.5_rk, 0.4_rk, 3.0_rk, 0.0_rk , &
      3.0_rk, 1.5_rk, auto, 3.0_rk, 0.0_rk, &
      3.0_rk, 0.5_rk, 1.5_rk, 3.0_rk, 0.0_rk, &
      3.0_rk, 0.5_rk, auto, -1.0_rk, 0.0_rk, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, 1.0e305_rk /), (/ 5, 23 /))
    real(kind=rk) :: v(4), a(4), bound
    integer :: unit, i, which, status
    open(newunit=unit, file='tests/data/cdfpoi_calls.csv', status='replace', action='write')
    write(unit, '(a)') '# which, p, q, s, xlam, status, bound, p, q, s, xlam (after)'
    do i = 1, size(r, 2)
      which = int(r(1, i))
      v = r(2:5, i)
      call resolve_pq(v(1), v(2))
      a = v
      status = 0
      bound = 0.0_rk
      call cdfpoi(which, a(1), a(2), a(3), a(4), status, bound)
      call emit(unit, which, v, status, bound, a)
    end do
    close(unit)
  end subroutine log_cdfpoi

  subroutine log_cdft()
    ! which, p, q, t, df
    real(kind=rk), parameter :: r(5, 24) = reshape((/ &
      1.0_rk, 0.0_rk, 0.0_rk, 2.228_rk, 10.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, -40.0_rk, 3.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 3.0_rk, &
      1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, &
      2.0_rk, 0.025_rk, auto, 0.0_rk, 10.0_rk, &
      2.0_rk, 0.975_rk, auto, 0.0_rk, 10.0_rk, &
      2.0_rk, 1.0e-10_rk, auto, 0.0_rk, 1.0_rk, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, 0.5_rk, &
      2.0_rk, 0.5_rk, auto, 0.0_rk, -1.0_rk, &
      2.0_rk, -0.5_rk, auto, 0.0_rk, 5.0_rk, &
      2.0_rk, 0.5_rk, 0.6_rk, 0.0_rk, 5.0_rk, &
      3.0_rk, 0.025_rk, auto, -2.228_rk, 0.0_rk, &
      3.0_rk, 0.975_rk, auto, 2.228_rk, 0.0_rk, &
      3.0_rk, 0.9_rk, auto, 1.5_rk, 0.0_rk, &
      3.0_rk, 0.6_rk, auto, 0.01_rk, 0.0_rk, &
      3.0_rk, 0.999_rk, auto, 1.0_rk, 0.0_rk, &
      3.0_rk, 0.6_rk, auto, 100.0_rk, 0.0_rk, &
      3.0_rk, 0.5_rk, auto, 1.0_rk, 0.0_rk, &
      3.0_rk, 1.5_rk, auto, 1.0_rk, 0.0_rk, &
      3.0_rk, 0.5_rk, 0.6_rk, 1.0_rk, 0.0_rk, &
      3.0_rk, 0.5_rk, -0.5_rk, 1.0_rk, 0.0_rk , &
      3.0_rk, 0.5_rk, 1.5_rk, 1.0_rk, 0.0_rk , &
      2.0_rk, 1.0e-100_rk, 1.0_rk, 0.0_rk, 1.0_rk , &
      2.0_rk, 1.0_rk, 1.0e-100_rk, 0.0_rk, 1.0_rk /), (/ 5, 24 /))
    real(kind=rk) :: v(4), a(4), bound
    integer :: unit, i, which, status
    open(newunit=unit, file='tests/data/cdft_calls.csv', status='replace', action='write')
    write(unit, '(a)') '# which, p, q, t, df, status, bound, p, q, t, df (after)'
    do i = 1, size(r, 2)
      which = int(r(1, i))
      v = r(2:5, i)
      call resolve_pq(v(1), v(2))
      a = v
      status = 0
      bound = 0.0_rk
      call cdft(which, a(1), a(2), a(3), a(4), status, bound)
      call emit(unit, which, v, status, bound, a)
    end do
    close(unit)
  end subroutine log_cdft

end program gen_dispatcher_calls
