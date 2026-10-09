! Reverse-communication traces of dinvr and dzror in cdflib.f90.
!
! Each case drives dstinv/dinvr or dstzr/dzror with a fixed objective
! function, chosen to reach every branch of the two state machines:
! monotone functions with the root inside and outside the bracket, step
! functions, plateaus, NaN and infinite values at chosen points, and
! non-monotone functions. Every evaluation request is written as
! "need <x> <fx>" with both values as IEEE binary64 bit patterns in
! hexadecimal, followed by the outcome: "conv <x>", "fail <qleft> <qhi> <x>",
! or "cap" if the case was cut after 5000 evaluations. Each case starts
! with "dinvr <small> <big> <start> <abstol> <reltol>" or "dzror <xlo>
! <xhi> <abstol> <reltol>"; the dinvr step parameters are 0.5, 0.5 and 5,
! as in every cdf* routine.
!
! Companion Rust test: the solver_traces unit test in src/search/mod.rs,
! which replays the recorded fx values and checks every requested x and
! the outcome bit for bit.

module objf
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer, parameter :: rk = kind(1.0d0)
contains
  function f(id, x) result(y)
    integer, intent(in) :: id
    real(rk), intent(in) :: x
    real(rk) :: y
    select case (id)
    case (1); y = x * x * x - 8.0d0
    case (2); y = 1.0d0 / x - 0.25d0
    case (3); y = x - 3.0d0
    case (4); y = x + 1.0d0
    case (5); y = x - 100.0d0
    case (6); y = -x - 1.0d0
    case (7); y = 100.0d0 - x
    case (8); if (x < 0.7d0) then; y = -1.0d0; else; y = 2.0d0; end if
    case (9); if (x < 0.7d0) then; y = 1.0d0; else; y = -2.0d0; end if
    case (10); y = 0.0d0
    case (11); y = 1.0d0
    case (12); y = -1.0d0
    case (13); y = ieee_value(y, ieee_quiet_nan)
    case (14); if (x <= 0.001d0) then; y = ieee_value(y, ieee_quiet_nan); else; y = x - 5.0d0; end if
    case (15); if (x >= 99.9d0) then; y = ieee_value(y, ieee_quiet_nan); else; y = x - 5.0d0; end if
    case (16); if (2.0d0 < x .and. x < 3.0d0) then; y = ieee_value(y, ieee_quiet_nan); else; y = x - 2.5d0; end if
    case (17); if (x > 10.0d0) then; y = ieee_value(y, ieee_quiet_nan); else; y = x - 50.0d0; end if
    case (18)
      if (x < 2.0d0) then; y = x - 2.0d0
      else if (x > 4.0d0) then; y = x - 4.0d0
      else; y = 0.0d0; end if
    case (19); y = (x - 0.3d0) * (x - 0.3d0) * (x - 0.3d0)
    case (20); y = x / (1.0d0 + x) - 0.999999d0
    case (21); y = 1.0d0 / (1.0d0 + x) - 1.0d-8
    case (22); if (x > 50.0d0) then; y = ieee_value(y, ieee_positive_inf); else; y = x - 60.0d0; end if
    case (23); if (x < 1.0d0) then; y = ieee_value(y, ieee_negative_inf); else; y = x - 2.0d0; end if
    case (24); if (x == 3.0d0) then; y = sign(0.0d0, -1.0d0); else; y = x - 3.0d0; end if
    case (25); y = (x - 1.0d0) * (x - 5.0d0) * (x - 9.0d0)
    case (26); y = x - 1.0d-12
    case (27); y = x - 99.999d0
    case (28); y = x - 0.3d0
    case (29); y = x * x - 0.5d0
    case (30); if (0.4d0 < x .and. x < 0.6d0) then; y = ieee_value(y, ieee_quiet_nan); else; y = x - 0.5d0; end if
    case (31); y = -(x * x) + 0.5d0
    case (32)
      if (x < 0.25d0) then; y = -1.0d0
      else if (x < 0.75d0) then; y = ieee_value(y, ieee_quiet_nan)
      else; y = 1.0d0; end if
    case (33); if (x > 0.0d0) then; y = 1.0d0; else; y = -1.0d0; end if
    case (34); if (x >= 99.9d0) then; y = ieee_value(y, ieee_quiet_nan); else; y = 5.0d0 + x; end if
    case default; stop 2
    end select
  end function
  function hx(x) result(s)
    real(rk), intent(in) :: x
    character(len=16) :: s
    write (s, '(z16.16)') transfer(x, 0_8)
  end function
  function tf(b) result(s)
    logical, intent(in) :: b
    character(len=1) :: s
    if (b) then; s = 'T'; else; s = 'F'; end if
  end function
end module

program gen_solver_traces
  use objf
  implicit none
  integer :: unit
  integer, parameter :: nd = 46, nz = 30
  real(rk) :: dsmall(nd), dbig(nd), dstart(nd)
  integer :: dfid(nd)
  real(rk) :: zlo(nz), zhi(nz)
  integer :: zfid(nz)
  integer, parameter :: ztol(6) = [19, 28, 29, 8, 26, 31]
  real(rk) :: atols(2), x, fx, xlo, xhi
  integer :: i, j, status, it
  logical :: qleft, qhi, done

  open(newunit=unit, file='tests/data/solver_traces.txt', status='replace', action='write')
  atols = [1.0d-10, 1.0d-50]
  call setd(1, 0.0d0, 100.0d0, 1.0d0, 1)
  call setd(2, 0.0d0, 100.0d0, 0.0d0, 1)
  call setd(3, 0.0d0, 100.0d0, 100.0d0, 1)
  call setd(4, -1.0d300, 1.0d300, 0.0d0, 1)
  call setd(5, 0.01d0, 1000.0d0, 10.0d0, 2)
  call setd(6, 0.0d0, 10.0d0, 3.0d0, 3)
  call setd(7, 1.0d0, 10.0d0, 5.0d0, 4)
  call setd(8, 1.0d0, 10.0d0, 5.0d0, 5)
  call setd(9, 1.0d0, 10.0d0, 5.0d0, 6)
  call setd(10, 1.0d0, 10.0d0, 5.0d0, 7)
  call setd(11, 0.0d0, 1.0d0, 0.5d0, 8)
  call setd(12, 0.0d0, 1.0d0, 0.9d0, 8)
  call setd(13, 0.0d0, 1.0d0, 0.1d0, 8)
  call setd(14, 0.0d0, 1.0d0, 0.0d0, 8)
  call setd(15, 0.0d0, 1.0d0, 1.0d0, 8)
  call setd(16, 0.0d0, 1.0d0, 0.5d0, 9)
  call setd(17, 0.0d0, 1.0d0, 0.5d0, 10)
  call setd(18, 0.0d0, 1.0d0, 0.5d0, 11)
  call setd(19, 0.0d0, 1.0d0, 0.5d0, 12)
  call setd(20, 0.0d0, 1.0d0, 0.5d0, 13)
  call setd(21, 0.0d0, 100.0d0, 1.0d0, 14)
  call setd(22, 0.0d0, 100.0d0, 1.0d0, 15)
  call setd(23, 0.0d0, 10.0d0, 1.0d0, 16)
  call setd(24, 0.0d0, 10.0d0, 5.0d0, 16)
  call setd(25, 0.0d0, 100.0d0, 1.0d0, 17)
  call setd(26, 0.0d0, 10.0d0, 1.0d0, 18)
  call setd(27, 0.0d0, 10.0d0, 8.0d0, 18)
  call setd(28, 0.0d0, 10.0d0, 3.0d0, 18)
  call setd(29, 0.0d0, 1.0d0, 0.9d0, 19)
  call setd(30, 0.0d0, 1.0d300, 5.0d0, 20)
  call setd(31, 0.0d0, 1.0d300, 5.0d0, 21)
  call setd(32, 0.0d0, 100.0d0, 1.0d0, 22)
  call setd(33, 0.0d0, 100.0d0, 5.0d0, 23)
  call setd(34, 0.0d0, 10.0d0, 3.0d0, 24)
  call setd(35, 0.0d0, 10.0d0, 4.5d0, 25)
  call setd(36, 0.0d0, 10.0d0, 0.5d0, 25)
  call setd(37, 0.0d0, 1.0d0, 0.5d0, 26)
  call setd(38, 0.0d0, 100.0d0, 1.0d0, 27)
  call setd(39, 0.0d0, 100.0d0, 100.0d0, 27)
  call setd(40, 0.0d0, 100.0d0, 0.0d0, 27)
  call setd(41, 0.0d0, 1.0d0, 0.5d0, 30)
  call setd(42, 0.0d0, 1.0d0, 0.1d0, 30)
  call setd(43, 0.0d0, 1.0d0, 0.9d0, 31)
  call setd(44, 0.0d0, 1.0d0, 0.5d0, 32)
  call setd(45, 0.0d0, 1.0d0, 0.1d0, 32)
  call setd(46, 0.0d0, 100.0d0, 1.0d0, 34)

  call setz(1, 0.0d0, 1.0d0, 8)
  call setz(2, 0.0d0, 1.0d0, 9)
  call setz(3, 0.0d0, 1.0d0, 19)
  call setz(4, 0.0d0, 1.0d0, 26)
  call setz(5, 0.0d0, 1.0d0, 10)
  call setz(6, 0.0d0, 1.0d0, 11)
  call setz(7, 0.0d0, 1.0d0, 12)
  call setz(8, 0.0d0, 1.0d0, 13)
  call setz(9, 1.0d0, 0.0d0, 8)
  call setz(10, 0.0d0, 10.0d0, 25)
  call setz(11, 0.0d0, 10.0d0, 16)
  call setz(12, 0.0d0, 10.0d0, 18)
  call setz(13, 0.0d0, 1.0d0, 28)
  call setz(14, 0.0d0, 1.0d0, 29)
  call setz(15, 0.0d0, 10.0d0, 3)
  call setz(16, 3.0d0, 10.0d0, 3)
  call setz(17, 0.0d0, 3.0d0, 3)
  call setz(18, 0.0d0, 100.0d0, 22)
  call setz(19, 0.0d0, 100.0d0, 23)
  call setz(20, 0.0d0, 100.0d0, 14)
  call setz(21, 0.0d0, 100.0d0, 15)
  call setz(22, 0.0d0, 1.0d0, 30)
  call setz(23, 0.0d0, 1.0d0, 31)
  call setz(24, 0.0d0, 1.0d0, 32)
  call setz(25, 1.0d0, 0.0d0, 29)
  call setz(26, 0.0d0, 10.0d0, 4)
  call setz(27, 0.0d0, 10.0d0, 6)
  call setz(28, 0.0d0, 10.0d0, 7)
  call setz(29, 0.0d0, 10.0d0, 5)
  call setz(30, ieee_value(x, ieee_negative_inf), ieee_value(x, ieee_positive_inf), 33)

  do i = 1, nd
    do j = 1, 2
      write (unit, '(a,5(1x,a))') 'dinvr', hx(dsmall(i)), hx(dbig(i)), hx(dstart(i)), &
        hx(atols(j)), hx(1.0d-8)
      call dstinv(dsmall(i), dbig(i), 0.5d0, 0.5d0, 5.0d0, atols(j), 1.0d-8)
      x = dstart(i)
      status = 0
      fx = 0.0d0
      qleft = .false.
      qhi = .false.
      call dinvr(status, x, fx, qleft, qhi)
      done = .false.
      do it = 1, 5000
        if (status /= 1) then
          done = .true.
          exit
        end if
        fx = f(dfid(i), x)
        write (unit, '(a,1x,a,1x,a)') 'need', hx(x), hx(fx)
        call dinvr(status, x, fx, qleft, qhi)
      end do
      if (.not. done) then
        if (status /= 1) then
          done = .true.
        end if
      end if
      if (.not. done) then
        write (unit, '(a)') 'cap'
      else if (status == 0) then
        write (unit, '(a,1x,a)') 'conv', hx(x)
      else
        write (unit, '(a,1x,a,1x,a,1x,a)') 'fail', tf(qleft), tf(qhi), hx(x)
      end if
    end do
  end do

  do i = 1, nz
    do j = 1, 2
      write (unit, '(a,4(1x,a))') 'dzror', hx(zlo(i)), hx(zhi(i)), hx(atols(j)), hx(1.0d-8)
      call dstzr(zlo(i), zhi(i), atols(j), 1.0d-8)
      status = 0
      fx = 0.0d0
      qleft = .false.
      qhi = .false.
      call dzror(status, x, fx, xlo, xhi, qleft, qhi)
      done = .false.
      do it = 1, 5000
        if (status /= 1) then
          done = .true.
          exit
        end if
        fx = f(zfid(i), x)
        write (unit, '(a,1x,a,1x,a)') 'need', hx(x), hx(fx)
        call dzror(status, x, fx, xlo, xhi, qleft, qhi)
      end do
      if (.not. done) then
        if (status /= 1) then
          done = .true.
        end if
      end if
      if (.not. done) then
        write (unit, '(a)') 'cap'
      else if (status == 0) then
        write (unit, '(a,1x,a)') 'conv', hx(x)
      else
        write (unit, '(a,1x,a,1x,a,1x,a)') 'fail', tf(qleft), tf(qhi), hx(x)
      end if
    end do
  end do

  ! dzror with zero tolerances, where a step can leave b unchanged and the
  ! divided difference guard d == a is reached. Each case is cut after 200
  ! evaluations.
  do i = 1, size(ztol)
    write (unit, '(a,4(1x,a))') 'dzror', hx(0.0d0), hx(1.0d0), hx(0.0d0), hx(0.0d0)
    call dstzr(0.0d0, 1.0d0, 0.0d0, 0.0d0)
    status = 0
    fx = 0.0d0
    qleft = .false.
    qhi = .false.
    call dzror(status, x, fx, xlo, xhi, qleft, qhi)
    done = .false.
    do it = 1, 200
      if (status /= 1) then
        done = .true.
        exit
      end if
      fx = f(ztol(i), x)
      write (unit, '(a,1x,a,1x,a)') 'need', hx(x), hx(fx)
      call dzror(status, x, fx, xlo, xhi, qleft, qhi)
    end do
    if (.not. done .and. status /= 1) done = .true.
    if (.not. done) then
      write (unit, '(a)') 'cap'
    else if (status == 0) then
      write (unit, '(a,1x,a)') 'conv', hx(x)
    else
      write (unit, '(a,1x,a,1x,a,1x,a)') 'fail', tf(qleft), tf(qhi), hx(x)
    end if
  end do

  close(unit)

contains
  subroutine setd(k, s, b, st, id)
    integer, intent(in) :: k, id
    real(rk), intent(in) :: s, b, st
    dsmall(k) = s; dbig(k) = b; dstart(k) = st; dfid(k) = id
  end subroutine
  subroutine setz(k, a, b, id)
    integer, intent(in) :: k, id
    real(rk), intent(in) :: a, b
    zlo(k) = a; zhi(k) = b; zfid(k) = id
  end subroutine
end program gen_solver_traces
