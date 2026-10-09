! Reference tables for beta_log, beta_inc, and dbetrm.

program gen_beta_kernels
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), external :: beta_log, dbetrm
  external :: beta_inc

  real(kind=rk), parameter :: ab(15) = (/ &
    0.1_rk, 0.5_rk, 0.999999999999_rk, 1.0_rk, 1.000000000001_rk, &
    1.999999999999_rk, 2.0_rk, 2.000000000001_rk, 5.0_rk, 7.999999999999_rk, &
    8.0_rk, 8.000000000001_rk, 10.0_rk, 30.0_rk, 100.0_rk /)
  real(kind=rk), parameter :: ab_inc(15) = (/ &
    0.5_rk, 0.999999999999_rk, 1.0_rk, 1.000000000001_rk, &
    1.999999999999_rk, 2.0_rk, 2.000000000001_rk, 5.0_rk, 7.999999999999_rk, &
    8.0_rk, 8.000000000001_rk, 10.0_rk, 15.0_rk, 30.0_rk, 100.0_rk /)
  real(kind=rk), parameter :: x_special(7) = (/ &
    1.0e-12_rk, 1.0e-6_rk, 0.05_rk, 0.5_rk, 0.95_rk, 0.999999_rk, 0.999999999999_rk /)
  integer :: i, j, k, unit, ierr
  real(kind=rk) :: a, b, x, y, w, w1

  ! beta_log
  open(newunit=unit, file='tests/data/beta_log.csv', status='replace', action='write')
  write(unit, '(a)') '# a, b, ln_beta(a, b)'
  do i = 1, size(ab)
    do j = 1, size(ab)
      a = ab(i)
      b = ab(j)
      call putval(unit, a, .false.)
      call putval(unit, b, .false.)
      call putval(unit, beta_log(a, b), .true.)
    end do
  end do
  close(unit)

  ! beta_inc
  open(newunit=unit, file='tests/data/beta_inc.csv', status='replace', action='write')
  write(unit, '(a)') '# a, b, x, P, Q'
  do i = 1, size(ab_inc)
    do j = 1, size(ab_inc)
      a = ab_inc(i)
      b = ab_inc(j)
      x = 0.05_rk
      do while (x < 1.0_rk)
        y = 1.0_rk - x
        w = 0.0_rk
        w1 = 0.0_rk
        ierr = 0
        call beta_inc(a, b, x, y, w, w1, ierr)
        if (ierr == 0) then
          call putval(unit, a, .false.)
          call putval(unit, b, .false.)
          call putval(unit, x, .false.)
          call putval(unit, w, .false.)
          call putval(unit, w1, .true.)
        end if
        x = x + 0.05_rk
      end do
      do k = 1, size(x_special)
        x = x_special(k)
        y = 1.0_rk - x
        w = 0.0_rk
        w1 = 0.0_rk
        ierr = 0
        call beta_inc(a, b, x, y, w, w1, ierr)
        if (ierr == 0) then
          call putval(unit, a, .false.)
          call putval(unit, b, .false.)
          call putval(unit, x, .false.)
          call putval(unit, w, .false.)
          call putval(unit, w1, .true.)
        end if
      end do
    end do
  end do
  close(unit)

  ! dbetrm
  call write_dbetrm()

  write(0, '(a)') 'wrote 3 tables under tests/data/'

contains
  subroutine write_dbetrm()
    integer :: u, i, j
    real(kind=rk) :: a, b
    real(kind=rk), parameter :: vals(*) = (/ &
      0.5_rk, 1.0_rk, 2.0_rk, 5.0_rk, 10.0_rk, 20.0_rk, 50.0_rk, 100.0_rk /)
    open(newunit=u, file='tests/data/dbetrm.csv', &
         status='replace', action='write')
    write(u, '(a)') '# a, b, dbetrm(a, b)'
    do i = 1, size(vals)
      do j = 1, size(vals)
        a = vals(i)
        b = vals(j)
        call putval(u, a, .false.)
        call putval(u, b, .false.)
        call putval(u, dbetrm(a, b), .true.)
      end do
    end do
    close(u)
  end subroutine write_dbetrm

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

end program gen_beta_kernels
