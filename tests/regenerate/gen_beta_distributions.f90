! Beta, Student's t, F distribution CDF reference tables.

program gen_beta_distributions
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  external :: cumbet, cumt, cumf

  real(kind=rk), parameter :: ab(5) = (/ 0.5_rk, 1.0_rk, 2.0_rk, 5.0_rk, 30.0_rk /)
  real(kind=rk), parameter :: dfs_t(6) = (/ 1.0_rk, 2.0_rk, 5.0_rk, 10.0_rk, 30.0_rk, 100.0_rk /)
  real(kind=rk), parameter :: dfs_f(5) = (/ 1.0_rk, 2.0_rk, 5.0_rk, 10.0_rk, 30.0_rk /)
  integer :: i, j, unit
  real(kind=rk) :: a, b, x, y, cum, ccum, t, df, dfn, dfd, fx

  ! Beta CDF
  open(newunit=unit, file='tests/data/beta_cdf.csv', status='replace', action='write')
  write(unit, '(a)') '# a, b, x, cdf, sf'
  do i = 1, size(ab)
    do j = 1, size(ab)
      a = ab(i)
      b = ab(j)
      x = 0.05_rk
      do while (x < 1.0_rk)
        y = 1.0_rk - x
        call cumbet(x, y, a, b, cum, ccum)
        call putval(unit, a, .false.)
        call putval(unit, b, .false.)
        call putval(unit, x, .false.)
        call putval(unit, cum, .false.)
        call putval(unit, ccum, .true.)
        x = x + 0.05_rk
      end do
    end do
  end do
  close(unit)

  ! Student's t CDF
  open(newunit=unit, file='tests/data/students_t_cdf.csv', status='replace', action='write')
  write(unit, '(a)') '# df, t, cdf, sf'
  do i = 1, size(dfs_t)
    df = dfs_t(i)
    t = -8.0_rk
    do while (t <= 8.0_rk + 1.0e-12_rk)
      call cumt(t, df, cum, ccum)
      call putval(unit, df, .false.)
      call putval(unit, t, .false.)
      call putval(unit, cum, .false.)
      call putval(unit, ccum, .true.)
      t = t + 0.25_rk
    end do
  end do
  close(unit)

  ! F CDF
  open(newunit=unit, file='tests/data/f_cdf.csv', status='replace', action='write')
  write(unit, '(a)') '# dfn, dfd, f, cdf, sf'
  do i = 1, size(dfs_f)
    do j = 1, size(dfs_f)
      dfn = dfs_f(i)
      dfd = dfs_f(j)
      fx = 0.05_rk
      do while (fx <= 10.0_rk + 1.0e-12_rk)
        call cumf(fx, dfn, dfd, cum, ccum)
        call putval(unit, dfn, .false.)
        call putval(unit, dfd, .false.)
        call putval(unit, fx, .false.)
        call putval(unit, cum, .false.)
        call putval(unit, ccum, .true.)
        fx = fx + 0.25_rk
      end do
    end do
  end do
  close(unit)

  write(0, '(a)') 'wrote 3 tables under tests/data/'

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
end program gen_beta_distributions
