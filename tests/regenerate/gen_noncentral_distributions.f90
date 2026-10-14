! Noncentral chi-squared and noncentral F CDF tables.

program gen_noncentral_distributions
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  external :: cumchn, cumfnc

  real(kind=rk), parameter :: dfs(5)  = (/ 1.0_rk, 2.0_rk, 5.0_rk, 10.0_rk, 30.0_rk /)
  real(kind=rk), parameter :: ncps(8) = (/ &
    0.0_rk, 1.0e-12_rk, 1.0e-10_rk, 1.0e-9_rk, 0.5_rk, 2.0_rk, 10.0_rk, 50.0_rk /)
  real(kind=rk), parameter :: dfns(4) = (/ 2.0_rk, 5.0_rk, 10.0_rk, 30.0_rk /)
  real(kind=rk), parameter :: dfds(4) = (/ 2.0_rk, 5.0_rk, 10.0_rk, 30.0_rk /)
  real(kind=rk), parameter :: ncps_f(7) = (/ &
    0.0_rk, 1.0e-12_rk, 1.0e-10_rk, 1.0e-9_rk, 1.0_rk, 5.0_rk, 20.0_rk /)

  integer :: i, j, k, unit
  real(kind=rk) :: df, ncp, mean, sd, x_max, step, x, cum, ccum
  real(kind=rk) :: dfn, dfd, fx, f_max

  ! noncentral chi-squared CDF
  open(newunit=unit, file='tests/data/chi_squared_noncentral_cdf.csv', &
       status='replace', action='write')
  write(unit, '(a)') '# df, ncp, x, cdf, sf'
  do i = 1, size(dfs)
    do j = 1, size(ncps)
      df = dfs(i)
      ncp = ncps(j)
      mean = df + ncp
      if (df + 2.0_rk * ncp <= 0.0_rk) then
        sd = 1.0_rk
      else
        sd = 2.0_rk * (df + 2.0_rk * ncp)
      end if
      if (sd <= 0.0_rk) sd = 1.0_rk
      if (sd < 1.0_rk)  sd = 1.0_rk
      x_max = mean + 5.0_rk * sd
      step = x_max / 30.0_rk
      if (step <= 0.0_rk) step = 0.1_rk
      x = step / 4.0_rk
      do while (x <= x_max + 1.0e-12_rk)
        call cumchn(x, df, ncp, cum, ccum)
        call putval(unit, df, .false.)
        call putval(unit, ncp, .false.)
        call putval(unit, x, .false.)
        call putval(unit, cum, .false.)
        call putval(unit, ccum, .true.)
        x = x + step
      end do
    end do
  end do
  ! A large noncentrality, where the backward sum stops after ntired terms.
  call cumchn(1.0e6_rk, 10.0_rk, 1.0e6_rk, cum, ccum)
  call putval(unit, 10.0_rk, .false.)
  call putval(unit, 1.0e6_rk, .false.)
  call putval(unit, 1.0e6_rk, .false.)
  call putval(unit, cum, .false.)
  call putval(unit, ccum, .true.)
  close(unit)

  ! noncentral F CDF
  open(newunit=unit, file='tests/data/fisher_snedecor_noncentral_cdf.csv', &
       status='replace', action='write')
  write(unit, '(a)') '# dfn, dfd, ncp, f, cdf, sf'
  do i = 1, size(dfns)
    do j = 1, size(dfds)
      do k = 1, size(ncps_f)
        dfn = dfns(i)
        dfd = dfds(j)
        ncp = ncps_f(k)
        f_max = 10.0_rk
        step = 0.25_rk
        fx = step
        do while (fx <= f_max + 1.0e-12_rk)
          call cumfnc(fx, dfn, dfd, ncp, cum, ccum)
          call putval(unit, dfn, .false.)
          call putval(unit, dfd, .false.)
          call putval(unit, ncp, .false.)
          call putval(unit, fx, .false.)
          call putval(unit, cum, .false.)
          call putval(unit, ccum, .true.)
          fx = fx + step
        end do
      end do
    end do
  end do
  close(unit)

  write(0, '(a)') 'wrote 2 tables under tests/data/'

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
end program gen_noncentral_distributions
