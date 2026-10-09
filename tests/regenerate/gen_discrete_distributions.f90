! Binomial, Poisson, Negative Binomial CDF reference tables.

program gen_discrete_distributions
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  external :: cumbin, cumpoi, cumnbn

  integer, parameter :: ns(4) = (/ 5, 10, 50, 200 /)
  real(kind=rk), parameter :: prs(5) = (/ 0.05_rk, 0.25_rk, 0.5_rk, 0.75_rk, 0.95_rk /)
  real(kind=rk), parameter :: lambdas(6) = (/ 0.5_rk, 1.0_rk, 5.0_rk, 12.0_rk, 50.0_rk, 200.0_rk /)
  integer, parameter :: rs(4) = (/ 1, 3, 10, 50 /)
  real(kind=rk), parameter :: prs_nbn(4) = (/ 0.1_rk, 0.25_rk, 0.5_rk, 0.9_rk /)

  integer :: i, j, unit, s_int, s_max
  real(kind=rk) :: s, n, pr, ompr, cum, ccum, lam, r_d

  ! Binomial CDF
  open(newunit=unit, file='tests/data/binomial_cdf.csv', status='replace', action='write')
  write(unit, '(a)') '# n, pr, s, cdf, sf'
  do i = 1, size(ns)
    n = real(ns(i), rk)
    do j = 1, size(prs)
      pr = prs(j)
      ompr = 1.0_rk - pr
      do s_int = 0, ns(i)
        s = real(s_int, rk)
        call cumbin(s, n, pr, ompr, cum, ccum)
        call putint(unit, ns(i), .false.)
        call putval(unit, pr, .false.)
        call putint(unit, s_int, .false.)
        call putval(unit, cum, .false.)
        call putval(unit, ccum, .true.)
      end do
    end do
  end do
  close(unit)

  ! Poisson CDF
  open(newunit=unit, file='tests/data/poisson_cdf.csv', status='replace', action='write')
  write(unit, '(a)') '# lambda, s, cdf, sf'
  do i = 1, size(lambdas)
    lam = lambdas(i)
    if (lam < 1.0_rk) then
      s_max = int(lam + 10.0_rk * 1.0_rk + 5)
    else
      s_max = int(lam + 10.0_rk * lam + 5)
    end if
    do s_int = 0, s_max
      s = real(s_int, rk)
      call cumpoi(s, lam, cum, ccum)
      call putval(unit, lam, .false.)
      call putint(unit, s_int, .false.)
      call putval(unit, cum, .false.)
      call putval(unit, ccum, .true.)
    end do
  end do
  close(unit)

  ! Negative Binomial CDF
  open(newunit=unit, file='tests/data/negative_binomial_cdf.csv', status='replace', action='write')
  write(unit, '(a)') '# n, pr, s, cdf, sf'
  do i = 1, size(rs)
    r_d = real(rs(i), rk)
    do j = 1, size(prs_nbn)
      pr = prs_nbn(j)
      ompr = 1.0_rk - pr
      s_max = int( real(rs(i), rk) * (1.0_rk - pr) / pr &
        + 10.0_rk * (real(rs(i), rk) * (1.0_rk - pr) / (pr * pr)) + 5.0_rk )
      if (s_max < 30) s_max = 30
      if (s_max > 1000) s_max = 1000
      do s_int = 0, s_max
        s = real(s_int, rk)
        call cumnbn(s, r_d, pr, ompr, cum, ccum)
        call putint(unit, rs(i), .false.)
        call putval(unit, pr, .false.)
        call putint(unit, s_int, .false.)
        call putval(unit, cum, .false.)
        call putval(unit, ccum, .true.)
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
end program gen_discrete_distributions
