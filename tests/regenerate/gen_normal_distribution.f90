! Reference tables for the Normal distribution via cdfnor.
!
! Build and run via tests/regenerate/regenerate.sh.

program gen_normal_distribution
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  external :: cdfnor

  real(kind=rk), dimension(3) :: means = (/ 0.0_rk, -1.5_rk, 3.25_rk /)
  real(kind=rk), dimension(3) :: sds   = (/ 1.0_rk,  0.5_rk, 2.75_rk /)
  integer :: i, unit, which, status
  real(kind=rk) :: x, mean, sd, p, q, bound
  real(kind=rk) :: xx, mm, ss

  ! normal_cdf
  open(newunit=unit, file='tests/data/normal_cdf.csv', status='replace', action='write')
  write(unit, '(a)') '# mean, sd, x, cdf, sf'
  do i = 1, 3
    mean = means(i)
    sd = sds(i)
    x = mean - 12.0_rk * sd
    do while (x <= mean + 12.0_rk * sd + 1.0e-12_rk)
      which = 1; status = 0; p = 0.0_rk; q = 0.0_rk; bound = 0.0_rk
      xx = x; mm = mean; ss = sd
      call cdfnor(which, p, q, xx, mm, ss, status, bound)
      if (status == 0) then
        call putval(unit, mean, .false.)
        call putval(unit, sd, .false.)
        call putval(unit, x, .false.)
        call putval(unit, p, .false.)
        call putval(unit, q, .true.)
      end if
      x = x + sd * 0.125_rk
    end do
  end do
  close(unit)

  ! normal_inverse_cdf: generate by x, store (mean, sd, p, q, x).
  open(newunit=unit, file='tests/data/normal_inverse_cdf.csv', status='replace', action='write')
  write(unit, '(a)') '# mean, sd, p, q, x'
  do i = 1, 3
    mean = means(i)
    sd = sds(i)
    x = mean - 7.0_rk * sd
    do while (x <= mean + 7.0_rk * sd + 1.0e-12_rk)
      which = 1; status = 0; p = 0.0_rk; q = 0.0_rk; bound = 0.0_rk
      xx = x; mm = mean; ss = sd
      call cdfnor(which, p, q, xx, mm, ss, status, bound)
      if (status == 0 .and. p > 0.0_rk .and. p < 1.0_rk &
          .and. q > 0.0_rk .and. q < 1.0_rk) then
        ! Make p and q exact complements of each other, so that the Rust
        ! inverse_cdf(p) and inverse_ccdf(q), which derive the other one,
        ! receive the same pair as cdfnor.
        q = 1.0_rk - p
        p = 1.0_rk - q
        which = 2; status = 0; bound = 0.0_rk
        xx = 0.0_rk; mm = mean; ss = sd
        call cdfnor(which, p, q, xx, mm, ss, status, bound)
        call putval(unit, mean, .false.)
        call putval(unit, sd, .false.)
        call putval(unit, p, .false.)
        call putval(unit, q, .false.)
        call putval(unit, xx, .true.)
      end if
      x = x + sd * 0.125_rk
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

end program gen_normal_distribution
