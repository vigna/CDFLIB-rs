! Reference tables for the erf and standard-normal kernels:
! error_f, error_fc, error_fc_scaled, cumnor, dinvnr.
!
! Build and run via tests/regenerate/regenerate.sh.

program gen_erf_normal_kernels
  implicit none
  integer, parameter :: rk = kind(1.0d0)

  real(kind=rk), external :: error_f, error_fc, dinvnr
  external :: cumnor

  call gen_error_f()
  call gen_error_fc()
  call gen_error_fc_scaled()
  call gen_cumnor()
  call gen_dinvnr()

  write(0, '(a)') 'wrote 5 tables under tests/data/'

contains

  ! Write a single double-precision value to `unit` followed by a comma,
  ! or by a newline if `last` is true. Uses 18 significant digits so each
  ! distinct f64 round-trips through decimal.
  subroutine putval(unit, v, last)
    integer, intent(in) :: unit
    real(kind=rk), intent(in) :: v
    logical, intent(in) :: last
    character(len=32) :: buf
    ! es25.17e3 = sign + digit + '.' + 17 frac + 'E' + signed-3 = 25 chars
    write(buf, '(es25.17e3)') v
    if (last) then
      write(unit, '(a)') trim(adjustl(buf))
    else
      write(unit, '(a, a)', advance='no') trim(adjustl(buf)), ','
    end if
  end subroutine putval

  subroutine gen_error_f()
    integer :: unit
    real(kind=rk) :: x
    open(newunit=unit, file='tests/data/error_f.csv', status='replace', action='write')
    write(unit, '(a)') '# x, erf(x)'
    x = -6.0_rk
    do while (x <= 6.0_rk + 1.0e-12_rk)
      call putval(unit, x, .false.)
      call putval(unit, error_f(x), .true.)
      x = x + 0.0625_rk
    end do
    close(unit)
  end subroutine gen_error_f

  subroutine gen_error_fc()
    integer :: unit, ind0
    real(kind=rk) :: x
    ind0 = 0
    open(newunit=unit, file='tests/data/error_fc.csv', status='replace', action='write')
    write(unit, '(a)') '# x, erfc(x)'
    x = -6.0_rk
    do while (x <= 30.0_rk + 1.0e-12_rk)
      call putval(unit, x, .false.)
      call putval(unit, error_fc(ind0, x), .true.)
      x = x + 0.0625_rk
    end do
    close(unit)
  end subroutine gen_error_fc

  subroutine gen_error_fc_scaled()
    integer :: unit, ind1
    real(kind=rk) :: x
    ind1 = 1
    open(newunit=unit, file='tests/data/error_fc_scaled.csv', status='replace', action='write')
    write(unit, '(a)') '# x, erfc(x)*exp(x^2)'
    x = -6.0_rk
    do while (x <= 60.0_rk + 1.0e-12_rk)
      call putval(unit, x, .false.)
      call putval(unit, error_fc(ind1, x), .true.)
      x = x + 0.0625_rk
    end do
    close(unit)
  end subroutine gen_error_fc_scaled

  subroutine gen_cumnor()
    integer :: unit
    real(kind=rk) :: x, cum, ccum
    open(newunit=unit, file='tests/data/cumnor.csv', status='replace', action='write')
    write(unit, '(a)') '# x, cum, ccum'
    x = -38.0_rk
    do while (x <= 38.0_rk + 1.0e-12_rk)
      call cumnor(x, cum, ccum)
      call putval(unit, x, .false.)
      call putval(unit, cum, .false.)
      call putval(unit, ccum, .true.)
      x = x + 0.0625_rk
    end do
    close(unit)
  end subroutine gen_cumnor

  subroutine gen_dinvnr()
    ! Parametrize by x rather than p; store (cum, ccum, back). Skips rows
    ! where cum/ccum has saturated to 0 or 1.
    integer :: unit
    real(kind=rk) :: x, cum, ccum, back
    open(newunit=unit, file='tests/data/dinvnr.csv', status='replace', action='write')
    write(unit, '(a)') '# p, q, x'
    x = -7.0_rk
    do while (x <= 7.0_rk + 1.0e-12_rk)
      call cumnor(x, cum, ccum)
      if (cum > 0.0_rk .and. cum < 1.0_rk .and. ccum > 0.0_rk .and. ccum < 1.0_rk) then
        back = dinvnr(cum, ccum)
        call putval(unit, cum, .false.)
        call putval(unit, ccum, .false.)
        call putval(unit, back, .true.)
      end if
      x = x + 0.0625_rk
    end do
    ! The endpoints, where the Newton steps never converge and the result
    ! is the NaN start value (cdflib.f90:8128, :8130).
    call putval(unit, 0.0_rk, .false.)
    call putval(unit, 1.0_rk, .false.)
    call putval(unit, dinvnr(0.0_rk, 1.0_rk), .true.)
    call putval(unit, 1.0_rk, .false.)
    call putval(unit, 0.0_rk, .false.)
    call putval(unit, dinvnr(1.0_rk, 0.0_rk), .true.)
    close(unit)
  end subroutine gen_dinvnr

end program gen_erf_normal_kernels
