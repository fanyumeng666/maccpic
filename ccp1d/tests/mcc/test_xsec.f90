! T4: cross-section lookup against direct interpolation of the original
!     tables, behaviour beyond the table ends and below thresholds, the
!     thresholds themselves and the null-collision frequencies.
program test_xsec
  use global
  use xsec
  implicit none
  integer :: nfail, k, i, n
  real(8), allocatable :: e(:), s(:)
  real(8) :: eps, ref, sig(nproc_e), sgi(nproc_i), maxrel, smax, u
  character(len=24) :: fe(4) = (/ 'electron_scattering.dat ', 'excitation_1.dat        ', 'excitation_2.dat        ', 'ionization.dat          ' /)
  character(len=24) :: fi(2) = (/ 'ion_scattering.dat      ', 'ion_back_scatter.dat    ' /)
  nfail = 0
  nx = 128; L = 0.067d0; n_gas = 9.64d20; xsec_dir = 'data/he'
  numax_emax_e = 1000d0; numax_emax_i = 1000d0
  call derive_parameters(1)
  call xsec_load()
  call random_seed()
  ! (a) electrons: grid knots coincide with the table knots -> agreement to round-off
  do k = 1, nproc_e
    call xsec_read_table('data/he/'//trim(fe(k)), e, s, n)
    smax = maxval(s)
    maxrel = 0d0
    do i = 1, 20000
      call random_number(u)
      eps = u*800d0
      ref = interp(e, s, n, eps, k > 1)
      call xsec_e(eps, sig)
      maxrel = max(maxrel, abs(sig(k) - ref)/smax)
    end do
    call report('T4a electron table '//trim(fe(k))//' lookup vs table (rel. to max sigma)', maxrel, 1d-12)
  end do
  ! (b) ions: grid 0.001 eV vs table knots at 1e-4 + n*0.01 -> small interpolation differences
  do k = 1, nproc_i
    call xsec_read_table('data/he/'//trim(fi(k)), e, s, n)
    maxrel = 0d0
    do i = 1, 20000
      call random_number(u)
      eps = 10d0**(-3d0 + 6d0*u)       ! 1e-3 .. 1e3 eV, log-uniform
      ref = interp(e, s, n, eps, .false.)
      call xsec_i(eps, sgi)
      maxrel = max(maxrel, abs(sgi(k) - ref)/ref)
    end do
    call report('T4b ion table '//trim(fi(k))//' lookup vs table (relative)', maxrel, 1d-6)
  end do
  ! (c) beyond the table: last value; below thresholds: zero; below the first ion energy: first value
  call xsec_read_table('data/he/ionization.dat', e, s, n)
  call xsec_e(900d0, sig)
  call report('T4c ionization sigma at 900 eV = last tabulated value', abs(sig(4) - s(n))/s(n), 1d-12)
  call xsec_e(19d0, sig)
  call report('T4c sigma(19 eV) of exc1, exc2, ion = 0', abs(sig(2)) + abs(sig(3)) + abs(sig(4)), 0d0)
  call xsec_e(24.5d0, sig)
  call report('T4c sigma_ion(24.5 eV) = 0', abs(sig(4)), 0d0)
  call xsec_read_table('data/he/ion_scattering.dat', e, s, n)
  call xsec_i(1d-5, sgi)
  call report('T4c sigma_iso(1e-5 eV) = first tabulated value', abs(sgi(1) - s(1))/s(1), 1d-12)
  ! (d) thresholds and nu_max (reference values from the plan)
  call report('T4d thresholds 19.82, 20.61, 24.59 eV', abs(eth_e(2) - 19.82d0) + abs(eth_e(3) - 20.61d0) + abs(eth_e(4) - 24.59d0), 1d-9)
  write(*,'(a,es13.5,a,f7.2,a,es13.5,a,f8.1,a)') '    T4d numax_e = ', numax_e, ' s^-1 at ', enumax_e, ' eV;  numax_i = ', numax_i, ' at ', enumax_i, ' eV'
  call report('T4d numax_e vs 8.654e7 s^-1 (relative)', abs(numax_e/8.654d7 - 1d0), 5d-3)
  call report('T4d numax_e energy vs 8.2 eV (eV)', abs(enumax_e - 8.2d0), 0.5d0)
  call report('T4d numax_i vs 2.793e7 s^-1 (relative)', abs(numax_i/2.793d7 - 1d0), 5d-3)
  ! (e) nu_tot never exceeds numax on a fine scan (including between the knots)
  maxrel = 0d0
  do i = 0, 2000000
    eps = i*5d-4
    call xsec_e(eps, sig)
    maxrel = max(maxrel, n_gas*sum(sig)*sqrt(2d0*eps*qe/me)/numax_e)
  end do
  call report('T4e max nu_tot_e/numax_e over a 5e-4 eV scan to 1000 eV (must be <= 1)', maxrel - 1d0, 0d0)
  maxrel = 0d0
  do i = 0, 1000000
    eps = i*1d-3
    call xsec_i(eps, sgi)
    maxrel = max(maxrel, n_gas*sum(sgi)*sqrt(2d0*eps*qe/(0.5d0*mion))/numax_i)
  end do
  call report('T4e max nu_tot_i/numax_i over a 1e-3 eV scan (must be <= 1)', maxrel - 1d0, 0d0)
  if (nfail == 0) then
    write(*,'(a)') 'test_xsec: PASS'
  else
    write(*,'(a,i0,a)') 'test_xsec: FAIL (', nfail, ' checks)'
    stop 1
  end if
contains
  function interp(e, s, n, eps, threshold) result(r)
    real(8), intent(in) :: e(:), s(:), eps
    integer, intent(in) :: n
    logical, intent(in) :: threshold
    real(8) :: r, f
    integer :: lo, hi, mid
    if (eps <= e(1)) then
      r = s(1)
      if (threshold) r = 0d0
    else if (eps >= e(n)) then
      r = s(n)
    else
      lo = 1; hi = n
      do while (hi - lo > 1)
        mid = (lo + hi)/2
        if (e(mid) <= eps) then
          lo = mid
        else
          hi = mid
        end if
      end do
      f = (eps - e(lo))/(e(hi) - e(lo))
      r = s(lo)*(1d0 - f) + s(hi)*f
    end if
  end function interp
  subroutine report(name, val, tol)
    character(len=*), intent(in) :: name
    real(8), intent(in) :: val, tol
    if (val <= tol) then
      write(*,'(a,a,es11.3,a,es9.1,a)') '    PASS ', name, val, '  (tol ', tol, ')'
    else
      write(*,'(a,a,es11.3,a,es9.1,a)') '    FAIL ', name, val, '  (tol ', tol, ')'
      nfail = nfail + 1
    end if
  end subroutine report
end program test_xsec
