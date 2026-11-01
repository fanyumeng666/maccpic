!=====================================================================
!  xsec: collision cross sections of e + He (elastic, two excitations,
!  ionization) and He+ + He (isotropic, backward) from the two-column
!  tables in xsec_dir (energy [eV], sigma [m^2]), put on uniform energy
!  grids and read back with linear interpolation in O(1):
!    electrons  de_e = 0.01  eV from e0_e = 0 eV      (same knots as the tables)
!    ions       de_i = 0.001 eV from e0_i = 1e-4 eV   (relative energy 1/2 mu g^2)
!  The grid origin is the smallest first energy of the species' tables, so
!  that the table knots (first energy + n*0.01 eV) are grid knots and the
!  piecewise-linear tables are reproduced exactly.  Beyond the last
!  tabulated energy the last value is used, below the first one the first
!  value (threshold processes: zero).  xsec_load
!  also finds the null-collision frequencies
!    numax_e = max_eps n_g sigma_tot(eps) v(eps),   v = sqrt(2 eps/m_e)
!    numax_i = max_eps n_g sigma_tot(eps) g(eps),   g = sqrt(2 eps/mu), mu = M/2
!=====================================================================
module xsec
  use global
  implicit none
  private
  public :: xsec_load, xsec_e, xsec_i, xsec_read_table, xsec_loaded, &
            nproc_e, nproc_i, de_e, de_i, eth_e, numax_e, numax_i, enumax_e, enumax_i, ngrid_e, ngrid_i
  integer, parameter :: nproc_e = 4, nproc_i = 2
  real(8), parameter :: de_e = 0.01d0, de_i = 0.001d0
  character(len=24), parameter :: fname_e(nproc_e) = (/ 'electron_scattering.dat ', 'excitation_1.dat        ', &
                                                        'excitation_2.dat        ', 'ionization.dat          ' /)
  character(len=24), parameter :: fname_i(nproc_i) = (/ 'ion_scattering.dat      ', 'ion_back_scatter.dat    ' /)
  integer :: ngrid_e = 0, ngrid_i = 0
  real(8) :: e0_e = 0d0, e0_i = 0d0                   ! grid origins [eV]
  real(8), allocatable :: sig_e(:,:), sig_i(:,:)      ! (process, 0:ngrid)
  type table
    real(8), allocatable :: e(:), s(:)
    integer :: n = 0
  end type table
  real(8) :: eth_e(nproc_e) = 0d0                     ! thresholds [eV] (first tabulated energy)
  real(8) :: numax_e = 0d0, numax_i = 0d0, enumax_e = 0d0, enumax_i = 0d0
  logical :: xsec_loaded = .false.
contains
  !
  subroutine xsec_load()
    integer :: k, j
    type(table) :: te(nproc_e), ti(nproc_i)
    real(8), allocatable :: tmp(:)
    real(8) :: v, nu, mu
    do k = 1, nproc_e
      call xsec_read_table(trim(xsec_dir)//'/'//trim(fname_e(k)), te(k)%e, te(k)%s, te(k)%n)
      eth_e(k) = 0d0
      if (k > 1) eth_e(k) = te(k)%e(1)
    end do
    do k = 1, nproc_i
      call xsec_read_table(trim(xsec_dir)//'/'//trim(fname_i(k)), ti(k)%e, ti(k)%s, ti(k)%n)
    end do
    e0_e = te(1)%e(1)
    do k = 2, nproc_e
      e0_e = min(e0_e, te(k)%e(1))
    end do
    e0_i = ti(1)%e(1)
    do k = 2, nproc_i
      e0_i = min(e0_i, ti(k)%e(1))
    end do
    ngrid_e = nint((numax_emax_e - e0_e)/de_e)
    ngrid_i = nint((numax_emax_i - e0_i)/de_i)
    if (allocated(sig_e)) deallocate(sig_e, sig_i)
    allocate(sig_e(nproc_e, 0:ngrid_e), sig_i(nproc_i, 0:ngrid_i), tmp(0:max(ngrid_e, ngrid_i)))
    do k = 1, nproc_e
      call to_grid(te(k)%e, te(k)%s, te(k)%n, e0_e, de_e, ngrid_e, tmp(0:ngrid_e), k > 1)
      sig_e(k,:) = tmp(0:ngrid_e)
    end do
    do k = 1, nproc_i
      call to_grid(ti(k)%e, ti(k)%s, ti(k)%n, e0_i, de_i, ngrid_i, tmp(0:ngrid_i), .false.)
      sig_i(k,:) = tmp(0:ngrid_i)
    end do
    deallocate(tmp)
    ! null-collision frequencies (maximum of the gridded piecewise-linear
    ! functions; the 1e-6 margin keeps the exceed counter quiet at the knots)
    numax_e = 0d0
    do j = 0, ngrid_e
      v = sqrt(2d0*(e0_e + j*de_e)*qe/me)
      nu = n_gas*sum(sig_e(:,j))*v
      if (nu > numax_e) then
        numax_e = nu; enumax_e = e0_e + j*de_e
      end if
    end do
    mu = 0.5d0*mion
    numax_i = 0d0
    do j = 0, ngrid_i
      v = sqrt(2d0*(e0_i + j*de_i)*qe/mu)
      nu = n_gas*(sig_i(1,j) + sig_i(2,j))*v
      if (nu > numax_i) then
        numax_i = nu; enumax_i = e0_i + j*de_i
      end if
    end do
    numax_e = numax_e*(1d0 + 1d-6)
    numax_i = numax_i*(1d0 + 1d-6)
    xsec_loaded = .true.
    if (myrank == 0) then
      write(*,'(a,a)') 'xsec: tables from ', trim(xsec_dir)
      write(*,'(a,3f8.3,a)') 'xsec: thresholds exc1, exc2, ion = ', eth_e(2:4), ' eV'
      write(*,'(a,es12.4,a,f7.2,a,es12.4,a,f8.1,a)') 'xsec: numax_e = ', numax_e, ' s^-1 at ', enumax_e, &
           ' eV;  numax_i = ', numax_i, ' s^-1 at ', enumax_i, ' eV'
      flush(6)
    end if
  end subroutine xsec_load
  !
  ! two-column table (energy, sigma), lines starting with # ignored
  subroutine xsec_read_table(fname, e, s, n)
    character(len=*), intent(in) :: fname
    real(8), allocatable, intent(out) :: e(:), s(:)
    integer, intent(out) :: n
    integer :: u, ios, i
    character(len=200) :: line
    open(newunit=u, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) call abort_run('xsec: cannot open '//trim(fname))
    n = 0
    do
      read(u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      n = n + 1
    end do
    if (n < 2) call abort_run('xsec: table too short: '//trim(fname))
    allocate(e(n), s(n))
    rewind(u)
    i = 0
    do
      read(u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      i = i + 1
      read(line, *, iostat=ios) e(i), s(i)
      if (ios /= 0) call abort_run('xsec: bad line in '//trim(fname))
      if (i > 1) then
        if (e(i) <= e(i-1)) call abort_run('xsec: energies not increasing in '//trim(fname))
      end if
    end do
    close(u)
  end subroutine xsec_read_table
  !
  ! linear interpolation of the table onto eps_j = e0 + j*de (constant beyond the ends)
  subroutine to_grid(e, s, n, e0, de, ng, sg, threshold)
    real(8), intent(in) :: e(:), s(:), e0, de
    integer, intent(in) :: n, ng
    real(8), intent(out) :: sg(0:ng)
    logical, intent(in) :: threshold
    integer :: j, m
    real(8) :: eps, f
    m = 1
    do j = 0, ng
      eps = e0 + j*de
      if (eps <= e(1)) then
        sg(j) = s(1)
        if (threshold) sg(j) = 0d0
      else if (eps >= e(n)) then
        sg(j) = s(n)
      else
        do while (e(m+1) < eps)
          m = m + 1
        end do
        f = (eps - e(m))/(e(m+1) - e(m))
        sg(j) = s(m)*(1d0 - f) + s(m+1)*f
      end if
    end do
  end subroutine to_grid
  !
  ! electron cross sections [m^2] at the kinetic energy eps [eV]
  subroutine xsec_e(eps, sig)
    real(8), intent(in) :: eps
    real(8), intent(out) :: sig(nproc_e)
    integer :: k
    real(8) :: f
    f = max(eps - e0_e, 0d0)/de_e
    k = int(f)
    if (k >= ngrid_e) then
      sig = sig_e(:, ngrid_e)
    else
      f = f - k
      sig = sig_e(:,k)*(1d0 - f) + sig_e(:,k+1)*f
    end if
  end subroutine xsec_e
  !
  ! ion cross sections [m^2] at the relative energy eps_rel [eV]
  subroutine xsec_i(eps, sig)
    real(8), intent(in) :: eps
    real(8), intent(out) :: sig(nproc_i)
    integer :: k
    real(8) :: f
    f = max(eps - e0_i, 0d0)/de_i
    k = int(f)
    if (k >= ngrid_i) then
      sig = sig_i(:, ngrid_i)
    else
      f = f - k
      sig = sig_i(:,k)*(1d0 - f) + sig_i(:,k+1)*f
    end if
  end subroutine xsec_i
  !
end module xsec
