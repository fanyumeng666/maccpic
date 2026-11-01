! T2: CIC deposition: single-particle weights and charge conservation,
!     a particle on a node, uniform random loading (interior density and
!     the half-cell electrode nodes), the species moments, and the sum
!     over ranks (run with -np 1 and -np 2: the 2-rank run splits the
!     same particles between the ranks and must give the same total).
program test_deposit
  use global
  use random
  use deposit
  implicit none
  integer :: nfail, j, isp
  integer(8) :: i, ntot
  real(8) :: w, tot, mean_in, sd, dev, rel
  real(8), allocatable :: dens(:), flux(:), ekd(:), ekxd(:)
  nfail = 0
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  nx = 64; L = 0.067d0; nppc = 400; n0 = 2.56d14; maxnp_factor = 2
  call derive_parameters(nprocs)
  call allocate_arrays()
  allocate(dens(0:nx), flux(0:nx), ekd(0:nx), ekxd(0:nx))
  !
  ! (a) one ion at x = 10.3 dx (on rank 0 only): weights 0.7 and 0.3, total charge q w
  if (myrank == 0) call add_particle(2, 10.3d0*dx, 0d0, 0d0, 0d0)
  call deposit_rho()
  w = qe*weight/dx
  dev = max(abs(rho(10) - 0.7d0*w), abs(rho(11) - 0.3d0*w))/w
  do j = 0, nx
    if (j /= 10 .and. j /= 11) dev = max(dev, abs(rho(j))/w)
  end do
  call report('T2a single particle weights (1-s, s)', dev, 1d-13)
  call report('T2a charge conservation sum(rho) dx - q w', abs(sum(rho)*dx - qe*weight)/(qe*weight), 1d-13)
  ! (b) particle exactly on node 20: everything on node 20
  np = 0
  if (myrank == 0) call add_particle(2, 20d0*dx, 0d0, 0d0, 0d0)
  call deposit_rho()
  call report('T2b particle on a node', abs(rho(20) - w)/w + abs(sum(rho)*dx - qe*weight)/(qe*weight), 1d-13)
  ! particle just below L lands in the last cell
  np = 0
  if (myrank == 0) call add_particle(2, L*(1d0 - 1d-15), 0d0, 0d0, 0d0)
  call deposit_rho()
  call report('T2b particle at L(1-1e-15) on node nx', abs(rho(nx) - w)/w, 1d-9)
  !
  ! (c) uniform random loading of both species, split over the ranks:
  !     rho ~ 0 by charge, n_i interior = n0 within 3 sigma, electrode nodes half
  np = 0
  call rng_init(777 + myrank)
  ntot = np_init_total/nprocs
  do isp = 1, nsp
    do i = 1, ntot
      call add_particle(isp, L*rand_u(), 1d3*isp, 0d0, 0d0)
    end do
  end do
  call deposit_number_density(2, dens)
  mean_in = sum(dens(1:nx-1))/(nx - 1)
  sd = n0/sqrt(dble(nppc))               ! 1 sigma of a node value (about nppc particles per cell)
  dev = 0d0
  do j = 1, nx - 1
    dev = max(dev, abs(dens(j) - n0)/sd)
  end do
  write(*,'(a,es12.4,a,es12.4)') '    T2c interior mean n_i = ', mean_in, '  n0 = ', n0
  call report('T2c interior mean density, deviation in sigma/sqrt(nx)', abs(mean_in - n0)/(sd/sqrt(dble(nx-1))), 3.5d0)
  call report('T2c max node deviation in sigma (nx-1 nodes)', dev, 4.5d0)
  call report('T2c electrode nodes = half density (sigma units)', &
       max(abs(dens(0) - 0.5d0*n0), abs(dens(nx) - 0.5d0*n0))/(sd/sqrt(2d0)), 4d0)
  ! moments: flux = n vx exactly (all vx equal), energy density = n (1/2 m vx^2)
  call deposit_species_moments(2, dens, flux, ekd, ekxd)
  rel = 0d0
  do j = 0, nx
    if (dens(j) > 0d0) then
      rel = max(rel, abs(flux(j) - 2d3*dens(j))/(2d3*dens(j)))
      rel = max(rel, abs(ekd(j) - 0.5d0*mion*4d6*dens(j))/(0.5d0*mion*4d6*dens(j)))
      rel = max(rel, abs(ekxd(j) - ekd(j))/ekd(j))
    end if
  end do
  call report('T2d moments flux = n vx, ekd = n m vx^2/2 (relative)', rel, 1d-13)
  ! total charge of both species: sum rho dx = weight (q_e N_e + q_i N_i) = 0 here
  call deposit_rho()
  tot = sum(rho)*dx
  call report('T2e total charge of a neutral loading / (q w N)', abs(tot)/(qe*weight*ntot*nprocs), 1d-12)
  !
  if (myrank == 0) then
    if (nfail == 0) then
      write(*,'(a,i0,a)') 'test_deposit (', nprocs, ' ranks): PASS'
    else
      write(*,'(a,i0,a,i0,a)') 'test_deposit (', nprocs, ' ranks): FAIL (', nfail, ' checks)'
    end if
  end if
  call MPI_FINALIZE(ierr)
  if (nfail /= 0) stop 1
contains
  subroutine report(name, val, tol)
    character(len=*), intent(in) :: name
    real(8), intent(in) :: val, tol
    if (myrank /= 0) return
    if (val <= tol) then
      write(*,'(a,a,es11.3,a,es9.1,a)') '    PASS ', name, val, '  (tol ', tol, ')'
    else
      write(*,'(a,a,es11.3,a,es9.1,a)') '    FAIL ', name, val, '  (tol ', tol, ')'
      nfail = nfail + 1
    end if
  end subroutine report
end program test_deposit
