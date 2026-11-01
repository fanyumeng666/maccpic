! T10: cold Langmuir oscillation in the bounded domain.  Electrons and
!      immobile ions on the same lattice (quiet start, rho = 0), electrons
!      displaced by A sin(2 pi x/L); the mode amplitude sum_j E_j sin(2 pi x_j/L)
!      must oscillate at omega_pe (leapfrog: (2/dt) asin(omega_pe dt/2)) with
!      constant amplitude.  Uses the production modules global, deposit,
!      poisson and push; only the electrons are pushed.
program langmuir
  use global
  use random
  use poisson
  use deposit
  use push
  implicit none
  integer, parameter :: nsteps = 1500
  integer :: nfail, j, k, ncross
  integer(8) :: i, ntot
  real(8) :: x0, amp, wpe, wnum, a(nsteps), tcross(nsteps), tfirst, tlast, omega, s, xd, ep, ke0, w0, e, emax, emin
  real(8) :: a1max, a2max, drift
  nfail = 0
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  nx = 128; L = 0.067d0; n0 = 2.56d14; nppc = 64; Te0 = 0d0; Ti0 = 0d0
  freq = 13.56d6; steps_per_period = 400; V0 = 0d0; maxnp_factor = 1
  call derive_parameters(nprocs)
  call allocate_arrays()
  call poisson_init(nx, dx)
  call push_init()
  amp = 0.1d0*dx
  ntot = np_init_total
  do i = myrank + 1, ntot, nprocs
    x0 = (i - 0.5d0)*L/ntot
    call add_particle(2, x0, 0d0, 0d0, 0d0)
    call add_particle(1, x0 + amp*sin(twopi*x0/L), 0d0, 0d0, 0d0)
  end do
  call deposit_rho()
  call poisson_solve(rho, 0d0, 0d0, phi, ex)
  do i = 1, np(1)
    xd = x(i,1)/dx; j = min(int(xd), nx - 1); s = xd - j
    ep = ex(j)*(1d0 - s) + ex(j+1)*s
    vx(i,1) = vx(i,1) - 0.5d0*qm(1)*ep*dt
  end do
  wpe = sqrt(n0*qe*qe/(eps0*me))
  wnum = 2d0/dt*asin(0.5d0*wpe*dt)
  ke0 = kinetic_energy(1)
  w0 = field_energy()
  emax = -1d300; emin = 1d300
  do k = 1, nsteps
    call push_species(1)
    call deposit_rho()
    it = it + 1; time = it*dt
    call poisson_solve(rho, 0d0, 0d0, phi, ex)
    a(k) = 0d0
    do j = 0, nx
      a(k) = a(k) + ex(j)*sin(twopi*j*dx/L)
    end do
    e = kinetic_energy(1)*weight + field_energy()
    emax = max(emax, e); emin = min(emin, e)
  end do
  ! zero crossings of the mode amplitude (linear interpolation)
  ncross = 0
  do k = 2, nsteps
    if (a(k-1)*a(k) < 0d0) then
      ncross = ncross + 1
      tcross(ncross) = ((k-1) - a(k-1)/(a(k) - a(k-1)))*dt
    end if
  end do
  tfirst = tcross(1); tlast = tcross(ncross)
  omega = pi*(ncross - 1)/(tlast - tfirst)
  a1max = maxval(abs(a(1:200)))
  a2max = maxval(abs(a(nsteps-199:nsteps)))
  drift = (emax - emin)/max(ke0*weight + w0, 1d-300)
  if (myrank == 0) then
    write(*,'(a,i0,a)') '    T10 ', ncross, ' zero crossings'
    write(*,'(a,es14.6,a,es14.6,a,es14.6)') '    T10 omega measured = ', omega, '  omega_pe = ', wpe, '  leapfrog omega = ', wnum
    write(*,'(a,es11.3)') '    T10 omega/omega_pe - 1 = ', omega/wpe - 1d0
    write(*,'(a,es11.3)') '    T10 omega/omega_leapfrog - 1 = ', omega/wnum - 1d0
    write(*,'(a,es11.3)') '    T10 amplitude change (last/first 200 steps - 1) = ', a2max/a1max - 1d0
    write(*,'(a,es11.3)') '    T10 total energy (KE + W_field) spread / initial = ', drift
    write(*,'(a,i0,a,i0)') '    T10 particles absorbed: ', count_abs(), '  electrons left: ', np(1)
  end if
  call check('T10 frequency vs omega_pe', abs(omega/wpe - 1d0), 1d-2)
  call check('T10 frequency vs leapfrog omega', abs(omega/wnum - 1d0), 2d-3)
  call check('T10 amplitude growth', abs(a2max/a1max - 1d0), 5d-2)
  ! the cold plasma (lambda_D = 0) heats numerically on the grid, so the total
  ! energy is reported but not checked; the sheath test (T11) checks it warm
  if (myrank == 0) then
    if (nfail == 0) then
      write(*,'(a)') 'langmuir: PASS'
    else
      write(*,'(a,i0,a)') 'langmuir: FAIL (', nfail, ' checks)'
    end if
  end if
  call MPI_FINALIZE(ierr)
  if (nfail /= 0) stop 1
contains
  function field_energy() result(w)
    real(8) :: w, ef(0:nx-1)
    call efield_faces(phi, ef)
    w = 0.5d0*eps0*dx*sum(ef*ef)
  end function field_energy
  function count_abs() result(n)
    integer(8) :: n
    n = ntot/nprocs - np(1)
    if (myrank < mod(ntot, int(nprocs,8))) n = n + 1
  end function count_abs
  subroutine check(name, val, tol)
    character(len=*), intent(in) :: name
    real(8), intent(in) :: val, tol
    if (myrank /= 0) return
    if (val <= tol) then
      write(*,'(a,a,es11.3,a,es9.1,a)') '    PASS ', name, val, '  (tol ', tol, ')'
    else
      write(*,'(a,a,es11.3,a,es9.1,a)') '    FAIL ', name, val, '  (tol ', tol, ')'
      nfail = nfail + 1
    end if
  end subroutine check
end program langmuir
