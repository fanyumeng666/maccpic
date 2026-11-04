! T12: electron swarm in a uniform field.  Electrons only, periodic in x,
!      no field solve, the production modules push the velocities and
!      mcc_electrons collides them with the helium background.  Pulsed
!      Townsend model: the electrons created by ionization are kept and the
!      population is renormalised after every step by deleting randomly
!      chosen electrons (the ions are discarded), so the normalised energy
!      distribution reaches the steady state of the temporal-growth
!      Boltzmann equation solved by tools/boltzmann2term.py with the same
!      cross sections.  For each E/N the mean energy, drift velocity <vx>,
!      effective ionization frequency and alpha/N, collision frequencies,
!      the energy distribution and two energy-balance checks are written to
!      <outdir>/swarm_pic.dat and swarm_eedf_<E/N>Td.dat.
!      Usage: swarm [outdir] [xsec_dir]
program swarm
  use global
  use random
  use xsec
  use mcc
  implicit none
  integer, parameter :: nen = 5
  real(8), parameter :: en_td(nen) = (/ 20d0, 50d0, 100d0, 150d0, 200d0 /)   ! above ~150 Td electrons run away in this cross-section model
  integer(8), parameter :: nel = 20000
  integer, parameter :: nsteps = 150000, navg = 100000, nbin = 600, hist_every = 20
  real(8), parameter :: debin = 0.5d0
  character(len=256) :: outdir, arg
  character(len=300) :: fn
  integer :: k, istep, u, uf, b
  integer(8) :: i, c_el, c_x1, c_x2, c_iz, nhist
  real(8) :: efield, f, sig, vold, vnew, work, sume, sumvx, eloss, einel, erem, ke0, ke1, tavg, e, vd
  integer(8) :: r, nrem
  real(8) :: nu_el, nu_x1, nu_x2, nu_iz, pin, ploss, ident, alpha_n, fb, hist(nbin), emean
  real(8) :: t0, t1
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  outdir = '.'
  if (command_argument_count() >= 1) call get_command_argument(1, outdir)
  xsec_dir = 'data/he'
  if (command_argument_count() >= 2) then
    call get_command_argument(2, arg); xsec_dir = arg
  end if
  nx = 100; L = 1d0; n_gas = 9.64d20; T_gas = 300d0
  nppc = 400; maxnp_factor = 2; freq = 13.56d6; steps_per_period = 400     ! dt of Case 1
  numax_emax_e = 1000d0; numax_emax_i = 1000d0
  call derive_parameters(1)
  call allocate_arrays()
  call rng_init(777)
  call xsec_load()
  call mcc_init()
  fb = (1d0 - exp(-numax_e*dt))/(numax_e*dt)
  write(fn,'(a,a)') trim(outdir), '/swarm_pic.dat'
  open(newunit=u, file=fn, status='replace', action='write')
  write(u,'(a)') '# ccp1d electron swarm (tests/swarm/swarm.f90): uniform field, periodic, He background, pulsed-Townsend (secondaries kept, random renormalisation)'
  write(u,'(a,es12.4,a,es12.4,a,i0,a,i0,a,i0)') '# n_gas = ', n_gas, ' m^-3, dt = ', dt, ' s, electrons = ', nel, ', steps = ', nsteps, ', averaged over the last ', navg
  write(u,'(a,f8.5,a)') '# one-collision-per-step factor (1-exp(-numax dt))/(numax dt) = ', fb, ' (all collision rates are lower by this factor)'
  write(u,'(a)') '# 1 E/N[Td]  2 E[V/m]  3 <eps>[eV]  4 v_d[m/s]  5 alpha/N = nu_iz/(N v_d)[m^2]  6 P_in=eEv_d[eV/s]'// &
       '  7 P_loss = inelastic + elastic + removed electrons [eV/s]  8 (P_in-P_loss)/P_in'// &
       '  9 nu_el[s^-1] 10 nu_exc1 11 nu_exc2 12 nu_iz (= net growth rate)  13 exact identity residual  14 wall[s]'
  do k = 1, nen
    efield = -en_td(k)*1d-21*n_gas            ! electrons drift towards +x
    f = qm(1)*efield*dt
    np = 0
    sig = sqrt(2d0*qe/me)                      ! 2 eV start
    do i = 1, nel
      call add_particle(1, L*rand_u(), sig*rand_gauss(), sig*rand_gauss(), sig*rand_gauss())
    end do
    c_el = 0; c_x1 = 0; c_x2 = 0; c_iz = 0; nhist = 0; nrem = 0
    work = 0d0; sume = 0d0; sumvx = 0d0; eloss = 0d0; einel = 0d0; erem = 0d0; hist = 0d0; ke0 = 0d0
    t0 = MPI_WTIME()
    do istep = 1, nsteps
      if (istep == nsteps - navg + 1) ke0 = kinetic(nel)
      do i = 1, nel
        vold = vx(i,1)
        vnew = vold + f
        vx(i,1) = vnew
        x(i,1) = modulo(x(i,1) + vnew*dt, L)
        if (istep > nsteps - navg) work = work + 0.5d0*me*(vnew*vnew - vold*vold)
      end do
      call mcc_step()
      ! pulsed Townsend: keep the new electrons, drop the ions, and delete randomly
      ! chosen electrons until the population is back to nel (their energy is booked)
      np(2) = 0
      e = 0d0
      do while (np(1) > nel)
        r = 1 + int(rand_u()*np(1))
        if (r > np(1)) r = np(1)
        e = e + 0.5d0*me*(vx(r,1)**2 + vy(r,1)**2 + vz(r,1)**2)
        call remove_particle(1, r)
        nrem = nrem + 1
      end do
      if (istep > nsteps - navg) then
        erem = erem + e
        eloss = eloss + eloss_el_step; einel = einel + einel_step
        c_el = c_el + n_el_step; c_x1 = c_x1 + n_exc_step(1); c_x2 = c_x2 + n_exc_step(2); c_iz = c_iz + n_ion_step
        do i = 1, nel
          sume = sume + 0.5d0*me*(vx(i,1)**2 + vy(i,1)**2 + vz(i,1)**2)
          sumvx = sumvx + vx(i,1)
        end do
        if (mod(istep, hist_every) == 0) then
          do i = 1, nel
            b = int(0.5d0*me*(vx(i,1)**2 + vy(i,1)**2 + vz(i,1)**2)/qe/debin) + 1
            if (b <= nbin) hist(b) = hist(b) + 1d0
          end do
          nhist = nhist + nel
        end if
      end if
    end do
    ke1 = kinetic(nel)
    t1 = MPI_WTIME()
    tavg = navg*dt
    emean = sume/(dble(nel)*navg)/qe
    vd = sumvx/(dble(nel)*navg)
    nu_el = c_el/(nel*tavg); nu_x1 = c_x1/(nel*tavg); nu_x2 = c_x2/(nel*tavg); nu_iz = c_iz/(nel*tavg)
    pin = work/(nel*tavg)/qe
    ploss = (einel + eloss + erem)/(nel*tavg)/qe
    ident = (work - (ke1 - ke0) - (einel + eloss + erem))/work
    alpha_n = nu_iz/(n_gas*vd)
    write(*,'(a,f6.0,a,f8.3,a,es11.4,a,es10.3,a,es10.3,a,es10.3,a,f7.4,a,es9.1,a,f6.1,a)') 'E/N = ', en_td(k), ' Td: <eps> = ', emean, &
         ' eV, v_d = ', vd, ' m/s, alpha/N = ', alpha_n, ' m^2, P_in = ', pin, ' eV/s, P_loss = ', ploss, &
         ' eV/s, (P_in-P_loss)/P_in = ', (pin - ploss)/pin, ', identity ', ident, '  (', t1 - t0, ' s)'
    write(u,'(f8.1,es14.6,f12.5,es14.6,es14.6,es14.6,es14.6,es12.3,4es13.5,es12.3,f10.1)') en_td(k), efield, emean, vd, alpha_n, &
         pin, ploss, (pin - ploss)/pin, nu_el, nu_x1, nu_x2, nu_iz, ident, t1 - t0
    flush(u)
    write(fn,'(a,a,i0,a)') trim(outdir), '/swarm_eedf_', nint(en_td(k)), 'Td.dat'
    open(newunit=uf, file=fn, status='replace', action='write')
    write(uf,'(a,f6.0,a)') '# electron energy distribution of the swarm at E/N = ', en_td(k), ' Td (pulsed Townsend), normalised to integral n deps = 1'
    write(uf,'(a)') '# 1 energy_bin_centre[eV]  2 n[eV^-1]'
    do b = 1, nbin
      write(uf,'(2es15.7)') (b - 0.5d0)*debin, hist(b)/(dble(nhist)*debin)
    end do
    close(uf)
  end do
  close(u)
  call MPI_FINALIZE(ierr)
contains
  function kinetic(n) result(ke)
    integer(8), intent(in) :: n
    real(8) :: ke
    integer(8) :: j
    ke = 0d0
    do j = 1, n
      ke = ke + vx(j,1)**2 + vy(j,1)**2 + vz(j,1)**2
    end do
    ke = 0.5d0*me*ke
  end function kinetic
end program swarm
