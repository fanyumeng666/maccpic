! T5-T7: the Monte Carlo collision module without fields or walls.
!   T5 collision rates: monoenergetic electrons (10, 30, 100 eV) and ions
!      (10 eV lab, cold gas): event counts per process versus
!      n_g sigma_k v * (1 - exp(-numax dt))/(numax dt); the P test count.
!   T6 kinematics: elastic energy loss formula, excitation and ionization
!      energies, isotropy of the new directions (chi^2), thermal velocities
!      of the created ions, centre-of-mass invariants of ion collisions.
!   T7 relaxation: 10 eV ions thermalize to the gas temperature; electrons
!      below the thresholds cool at the elastic rate 2 m/M nu_el eps.
program test_mcc
  use global
  use random
  use xsec
  use mcc
  implicit none
  integer :: nfail
  nfail = 0
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  nx = 32; L = 0.067d0; n_gas = 9.64d20; T_gas = 300d0; xsec_dir = 'data/he'
  numax_emax_e = 1000d0; numax_emax_i = 1000d0
  nppc = 8000; maxnp_factor = 2; freq = 13.56d6; steps_per_period = 400    ! Case 1 dt, maxnp = 512000
  call derive_parameters(1)
  call allocate_arrays()
  call rng_init(2024)
  call xsec_load()
  !
  call rates_electrons(10d0, 100000_8, 600)
  call rates_electrons(30d0, 100000_8, 600)
  call rates_electrons(100d0, 100000_8, 600)
  call rates_ions(10d0, 100000_8, 1000)
  call kinematics_electrons()
  call kinematics_ions()
  call isotropy()
  call relax_ions()
  call relax_electrons()
  !
  if (nfail == 0) then
    write(*,'(a)') 'test_mcc: PASS'
  else
    write(*,'(a,i0,a)') 'test_mcc: FAIL (', nfail, ' checks)'
  end if
  call MPI_FINALIZE(ierr)
  if (nfail /= 0) stop 1
contains
  !
  subroutine load_mono(isp, eps0, n)
    integer, intent(in) :: isp
    real(8), intent(in) :: eps0
    integer(8), intent(in) :: n
    integer(8) :: i
    real(8) :: v0, ux, uy, uz
    np = 0
    v0 = sqrt(2d0*eps0*qe/m(isp))
    do i = 1, n
      call isotropic_direction(ux, uy, uz)
      call add_particle(isp, L*rand_u(), v0*ux, v0*uy, v0*uz)
    end do
  end subroutine load_mono
  !
  ! rescale all speeds of species isp to eps0 (random direction if at rest), drop products
  subroutine reset_mono(isp, eps0, n)
    integer, intent(in) :: isp
    real(8), intent(in) :: eps0
    integer(8), intent(in) :: n
    integer(8) :: i
    real(8) :: v0, v, ux, uy, uz
    v0 = sqrt(2d0*eps0*qe/m(isp))
    np(isp) = n
    np(3 - isp) = 0
    do i = 1, n
      v = sqrt(vx(i,isp)**2 + vy(i,isp)**2 + vz(i,isp)**2)
      if (v > 0d0) then
        vx(i,isp) = vx(i,isp)*v0/v; vy(i,isp) = vy(i,isp)*v0/v; vz(i,isp) = vz(i,isp)*v0/v
      else
        call isotropic_direction(ux, uy, uz)
        vx(i,isp) = v0*ux; vy(i,isp) = v0*uy; vz(i,isp) = v0*uz
      end if
    end do
  end subroutine reset_mono
  !
  subroutine rates_electrons(eps0, n, nsteps)
    real(8), intent(in) :: eps0
    integer(8), intent(in) :: n
    integer, intent(in) :: nsteps
    integer :: it, k
    integer(8) :: cnt(4), ntest
    real(8) :: sig(nproc_e), nu(nproc_e), v0, fb, expc, z, zmax
    character(len=10) :: pname(4) = (/ 'elastic   ', 'exc 19.82 ', 'exc 20.61 ', 'ionization' /)
    T_gas = 300d0
    call mcc_init()
    call load_mono(1, eps0, n)
    cnt = 0; ntest = 0
    do it = 1, nsteps
      call mcc_electrons_only()
      cnt(1) = cnt(1) + n_el_step; cnt(2) = cnt(2) + n_exc_step(1); cnt(3) = cnt(3) + n_exc_step(2); cnt(4) = cnt(4) + n_ion_step
      ntest = ntest + n_test_step(1)
      call reset_mono(1, eps0, n)
    end do
    v0 = sqrt(2d0*eps0*qe/me)
    call xsec_e(eps0, sig)
    nu = n_gas*sig*v0
    fb = (1d0 - exp(-numax_e*dt))/(numax_e*dt)
    write(*,'(a,f6.1,a,es10.3,a,f7.4)') '    T5 electrons at ', eps0, ' eV: numax dt = ', numax_e*dt, '  one-collision-per-step factor = ', fb
    expc = dble(n)*nsteps*p_coll_e
    z = (ntest - expc)/sqrt(expc)
    call report('T5 electrons: number of P tests passed, z-score', abs(z), 4d0)
    zmax = 0d0
    do k = 1, 4
      expc = dble(n)*nsteps*dt*nu(k)*fb
      if (expc > 20d0) then
        z = (cnt(k) - expc)/sqrt(expc)
        write(*,'(a,a,a,i9,a,f12.1,a,f7.2,a,f10.4,a)') '       ', pname(k), ': counted ', cnt(k), '  expected ', expc, '  z = ', z, &
             '  (nu_k/numax = ', nu(k)/numax_e, ')'
        zmax = max(zmax, abs(z))
      else
        write(*,'(a,a,a,i9,a,f12.3)') '       ', pname(k), ': counted ', cnt(k), '  expected ', expc
        if (cnt(k) > 0 .and. expc == 0d0) zmax = 100d0
      end if
    end do
    call report('T5 electrons: max |z| of the process rates', zmax, 4d0)
    call report('T5 electrons: sum nu_k <= numax (nu_tot/numax)', sum(nu)/numax_e - 1d0, 0d0)
  end subroutine rates_electrons
  !
  subroutine mcc_electrons_only()
    ! mcc_step resets the counters and runs both species; here the ion population is empty
    call mcc_step()
  end subroutine mcc_electrons_only
  !
  subroutine rates_ions(eps0, n, nsteps)
    real(8), intent(in) :: eps0
    integer(8), intent(in) :: n
    integer, intent(in) :: nsteps
    integer :: it, k
    integer(8) :: cnt(2), ntest
    real(8) :: sig(nproc_i), nu(2), v0, fb, expc, z, zmax, eps_rel
    character(len=10) :: pname(2) = (/ 'isotropic ', 'backward  ' /)
    T_gas = 0d0                        ! cold gas: g = v_i exactly
    call mcc_init()
    call load_mono(2, eps0, n)
    cnt = 0; ntest = 0
    do it = 1, nsteps
      call mcc_step()
      cnt = cnt + n_icoll_step
      ntest = ntest + n_test_step(2)
      call reset_mono(2, eps0, n)
    end do
    v0 = sqrt(2d0*eps0*qe/mion)
    eps_rel = 0.5d0*(0.5d0*mion)*v0*v0/qe
    call xsec_i(eps_rel, sig)
    nu = n_gas*sig*v0
    fb = (1d0 - exp(-numax_i*dt))/(numax_i*dt)
    write(*,'(a,f6.1,a,f6.2,a,es10.3,a,f7.4)') '    T5 ions at ', eps0, ' eV lab (eps_rel = ', eps_rel, ' eV): numax dt = ', numax_i*dt, '  factor = ', fb
    expc = dble(n)*nsteps*p_coll_i
    call report('T5 ions: number of P tests passed, z-score', abs(ntest - expc)/sqrt(expc), 4d0)
    zmax = 0d0
    do k = 1, 2
      expc = dble(n)*nsteps*dt*nu(k)*fb
      z = (cnt(k) - expc)/sqrt(expc)
      write(*,'(a,a,a,i9,a,f12.1,a,f7.2)') '       ', pname(k), ': counted ', cnt(k), '  expected ', expc, '  z = ', z
      zmax = max(zmax, abs(z))
    end do
    call report('T5 ions: max |z| of the process rates', zmax, 4d0)
    T_gas = 300d0
  end subroutine rates_ions
  !
  subroutine kinematics_electrons()
    integer(8), parameter :: n = 200000
    integer(8) :: i, nchanged, nbad, nel, nx1, nx2, nio
    real(8), allocatable :: v0x(:), v0y(:), v0z(:)
    real(8) :: eps0, e0, e1, cchi, pred, err, maxerr, vn2, v1
    real(8), parameter :: floss = 2d0*me/mion
    T_gas = 300d0
    call mcc_init()
    allocate(v0x(n), v0y(n), v0z(n))
    ! (a) 10 eV: only elastic collisions possible, check eps' = eps [1 - 2m/M (1 - cos chi)]
    eps0 = 10d0
    call load_mono(1, eps0, n)
    v0x = vx(1:n,1); v0y = vy(1:n,1); v0z = vz(1:n,1)
    call mcc_step()
    maxerr = 0d0; nchanged = 0
    do i = 1, n
      if (vx(i,1) == v0x(i) .and. vy(i,1) == v0y(i) .and. vz(i,1) == v0z(i)) cycle
      nchanged = nchanged + 1
      v1 = sqrt(vx(i,1)**2 + vy(i,1)**2 + vz(i,1)**2)
      cchi = (vx(i,1)*v0x(i) + vy(i,1)*v0y(i) + vz(i,1)*v0z(i))/(v1*sqrt(v0x(i)**2 + v0y(i)**2 + v0z(i)**2))
      e1 = 0.5d0*me*v1*v1/qe
      pred = eps0*(1d0 - floss*(1d0 - cchi))
      maxerr = max(maxerr, abs(e1 - pred)/eps0)
    end do
    write(*,'(a,i0,a,i0,a)') '    T6a 10 eV electrons: ', nchanged, ' collided in one step (counter ', n_el_step, ')'
    call report('T6a elastic energy loss formula, max relative error', maxerr, 1d-12)
    call report('T6a no inelastic events below threshold', dble(n_exc_step(1) + n_exc_step(2) + n_ion_step + (np(1) - n)), 0d0)
    call report('T6a changed particles = elastic counter', dble(abs(nchanged - n_el_step)), 0d0)
    ! (b) 100 eV: classify every changed electron by its energy
    eps0 = 100d0
    call load_mono(1, eps0, n)
    v0x = vx(1:n,1); v0y = vy(1:n,1); v0z = vz(1:n,1)
    call mcc_step()
    nel = 0; nx1 = 0; nx2 = 0; nio = 0; nbad = 0
    do i = 1, n
      if (vx(i,1) == v0x(i) .and. vy(i,1) == v0y(i) .and. vz(i,1) == v0z(i)) cycle
      v1 = sqrt(vx(i,1)**2 + vy(i,1)**2 + vz(i,1)**2)
      e1 = 0.5d0*me*v1*v1/qe
      cchi = (vx(i,1)*v0x(i) + vy(i,1)*v0y(i) + vz(i,1)*v0z(i))/(v1*sqrt(v0x(i)**2 + v0y(i)**2 + v0z(i)**2))
      if (abs(e1 - eps0*(1d0 - floss*(1d0 - cchi))) < 1d-9*eps0) then
        nel = nel + 1
      else if (abs(e1 - (eps0 - eth_e(2))) < 1d-9*eps0) then
        nx1 = nx1 + 1
      else if (abs(e1 - (eps0 - eth_e(3))) < 1d-9*eps0) then
        nx2 = nx2 + 1
      else if (abs(e1 - 0.5d0*(eps0 - eth_e(4))) < 1d-9*eps0) then
        nio = nio + 1
      else
        nbad = nbad + 1
      end if
    end do
    err = 0d0
    do i = n + 1, np(1)                 ! created electrons: (eps0 - 24.59)/2 each
      e0 = 0.5d0*me*(vx(i,1)**2 + vy(i,1)**2 + vz(i,1)**2)/qe
      err = max(err, abs(e0 - 0.5d0*(eps0 - eth_e(4)))/eps0)
    end do
    write(*,'(a,4i8,a,4i8)') '    T6b 100 eV electrons classified el/exc1/exc2/ion: ', nel, nx1, nx2, nio, '   counters: ', &
         n_el_step, n_exc_step(1), n_exc_step(2), n_ion_step
    call report('T6b every collided electron matches one process energy (unclassified)', dble(nbad), 0d0)
    call report('T6b classification = counters', dble(abs(nel - n_el_step) + abs(nx1 - n_exc_step(1)) + abs(nx2 - n_exc_step(2)) + abs(nio - n_ion_step)), 0d0)
    call report('T6b new electrons = ionizations; new ions = ionizations', dble(abs(np(1) - n - n_ion_step) + abs(np(2) - n_ion_step)), 0d0)
    call report('T6b energy of the created electrons = (eps - 24.59)/2 (relative)', err, 1d-12)
    call report('T6b inelastic ledger = 19.82 n1 + 20.61 n2 + 24.59 n_ion (relative)', &
         abs(einel_step/qe - (eth_e(2)*n_exc_step(1) + eth_e(3)*n_exc_step(2) + eth_e(4)*n_ion_step))/(einel_step/qe), 1d-12)
    ! (c) velocity variance of the created ions over several steps: kT_g/M per component
    vn2 = 0d0; nio = 0
    do i = 1, 40
      call load_mono(1, eps0, n)
      call mcc_step()
      vn2 = vn2 + sum(vx(1:np(2),2)**2 + vy(1:np(2),2)**2 + vz(1:np(2),2)**2)
      nio = nio + np(2)
    end do
    vn2 = vn2/(3d0*nio)
    write(*,'(a,i0,a,es12.4,a,es12.4)') '    T6c ', nio, ' created ions: <v^2>/3 = ', vn2, '  kT_g/M = ', kb*T_gas/mion
    call report('T6c created ions thermal at T_g (relative, 3 sigma = ', abs(vn2/(kb*T_gas/mion) - 1d0), 3d0*sqrt(2d0/(3d0*nio)))
    deallocate(v0x, v0y, v0z)
  end subroutine kinematics_electrons
  !
  subroutine kinematics_ions()
    integer(8), parameter :: n = 200000
    integer(8) :: i, niso, nback, nbad
    real(8), allocatable :: v0x(:), v0y(:), v0z(:)
    real(8) :: dxx, dyy, dzz, g, v1, eexp
    T_gas = 0d0                          ! v_n = 0: isotropic -> |v' - v/2| = |v|/2, backward -> v' = 0
    call mcc_init()
    allocate(v0x(n), v0y(n), v0z(n))
    call load_mono(2, 10d0, n)
    v0x = vx(1:n,2); v0y = vy(1:n,2); v0z = vz(1:n,2)
    call mcc_step()
    niso = 0; nback = 0; nbad = 0; eexp = 0d0
    do i = 1, n
      if (vx(i,2) == v0x(i) .and. vy(i,2) == v0y(i) .and. vz(i,2) == v0z(i)) cycle
      v1 = sqrt(vx(i,2)**2 + vy(i,2)**2 + vz(i,2)**2)
      g = sqrt(v0x(i)**2 + v0y(i)**2 + v0z(i)**2)
      dxx = vx(i,2) - 0.5d0*v0x(i); dyy = vy(i,2) - 0.5d0*v0y(i); dzz = vz(i,2) - 0.5d0*v0z(i)
      eexp = eexp + 0.5d0*mion*(g*g - v1*v1)              ! energy given to the gas by this collision
      if (v1 == 0d0) then
        nback = nback + 1
      else if (abs(sqrt(dxx*dxx + dyy*dyy + dzz*dzz) - 0.5d0*g) < 1d-12*g) then
        niso = niso + 1
      else
        nbad = nbad + 1
      end if
    end do
    write(*,'(a,2i8,a,2i8)') '    T6d cold-gas ion collisions iso/back: ', niso, nback, '   counters: ', n_icoll_step
    call report("T6d every collided ion satisfies the CM invariant or v' = v_n (unclassified)", dble(nbad), 0d0)
    call report('T6d classification = counters', dble(abs(niso - n_icoll_step(1)) + abs(nback - n_icoll_step(2))), 0d0)
    call report('T6d ion energy ledger = sum(E0 - E1) (relative)', abs(egas_i_step - eexp)/eexp, 1d-10)
    T_gas = 300d0
    deallocate(v0x, v0y, v0z)
  end subroutine kinematics_ions
  !
  subroutine isotropy()
    integer, parameter :: nb = 20
    integer(8), parameter :: n = 200000
    integer(8) :: i, nc, h1(nb), h2(nb), h3(nb)
    integer :: it, k
    real(8), allocatable :: v0x(:), v0y(:), v0z(:)
    real(8) :: v1, v0, cz, cchi, ph, chi2(3), expc
    T_gas = 300d0
    call mcc_init()
    allocate(v0x(n), v0y(n), v0z(n))
    h1 = 0; h2 = 0; h3 = 0; nc = 0
    do it = 1, 30
      call load_mono(1, 10d0, n)
      v0x = vx(1:n,1); v0y = vy(1:n,1); v0z = vz(1:n,1)
      call mcc_step()
      do i = 1, n
        if (vx(i,1) == v0x(i) .and. vy(i,1) == v0y(i) .and. vz(i,1) == v0z(i)) cycle
        nc = nc + 1
        v1 = sqrt(vx(i,1)**2 + vy(i,1)**2 + vz(i,1)**2)
        v0 = sqrt(v0x(i)**2 + v0y(i)**2 + v0z(i)**2)
        cz = vz(i,1)/v1                                                ! polar cosine of the new direction
        cchi = (vx(i,1)*v0x(i) + vy(i,1)*v0y(i) + vz(i,1)*v0z(i))/(v1*v0)  ! cosine of the scattering angle
        ph = atan2(vy(i,1), vx(i,1))                                   ! azimuth of the new direction
        k = min(nb, int((cz + 1d0)/2d0*nb) + 1);  h1(k) = h1(k) + 1
        k = min(nb, int((cchi + 1d0)/2d0*nb) + 1); h2(k) = h2(k) + 1
        k = min(nb, int((ph + pi)/(2d0*pi)*nb) + 1); h3(k) = h3(k) + 1
      end do
    end do
    expc = dble(nc)/nb
    chi2 = 0d0
    do k = 1, nb
      chi2(1) = chi2(1) + (h1(k) - expc)**2/expc
      chi2(2) = chi2(2) + (h2(k) - expc)**2/expc
      chi2(3) = chi2(3) + (h3(k) - expc)**2/expc
    end do
    write(*,'(a,i0,a,3f8.1,a)') '    T6e isotropy of ', nc, ' elastic scatterings: chi^2 (19 dof) for cos(theta_z), cos(chi), phi = ', chi2, '  (mean 19, p=0.001 at 43.8)'
    call report('T6e max chi^2 of the three uniformity tests', maxval(chi2), 43.8d0)
    deallocate(v0x, v0y, v0z)
  end subroutine isotropy
  !
  subroutine relax_ions()
    integer(8), parameter :: n = 20000
    integer :: it, nsteps
    real(8) :: emean, vv(3), tgt
    T_gas = 300d0
    freq = 1d6; steps_per_period = 500           ! dt = 2e-9 s
    call derive_parameters(1)
    call mcc_init()
    call load_mono(2, 10d0, n)
    nsteps = 10000                                ! 20 us, about 20 thermal collisions per ion
    write(*,'(a,es10.3,a,f6.3)') '    T7a ions 10 eV -> 300 K gas, dt = ', dt, ' s, numax dt = ', numax_i*dt
    do it = 1, nsteps
      call mcc_step()
      if (mod(it, 2000) == 0) then
        emean = 0.5d0*mion*sum(vx(1:n,2)**2 + vy(1:n,2)**2 + vz(1:n,2)**2)/n/qe
        write(*,'(a,f7.2,a,f9.5,a)') '        t = ', it*dt*1d6, ' us: <eps_i> = ', emean, ' eV'
      end if
    end do
    emean = 0.5d0*mion*sum(vx(1:n,2)**2 + vy(1:n,2)**2 + vz(1:n,2)**2)/n/qe
    vv(1) = sum(vx(1:n,2)**2)/n; vv(2) = sum(vy(1:n,2)**2)/n; vv(3) = sum(vz(1:n,2)**2)/n
    tgt = kb*T_gas/mion
    write(*,'(a,f8.5,a,f8.5,a,3f8.4)') '    T7a final <eps_i> = ', emean, ' eV, 3/2 kT_g = ', 1.5d0*kb*T_gas/qe, ' eV; <v_k^2>/(kT/M) = ', vv/tgt
    call report('T7a ion mean energy / (3/2 kT_g) - 1', abs(emean/(1.5d0*kb*T_gas/qe) - 1d0), 0.05d0)
    call report('T7a max |<v_k^2>/(kT_g/M) - 1|', maxval(abs(vv/tgt - 1d0)), 0.06d0)
    call report('T7a ions neither created nor lost', dble(abs(np(2) - n) + np(1)), 0d0)
  end subroutine relax_ions
  !
  subroutine relax_electrons()
    integer(8), parameter :: n = 20000
    integer :: it, nsteps, k
    real(8) :: emean, eode, fb, h, k1, k2, k3, k4
    freq = 1d6; steps_per_period = 5000           ! dt = 2e-10 s
    call derive_parameters(1)
    call mcc_init()
    call load_mono(1, 10d0, n)
    nsteps = 10000                                ! 2 us
    fb = (1d0 - exp(-numax_e*dt))/(numax_e*dt)
    do it = 1, nsteps
      call mcc_step()
    end do
    emean = 0.5d0*me*sum(vx(1:n,1)**2 + vy(1:n,1)**2 + vz(1:n,1)**2)/n/qe
    ! reference: d eps/dt = -(2 m/M) nu_el(eps) eps fb (isotropic: <1 - cos chi> = 1), RK4
    eode = 10d0; h = 1d-9
    do k = 1, nint(nsteps*dt/h)
      k1 = cool_rate(eode, fb); k2 = cool_rate(eode + 0.5d0*h*k1, fb); k3 = cool_rate(eode + 0.5d0*h*k2, fb); k4 = cool_rate(eode + h*k3, fb)
      eode = eode + h*(k1 + 2d0*k2 + 2d0*k3 + k4)/6d0
    end do
    write(*,'(a,f9.5,a,f9.5,a,f8.5,a)') '    T7b electrons 10 eV, elastic cooling 2 us: <eps> = ', emean, ' eV, ODE = ', eode, &
         ' eV (lost ', 10d0 - eode, ' eV)'
    call report('T7b |<eps> - ODE| / energy lost', abs(emean - eode)/(10d0 - eode), 0.03d0)
    call report('T7b no inelastic events, no new particles', dble(abs(np(1) - n) + np(2)), 0d0)
  end subroutine relax_electrons
  !
  ! d eps/dt of the elastic cooling, eps in eV
  function cool_rate(e, fb) result(r)
    real(8), intent(in) :: e, fb
    real(8) :: r, sg(nproc_e)
    call xsec_e(e, sg)
    r = -(2d0*me/mion)*n_gas*sg(1)*sqrt(2d0*e*qe/me)*e*fb
  end function cool_rate
  !
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
end program test_mcc
