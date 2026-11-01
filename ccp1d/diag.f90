!=====================================================================
!  diag: time histories, time averages over the measurement window
!  (the last n_avg periods), wall-energy histograms, the x-t maps,
!  the energy balance, info.txt and the restart files.
!
!  Accumulators come in two kinds.  Particle-based ones are summed over
!  this rank's particles only ("per unit weight") and are reduced over
!  the ranks once at the end of the window (prefix w_, rank files on
!  restart).  Grid-based ones are identical on every rank and are never
!  reduced (restart_grid.bin).  All files are written by rank 0.
!=====================================================================
module diag
  use global
  use random
  use poisson
  use deposit
  use push
  use xsec
  use mcc
  use version
  implicit none
  private
  public :: diag_init, diag_step, diag_period, diag_final, diag_window_check, &
            restart_write, restart_read, write_info_header, dump_initial_state, &
            history_open, in_window, numax_exceed, wall_clock
  !
  ! ---- per-period tallies (rank-local, reset in diag_period) ----
  real(8)    :: p_work(nsp) = 0d0, p_ewall(2,nsp) = 0d0, p_vi = 0d0
  integer(8) :: p_nwall(2,nsp) = 0
  real(8)    :: p_gauss_max = 0d0            ! grid-based
  integer(8) :: p_nion = 0, p_nexc = 0        ! filled by mcc (M2)
  ! ---- cumulative since t = 0 (rank-local) ----
  real(8)    :: c_work(nsp) = 0d0, c_ewall(nsp) = 0d0
  integer(8) :: numax_exceed = 0
  !
  ! ---- measurement window ----
  logical    :: in_window = .false.
  integer(8) :: nstep_win = 0                                   ! grid-based
  real(8), allocatable :: w_dens(:,:), w_flux(:,:), w_ekd(:,:), w_ekxd(:,:), w_pe(:,:), w_pw(:,:)   ! (0:nx,nsp)
  real(8), allocatable :: w_ioniz(:)                            ! (0:nx), M2
  real(8), allocatable :: w_phi(:), w_ex(:)                     ! grid-based
  integer, parameter   :: nxt = 7                               ! ne, ni, Je, Ji, Je.E, ionization, electron work density
  real(8), allocatable :: w_xt(:,:,:)                           ! (0:nx, 0:n_phase-1, nxt)
  real(8), allocatable :: w_xt_phi(:,:)                         ! grid-based
  integer(8), allocatable :: w_xt_count(:)                      ! grid-based
  real(8), allocatable :: w_hwall(:,:,:)                        ! (nbin_iedf, side, isp)
  real(8), allocatable :: w_eepf(:)                             ! (nbin_eepf)
  real(8) :: w_work(nsp), w_ewall(2,nsp), w_nwall(2,nsp), w_pje(nsp), w_vi_cond
  real(8) :: w_ke_start(nsp), w_ke_end(nsp)
  real(8) :: w_nion, w_nexc(2), w_eloss_el, w_einel, w_egas_i, w_enewion   ! collision ledgers
  real(8) :: w_vi_disp, w_wf_start, w_wf_end, w_time_start, w_time_end ! grid-based
  integer :: w_period_start = 0                                        ! grid-based
  !
  integer :: u_hist = -1
  real(8) :: wall0 = 0d0
  real(8), allocatable :: dens(:), flux(:), ekd(:), ekxd(:), emid(:), ef(:)
  character(len=8), parameter :: restart_magic = 'CCP1DR03'
  !
contains
  !
  !--------------------------------------------------------------------
  subroutine diag_init()
    allocate(w_dens(0:nx,nsp), w_flux(0:nx,nsp), w_ekd(0:nx,nsp), w_ekxd(0:nx,nsp), w_pe(0:nx,nsp), w_pw(0:nx,nsp), w_ioniz(0:nx))
    allocate(w_phi(0:nx), w_ex(0:nx))
    allocate(w_xt(0:nx,0:n_phase-1,nxt), w_xt_phi(0:nx,0:n_phase-1), w_xt_count(0:n_phase-1))
    allocate(w_hwall(nbin_iedf,2,nsp), w_eepf(nbin_eepf))
    allocate(dens(0:nx), flux(0:nx), ekd(0:nx), ekxd(0:nx), emid(0:nx), ef(0:nx-1))
    call window_reset()
    in_window = .false.
    wall0 = MPI_WTIME()
  end subroutine diag_init
  !
  function wall_clock() result(t)
    real(8) :: t
    t = MPI_WTIME() - wall0
  end function wall_clock
  !
  subroutine window_reset()
    w_dens = 0d0; w_flux = 0d0; w_ekd = 0d0; w_ekxd = 0d0; w_pe = 0d0; w_pw = 0d0; w_ioniz = 0d0
    w_phi = 0d0; w_ex = 0d0; w_xt = 0d0; w_xt_phi = 0d0; w_xt_count = 0
    w_hwall = 0d0; w_eepf = 0d0
    w_work = 0d0; w_ewall = 0d0; w_nwall = 0d0; w_pje = 0d0; w_vi_cond = 0d0
    w_ke_start = 0d0; w_ke_end = 0d0
    w_nion = 0d0; w_nexc = 0d0; w_eloss_el = 0d0; w_einel = 0d0; w_egas_i = 0d0; w_enewion = 0d0
    w_vi_disp = 0d0; w_wf_start = 0d0; w_wf_end = 0d0; w_time_start = 0d0; w_time_end = 0d0
    nstep_win = 0
  end subroutine window_reset
  !
  ! trapezoidal integral over the nodes
  function trapz(f) result(s)
    real(8), intent(in) :: f(0:nx)
    real(8) :: s
    s = dx*(sum(f(1:nx-1)) + 0.5d0*(f(0) + f(nx)))
  end function trapz
  !
  ! electrostatic energy per unit area [J/m^2] from the face fields
  function wfield() result(w)
    real(8) :: w
    call efield_faces(phi, ef)
    w = 0.5d0*eps0*dx*sum(ef*ef)
  end function wfield
  !
  ! sum a rank-local array onto rank 0
  subroutine reduce0(a, n)
    integer, intent(in) :: n
    real(8), intent(inout) :: a(*)
    real(8) :: dummy
    if (myrank == 0) then
      call MPI_REDUCE(MPI_IN_PLACE, a, n, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    else
      call MPI_REDUCE(a, dummy, n, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    end if
  end subroutine reduce0
  !
  !--------------------------------------------------------------------
  ! called at the start of every period: open or (re)set the window
  subroutine diag_window_check()
    if (in_window .and. iperiod < n_periods - n_avg) then
      if (myrank == 0) write(*,'(a)') 'ccp1d: n_periods was increased after the window had started;'// &
           ' the window accumulators are reset'
      in_window = .false.
    end if
    if (.not. in_window .and. iperiod >= n_periods - n_avg) then
      call window_reset()
      in_window = .true.
      w_period_start = iperiod + 1
      w_ke_start(1) = kinetic_energy(1)
      w_ke_start(2) = kinetic_energy(2)
      w_wf_start = wfield()
      w_time_start = time
      if (myrank == 0) write(*,'(a,i0)') 'ccp1d: measurement window starts with period ', iperiod + 1
    end if
    tally_work = in_window
  end subroutine diag_window_check
  !
  !--------------------------------------------------------------------
  ! called after the field solve of every step (it, time already advanced)
  subroutine diag_step()
    integer :: isp, j, p, k
    integer(8) :: i
    real(8) :: g1, g2, res, vmid, iramo, idisp, e, scale
    !
    ! external circuit (Shockley-Ramo): I = (weight/L) sum q vx + (eps0/L) dV/dt, power V(t_{n+1/2}) I
    vmid = v_drive(time - 0.5d0*dt)
    iramo = 0d0
    do isp = 1, nsp
      iramo = iramo + q(isp)*vxsum_step(isp)/L          ! per unit weight
    end do
    idisp = eps0*(v_drive(time) - v_drive(time - dt))/(dt*L)
    p_vi = p_vi + vmid*(iramo*weight + idisp)
    p_work = p_work + work_step
    p_ewall = p_ewall + ewall_step
    p_nwall = p_nwall + nwall_step
    p_nion = p_nion + n_ion_step
    p_nexc = p_nexc + n_exc_step(1) + n_exc_step(2)
    numax_exceed = numax_exceed + numax_exceed_step
    c_work = c_work + work_step
    do isp = 1, nsp
      c_ewall(isp) = c_ewall(isp) + ewall_step(1,isp) + ewall_step(2,isp)
    end do
    ! discrete Gauss law: eps0 (E_{nx-1/2} - E_{1/2}) = dx sum_{j=1}^{nx-1} rho_j (exact for the solver)
    call efield_faces(phi, ef)
    g1 = eps0*(ef(nx-1) - ef(0))
    g2 = dx*sum(rho(1:nx-1))
    scale = max(abs(g1), abs(g2), qe*weight)
    res = abs(g1 - g2)/scale
    p_gauss_max = max(p_gauss_max, res)
    !
    if (in_window) then
      nstep_win = nstep_win + 1
      emid = 0.5d0*(ex_old + ex)
      p = int(mod(it, int(steps_per_period,8))*int(n_phase,8)/int(steps_per_period,8))
      w_xt_count(p) = w_xt_count(p) + 1
      do isp = 1, nsp
        call deposit_species_moments(isp, dens, flux, ekd, ekxd)
        w_dens(:,isp) = w_dens(:,isp) + dens
        w_flux(:,isp) = w_flux(:,isp) + flux
        w_ekd(:,isp)  = w_ekd(:,isp)  + ekd
        w_ekxd(:,isp) = w_ekxd(:,isp) + ekxd
        flux = q(isp)*flux                       ! current density per unit weight
        w_pe(:,isp) = w_pe(:,isp) + flux*emid
        w_pje(isp) = w_pje(isp) + trapz(flux*emid)
        w_xt(:,p,isp)   = w_xt(:,p,isp)   + dens
        w_xt(:,p,2+isp) = w_xt(:,p,2+isp) + flux
        if (isp == 1) w_xt(:,p,5) = w_xt(:,p,5) + flux*emid
      end do
      w_pw = w_pw + work_dep
      w_xt(:,p,7) = w_xt(:,p,7) + work_dep(:,1)
      w_phi = w_phi + phi
      w_ex = w_ex + ex
      w_xt_phi(:,p) = w_xt_phi(:,p) + phi
      w_vi_cond = w_vi_cond + vmid*iramo
      w_vi_disp = w_vi_disp + vmid*idisp
      w_work = w_work + work_step
      w_ewall = w_ewall + ewall_step
      w_nwall = w_nwall + dble(nwall_step)
      w_nion = w_nion + dble(n_ion_step)
      w_nexc = w_nexc + dble(n_exc_step)
      w_eloss_el = w_eloss_el + eloss_el_step
      w_einel = w_einel + einel_step
      w_egas_i = w_egas_i + egas_i_step
      w_enewion = w_enewion + enewion_step
      ! ionization positions -> rate profile and x-t map (CIC, per unit weight)
      do k = 1, nion_buf
        g1 = ion_x(k)/dx
        j = min(int(g1), nx - 1)
        g2 = g1 - j
        w_ioniz(j)   = w_ioniz(j)   + (1d0 - g2)/dx
        w_ioniz(j+1) = w_ioniz(j+1) + g2/dx
        w_xt(j,p,6)   = w_xt(j,p,6)   + (1d0 - g2)/dx
        w_xt(j+1,p,6) = w_xt(j+1,p,6) + g2/dx
      end do
      ! wall energy histograms (bin k covers [(k-1) de, k de))
      do k = 1, nabs
        j = int(abs_energy(k)/(qe*iedf_de)) + 1
        if (j >= 1 .and. j <= nbin_iedf) w_hwall(j,abs_side(k),abs_isp(k)) = w_hwall(j,abs_side(k),abs_isp(k)) + 1d0
      end do
      ! electron energy histogram in the central 10 % of the gap
      do i = 1, np(1)
        if (abs(x(i,1) - 0.5d0*L) < 0.05d0*L) then
          e = 0.5d0*me*(vx(i,1)**2 + vy(i,1)**2 + vz(i,1)**2)/qe
          j = int(e/eepf_de) + 1
          if (j <= nbin_eepf) w_eepf(j) = w_eepf(j) + 1d0
        end if
      end do
    end if
    nabs = 0
  end subroutine diag_step
  !
  !--------------------------------------------------------------------
  subroutine history_open(append)
    logical, intent(in) :: append
    if (myrank /= 0) return
    if (append) then
      open(newunit=u_hist, file='history.dat', status='old', position='append', action='write')
    else
      open(newunit=u_hist, file='history.dat', status='replace', action='write')
      write(u_hist,'(a)') '# ccp1d history: one line per RF period (period 0 = initial state)'
      write(u_hist,'(a)') '# 1 period  2 t[s]  3 N_e  4 N_i  5 <n_e>[m^-3]  6 <n_i>[m^-3]  7 n_e_center[m^-3]  8 n_i_center[m^-3]'
      write(u_hist,'(a)') '# 9 Gi_L[m^-2s^-1] 10 Gi_R 11 Ge_L 12 Ge_R  13 Ei_L[eV] 14 Ei_R[eV]  15 P_e[W/m^2] 16 P_i[W/m^2] 17 P_VI[W/m^2]'
      write(u_hist,'(a)') '# 18 n_ioniz 19 n_exc  20 KE_e[J/m^2] 21 KE_i[J/m^2] 22 W_field[J/m^2] 23 E_wall_cum[J/m^2] 24 Work_cum[J/m^2]'
      write(u_hist,'(a)') '# 25 gauss_res_max  26 wall[s]  27 Tx_e[eV] 28 Tperp_e[eV] 29 Tx_i[eV] 30 Tperp_i[eV]  31 numax_exceed'
      write(u_hist,'(a)') '# Tx = m<vx^2>/e, Tperp = m<vy^2+vz^2>/(2e): temperatures of a Maxwellian with the same mean energies'
      write(u_hist,'(a)') '# P_e, P_i = work of the field on the species (exact); P_VI = <V I> of the external circuit (Shockley-Ramo current)'
      write(u_hist,'(a)') '# centre density = mean of the 5 central nodes at the end of the period; fluxes, energies, power = period averages;'
      write(u_hist,'(a)') '# KE with vx at t-dt/2 (after the push), exact identity: KE_e+KE_i - KE(0) = Work_cum - E_wall_cum'
    end if
  end subroutine history_open
  !
  ! history line at the end of a period (also called with iperiod = 0 after init)
  subroutine diag_period()
    real(8) :: buf(25), dc(nsp), ke(nsp), wf, gam(2,nsp), ei(2,nsp), pw(nsp), nn(nsp), ekx(nsp), ekp(nsp), tx(nsp), tp(nsp)
    integer :: isp, j0, j1
    j0 = max(nx/2 - 2, 1); j1 = min(nx/2 + 2, nx - 1)
    do isp = 1, nsp
      call deposit_number_density(isp, dens)
      dc(isp) = sum(dens(j0:j1))/dble(j1 - j0 + 1)
      call kinetic_sums(isp, ekx(isp), ekp(isp))
      ke(isp) = ekx(isp) + ekp(isp)
    end do
    buf(1:2)   = dble(np)
    buf(3:4)   = p_work
    buf(5:6)   = p_ewall(1:2,1)
    buf(7:8)   = p_ewall(1:2,2)
    buf(9:10)  = dble(p_nwall(1:2,1))
    buf(11:12) = dble(p_nwall(1:2,2))
    buf(13:14) = c_work
    buf(15:16) = c_ewall
    buf(17:18) = ke
    buf(19)    = dble(p_nion)
    buf(20)    = dble(p_nexc)
    buf(21)    = p_vi
    buf(22:23) = ekx
    buf(24:25) = ekp
    call reduce0(buf, 25)
    wf = wfield()
    if (myrank == 0) then
      do isp = 1, nsp
        nn(isp) = buf(isp)
        gam(1,isp) = buf(9 + 2*(isp-1))*weight/period
        gam(2,isp) = buf(10 + 2*(isp-1))*weight/period
        ei(:,isp) = 0d0
        if (buf(9 + 2*(isp-1))  > 0d0) ei(1,isp) = buf(5 + 2*(isp-1))/buf(9 + 2*(isp-1))/qe
        if (buf(10 + 2*(isp-1)) > 0d0) ei(2,isp) = buf(6 + 2*(isp-1))/buf(10 + 2*(isp-1))/qe
        pw(isp) = buf(2 + isp)*weight/period
      end do
      do isp = 1, nsp
        tx(isp) = 0d0; tp(isp) = 0d0
        if (nn(isp) > 0d0) then
          tx(isp) = 2d0*buf(21 + isp)/nn(isp)/qe
          tp(isp) = buf(23 + isp)/nn(isp)/qe
        end if
      end do
      write(u_hist,'(i8,es16.8,2f14.0,21es15.7,es13.5,4es13.5,i10)') iperiod, time, nn(1), nn(2), &
           nn(1)*weight/L, nn(2)*weight/L, dc(1), dc(2), &
           gam(1,2), gam(2,2), gam(1,1), gam(2,1), ei(1,2), ei(2,2), pw(1), pw(2), buf(21)/period, &
           buf(19), buf(20), buf(17)*weight, buf(18)*weight, wf, (buf(15) + buf(16))*weight, &
           (buf(13) + buf(14))*weight, p_gauss_max, wall_clock(), tx(1), tp(1), tx(2), tp(2), numax_exceed
      flush(u_hist)
      if (iperiod <= 5 .or. mod(iperiod, max(1, n_periods/50)) == 0 .or. iperiod == n_periods) then
        write(*,'(a,i6,a,i6,a,es11.3,a,es11.3,a,f9.1,a)') 'period ', iperiod, '/', n_periods, &
             '  N_e=', nn(1), '  N_i=', nn(2), '  wall=', wall_clock(), ' s'
        flush(6)
      end if
    end if
    p_work = 0d0; p_ewall = 0d0; p_nwall = 0; p_gauss_max = 0d0; p_nion = 0; p_nexc = 0; p_vi = 0d0
  end subroutine diag_period
  !
  !--------------------------------------------------------------------
  subroutine write_info_header()
    integer :: u
    real(8) :: wpe, ld, te_ev
    character(len=8) :: d
    character(len=10) :: t
    if (myrank /= 0) return
    call date_and_time(date=d, time=t)
    open(newunit=u, file='info.txt', status='replace', action='write')
    write(u,'(a)') '# ccp1d run information (key = value  # unit)'
    write(u,'(a,a)') 'git_commit = ', trim(git_commit)
    write(u,'(a,a,a,a,a,a,a,a,a,a,a)') 'date_start = ', d(1:4), '-', d(5:6), '-', d(7:8), ' ', t(1:2), ':', t(3:4), ':'//t(5:6)
    write(u,'(a,i0)') 'nprocs = ', nprocs
    call print_input(u)
    te_ev = kb*Te0/qe
    wpe = sqrt(n0*qe*qe/(eps0*me))
    ld = sqrt(eps0*kb*Te0/(n0*qe*qe))
    write(u,'(a)') '# ---- initial plasma parameters ----'
    write(u,'(a,es14.6,a)') 'Te0_eV = ', te_ev, '  # eV'
    write(u,'(a,es14.6,a)') 'wpe_initial = ', wpe, '  # rad/s'
    write(u,'(a,es14.6)')   'wpe_dt_initial = ', wpe*dt
    write(u,'(a,es14.6,a)') 'lambdaD_initial = ', ld, '  # m'
    write(u,'(a,es14.6)')   'lambdaD_dx_initial = ', ld/dx
    write(u,'(a,es14.6,a)') 'vth_e_initial = ', sqrt(kb*Te0/me), '  # m/s  sqrt(kT/m)'
    write(u,'(a,es14.6)')   'vth_e_dt_dx_initial = ', sqrt(kb*Te0/me)*dt/dx
    write(u,'(a,es14.6,a)') 'bohm_speed_initial = ', sqrt(kb*Te0/mion), '  # m/s'
    write(u,'(a,es14.6,a)') 'ion_transit_half_gap = ', 0.5d0*L/sqrt(kb*Te0/mion), '  # s'
    write(u,'(a,l1)')       'mcc_on = ', mcc_on
    if (mcc_on) then
      write(u,'(a)') '# ---- collisions (null-collision method) ----'
      write(u,'(a,es14.6,a)') 'numax_e = ', numax_e, '  # s^-1'
      write(u,'(a,es14.6,a)') 'numax_e_energy = ', enumax_e, '  # eV'
      write(u,'(a,es14.6,a)') 'P_coll_e = ', p_coll_e, '  # per step, 1 - exp(-numax_e dt)'
      write(u,'(a,es14.6,a)') 'numax_i = ', numax_i, '  # s^-1'
      write(u,'(a,es14.6,a)') 'numax_i_energy = ', enumax_i, '  # eV (relative energy)'
      write(u,'(a,es14.6,a)') 'P_coll_i = ', p_coll_i, '  # per step'
      write(u,'(a,3f8.3,a)')  'thresholds = ', eth_e(2:4), '  # eV (exc1, exc2, ionization)'
      write(u,'(a,es14.6,a)') 'gas_pressure = ', n_gas*kb*T_gas, '  # Pa'
    end if
    close(u)
  end subroutine write_info_header
  !
  !--------------------------------------------------------------------
  ! t = 0 particles and grid (dump_initial = .true.), before the half step
  subroutine dump_initial_state()
    integer :: u, isp, j
    integer(8) :: i
    character(len=40) :: fn
    real(8) :: ne(0:nx), ni(0:nx)
    write(fn,'(a,i3.3,a)') 'init_particles_rank', myrank, '.dat'
    open(newunit=u, file=fn, status='replace', action='write')
    write(u,'(a)') '# initial particles of this rank at t = 0 (before the leapfrog half step)'
    write(u,'(a)') '# 1 species(1=e,2=He+)  2 x[m]  3 vx[m/s]  4 vy[m/s]  5 vz[m/s]'
    do isp = 1, nsp
      do i = 1, np(isp)
        write(u,'(i2,4es16.8)') isp, x(i,isp), vx(i,isp), vy(i,isp), vz(i,isp)
      end do
    end do
    close(u)
    call deposit_number_density(1, ne)
    call deposit_number_density(2, ni)
    if (myrank == 0) then
      open(newunit=u, file='init_grid.dat', status='replace', action='write')
      write(u,'(a)') '# initial grid quantities at t = 0 (rho deposited from the loaded particles, Poisson solved with V(0))'
      write(u,'(a)') '# electrode nodes 0 and nx: densities multiplied by 2 (half-cell volume)'
      write(u,'(a)') '# 1 x[m]  2 rho[C/m^3]  3 phi[V]  4 E[V/m]  5 n_e[m^-3]  6 n_i[m^-3]'
      do j = 0, nx
        if (j == 0 .or. j == nx) then
          write(u,'(6es16.8)') j*dx, rho(j), phi(j), ex(j), 2d0*ne(j), 2d0*ni(j)
        else
          write(u,'(6es16.8)') j*dx, rho(j), phi(j), ex(j), ne(j), ni(j)
        end if
      end do
      close(u)
    end if
  end subroutine dump_initial_state
  !
  !--------------------------------------------------------------------
  ! end of the run: reduce the window accumulators and write everything
  subroutine diag_final()
    integer :: u, j, k, isp
    real(8) :: twin, sc(30), fac, e, ppart(nsp), pwall(nsp), pke(nsp), pinel, pel, pje(nsp), pvi, pgas, pnew
    real(8) :: resid, ptot, dcen(nsp), eps_c, wpe, ld, te_eff
    real(8), allocatable :: prof(:,:)
    character(len=20) :: fn
    character(len=8) :: d
    character(len=10) :: t
    if (.not. in_window .or. nstep_win == 0) then
      if (myrank == 0) write(*,'(a)') 'ccp1d: no measurement window data; profiles not written'
      return
    end if
    w_ke_end(1) = kinetic_energy(1)
    w_ke_end(2) = kinetic_energy(2)
    w_wf_end = wfield()
    w_time_end = time
    twin = w_time_end - w_time_start
    ! ---- reduce the particle-based accumulators ----
    call reduce0(w_dens, (nx+1)*nsp)
    call reduce0(w_flux, (nx+1)*nsp)
    call reduce0(w_ekd,  (nx+1)*nsp)
    call reduce0(w_ekxd, (nx+1)*nsp)
    call reduce0(w_pe,   (nx+1)*nsp)
    call reduce0(w_pw,   (nx+1)*nsp)
    call reduce0(w_ioniz, nx+1)
    call reduce0(w_xt, (nx+1)*n_phase*nxt)
    call reduce0(w_hwall, nbin_iedf*2*nsp)
    call reduce0(w_eepf, nbin_eepf)
    sc = 0d0
    sc(1:2) = w_work; sc(3:4) = w_ewall(1:2,1); sc(5:6) = w_ewall(1:2,2)
    sc(7:8) = w_nwall(1:2,1); sc(9:10) = w_nwall(1:2,2); sc(11:12) = w_pje; sc(13) = w_vi_cond
    sc(14:15) = w_ke_start; sc(16:17) = w_ke_end; sc(18) = w_nion; sc(19:20) = w_nexc
    sc(21) = w_eloss_el; sc(22) = w_einel; sc(23) = w_egas_i; sc(24) = w_enewion
    call reduce0(sc, 30)
    if (myrank /= 0) return
    !
    fac = weight/dble(nstep_win)
    ! ---- profiles.dat ----
    allocate(prof(0:nx,13))
    do j = 0, nx
      prof(j,12) = w_pw(j,1)*weight/(dble(nstep_win)*dt*dx)
      prof(j,13) = w_pw(j,2)*weight/(dble(nstep_win)*dt*dx)
      prof(j,1) = w_dens(j,1)*fac
      prof(j,2) = w_dens(j,2)*fac
      prof(j,3) = w_phi(j)/dble(nstep_win)
      prof(j,4) = w_ex(j)/dble(nstep_win)
      prof(j,5) = 0d0; prof(j,10) = 0d0; prof(j,11) = 0d0
      if (w_dens(j,1) > 0d0) prof(j,5)  = w_ekd(j,1)/w_dens(j,1)/qe
      if (w_dens(j,1) > 0d0) prof(j,10) = w_ekxd(j,1)/w_dens(j,1)/qe
      if (w_dens(j,2) > 0d0) prof(j,11) = w_ekd(j,2)/w_dens(j,2)/qe
      prof(j,6) = w_ioniz(j)*weight/(dble(nstep_win)*dt)
      prof(j,7) = w_pe(j,1)*fac
      prof(j,8) = w_pe(j,2)*fac
      prof(j,9) = w_flux(j,2)*fac
      if (j == 0 .or. j == nx) then
        prof(j,1:2) = 2d0*prof(j,1:2)
        prof(j,6:9) = 2d0*prof(j,6:9)
        prof(j,12:13) = 2d0*prof(j,12:13)
      end if
    end do
    open(newunit=u, file='profiles.dat', status='replace', action='write')
    write(u,'(a,i0,a,i0,a,i0,a)') '# ccp1d time-averaged profiles over periods ', w_period_start, '..', iperiod, &
         ' (', nstep_win, ' steps)'
    write(u,'(a)') '# electrode nodes 0 and nx: densities, rates, J.E and flux multiplied by 2 (half-cell volume)'
    write(u,'(a)') '# 1 x[m]  2 n_e[m^-3]  3 n_i[m^-3]  4 phi[V]  5 E[V/m]  6 <eps_e>[eV]  7 ionization_rate[m^-3s^-1]'// &
         '  8 <J_e.E>[W/m^3]  9 <J_i.E>[W/m^3]  10 <n_i v_i>[m^-2s^-1]  11 <eps_x,e>[eV]  12 <eps_i>[eV]'// &
         '  13 P_work_e[W/m^3]  14 P_work_i[W/m^3]'
    write(u,'(a)') '# P_work = work of the field on the particles deposited at the old positions (integrates exactly to P_part); <J.E> = grid product'
    do j = 0, nx
      write(u,'(14es15.7)') j*dx, prof(j,1:13)
    end do
    close(u)
    ! ---- wall energy distributions: flux per energy interval [m^-2 s^-1 eV^-1] ----
    do isp = 1, nsp
      do k = 1, 2
        if (isp == 2 .and. k == 1) fn = 'iedf_left.dat'
        if (isp == 2 .and. k == 2) fn = 'iedf_right.dat'
        if (isp == 1 .and. k == 1) fn = 'eedf_left.dat'
        if (isp == 1 .and. k == 2) fn = 'eedf_right.dat'
        open(newunit=u, file=trim(fn), status='replace', action='write')
        write(u,'(a,a,a,i0,a)') '# ccp1d energy distribution of the ', trim(spname(isp)), &
             ' flux to the electrode (side ', k, ': 1 = x=0 driven, 2 = x=L grounded), window average'
        write(u,'(a,es14.6,a)') '# total flux = ', sc(6 + 2*isp + k)*weight/twin, ' m^-2 s^-1 (sum of column 2 times de)'
        write(u,'(a)') '# 1 energy_bin_centre[eV]  2 flux_per_eV[m^-2 s^-1 eV^-1]'
        do j = 1, nbin_iedf
          write(u,'(2es15.7)') (j - 0.5d0)*iedf_de, w_hwall(j,k,isp)*weight/(iedf_de*twin)
        end do
        close(u)
      end do
    end do
    ! ---- EEPF in the central 10 % ----
    open(newunit=u, file='eepf.dat', status='replace', action='write')
    write(u,'(a)') '# ccp1d electron energy distribution in the central 10 % of the gap, window average'
    write(u,'(a)') '# F(eps) normalised so that integral F d(eps) = n_e; f_p = F/sqrt(eps) (EEPF, straight line for a Maxwellian)'
    write(u,'(a)') '# 1 energy_bin_centre[eV]  2 F[m^-3 eV^-1]  3 f_p[m^-3 eV^-3/2]'
    do j = 1, nbin_eepf
      e = (j - 0.5d0)*eepf_de
      write(u,'(3es15.7)') e, w_eepf(j)*weight/(eepf_de*0.1d0*L*dble(nstep_win)), &
           w_eepf(j)*weight/(eepf_de*0.1d0*L*dble(nstep_win))/sqrt(e)
    end do
    close(u)
    ! ---- x-t maps ----
    call write_xt('xt_phi.dat',   'phi[V]',             w_xt_phi, 1d0)
    call write_xt('xt_ne.dat',    'n_e[m^-3]',          w_xt(:,:,1), weight)
    call write_xt('xt_ni.dat',    'n_i[m^-3]',          w_xt(:,:,2), weight)
    call write_xt('xt_Je.dat',    'J_e[A/m^2]',         w_xt(:,:,3), weight)
    call write_xt('xt_Ji.dat',    'J_i[A/m^2]',         w_xt(:,:,4), weight)
    call write_xt('xt_PeE.dat',   'J_e.E[W/m^3]',       w_xt(:,:,5), weight)
    call write_xt('xt_ioniz.dat', 'ionization[m^-3s^-1]', w_xt(:,:,6), weight/dt)
    call write_xt('xt_Pe_work.dat', 'electron heating P_work_e[W/m^3] (exact work deposit)', w_xt(:,:,7), weight/(dt*dx))
    ! ---- energy balance and summary -> info.txt ----
    do isp = 1, nsp
      ppart(isp) = sc(isp)*weight/twin
      pwall(isp) = (sc(1 + 2*isp) + sc(2 + 2*isp))*weight/twin
      pke(isp)   = (sc(15 + isp) - sc(13 + isp))*weight/twin
      pje(isp)   = sc(10 + isp)*weight/dble(nstep_win)
    end do
    pinel = sc(22)*weight/twin
    pel   = sc(21)*weight/twin
    pgas  = sc(23)*weight/twin
    pnew  = sc(24)*weight/twin
    pvi   = (sc(13)*weight + w_vi_disp)/dble(nstep_win)
    ptot  = ppart(1) + ppart(2)
    resid = ptot + pnew - pwall(1) - pwall(2) - pke(1) - pke(2) - pinel - pel - pgas
    do isp = 1, nsp
      dcen(isp) = sum(prof(nx/2-2:nx/2+2, isp))/5d0
    end do
    eps_c = sum(prof(nx/2-2:nx/2+2, 5))/5d0
    te_eff = 2d0/3d0*eps_c
    wpe = sqrt(max(dcen(1),1d0)*qe*qe/(eps0*me))
    ld = sqrt(eps0*max(te_eff,1d-30)*qe/(max(dcen(1),1d0)*qe*qe))
    call date_and_time(date=d, time=t)
    open(newunit=u, file='info.txt', status='old', position='append', action='write')
    write(u,'(a)') '# ---- results (time averages over the measurement window) ----'
    write(u,'(a,i0,a,i0)') 'window_periods = ', w_period_start, '..', iperiod
    write(u,'(a,es14.6,a)') 'window_duration = ', twin, '  # s'
    write(u,'(a,i0)')       'window_steps = ', nstep_win
    write(u,'(a,es14.6,a)') 'n_e_center = ', dcen(1), '  # m^-3 (mean of 5 central nodes)'
    write(u,'(a,es14.6,a)') 'n_i_center = ', dcen(2), '  # m^-3'
    write(u,'(a,es14.6,a)') 'eps_e_center = ', eps_c, '  # eV mean electron energy'
    write(u,'(a,es14.6,a)') 'T_eff_center = ', te_eff, '  # eV (2/3 of the mean energy)'
    write(u,'(a,es14.6,a)') 'Tx_e_center = ', 2d0*sum(prof(nx/2-2:nx/2+2, 10))/5d0, '  # eV (m<vx^2>/e)'
    write(u,'(a,es14.6,a)') 'eps_i_center = ', sum(prof(nx/2-2:nx/2+2, 11))/5d0, '  # eV mean ion energy'
    write(u,'(a,es14.6)')   'wpe_dt_steady = ', wpe*dt
    write(u,'(a,es14.6)')   'lambdaD_dx_steady = ', ld/dx
    write(u,'(a,es14.6,a)') 'Gamma_i_left = ',  sc(9)*weight/twin,  '  # m^-2 s^-1'
    write(u,'(a,es14.6,a)') 'Gamma_i_right = ', sc(10)*weight/twin, '  # m^-2 s^-1'
    write(u,'(a,es14.6,a)') 'Gamma_e_left = ',  sc(7)*weight/twin,  '  # m^-2 s^-1'
    write(u,'(a,es14.6,a)') 'Gamma_e_right = ', sc(8)*weight/twin,  '  # m^-2 s^-1'
    if (sc(9) > 0d0)  write(u,'(a,es14.6,a)') 'E_i_mean_left = ',  sc(5)/sc(9)/qe,  '  # eV mean ion impact energy'
    if (sc(10) > 0d0) write(u,'(a,es14.6,a)') 'E_i_mean_right = ', sc(6)/sc(10)/qe, '  # eV'
    if (sc(7) > 0d0)  write(u,'(a,es14.6,a)') 'E_e_mean_left = ',  sc(3)/sc(7)/qe,  '  # eV mean electron impact energy'
    if (sc(8) > 0d0)  write(u,'(a,es14.6,a)') 'E_e_mean_right = ', sc(4)/sc(8)/qe,  '  # eV'
    write(u,'(a)') '# energy balance per unit electrode area [W/m^2]; P_part = exact work of the field on the particles'
    write(u,'(a,es14.6,a)') 'P_part_e = ', ppart(1), '  # W/m^2'
    write(u,'(a,es14.6,a)') 'P_part_i = ', ppart(2), '  # W/m^2'
    write(u,'(a,es14.6,a)') 'P_part_total = ', ptot, '  # W/m^2'
    write(u,'(a,es14.6,a)') 'P_wall_e = ', pwall(1), '  # W/m^2 electron kinetic energy flux to both electrodes'
    write(u,'(a,es14.6,a)') 'P_wall_i = ', pwall(2), '  # W/m^2'
    write(u,'(a,es14.6,a)') 'P_dKE_e = ', pke(1), '  # W/m^2 change of stored kinetic energy / window'
    write(u,'(a,es14.6,a)') 'P_dKE_i = ', pke(2), '  # W/m^2'
    write(u,'(a,es14.6,a)') 'P_inelastic = ', pinel, '  # W/m^2 excitation + ionization threshold energy'
    write(u,'(a,es14.6,a)') 'P_elastic_e = ', pel, '  # W/m^2 electron elastic energy loss to the gas'
    write(u,'(a,es14.6,a)') 'P_gas_i = ', pgas, '  # W/m^2 ion kinetic energy lost to the gas in collisions'
    write(u,'(a,es14.6,a)') 'P_newion = ', pnew, '  # W/m^2 thermal energy of the created ions (from the gas)'
    write(u,'(a,es14.6,a)') 'closure_residual = ', resid, '  # W/m^2  P_part + P_newion - walls - dKE - inelastic - elastic - P_gas_i'
    write(u,'(a,es14.6,a)') 'closure_relative = ', abs(resid)/max(abs(ptot), 1d-300), '  # should be round-off'
    write(u,'(a,es14.6,a)') 'P_work_e_integral = ', trapz(prof(:,12)), '  # W/m^2  integral of the deposited electron work (= P_part_e to round-off)'
    write(u,'(a,es14.6,a)') 'P_work_i_integral = ', trapz(prof(:,13)), '  # W/m^2'
    write(u,'(a,es14.6,a)') 'P_JE_e = ', pje(1), '  # W/m^2  integral <J_e.E> dx'
    write(u,'(a,es14.6,a)') 'P_JE_i = ', pje(2), '  # W/m^2'
    write(u,'(a,es14.6,a)') 'P_JE_total = ', pje(1) + pje(2), '  # W/m^2'
    write(u,'(a,es14.6,a)') 'P_VI = ', pvi, '  # W/m^2  <V(t) I(t)> delivered by the source, I = (w/L) sum q vx + (eps0/L) dV/dt (Shockley-Ramo)'
    write(u,'(a,es14.6,a)') 'P_VI_particles = ', sc(13)*weight/dble(nstep_win), '  # W/m^2'
    write(u,'(a,es14.6,a)') 'P_VI_vacuum_displacement = ', w_vi_disp/dble(nstep_win), '  # W/m^2 (reactive, ~0 over full periods)'
    write(u,'(a,es14.6,a)') 'dW_field_dt = ', (w_wf_end - w_wf_start)/twin, '  # W/m^2 change of field energy / window'
    write(u,'(a,es14.6,a)') 'P_numerical = ', pvi - (w_wf_end - w_wf_start)/twin - ptot, '  # W/m^2  P_VI - dW_field_dt - P_part: energy error of the scheme (numerical heating if < 0)'
    write(u,'(a,es14.6,a)') 'P_numerical_relative = ', abs(pvi - (w_wf_end - w_wf_start)/twin - ptot)/max(abs(ptot), abs(pvi), 1d-300), '  # relative to the absorbed power'
    write(u,'(a,es14.6)')   'ionizations_per_period = ', sc(18)/dble(n_avg)
    write(u,'(a,es14.6)')   'excitations_per_period = ', (sc(19) + sc(20))/dble(n_avg)
    write(u,'(a,i0)')       'numax_exceeded_count = ', numax_exceed
    write(u,'(a,i0)')       'periods_done = ', iperiod
    write(u,'(a,es14.6,a)') 'wall_time = ', wall_clock(), '  # s (this segment)'
    write(u,'(a,a,a,a,a,a,a,a,a,a,a)') 'date_end = ', d(1:4), '-', d(5:6), '-', d(7:8), ' ', t(1:2), ':', t(3:4), ':'//t(5:6)
    close(u)
    deallocate(prof)
  end subroutine diag_final
  !
  subroutine write_xt(fn, what, a, factor)
    character(len=*), intent(in) :: fn, what
    real(8), intent(in) :: a(0:nx,0:n_phase-1), factor
    integer :: u, p
    open(newunit=u, file=fn, status='replace', action='write')
    write(u,'(a,a,a)') '# ccp1d x-t map of ', trim(what), ': row p = phase bin p (t/T in [p/n_phase, (p+1)/n_phase)), column j = node x = j*dx'
    write(u,'(a,i0,a,i0,a,es14.6,a,es14.6,a,i0,a,i0)') '# nx = ', nx, '  n_phase = ', n_phase, '  dx = ', dx, &
         '  period = ', period, '  periods ', w_period_start, '..', iperiod
    do p = 0, n_phase - 1
      if (w_xt_count(p) > 0) then
        write(u,'(*(es15.7))') a(:,p)*factor/dble(w_xt_count(p))
      else
        write(u,'(*(es15.7))') a(:,p)*0d0
      end if
    end do
    close(u)
  end subroutine write_xt
  !
  !--------------------------------------------------------------------
  ! restart files: this rank's particles, random state and particle-based
  ! accumulators; rank 0 also writes the grid file
  subroutine restart_write()
    integer :: u, isp, s(rng_state_size)
    character(len=40) :: fn
    write(fn,'(a,i3.3,a)') 'restart_rank', myrank, '.bin'
    open(newunit=u, file=fn, status='replace', access='stream', form='unformatted', action='write')
    write(u) restart_magic, nx, nprocs, it, iperiod
    do isp = 1, nsp
      write(u) np(isp)
      write(u) x(1:np(isp),isp), vx(1:np(isp),isp), vy(1:np(isp),isp), vz(1:np(isp),isp)
    end do
    call rng_get_state(s)
    write(u) s
    write(u) c_work, c_ewall, numax_exceed
    write(u) in_window
    if (in_window) then
      write(u) w_dens, w_flux, w_ekd, w_ekxd, w_pe, w_pw, w_ioniz, w_xt, w_hwall, w_eepf
      write(u) w_work, w_ewall, w_nwall, w_pje, w_vi_cond, w_ke_start, w_nion, w_nexc, w_eloss_el, w_einel, w_egas_i, w_enewion
    end if
    close(u)
    if (myrank == 0) then
      open(newunit=u, file='restart_grid.bin', status='replace', access='stream', form='unformatted', action='write')
      write(u) restart_magic, nx, nprocs, it, iperiod, time
      write(u) phi, ex
      write(u) in_window
      if (in_window) then
        write(u) nstep_win, w_period_start, w_phi, w_ex, w_xt_phi, w_xt_count
        write(u) w_vi_disp, w_wf_start, w_time_start
      end if
      close(u)
      write(*,'(a,i0)') 'ccp1d: restart files written at period ', iperiod
      flush(6)
    end if
  end subroutine restart_write
  !
  subroutine restart_read()
    integer :: u, isp, s(rng_state_size), nx_f, np_f, ios
    integer(8) :: n
    character(len=8) :: magic
    character(len=40) :: fn
    write(fn,'(a,i3.3,a)') 'restart_rank', myrank, '.bin'
    open(newunit=u, file=fn, status='old', access='stream', form='unformatted', action='read', iostat=ios)
    if (ios /= 0) call abort_run('cannot open '//trim(fn))
    read(u) magic, nx_f, np_f, it, iperiod
    if (magic /= restart_magic) call abort_run('bad magic in '//trim(fn))
    if (nx_f /= nx) call abort_run('nx differs from the restart file')
    if (np_f /= nprocs) call abort_run('number of ranks differs from the restart file')
    do isp = 1, nsp
      read(u) n
      if (n > maxnp) call abort_run('restart file has more particles than maxnp: increase maxnp_factor')
      np(isp) = n
      read(u) x(1:n,isp), vx(1:n,isp), vy(1:n,isp), vz(1:n,isp)
    end do
    read(u) s
    call rng_set_state(s)
    read(u) c_work, c_ewall, numax_exceed
    read(u) in_window
    if (in_window) then
      read(u) w_dens, w_flux, w_ekd, w_ekxd, w_pe, w_pw, w_ioniz, w_xt, w_hwall, w_eepf
      read(u) w_work, w_ewall, w_nwall, w_pje, w_vi_cond, w_ke_start, w_nion, w_nexc, w_eloss_el, w_einel, w_egas_i, w_enewion
    end if
    close(u)
    open(newunit=u, file='restart_grid.bin', status='old', access='stream', form='unformatted', action='read', iostat=ios)
    if (ios /= 0) call abort_run('cannot open restart_grid.bin')
    read(u) magic, nx_f, np_f, it, iperiod, time
    if (magic /= restart_magic .or. nx_f /= nx) call abort_run('restart_grid.bin does not match')
    read(u) phi, ex
    read(u) in_window
    if (in_window) then
      read(u) nstep_win, w_period_start, w_phi, w_ex, w_xt_phi, w_xt_count
      read(u) w_vi_disp, w_wf_start, w_time_start
    end if
    close(u)
    tally_work = in_window
    if (myrank == 0) then
      write(*,'(a,i0,a,es12.4,a)') 'ccp1d: restart from period ', iperiod, ' (t = ', time, ' s)'
      open(newunit=u, file='info.txt', status='old', position='append', action='write')
      write(u,'(a,i0,a,i0)') 'restart_at_period = ', iperiod, '  # continuing to ', n_periods
      close(u)
    end if
  end subroutine restart_read
  !
end module diag
