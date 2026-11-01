!=====================================================================
!  mcc: Monte Carlo collisions with the fixed helium background by the
!  null-collision method (Skullerud 1968; Vahedi & Surendra 1995).
!  Every particle of a species collides with the same probability
!  P = 1 - exp(-numax dt) per step; the process is then chosen from the
!  real frequencies nu_k(eps) = n_g sigma_k v, the remainder up to numax
!  being a null collision.  Rules of the Turner et al. (2013) benchmark:
!    electrons: all scattering isotropic; elastic energy loss
!               eps' = eps [1 - 2 m_e/M (1 - cos chi)]; excitation loses
!               the threshold; ionization shares the remaining energy
!               equally between the two electrons and creates a He+ with
!               the gas Maxwellian at the same position, same weight.
!    ions:      target velocity drawn from the gas Maxwellian, relative
!               energy 1/2 mu g^2 (mu = M/2); isotropic scattering in the
!               centre-of-mass frame, v' = v_cm + |g|/2 n; backward
!               scattering n = -g/|g| gives v' = v_n (charge exchange).
!  Per-step ledgers (per unit weight, energies in J) are read by diag.
!=====================================================================
module mcc
  use global
  use random
  use xsec
  implicit none
  private
  public :: mcc_init, mcc_step, mcc_electrons, mcc_ions, isotropic_direction, &
            p_coll_e, p_coll_i, n_test_step, n_el_step, n_exc_step, n_ion_step, n_icoll_step, &
            numax_exceed_step, eloss_el_step, einel_step, enewion_step, egas_i_step, nion_buf, ion_x
  real(8) :: p_coll_e = 0d0, p_coll_i = 0d0
  integer(8) :: n_test_step(nsp) = 0          ! particles that passed the P test
  integer(8) :: n_el_step = 0, n_exc_step(2) = 0, n_ion_step = 0, n_icoll_step(2) = 0
  integer(8) :: numax_exceed_step = 0
  real(8) :: eloss_el_step = 0d0              ! electron energy lost to the gas in elastic collisions
  real(8) :: einel_step = 0d0                 ! threshold energy spent in excitation and ionization
  real(8) :: enewion_step = 0d0               ! kinetic energy of the created ions (from the gas)
  real(8) :: egas_i_step = 0d0                ! ion kinetic energy lost to the gas (negative: gained)
  integer :: nion_buf = 0                     ! ionization positions of this step
  real(8), allocatable :: ion_x(:)
  real(8) :: vn_sigma = 0d0, mu_i = 0d0
contains
  !
  subroutine mcc_init()
    if (.not. xsec_loaded) call abort_run('mcc_init called before xsec_load')
    p_coll_e = 1d0 - exp(-numax_e*dt)
    p_coll_i = 1d0 - exp(-numax_i*dt)
    vn_sigma = sqrt(kb*T_gas/mion)
    mu_i = 0.5d0*mion
    if (.not. allocated(ion_x)) allocate(ion_x(4096))
    nion_buf = 0
  end subroutine mcc_init
  !
  subroutine mcc_step()
    n_test_step = 0; n_el_step = 0; n_exc_step = 0; n_ion_step = 0; n_icoll_step = 0; numax_exceed_step = 0
    eloss_el_step = 0d0; einel_step = 0d0; enewion_step = 0d0; egas_i_step = 0d0
    nion_buf = 0
    call mcc_electrons()
    call mcc_ions()
  end subroutine mcc_step
  !
  ! unit vector uniformly distributed on the sphere
  subroutine isotropic_direction(ux, uy, uz)
    real(8), intent(out) :: ux, uy, uz
    real(8) :: ct, st, ph
    ct = 1d0 - 2d0*rand_u()
    st = sqrt(max(0d0, 1d0 - ct*ct))
    ph = twopi*rand_u()
    ux = st*cos(ph); uy = st*sin(ph); uz = ct
  end subroutine isotropic_direction
  !
  subroutine mcc_electrons()
    integer(8) :: i, n0
    real(8) :: v2, v, eps, sig(nproc_e), nu(nproc_e), nutot, r, c1, c2, c3
    real(8) :: ux, uy, uz, cchi, eps1, eps2, vnew, vix, viy, viz
    real(8), parameter :: floss = 2d0*me/mion
    n0 = np(1)                                  ! electrons created this step are not tested again
    do i = 1, n0
      if (rand_u() >= p_coll_e) cycle
      n_test_step(1) = n_test_step(1) + 1
      v2 = vx(i,1)**2 + vy(i,1)**2 + vz(i,1)**2
      if (v2 <= 0d0) cycle
      v = sqrt(v2)
      eps = 0.5d0*me*v2/qe                      ! eV
      call xsec_e(eps, sig)
      nu = n_gas*sig*v
      c1 = nu(1); c2 = c1 + nu(2); c3 = c2 + nu(3); nutot = c3 + nu(4)
      if (nutot > numax_e) numax_exceed_step = numax_exceed_step + 1
      r = rand_u()*numax_e
      if (r >= nutot) cycle                     ! null collision
      call isotropic_direction(ux, uy, uz)
      if (r < c1) then
        ! elastic: new direction, energy loss 2 m/M (1 - cos chi), chi = angle between old and new direction
        cchi = (vx(i,1)*ux + vy(i,1)*uy + vz(i,1)*uz)/v
        eps1 = eps*(1d0 - floss*(1d0 - cchi))
        vnew = sqrt(2d0*eps1*qe/me)
        vx(i,1) = vnew*ux; vy(i,1) = vnew*uy; vz(i,1) = vnew*uz
        eloss_el_step = eloss_el_step + (eps - eps1)*qe
        n_el_step = n_el_step + 1
      else if (r < c3) then
        ! excitation 1 or 2: threshold lost, new isotropic direction
        if (r < c2) then
          eps1 = eps - eth_e(2); einel_step = einel_step + eth_e(2)*qe; n_exc_step(1) = n_exc_step(1) + 1
        else
          eps1 = eps - eth_e(3); einel_step = einel_step + eth_e(3)*qe; n_exc_step(2) = n_exc_step(2) + 1
        end if
        vnew = sqrt(2d0*max(eps1, 0d0)*qe/me)
        vx(i,1) = vnew*ux; vy(i,1) = vnew*uy; vz(i,1) = vnew*uz
      else
        ! ionization: remaining energy shared equally, independent isotropic directions;
        ! new electron and a thermal He+ at the parent position
        eps1 = 0.5d0*max(eps - eth_e(4), 0d0)
        eps2 = eps1
        einel_step = einel_step + eth_e(4)*qe
        vnew = sqrt(2d0*eps1*qe/me)
        vx(i,1) = vnew*ux; vy(i,1) = vnew*uy; vz(i,1) = vnew*uz
        call isotropic_direction(ux, uy, uz)
        vnew = sqrt(2d0*eps2*qe/me)
        call add_particle(1, x(i,1), vnew*ux, vnew*uy, vnew*uz)
        vix = vn_sigma*rand_gauss(); viy = vn_sigma*rand_gauss(); viz = vn_sigma*rand_gauss()
        call add_particle(2, x(i,1), vix, viy, viz)
        enewion_step = enewion_step + 0.5d0*mion*(vix*vix + viy*viy + viz*viz)
        n_ion_step = n_ion_step + 1
        call record_ionization(x(i,1))
      end if
    end do
  end subroutine mcc_electrons
  !
  subroutine mcc_ions()
    integer(8) :: i, n0
    real(8) :: vnx, vny, vnz, gx, gy, gz, g2, g, eps, sig(nproc_i), nu1, nu2, r, ux, uy, uz, hg
    real(8) :: cx, cy, cz, e0, e1
    n0 = np(2)
    do i = 1, n0
      if (rand_u() >= p_coll_i) cycle
      n_test_step(2) = n_test_step(2) + 1
      vnx = vn_sigma*rand_gauss(); vny = vn_sigma*rand_gauss(); vnz = vn_sigma*rand_gauss()
      gx = vx(i,2) - vnx; gy = vy(i,2) - vny; gz = vz(i,2) - vnz
      g2 = gx*gx + gy*gy + gz*gz
      if (g2 <= 0d0) cycle
      g = sqrt(g2)
      eps = 0.5d0*mu_i*g2/qe                    ! relative energy [eV]
      call xsec_i(eps, sig)
      nu1 = n_gas*sig(1)*g
      nu2 = n_gas*sig(2)*g
      if (nu1 + nu2 > numax_i) numax_exceed_step = numax_exceed_step + 1
      r = rand_u()*numax_i
      if (r >= nu1 + nu2) cycle                 ! null collision
      e0 = 0.5d0*mion*(vx(i,2)**2 + vy(i,2)**2 + vz(i,2)**2)
      if (r < nu1) then
        ! isotropic scattering in the centre-of-mass frame
        cx = 0.5d0*(vx(i,2) + vnx); cy = 0.5d0*(vy(i,2) + vny); cz = 0.5d0*(vz(i,2) + vnz)
        call isotropic_direction(ux, uy, uz)
        hg = 0.5d0*g
        vx(i,2) = cx + hg*ux; vy(i,2) = cy + hg*uy; vz(i,2) = cz + hg*uz
        n_icoll_step(1) = n_icoll_step(1) + 1
      else
        ! backward scattering = charge exchange: the ion takes the neutral's velocity
        vx(i,2) = vnx; vy(i,2) = vny; vz(i,2) = vnz
        n_icoll_step(2) = n_icoll_step(2) + 1
      end if
      e1 = 0.5d0*mion*(vx(i,2)**2 + vy(i,2)**2 + vz(i,2)**2)
      egas_i_step = egas_i_step + (e0 - e1)
    end do
  end subroutine mcc_ions
  !
  subroutine record_ionization(xp)
    real(8), intent(in) :: xp
    real(8), allocatable :: tmp(:)
    if (nion_buf >= size(ion_x)) then
      allocate(tmp(2*size(ion_x)))
      tmp(1:nion_buf) = ion_x(1:nion_buf)
      call move_alloc(tmp, ion_x)
    end if
    nion_buf = nion_buf + 1
    ion_x(nion_buf) = xp
  end subroutine record_ionization
  !
end module mcc
