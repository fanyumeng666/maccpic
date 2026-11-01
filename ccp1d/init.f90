!=====================================================================
!  init: fresh start.  Loads np_init_total/nprocs macro-particles of
!  each species on this rank, uniformly random in [0, L) with
!  Maxwellian velocities (Box-Muller), deposits rho, solves the field
!  at t = 0 with V(0), and pushes the velocities back by dt/2 so that
!  the leapfrog starts with vx at t = -dt/2.
!=====================================================================
module init
  use global
  use random
  use poisson
  use deposit
  use diag
  implicit none
contains
  !
  subroutine init_fresh()
    integer :: isp, j
    integer(8) :: nr, i
    real(8) :: sig, xd, s, ep, temp(nsp), xp, v1, v2, v3
    temp(1) = Te0; temp(2) = Ti0
    nr = np_init_total/nprocs
    if (myrank < mod(np_init_total, int(nprocs,8))) nr = nr + 1
    if (nr > maxnp) call abort_run('initial particles exceed maxnp')
    do isp = 1, nsp
      sig = sqrt(kb*temp(isp)/m(isp))
      do i = 1, nr
        ! separate statements: Fortran does not define the evaluation order of
        ! function references inside one argument list, and the draw order
        ! (x, vx, vy, vz: 7 uniforms per particle) must be reproducible
        xp = L*rand_u()
        v1 = sig*rand_gauss()
        v2 = sig*rand_gauss()
        v3 = sig*rand_gauss()
        call add_particle(isp, xp, v1, v2, v3)
      end do
    end do
    it = 0; iperiod = 0; time = 0d0
    call deposit_rho()
    call poisson_solve(rho, v_drive(0d0), 0d0, phi, ex)
    ex_old = ex
    if (dump_initial) call dump_initial_state()
    ! leapfrog start: v(-dt/2) = v(0) - (q/m) E_p dt/2
    do isp = 1, nsp
      do i = 1, np(isp)
        xd = x(i,isp)/dx
        j = min(int(xd), nx - 1)
        s = xd - j
        ep = ex(j)*(1d0 - s) + ex(j+1)*s
        vx(i,isp) = vx(i,isp) - 0.5d0*qm(isp)*ep*dt
      end do
    end do
  end subroutine init_fresh
  !
end module init
