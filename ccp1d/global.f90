!=====================================================================
!  global: physical constants, species, MPI rank variables, the grid
!  arrays (nodes 0..nx), the particle arrays of this rank and the step
!  counters.  Particle slot i of species isp holds x(i,isp) [m] and
!  vx, vy, vz(i,isp) [m/s]; slots 1..np(isp) are in use.
!=====================================================================
module global
  use input
  use mpi
  implicit none
  !
  ! ---- constants (CODATA 2018; M_He+ as in Turner et al. 2013) ----
  real(8), parameter :: pi    = 3.141592653589793d0, twopi = 2d0*pi
  real(8), parameter :: qe    = 1.602176634d-19      ! elementary charge [C]
  real(8), parameter :: me    = 9.1093837015d-31     ! electron mass [kg]
  real(8), parameter :: mion  = 6.67d-27             ! He+ mass [kg]
  real(8), parameter :: eps0  = 8.8541878128d-12     ! vacuum permittivity [F/m]
  real(8), parameter :: kb    = 1.380649d-23         ! Boltzmann constant [J/K]
  !
  ! ---- species: 1 = electrons, 2 = He+ ----
  integer, parameter :: nsp = 2
  real(8), parameter :: q(nsp) = (/ -qe, qe /)       ! charge [C]
  real(8), parameter :: m(nsp) = (/ me, mion /)      ! mass [kg]
  real(8), parameter :: qm(nsp) = q/m                ! charge-to-mass ratio [C/kg]
  character(len=4), parameter :: spname(nsp) = (/ 'e   ', 'He+ ' /)
  !
  ! ---- MPI ----
  integer :: myrank = 0, nprocs = 1, ierr = 0
  !
  ! ---- grid (nodes j = 0..nx at x_j = j*dx) ----
  real(8), allocatable :: phi(:)      ! potential [V]
  real(8), allocatable :: ex(:)       ! electric field [V/m]
  real(8), allocatable :: ex_old(:)   ! E of the previous step (for J.E and the displacement current)
  real(8), allocatable :: rho(:)      ! charge density [C/m^3], summed over ranks
  !
  ! ---- particles of this rank ----
  real(8), allocatable :: x(:,:), vx(:,:), vy(:,:), vz(:,:)
  integer(8) :: np(nsp) = 0
  !
  ! ---- time ----
  integer(8) :: it = 0                ! steps done
  integer    :: iperiod = 0           ! periods done
  real(8)    :: time = 0d0            ! [s]
  logical    :: mcc_on = .false.      ! collisions active (n_gas > 0 and tables loaded)
  !
contains
  !
  subroutine allocate_arrays()
    allocate(phi(0:nx), ex(0:nx), ex_old(0:nx), rho(0:nx))
    phi = 0d0; ex = 0d0; ex_old = 0d0; rho = 0d0
    allocate(x(maxnp,nsp), vx(maxnp,nsp), vy(maxnp,nsp), vz(maxnp,nsp))
    x = 0d0; vx = 0d0; vy = 0d0; vz = 0d0
    np = 0
  end subroutine allocate_arrays
  !
  ! append a particle; abort when the slots are exhausted
  subroutine add_particle(isp, xp, vxp, vyp, vzp)
    integer, intent(in) :: isp
    real(8), intent(in) :: xp, vxp, vyp, vzp
    if (np(isp) >= maxnp) call abort_run('np > maxnp for species '//trim(spname(isp))// &
         ': increase maxnp_factor')
    np(isp) = np(isp) + 1
    x(np(isp),isp) = xp; vx(np(isp),isp) = vxp; vy(np(isp),isp) = vyp; vz(np(isp),isp) = vzp
  end subroutine add_particle
  !
  ! delete slot i by moving the last particle into it (order is irrelevant)
  subroutine remove_particle(isp, i)
    integer, intent(in) :: isp
    integer(8), intent(in) :: i
    integer(8) :: n
    n = np(isp)
    if (i /= n) then
      x(i,isp) = x(n,isp); vx(i,isp) = vx(n,isp); vy(i,isp) = vy(n,isp); vz(i,isp) = vz(n,isp)
    end if
    np(isp) = n - 1
  end subroutine remove_particle
  !
  ! fatal error on any rank: print and bring the whole job down
  subroutine abort_run(msg)
    character(len=*), intent(in) :: msg
    write(*,'(a,i0,a,a)') 'ccp1d ABORT (rank ', myrank, '): ', trim(msg)
    flush(6)
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    stop 1
  end subroutine abort_run
  !
end module global
