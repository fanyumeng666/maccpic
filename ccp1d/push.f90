!=====================================================================
!  push: leapfrog advance of one species with the field of the current
!  step.  On entry x is at t_n and vx at t_{n-1/2}; on exit x is at
!  t_{n+1} and vx at t_{n+1/2}.  Particles that leave [0, L) are
!  absorbed by the electrode: their energy and species go into a buffer
!  that diag bins into the wall histograms, and the slot is reused.
!  The routine also sums, for this step and rank,
!    work_step(isp)  = sum q E_p dt (v- + v+)/2 = exact kinetic energy
!                      gained from the field (per unit weight)
!    ekin_step(isp)  = sum 1/2 m v^2 of the particles still inside,
!                      with vx at t_{n+1/2}
!    nwall_step, ewall_step = absorbed count and kinetic energy per
!                      electrode (1 = x < 0, 2 = x >= L)
!    vxsum_step(isp) = sum of vx(t_{n+1/2}) over all pushed particles:
!                      by the Shockley-Ramo theorem the current in the
!                      external circuit is (weight/L) sum q vx.
!  The exact bookkeeping identity is then
!    ekin(after step) - ekin(before step) = work_step - sum ewall_step.
!=====================================================================
module push
  use global
  implicit none
  real(8)    :: work_step(nsp) = 0d0, ekin_step(nsp) = 0d0, vxsum_step(nsp) = 0d0
  integer(8) :: nwall_step(2,nsp) = 0
  real(8)    :: ewall_step(2,nsp) = 0d0
  ! work of the field on the particles deposited on the nodes with the CIC
  ! weights of the old position (per unit weight, J per step): integrates
  ! exactly to work_step, unlike the grid product J.E (diag uses it in the
  ! measurement window: tally_work = in_window)
  logical    :: tally_work = .false.
  real(8), allocatable :: work_dep(:,:)
  ! absorbed-particle buffer of this step
  integer :: nabs = 0
  integer, allocatable :: abs_side(:), abs_isp(:)
  real(8), allocatable :: abs_energy(:)
  real(8), allocatable, private :: a(:)     ! (q/m) E dt on the nodes
contains
  !
  subroutine push_init()
    allocate(a(0:nx), work_dep(0:nx,nsp))
    work_dep = 0d0
    allocate(abs_side(4096), abs_isp(4096), abs_energy(4096))
    nabs = 0
  end subroutine push_init
  !
  subroutine push_species(isp)
    integer, intent(in) :: isp
    integer(8) :: i
    integer :: j, side
    real(8) :: xd, s, ep, vold, vnew, xi, f, wsum, esum, vsum, e, hm, wp
    f = qm(isp)*dt
    a = f*ex
    hm = 0.5d0*m(isp)
    wsum = 0d0; esum = 0d0; vsum = 0d0
    nwall_step(:,isp) = 0; ewall_step(:,isp) = 0d0
    if (tally_work) work_dep(:,isp) = 0d0
    i = 1
    do while (i <= np(isp))
      xd = x(i,isp)/dx
      j = min(int(xd), nx - 1)
      s = xd - j
      ep = a(j)*(1d0 - s) + a(j+1)*s          ! (q/m) E_p dt
      vold = vx(i,isp)
      vnew = vold + ep
      vx(i,isp) = vnew
      xi = x(i,isp) + vnew*dt
      x(i,isp) = xi
      wsum = wsum + ep*(vold + vnew)           ! 2/m * work
      vsum = vsum + vnew
      if (tally_work) then
        wp = hm*ep*(vold + vnew)               ! q E_p dt (v- + v+)/2
        work_dep(j,isp)   = work_dep(j,isp)   + wp*(1d0 - s)
        work_dep(j+1,isp) = work_dep(j+1,isp) + wp*s
      end if
      if (xi < 0d0 .or. xi >= L) then
        side = 2
        if (xi < 0d0) side = 1
        e = hm*(vnew*vnew + vy(i,isp)**2 + vz(i,isp)**2)
        nwall_step(side,isp) = nwall_step(side,isp) + 1
        ewall_step(side,isp) = ewall_step(side,isp) + e
        call record_absorbed(side, isp, e)
        call remove_particle(isp, i)           ! the last particle now sits in slot i
      else
        esum = esum + vnew*vnew + vy(i,isp)**2 + vz(i,isp)**2
        i = i + 1
      end if
    end do
    work_step(isp) = hm*wsum
    ekin_step(isp) = hm*esum
    vxsum_step(isp) = vsum
  end subroutine push_species
  !
  subroutine record_absorbed(side, isp, e)
    integer, intent(in) :: side, isp
    real(8), intent(in) :: e
    integer, allocatable :: ti(:)
    real(8), allocatable :: tr(:)
    integer :: n
    if (nabs >= size(abs_side)) then
      n = 2*size(abs_side)
      allocate(ti(n)); ti(1:nabs) = abs_side(1:nabs); call move_alloc(ti, abs_side)
      allocate(ti(n)); ti(1:nabs) = abs_isp(1:nabs);  call move_alloc(ti, abs_isp)
      allocate(tr(n)); tr(1:nabs) = abs_energy(1:nabs); call move_alloc(tr, abs_energy)
    end if
    nabs = nabs + 1
    abs_side(nabs) = side; abs_isp(nabs) = isp; abs_energy(nabs) = e
  end subroutine record_absorbed
  !
  ! x-directed and perpendicular kinetic energy sums of this rank (per unit weight) [J]
  subroutine kinetic_sums(isp, ekx, ekperp)
    integer, intent(in) :: isp
    real(8), intent(out) :: ekx, ekperp
    integer(8) :: i
    ekx = 0d0; ekperp = 0d0
    do i = 1, np(isp)
      ekx = ekx + vx(i,isp)**2
      ekperp = ekperp + vy(i,isp)**2 + vz(i,isp)**2
    end do
    ekx = 0.5d0*m(isp)*ekx
    ekperp = 0.5d0*m(isp)*ekperp
  end subroutine kinetic_sums
  !
  ! kinetic energy of the particles of this rank (per unit weight) [J]
  function kinetic_energy(isp) result(ek)
    integer, intent(in) :: isp
    real(8) :: ek
    integer(8) :: i
    ek = 0d0
    do i = 1, np(isp)
      ek = ek + vx(i,isp)**2 + vy(i,isp)**2 + vz(i,isp)**2
    end do
    ek = 0.5d0*m(isp)*ek
  end function kinetic_energy
  !
end module push
