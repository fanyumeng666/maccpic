!=====================================================================
!  deposit: linear (CIC) weighting of the particles of this rank onto
!  the nodes j = 0..nx.  With j = floor(x/dx) and s = x/dx - j a
!  particle contributes (1-s)/dx to node j and s/dx to node j+1, so the
!  node values are densities (per m^3 per unit weight).  The electrode
!  nodes 0 and nx only see half a cell; the profile output corrects
!  that with a factor 2, the Poisson solver does not use them.
!  deposit_rho and deposit_number_density return sums over all ranks,
!  deposit_species_moments returns the rank-local sums (reduced once at
!  the end of the measurement window).
!=====================================================================
module deposit
  use global
  implicit none
contains
  !
  ! charge density [C/m^3] of all species, all ranks
  subroutine deposit_rho()
    integer :: isp, j
    integer(8) :: i
    real(8) :: s, xd, qw
    rho = 0d0
    do isp = 1, nsp
      qw = q(isp)*weight/dx
      do i = 1, np(isp)
        xd = x(i,isp)/dx
        j = min(int(xd), nx - 1)
        s = xd - j
        rho(j)   = rho(j)   + qw*(1d0 - s)
        rho(j+1) = rho(j+1) + qw*s
      end do
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE, rho, nx + 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  end subroutine deposit_rho
  !
  ! number density [m^-3] of one species, all ranks
  subroutine deposit_number_density(isp, dens)
    integer, intent(in)  :: isp
    real(8), intent(out) :: dens(0:nx)
    integer :: j
    integer(8) :: i
    real(8) :: s, xd, w
    dens = 0d0
    w = weight/dx
    do i = 1, np(isp)
      xd = x(i,isp)/dx
      j = min(int(xd), nx - 1)
      s = xd - j
      dens(j)   = dens(j)   + w*(1d0 - s)
      dens(j+1) = dens(j+1) + w*s
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE, dens, nx + 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  end subroutine deposit_number_density
  !
  ! rank-local moments of one species per unit weight: number density
  ! (1/m^3), number flux density n*vx (1/m^2/s), kinetic energy density
  ! (J/m^3) and its x-directed part.  Multiply by weight for physical values.
  subroutine deposit_species_moments(isp, dens, flux, ekd, ekxd)
    integer, intent(in)  :: isp
    real(8), intent(out) :: dens(0:nx), flux(0:nx), ekd(0:nx), ekxd(0:nx)
    integer :: j
    integer(8) :: i
    real(8) :: s, xd, wl, wr, v, ek, ekx, hm
    dens = 0d0; flux = 0d0; ekd = 0d0; ekxd = 0d0
    hm = 0.5d0*m(isp)
    do i = 1, np(isp)
      xd = x(i,isp)/dx
      j = min(int(xd), nx - 1)
      s = xd - j
      wl = (1d0 - s)/dx
      wr = s/dx
      v = vx(i,isp)
      ekx = hm*v*v
      ek = ekx + hm*(vy(i,isp)**2 + vz(i,isp)**2)
      dens(j)   = dens(j)   + wl
      dens(j+1) = dens(j+1) + wr
      flux(j)   = flux(j)   + wl*v
      flux(j+1) = flux(j+1) + wr*v
      ekd(j)    = ekd(j)    + wl*ek
      ekd(j+1)  = ekd(j+1)  + wr*ek
      ekxd(j)   = ekxd(j)   + wl*ekx
      ekxd(j+1) = ekxd(j+1) + wr*ekx
    end do
  end subroutine deposit_species_moments
  !
end module deposit
