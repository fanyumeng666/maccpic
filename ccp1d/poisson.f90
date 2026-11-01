!=====================================================================
!  poisson: one-dimensional Poisson equation on the nodes j = 0..nx,
!      (phi(j-1) - 2 phi(j) + phi(j+1))/dx^2 = -rho(j)/eps0,  j = 1..nx-1
!  with Dirichlet values phi(0) = v_left, phi(nx) = v_right, solved by
!  the Thomas algorithm.  The forward-elimination coefficients are
!  constant and computed once in poisson_init.  The field is
!  E = -dphi/dx: central differences inside, second-order one-sided
!  differences at the two electrodes.  Independent of the other
!  modules (used by tests/poisson).
!=====================================================================
module poisson
  implicit none
  private
  public :: poisson_init, poisson_solve, efield_faces
  real(8), parameter :: eps0 = 8.8541878128d-12
  integer :: nxp = 0
  real(8) :: dxp = 0d0
  real(8), allocatable :: cp(:), rb(:)   ! modified upper coefficients, 1/(b - a*cp(j-1))
contains
  !
  subroutine poisson_init(nx, dx)
    integer, intent(in) :: nx
    real(8), intent(in) :: dx
    integer :: j
    nxp = nx; dxp = dx
    if (allocated(cp)) deallocate(cp, rb)
    allocate(cp(1:nx-1), rb(1:nx-1))
    ! tridiagonal system for j = 1..nx-1: a = 1, b = -2, c = 1
    rb(1) = 1d0/(-2d0)
    cp(1) = 1d0*rb(1)
    do j = 2, nx - 1
      rb(j) = 1d0/(-2d0 - cp(j-1))
      cp(j) = rb(j)
    end do
  end subroutine poisson_init
  !
  subroutine poisson_solve(rho, v_left, v_right, phi, ex)
    real(8), intent(in)  :: rho(0:nxp), v_left, v_right
    real(8), intent(out) :: phi(0:nxp), ex(0:nxp)
    real(8) :: d(1:nxp-1), f
    integer :: j, nx
    nx = nxp
    f = -dxp*dxp/eps0
    do j = 1, nx - 1
      d(j) = f*rho(j)
    end do
    d(1) = d(1) - v_left
    d(nx-1) = d(nx-1) - v_right
    ! forward sweep
    d(1) = d(1)*rb(1)
    do j = 2, nx - 1
      d(j) = (d(j) - d(j-1))*rb(j)
    end do
    ! back substitution
    phi(0) = v_left
    phi(nx) = v_right
    phi(nx-1) = d(nx-1)
    do j = nx - 2, 1, -1
      phi(j) = d(j) - cp(j)*phi(j+1)
    end do
    ! E = -dphi/dx
    do j = 1, nx - 1
      ex(j) = -(phi(j+1) - phi(j-1))/(2d0*dxp)
    end do
    ex(0)  = -(-3d0*phi(0) + 4d0*phi(1) - phi(2))/(2d0*dxp)
    ex(nx) = -(3d0*phi(nx) - 4d0*phi(nx-1) + phi(nx-2))/(2d0*dxp)
  end subroutine poisson_solve
  !
  ! field on the cell faces, E(j+1/2) = -(phi(j+1) - phi(j))/dx, j = 0..nx-1:
  ! the quantity for which the discrete Gauss law and energy are exact
  subroutine efield_faces(phi, ef)
    real(8), intent(in)  :: phi(0:nxp)
    real(8), intent(out) :: ef(0:nxp-1)
    integer :: j
    do j = 0, nxp - 1
      ef(j) = -(phi(j+1) - phi(j))/dxp
    end do
  end subroutine efield_faces
  !
end module poisson
