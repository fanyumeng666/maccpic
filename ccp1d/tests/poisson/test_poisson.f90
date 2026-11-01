! T1: Poisson solver against analytic solutions, second-order convergence,
!     the discrete residual and the discrete Gauss law.
program test_poisson
  use poisson
  implicit none
  real(8), parameter :: eps0 = 8.8541878128d-12, pi = 3.141592653589793d0
  real(8), parameter :: L = 0.067d0
  integer :: nx, j, k, nfail
  real(8) :: dx, v0, rho0, err, errs(3), ratio, g1, g2, res, maxres
  real(8), allocatable :: rho(:), phi(:), ex(:), ef(:), exact(:)
  nfail = 0
  !
  ! (a) rho = 0, phi(0) = V0: linear potential, constant field -V0/L
  nx = 128; dx = L/nx; v0 = 450d0
  allocate(rho(0:nx), phi(0:nx), ex(0:nx), ef(0:nx-1), exact(0:nx))
  call poisson_init(nx, dx)
  rho = 0d0
  call poisson_solve(rho, v0, 0d0, phi, ex)
  err = 0d0
  do j = 0, nx
    err = max(err, abs(phi(j) - v0*(1d0 - j*dx/L))/v0, abs(ex(j) - v0/L)/(v0/L))   ! E = -dphi/dx = +V0/L
  end do
  call report('T1a linear potential, rel. error', err, 1d-10)
  !
  ! (b) rho = const, grounded: phi = rho0 x (L - x)/(2 eps0)
  rho0 = 1d-6
  rho = rho0
  call poisson_solve(rho, 0d0, 0d0, phi, ex)
  do j = 0, nx
    exact(j) = rho0*(j*dx)*(L - j*dx)/(2d0*eps0)
  end do
  err = 0d0
  do j = 0, nx
    err = max(err, abs(phi(j) - exact(j))/maxval(exact))
    err = max(err, abs(ex(j) - rho0*(j*dx - 0.5d0*L)/eps0)/(rho0*0.5d0*L/eps0))
  end do
  call report('T1b parabola (phi and E), rel. error', err, 1d-10)
  !
  ! (c) rho = sin(pi x/L): error ratio 4 when nx doubles (second order)
  do k = 1, 3
    nx = 32*2**(k-1); dx = L/nx
    deallocate(rho, phi, ex, ef, exact)
    allocate(rho(0:nx), phi(0:nx), ex(0:nx), ef(0:nx-1), exact(0:nx))
    call poisson_init(nx, dx)
    do j = 0, nx
      rho(j) = rho0*sin(pi*j*dx/L)
      exact(j) = rho0*(L/pi)**2*sin(pi*j*dx/L)/eps0
    end do
    call poisson_solve(rho, 0d0, 0d0, phi, ex)
    errs(k) = maxval(abs(phi - exact))/maxval(exact)
  end do
  ratio = errs(2)/errs(3)
  write(*,'(a,3es12.4)') '    T1c errors nx = 32, 64, 128:', errs
  call report('T1c convergence ratio err(64)/err(128) - 4', abs(ratio - 4d0), 0.1d0)
  !
  ! (d) random rho: residual of the discrete equation, (e) discrete Gauss law
  nx = 128; dx = L/nx
  deallocate(rho, phi, ex, ef, exact)
  allocate(rho(0:nx), phi(0:nx), ex(0:nx), ef(0:nx-1), exact(0:nx))
  call poisson_init(nx, dx)
  call random_seed()
  call random_number(rho)
  rho = (rho - 0.5d0)*1d-5
  call poisson_solve(rho, 123d0, 0d0, phi, ex)
  maxres = 0d0
  do j = 1, nx - 1
    res = abs((phi(j-1) - 2d0*phi(j) + phi(j+1))/dx**2 + rho(j)/eps0)
    maxres = max(maxres, res/(maxval(abs(rho))/eps0))
  end do
  call report('T1d discrete residual (relative)', maxres, 1d-12)
  call efield_faces(phi, ef)
  g1 = eps0*(ef(nx-1) - ef(0))
  g2 = dx*sum(rho(1:nx-1))
  ! normalise with sum |rho| dx: the signed sum of a random rho nearly cancels
  call report('T1e discrete Gauss law eps0 dE_face = int rho (rel. to int |rho|)', abs(g1 - g2)/(dx*sum(abs(rho(1:nx-1)))), 1d-13)
  !
  if (nfail == 0) then
    write(*,'(a)') 'test_poisson: PASS'
  else
    write(*,'(a,i0,a)') 'test_poisson: FAIL (', nfail, ' checks)'
    stop 1
  end if
contains
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
end program test_poisson
