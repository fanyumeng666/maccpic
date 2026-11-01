! Serial correlation of the uniform generator at lags 1..16 and the
! correlation between the Box-Muller normals of consecutive particles
! (7 draws apart) and between vx, vy, vz of the same particle.
program test_corr
  use random
  implicit none
  integer, parameter :: n = 20000000, maxlag = 16
  real(8), parameter :: twopi = 6.283185307179586d0
  real(8), allocatable :: u(:)
  real(8) :: c, m, v, g(3), gp(3), sgg(3), sgx(3), x
  integer :: lag, i, k, gen
  character(len=10) :: gname(2) = (/ 'ran2      ', 'intrinsic ' /)
  allocate(u(n))
  do gen = 1, 2
    call rng_init(12345)
    call random_seed()
    do i = 1, n
      if (gen == 1) then
        u(i) = rand_u()
      else
        call random_number(u(i))
      end if
    end do
    m = sum(u)/n; v = sum((u - m)**2)/n
    write(*,'(a,a,a,f9.6,a,f9.6)') '==== ', trim(gname(gen)), ': mean ', m, '  variance ', v
    write(*,'(a)') '  lag  autocorrelation (se = 2.2e-4)'
    do lag = 1, maxlag
      c = sum((u(1:n-lag) - m)*(u(1+lag:n) - m))/((n - lag)*v)
      write(*,'(i5,f12.5)') lag, c
    end do
    ! normals in the loading pattern: correlation of the same component between
    ! consecutive particles, and of x with vx
    sgg = 0d0; sgx = 0d0; gp = 0d0
    do i = 0, n/7 - 1
      x = u(7*i + 1) - 0.5d0
      do k = 1, 3
        g(k) = sqrt(-2d0*log(u(7*i + 2*k)))*cos(twopi*u(7*i + 2*k + 1))
        sgg(k) = sgg(k) + g(k)*gp(k)
        sgx(k) = sgx(k) + g(k)*x
      end do
      gp = g
    end do
    write(*,'(a,3f10.5,a)') '  corr(v_k of particle i, v_k of particle i-1):', sgg/(n/7), '  (se 6e-4)'
    write(*,'(a,3f10.5,a)') '  corr(x, v_k) of the same particle:          ', sgx/(n/7)/sqrt(1d0/12d0), '  (se 6e-4)'
  end do
end program test_corr
