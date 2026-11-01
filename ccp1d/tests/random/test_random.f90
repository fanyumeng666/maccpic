! T3: uniform and normal deviates (moments), state save/restore, seeds
program test_random
  use random
  implicit none
  integer, parameter :: n = 10000000
  integer :: i, nfail, s(rng_state_size)
  real(8) :: u, g, m1, m2, m4, a(1000), b(1000)
  nfail = 0
  call rng_init(12345)
  m1 = 0d0; m2 = 0d0
  do i = 1, n
    u = rand_u()
    m1 = m1 + u; m2 = m2 + u*u
    if (u <= 0d0 .or. u >= 1d0) nfail = nfail + 1
  end do
  m1 = m1/n; m2 = m2/n - m1*m1
  call report('T3a uniform mean - 1/2 (sigma units)', abs(m1 - 0.5d0)/(sqrt(1d0/12d0)/sqrt(dble(n))), 4d0)
  call report('T3a uniform variance - 1/12 (relative)', abs(m2 - 1d0/12d0)*12d0, 2d-3)
  m1 = 0d0; m2 = 0d0; m4 = 0d0
  do i = 1, n
    g = rand_gauss()
    m1 = m1 + g; m2 = m2 + g*g; m4 = m4 + g**4
  end do
  m1 = m1/n; m2 = m2/n; m4 = m4/n
  call report('T3b gauss mean (sigma units)', abs(m1)*sqrt(dble(n)), 4d0)
  call report('T3b gauss variance - 1 (relative)', abs(m2 - 1d0), 2d-3)
  call report('T3b gauss kurtosis - 3 (relative)', abs(m4/m2**2 - 3d0)/3d0, 5d-3)
  ! state save / restore: identical sequence
  call rng_get_state(s)
  do i = 1, 1000
    a(i) = rand_gauss()
  end do
  call rng_set_state(s)
  do i = 1, 1000
    b(i) = rand_gauss()
  end do
  call report('T3c state restore reproduces the sequence (max diff)', maxval(abs(a - b)), 0d0)
  ! different seeds: different sequences
  call rng_init(1)
  do i = 1, 1000
    a(i) = rand_u()
  end do
  call rng_init(2)
  do i = 1, 1000
    b(i) = rand_u()
  end do
  call report('T3d seeds 1 and 2: number of coincidences', dble(count(abs(a - b) <= 0d0)), 0d0)
  if (nfail == 0) then
    write(*,'(a)') 'test_random: PASS'
  else
    write(*,'(a,i0,a)') 'test_random: FAIL (', nfail, ' checks)'
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
end program test_random
