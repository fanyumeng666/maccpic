!=====================================================================
!  random: L'Ecuyer generator with Bays-Durham shuffle (ran2 of
!  Numerical Recipes, the algorithm of em2's dran2.f90) with the state
!  in module variables so that it can be saved in the restart files,
!  plus a Box-Muller normal deviate.  Not thread safe (MPI only).
!=====================================================================
module random
  implicit none
  private
  public :: rng_init, rand_u, rand_gauss, rng_get_state, rng_set_state, rng_state_size
  integer, parameter :: IM1 = 2147483563, IM2 = 2147483399, IMM1 = IM1 - 1
  integer, parameter :: IA1 = 40014, IA2 = 40692, IQ1 = 53668, IQ2 = 52774
  integer, parameter :: IR1 = 12211, IR2 = 3791
  integer, parameter :: NTAB = 32, NDIV = 1 + IMM1/NTAB
  real(8), parameter :: AM = 1d0/IM1, EPS = 1.2d-7, RNMX = 1d0 - EPS
  real(8), parameter :: twopi = 6.283185307179586d0
  integer, parameter :: rng_state_size = NTAB + 3
  integer, save :: idum = 0, idum2 = 123456789, iy = 0, iv(NTAB) = 0
contains
  !
  subroutine rng_init(seed)
    integer, intent(in) :: seed
    integer :: j, k
    idum = max(seed, 1)
    idum2 = idum
    do j = NTAB + 8, NTAB + 1, -1   ! 8 warm-up draws
      k = idum/IQ1
      idum = IA1*(idum - k*IQ1) - k*IR1
      if (idum < 0) idum = idum + IM1
    end do
    do j = NTAB, 1, -1              ! fill the shuffle table
      k = idum/IQ1
      idum = IA1*(idum - k*IQ1) - k*IR1
      if (idum < 0) idum = idum + IM1
      iv(j) = idum
    end do
    iy = iv(1)
  end subroutine rng_init
  !
  ! uniform deviate in (0, 1)
  function rand_u() result(r)
    real(8) :: r
    integer :: j, k
    k = idum/IQ1
    idum = IA1*(idum - k*IQ1) - k*IR1
    if (idum < 0) idum = idum + IM1
    k = idum2/IQ2
    idum2 = IA2*(idum2 - k*IQ2) - k*IR2
    if (idum2 < 0) idum2 = idum2 + IM2
    j = 1 + iy/NDIV
    iy = iv(j) - idum2
    iv(j) = idum
    if (iy < 1) iy = iy + IMM1
    r = min(AM*iy, RNMX)
  end function rand_u
  !
  ! standard normal deviate (Box-Muller, no caching so that the state is complete)
  function rand_gauss() result(g)
    real(8) :: g, u1, u2
    u1 = rand_u()
    u2 = rand_u()
    g = sqrt(-2d0*log(u1))*cos(twopi*u2)
  end function rand_gauss
  !
  subroutine rng_get_state(s)
    integer, intent(out) :: s(rng_state_size)
    s(1) = idum; s(2) = idum2; s(3) = iy; s(4:NTAB+3) = iv
  end subroutine rng_get_state
  !
  subroutine rng_set_state(s)
    integer, intent(in) :: s(rng_state_size)
    idum = s(1); idum2 = s(2); iy = s(3); iv = s(4:NTAB+3)
  end subroutine rng_set_state
  !
end module random
