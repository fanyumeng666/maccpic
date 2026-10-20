!#############
FUNCTION dran2(idum)
!#############
  !
  ! Uniform random number in (0,1): L'Ecuyer generator with Bays-Durham
  ! shuffle (ran2 of Numerical Recipes). Call with idum < 0 to initialize.
  !
  implicit real(kind=8)(a-h,o-z),integer(kind=4)(i-n)
  INTEGER(kind=4) :: idum
  real(kind=8) :: dran2
  integer(kind=4), parameter :: IM1=2147483563,IM2=2147483399,IMM1=IM1-1
  integer(kind=4), parameter :: IA1=40014,IA2=40692,IQ1=53668,IQ2=52774
  integer(kind=4), parameter :: IR1=12211,IR2=3791
  integer(kind=4), parameter :: NTAB=32,NDIV=1+IMM1/NTAB
  real(kind=8), parameter :: AM=1.d0/IM1,EPS=1.2d-7,RNMX=1.d0-EPS
  INTEGER(kind=4) :: idum2,j,k,iv(NTAB),iy
  SAVE :: iv,iy,idum2
  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
  if (idum <= 0) then
    idum=max0(-idum,1)
    idum2=idum
    do j=NTAB+8,1,-1
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum < 0) idum=idum+IM1
      if (j <= NTAB) iv(j)=idum
    enddo
    iy=iv(1)
  endif
  k=idum/IQ1
  idum=IA1*(idum-k*IQ1)-k*IR1
  if (idum < 0) idum=idum+IM1
  k=idum2/IQ2
  idum2=IA2*(idum2-k*IQ2)-k*IR2
  if (idum2 < 0) idum2=idum2+IM2
  j=1+iy/NDIV
  iy=iv(j)-idum2
  iv(j)=idum
  if(iy < 1)iy=iy+IMM1
  dran2=min(AM*iy,RNMX)
  return
END function dran2
