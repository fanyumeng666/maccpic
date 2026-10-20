!###############
subroutine smear(runtimes)
!###############
  !
  ! subroutine smear randomizes vector arrays by random pair exchange.
  ! It swaps vx, vy, vz (independently) among slots 1..np/2 of every
  ! species, runtimes times; the position swaps are commented out.

  use global, only : nsp,maxnp,x,y,vx,vy,vz,np,MYRANK
  implicit none
  ! ----- MPI -----
  include 'mpif.h'
  ! ---------------
  integer(8) i,ii,nph,mynp
  integer(4) isp,idum,rt,runtimes
  real(8) xx,vxx,vyx,vzx,dran2
  !
  do rt=1,runtimes
    do isp=1,nsp
      idum=-7146+MYRANK
      nph=np(isp)/2
      mynp=MYRANK*maxnp
      do i=1,nph
        !     ii=nph*dran2(idum)+1
        !     xx         = x(mynp+i ,isp)
        !     x(mynp+i ,isp)  = x(mynp+ii,isp)
        !     x(mynp+ii,isp)  = xx
        !!     ii=nph*dran2(idum)+1
        !     xx         = y(mynp+i ,isp)
        !     y(mynp+i ,isp)  = y(mynp+ii,isp)
        !     y(mynp+ii,isp)  = xx
        ii=nph*dran2(idum)+1
        vxx        = vx(mynp+i ,isp)
        vx(mynp+i ,isp) = vx(mynp+ii,isp)
        vx(mynp+ii,isp) = vxx
        ii=nph*dran2(idum)+1
        vyx        = vy(mynp+i ,isp)
        vy(mynp+i ,isp) = vy(mynp+ii,isp)
        vy(mynp+ii,isp) = vyx
        ii=nph*dran2(idum)+1
        vzx        = vz(mynp+i ,isp)
        vz(mynp+i ,isp) = vz(mynp+ii,isp)
        vz(mynp+ii,isp) = vzx
      enddo
    enddo
  enddo
end subroutine smear
