!###############
subroutine setv
!###############
  !
  ! subroutine setv converts particle velocities at t=0 to computer
  ! normalization at t=-dt/2.
  ! (one Boris push of -dt/2 with the initial E and B, so that the
  ! velocities start half a step behind the positions; linear weighting)
  !
  use global, only : com_den,nx,ny,nx1,ny1,nx2,ny2,nsp,maxnp,dt,dx,dy,c,&
                     bx,by,bz,ex,ey,ez,bxm,bym,bzm,time,np,tx,ty,am,&
                     x,y,vx,vy,vz,gamma,NPROCS,MYRANK,nthreads
  !$ use omp_lib
  implicit none
  ! ----- MPI -----
  include 'mpif.h'
  ! ---------------
  integer(8) j
  integer(4) isp,ix,iy,ix1,iy1,ix1p,iy1p,ixh,iyh,ixhp,iyhp
  real(8) ae,drx,dlx,duy,dly,drxh,dlxh,duyh,dlyh,sh11,sh12,sh21,sh22
  real(8) s1h1,s1h2,s2h1,s2h2,sh1h1,sh1h2,sh2h1,sh2h2
  real(8) aexpt,aeypt,aezpt,abxpt,abypt,abzpt,gvxs,gvys,gvzs,dcg,vvx,vvy,vvz,f
  !
  real(8), dimension(nsp) :: dn,q,qm
  !
  !     local symbols:
  !
  real(8), dimension(0:nx2,0:ny2) :: aex,aey,aez,abx,aby,abz
  !
  call com_den(dn,q,qm)
  !
  do isp=1,nsp
    !
    !  ..... set various constants.
    !
    ae=-0.25d0*qm(isp)*dt/c
    !
    ! first time step duties:
    ! renormalize the fields and zero arrays.
    !
    !$ call omp_set_dynamic(.false.)
    !$ call omp_set_num_threads(nthreads)
    !$OMP PARALLEL PRIVATE(ix,iy) &
    !$OMP SHARED(ae,ex,ey,ez,bx,by,bz,aex,aey,aez,abx,aby,abz,tx,ty,isp)
    !$OMP DO
    do ix=0,nx2
      do iy=0,ny2
        aex(ix,iy)=ae*ex(ix,iy)
        aey(ix,iy)=ae*ey(ix,iy)
        aez(ix,iy)=ae*ez(ix,iy)
        ! B at t = 0 used directly, as in move (bxm is not set yet here)
        !        abx(ix,iy)=ae*0.5d0*(bx(ix,iy)+bxm(ix,iy))-0.5d0*tx(isp)
        !        aby(ix,iy)=ae*0.5d0*(by(ix,iy)+bym(ix,iy))-0.5d0*ty(isp)
        !        abz(ix,iy)=ae*0.5d0*(bz(ix,iy)+bzm(ix,iy))
        abx(ix,iy)=ae*bx(ix,iy)-0.5d0*tx(isp)
        aby(ix,iy)=ae*by(ix,iy)-0.5d0*ty(isp)
        abz(ix,iy)=ae*bz(ix,iy)
      enddo
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL
    !
    !$OMP PARALLEL DEFAULT(PRIVATE) &
    !$OMP SHARED(aex,aey,aez,abx,aby,abz,isp,np,vx,vy,vz,x,y,gamma)
    !$OMP DO
    do j=MYRANK*maxnp+1,MYRANK*maxnp+np(isp)
      !
      ! interpolate fields to particle position.
      !
      ix1 =x(j,isp)+1;         iy1 =y(j,isp)+1
      drx =x(j,isp)+1.d0-ix1;  duy =y(j,isp)+1.d0-iy1
      dlx =1.d0-drx;           dly =1.d0-duy
      ixh =x(j,isp)+0.5d0;     iyh =y(j,isp)+0.5d0
      drxh=x(j,isp)+0.5d0-ixh; duyh=y(j,isp)+0.5d0-iyh
      dlxh=1.d0-drxh;          dlyh=1.d0-duyh
      !   s11=dlx*dly; s12=dlx*duy; s21=drx*dly; s22=drx*duy
      sh11 =dlxh*dly;  sh12 =dlxh*duy;  sh21 =drxh*dly;  sh22 =drxh*duy
      s1h1 =dlx*dlyh;  s1h2 =dlx*duyh;  s2h1 =drx*dlyh;  s2h2 =drx*duyh
      sh1h1=dlxh*dlyh; sh1h2=dlxh*duyh; sh2h1=drxh*dlyh; sh2h2=drxh*duyh
      ix1p=ix1+1; iy1p=iy1+1
      ixhp=ixh+1; iyhp=iyh+1
      aexpt=aex(ixh ,iy1 )*sh11  + aex(ixhp,iy1 )*sh21  +&
            aex(ixhp,iy1p)*sh22  + aex(ixh ,iy1p)*sh12
      aeypt=aey(ix1 ,iyh )*s1h1  + aey(ix1p,iyh )*s2h1  +&
            aey(ix1p,iyhp)*s2h2  + aey(ix1 ,iyhp)*s1h2
      aezpt=aez(ixh ,iyh )*sh1h1 + aez(ixhp,iyh )*sh2h1 +&
            aez(ixhp,iyhp)*sh2h2 + aez(ixh ,iyhp)*sh1h2
      abxpt=abx(ixh ,iy1 )*sh11  + abx(ixhp,iy1 )*sh21  +&
            abx(ixhp,iy1p)*sh22  + abx(ixh ,iy1p)*sh12
      abypt=aby(ix1 ,iyh )*s1h1  + aby(ix1p,iyh )*s2h1  +&
            aby(ix1p,iyhp)*s2h2  + aby(ix1 ,iyhp)*s1h2
      abzpt=abz(ixh ,iyh )*sh1h1 + abz(ixhp,iyh )*sh2h1 +&
            abz(ixhp,iyhp)*sh2h2 + abz(ixh ,iyhp)*sh1h2
      gamma(j,isp)=1.d0/sqrt(1.d0-vx(j,isp)*vx(j,isp)-vy(j,isp)*vy(j,isp)-vz(j,isp)*vz(j,isp))
      !
      ! half-acceleration along e.
      !
      gvxs=vx(j,isp)*gamma(j,isp)+aexpt
      gvys=vy(j,isp)*gamma(j,isp)+aeypt
      gvzs=vz(j,isp)*gamma(j,isp)+aezpt
      !
      ! first half-rotation, using v x b / gamma.
      !
      gamma(j,isp)=sqrt(1.d0+gvxs*gvxs+gvys*gvys+gvzs*gvzs)
      dcg=1.d0 / (c*gamma(j,isp))
      abxpt=dcg*abxpt
      abypt=dcg*abypt
      abzpt=dcg*abzpt
      vvx=gvxs+gvys*abzpt-gvzs*abypt
      vvy=gvys+gvzs*abxpt-gvxs*abzpt
      vvz=gvzs+gvxs*abypt-gvys*abxpt
      !
      ! second half-rotation and second half-acceleration.
      !
      f=2.d0/(1.d0+abxpt*abxpt+abypt*abypt+abzpt*abzpt)
      abxpt=f*abxpt
      abypt=f*abypt
      abzpt=f*abzpt
      gvxs=gvxs+vvy*abzpt-vvz*abypt+aexpt
      gvys=gvys+vvz*abxpt-vvx*abzpt+aeypt
      gvzs=gvzs+vvx*abypt-vvy*abxpt+aezpt
      !
      ! use v=p / (mc*gamma).
      !
      gamma(j,isp)=sqrt(1.d0+gvxs*gvxs+gvys*gvys+gvzs*gvzs)
      vx(j,isp)=gvxs/gamma(j,isp)
      vy(j,isp)=gvys/gamma(j,isp)
      vz(j,isp)=gvzs/gamma(j,isp)
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL
    !
  enddo
  !
  return
end subroutine setv
