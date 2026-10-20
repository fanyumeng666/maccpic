!###############
subroutine move
!###############
  !
  !  subroutine move advances the particle velocities by one time step
  !  (relativistic Boris scheme) and sums the kinetic energy ake.
  !  Fields are interpolated with second-order (quadratic spline) weights
  !  on the staggered grid. Positions are advanced in curnt.
  !
  use global, only : com_den,nx,nx1,nx2,ny,ny1,ny2,nsp,maxnp,dt,dx,dy,c,&
                     bx,by,bz,ex,ey,ez,bxm,bym,bzm,time,&
                     x,y,vx,vy,vz,gamma,tx,ty,am,np,ake,IERR,MYRANK,nthreads
  !$ use omp_lib
  implicit none
  ! ----- MPI -----
  include 'mpif.h'
  integer(4) ISTATUS(MPI_STATUS_SIZE)
  real(8), dimension(nsp) :: ake0
  ! ---------------
  integer(8) j
  integer(4) isp,ix,iy!,ix1,iy1,ix1p,iy1p,ixh,iyh,ixhp,iyhp
  integer(4) ixm,ixp,iym,iyp
  real(8) ae!,drx,dlx,duy,dly,drxh,dlxh,duyh,dlyh,sh11,sh12,sh21,sh22
  ! real(8) s1h1,s1h2,s2h1,s2h2,sh1h1,sh1h2,sh2h1,sh2h2
  real(8)    s2x,s2y,s2xp,s2yp,s2xm,s2ym
  real(8)    qxy,qxpy,qxmy,qxyp,qxpyp,qxmyp,qxym,qxpym,qxmym
  real(8) aexpt,aeypt,aezpt,abxpt,abypt,abzpt,gvxs,gvys,gvzs,dcg,vvx,vvy,vvz,f
  !
  real(8), dimension(nsp) :: q,dn,qm
  !
  ! local symbols:
  !
  real(8), dimension(0:nx2,0:ny2) :: aex,aey,aez,abx,aby,abz
  !
  call com_den(dn,q,qm)
  !
  do isp=1,nsp
    ake(isp)=0.d0
    !
    !  ..... set various constants.
    !
    ae=0.5d0*qm(isp)*dt/c
    !
    ! scale E and B by (q/m)(dt/2) once per grid point. After fields, B is
    ! at the same time level as E (two half steps), so B is used directly;
    ! the former average with bxm (B of the previous step) is kept commented.
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
        !        abx(ix,iy)=ae*0.5d0*(bx(ix,iy)+bxm(ix,iy))+tx(isp)
        !        aby(ix,iy)=ae*0.5d0*(by(ix,iy)+bym(ix,iy))+ty(isp)
        !        abz(ix,iy)=ae*0.5d0*(bz(ix,iy)+bzm(ix,iy))
        abx(ix,iy)=ae*bx(ix,iy)+tx(isp)
        aby(ix,iy)=ae*by(ix,iy)+ty(isp)
        abz(ix,iy)=ae*bz(ix,iy)
      enddo
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL

    !$OMP PARALLEL DEFAULT(PRIVATE) &
    !$OMP SHARED(aex,aey,aez,abx,aby,abz,isp,np,vx,vy,vz,x,y,gamma,ake,MYRANK)
    !$OMP DO REDUCTION(+:ake)
    do j=MYRANK*maxnp+1,MYRANK*maxnp+np(isp)
      !
      ! interpolate fields to particle position.
      ! Ex and Bx: x on half-integer nodes, y on integer nodes
      !
      ix = x(j,isp) + 1.0d0;  iy = y(j,isp) + 1.5d0
      ixp= ix + 1;              iyp= iy + 1
      ixm= ix - 1;              iym= iy - 1
      s2x = 0.75d0 - (x(j,isp) + 0.5d0 - ix)**2
      s2xp= 0.5d0 * ( 1.0d0 + x(j,isp) - ix)**2
      s2xm= 0.5d0 * (       - x(j,isp) + ix)**2
      s2y = 0.75d0 - (y(j,isp) + 1.d0 - iy)**2
      s2yp= 0.5d0 * ( 1.5d0 + y(j,isp) - iy)**2
      s2ym= 0.5d0 * (-0.5d0 - y(j,isp) + iy)**2
      !
      qxy   = s2x *s2y
      qxpy  = s2xp*s2y
      qxmy  = s2xm*s2y
      qxyp  = s2x *s2yp
      qxpyp = s2xp*s2yp
      qxmyp = s2xm*s2yp
      qxym  = s2x *s2ym
      qxpym = s2xp*s2ym
      qxmym = s2xm*s2ym
      aexpt = aex(ix ,iy)*qxy  + aex(ix ,iyp)*qxyp  + aex(ix ,iym)*qxym  + &
              aex(ixp,iy)*qxpy + aex(ixp,iyp)*qxpyp + aex(ixp,iym)*qxpym + &
              aex(ixm,iy)*qxmy + aex(ixm,iyp)*qxmyp + aex(ixm,iym)*qxmym
      abxpt = abx(ix ,iy)*qxy  + abx(ix ,iyp)*qxyp  + abx(ix ,iym)*qxym  + &
              abx(ixp,iy)*qxpy + abx(ixp,iyp)*qxpyp + abx(ixp,iym)*qxpym + &
              abx(ixm,iy)*qxmy + abx(ixm,iyp)*qxmyp + abx(ixm,iym)*qxmym
      !
      ! Ez and Bz: half-integer nodes in x and y
      iy = y(j,isp) + 1.0d0
      iyp= iy + 1
      iym= iy - 1
      s2y = 0.75d0 - (y(j,isp) + 0.5d0 - iy)**2
      s2yp= 0.5d0 * ( 1.0d0 + y(j,isp) - iy)**2
      s2ym= 0.5d0 * (       - y(j,isp) + iy)**2
      !
      qxy   = s2x *s2y
      qxpy  = s2xp*s2y
      qxmy  = s2xm*s2y
      qxyp  = s2x *s2yp
      qxpyp = s2xp*s2yp
      qxmyp = s2xm*s2yp
      qxym  = s2x *s2ym
      qxpym = s2xp*s2ym
      qxmym = s2xm*s2ym
      aezpt = aez(ix ,iy)*qxy  + aez(ix ,iyp)*qxyp  + aez(ix ,iym)*qxym  + &
              aez(ixp,iy)*qxpy + aez(ixp,iyp)*qxpyp + aez(ixp,iym)*qxpym + &
              aez(ixm,iy)*qxmy + aez(ixm,iyp)*qxmyp + aez(ixm,iym)*qxmym
      abzpt = abz(ix ,iy)*qxy  + abz(ix ,iyp)*qxyp  + abz(ix ,iym)*qxym  + &
              abz(ixp,iy)*qxpy + abz(ixp,iyp)*qxpyp + abz(ixp,iym)*qxpym + &
              abz(ixm,iy)*qxmy + abz(ixm,iyp)*qxmyp + abz(ixm,iym)*qxmym
      !
      ! Ey and By: integer nodes in x, half-integer nodes in y
      ix = x(j,isp) + 1.5d0
      ixp= ix + 1
      ixm= ix - 1
      s2x = 0.75d0 - (x(j,isp) + 1.0d0 - ix)**2
      s2xp= 0.5d0 * ( 1.5d0 + x(j,isp) - ix)**2
      s2xm= 0.5d0 * (-0.5d0 - x(j,isp) + ix)**2
      !
      qxy   = s2x *s2y
      qxpy  = s2xp*s2y
      qxmy  = s2xm*s2y
      qxyp  = s2x *s2yp
      qxpyp = s2xp*s2yp
      qxmyp = s2xm*s2yp
      qxym  = s2x *s2ym
      qxpym = s2xp*s2ym
      qxmym = s2xm*s2ym
      aeypt = aey(ix ,iy)*qxy  + aey(ix ,iyp)*qxyp  + aey(ix ,iym)*qxym  + &
              aey(ixp,iy)*qxpy + aey(ixp,iyp)*qxpyp + aey(ixp,iym)*qxpym + &
              aey(ixm,iy)*qxmy + aey(ixm,iyp)*qxmyp + aey(ixm,iym)*qxmym
      abypt = aby(ix ,iy)*qxy  + aby(ix ,iyp)*qxyp  + aby(ix ,iym)*qxym  + &
              aby(ixp,iy)*qxpy + aby(ixp,iyp)*qxpyp + aby(ixp,iym)*qxpym + &
              aby(ixm,iy)*qxmy + aby(ixm,iyp)*qxmyp + aby(ixm,iym)*qxmym

      !     ix1 =x(j,isp)+1;         iy1 =y(j,isp)+1
      !     drx =x(j,isp)+1.d0-ix1;  duy =y(j,isp)+1.d0-iy1
      !     dlx =1.d0-drx;           dly =1.d0-duy
      !     ixh =x(j,isp)+0.5d0;     iyh =y(j,isp)+0.5d0
      !     drxh=x(j,isp)+0.5d0-ixh; duyh=y(j,isp)+0.5d0-iyh
      !     dlxh=1.d0-drxh;          dlyh=1.d0-duyh
      !     sh11 =dlxh*dly;  sh12 =dlxh*duy;  sh21 =drxh*dly;  sh22 =drxh*duy
      !     s1h1 =dlx*dlyh;  s1h2 =dlx*duyh;  s2h1 =drx*dlyh;  s2h2 =drx*duyh
      !     sh1h1=dlxh*dlyh; sh1h2=dlxh*duyh; sh2h1=drxh*dlyh; sh2h2=drxh*duyh
      !     ix1p=ix1+1; iy1p=iy1+1
      !     ixhp=ixh+1; iyhp=iyh+1
      !     aexpt=aex(ixh ,iy1 )*sh11  + aex(ixhp,iy1 )*sh21  +&
      !           aex(ixhp,iy1p)*sh22  + aex(ixh ,iy1p)*sh12
      !     aeypt=aey(ix1 ,iyh )*s1h1  + aey(ix1p,iyh )*s2h1  +&
      !           aey(ix1p,iyhp)*s2h2  + aey(ix1 ,iyhp)*s1h2
      !     aezpt=aez(ixh ,iyh )*sh1h1 + aez(ixhp,iyh )*sh2h1 +&
      !           aez(ixhp,iyhp)*sh2h2 + aez(ixh ,iyhp)*sh1h2
      !     abxpt=abx(ixh ,iy1 )*sh11  + abx(ixhp,iy1 )*sh21  +&
      !           abx(ixhp,iy1p)*sh22  + abx(ixh ,iy1p)*sh12
      !     abypt=aby(ix1 ,iyh )*s1h1  + aby(ix1p,iyh )*s2h1  +&
      !           aby(ix1p,iyhp)*s2h2  + aby(ix1 ,iyhp)*s1h2
      !     abzpt=abz(ixh ,iyh )*sh1h1 + abz(ixhp,iyh )*sh2h1 +&
      !           abz(ixhp,iyhp)*sh2h2 + abz(ixh ,iyhp)*sh1h2
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
      ake(isp)=ake(isp)+gamma(j,isp)-1.d0
      dcg=1.d0 / (c*gamma(j,isp))
      abxpt=dcg*abxpt; abypt=dcg*abypt; abzpt=dcg*abzpt
      vvx=gvxs+gvys*abzpt-gvzs*abypt
      vvy=gvys+gvzs*abxpt-gvxs*abzpt
      vvz=gvzs+gvxs*abypt-gvys*abxpt
      !
      ! second half-rotation and second half-acceleration.
      !
      f=2.d0/(1.d0+abxpt*abxpt+abypt*abypt+abzpt*abzpt)
      abxpt=f*abxpt; abypt=f*abypt; abzpt=f*abzpt
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
    ake(isp)=ake(isp)*am(isp)*dn(isp)
    !
  enddo
  !
  ! ----- MPI -----
  ! total kinetic energy of each species over all ranks
  call MPI_ALLREDUCE(ake,ake0,nsp,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
  ake=ake0
  ! ---------------
end subroutine move
