!###############
subroutine curnt
!###############
  !
  !  subroutine curnt deposits the current at the old particle position
  !  (jxm, jym, jzm), advances the positions (periodic in x and y), then
  !  deposits rho and the current at the new position (jxp, jyp, jzp).
  !  fields uses the average of the two currents. Quadratic-spline
  !  weights with the same stencils as in move.
  !
  use global, only : com_den,nx,ny,nx1,ny1,nsp,maxnp,dx,dy,dt,&
                     jxm,jxp,jym,jyp,jzm,jzp,rho,time,xn,yn,&
                     x,y,vx,vy,vz,gamma,np,IERR,MYRANK,nthreads
  !!$ use omp_lib
  !
  implicit none
  ! ----- MPI -----
  include 'mpif.h'
  real(8), dimension(0:nx1+1,0:ny1+1) :: rho0,jxm0,jym0,jzm0,jxp0,jyp0,jzp0
  ! ---------------
  !
  integer(8) j
  integer(4) isp!,ix1,iy1,ix1p,iy1p,ixh,iyh,ixhp,iyhp
  real(8) qdxdn,qdtdn!,drx,dlx,duy,dly,drxh,dlxh,duyh,dlyh
  ! real(8) q11,q12,q21,q22
  ! real(8) qvx11,qvx12,qvx21,qvx22
  ! real(8) qvy11,qvy12,qvy21,qvy22
  ! real(8) qvz11,qvz12,qvz21,qvz22
  integer(4) ix,iy,ixm,ixp,iym,iyp
  real(8) s2x,s2y,s2xp,s2yp,s2xm,s2ym
  real(8) qxy,qxpy,qxmy,qxyp,qxpyp,qxmyp,qxym,qxpym,qxmym
  !
  real(8), dimension(nsp) :: dn,q,qm
  !
  call com_den(dn,q,qm)
  !
  do isp=1,nsp
    !
    !  ..... set various constants.
    !
    qdxdn=q(isp)*dn(isp)/(dx*dy)
    qdtdn=q(isp)*dn(isp)/(dx*dy)
    !
    !!$ call omp_set_dynamic(.false.)
    !!$ call omp_set_num_threads(nthreads)
    !!$OMP PARALLEL DEFAULT(PRIVATE) &
    !!$OMP SHARED(isp,np,vx,vy,vz,x,y,jxm,jym,jzm,jxp,jyp,jzp,rho,qdtdn,qdxdn,MYRANK)
    !!$OMP DO REDUCTION(+:rho,jxm,jym,jzm,jxp,jyp,jzp)
    do j=MYRANK*maxnp+1,MYRANK*maxnp+np(isp)
      !
      ! compute old particle position quantities.
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
      qxy   = s2x *s2y*qdtdn
      qxpy  = s2xp*s2y*qdtdn
      qxmy  = s2xm*s2y*qdtdn
      qxyp  = s2x *s2yp*qdtdn
      qxpyp = s2xp*s2yp*qdtdn
      qxmyp = s2xm*s2yp*qdtdn
      qxym  = s2x *s2ym*qdtdn
      qxpym = s2xp*s2ym*qdtdn
      qxmym = s2xm*s2ym*qdtdn
      jxm(ix ,iy ) = jxm(ix ,iy ) + qxy  * vx(j,isp)
      jxm(ixp,iy ) = jxm(ixp,iy ) + qxpy * vx(j,isp)
      jxm(ixm,iy ) = jxm(ixm,iy ) + qxmy * vx(j,isp)
      jxm(ix ,iyp) = jxm(ix ,iyp) + qxyp * vx(j,isp)
      jxm(ixp,iyp) = jxm(ixp,iyp) + qxpyp* vx(j,isp)
      jxm(ixm,iyp) = jxm(ixm,iyp) + qxmyp* vx(j,isp)
      jxm(ix ,iym) = jxm(ix ,iym) + qxym * vx(j,isp)
      jxm(ixp,iym) = jxm(ixp,iym) + qxpym* vx(j,isp)
      jxm(ixm,iym) = jxm(ixm,iym) + qxmym* vx(j,isp)
      !
      iy = y(j,isp) + 1.0d0
      iyp= iy + 1
      iym= iy - 1
      s2y = 0.75d0 - (y(j,isp) + 0.5d0 - iy)**2
      s2yp= 0.5d0 * ( 1.0d0 + y(j,isp) - iy)**2
      s2ym= 0.5d0 * (       - y(j,isp) + iy)**2
      !
      qxy   = s2x *s2y*qdtdn
      qxpy  = s2xp*s2y*qdtdn
      qxmy  = s2xm*s2y*qdtdn
      qxyp  = s2x *s2yp*qdtdn
      qxpyp = s2xp*s2yp*qdtdn
      qxmyp = s2xm*s2yp*qdtdn
      qxym  = s2x *s2ym*qdtdn
      qxpym = s2xp*s2ym*qdtdn
      qxmym = s2xm*s2ym*qdtdn
      jzm(ix ,iy ) = jzm(ix ,iy ) + qxy  * vz(j,isp)
      jzm(ixp,iy ) = jzm(ixp,iy ) + qxpy * vz(j,isp)
      jzm(ixm,iy ) = jzm(ixm,iy ) + qxmy * vz(j,isp)
      jzm(ix ,iyp) = jzm(ix ,iyp) + qxyp * vz(j,isp)
      jzm(ixp,iyp) = jzm(ixp,iyp) + qxpyp* vz(j,isp)
      jzm(ixm,iyp) = jzm(ixm,iyp) + qxmyp* vz(j,isp)
      jzm(ix ,iym) = jzm(ix ,iym) + qxym * vz(j,isp)
      jzm(ixp,iym) = jzm(ixp,iym) + qxpym* vz(j,isp)
      jzm(ixm,iym) = jzm(ixm,iym) + qxmym* vz(j,isp)
      !
      ix = x(j,isp) + 1.5d0
      ixp= ix + 1
      ixm= ix - 1
      s2x = 0.75d0 - (x(j,isp) + 1.0d0 - ix)**2
      s2xp= 0.5d0 * ( 1.5d0 + x(j,isp) - ix)**2
      s2xm= 0.5d0 * (-0.5d0 - x(j,isp) + ix)**2
      !
      qxy   = s2x *s2y*qdtdn
      qxpy  = s2xp*s2y*qdtdn
      qxmy  = s2xm*s2y*qdtdn
      qxyp  = s2x *s2yp*qdtdn
      qxpyp = s2xp*s2yp*qdtdn
      qxmyp = s2xm*s2yp*qdtdn
      qxym  = s2x *s2ym*qdtdn
      qxpym = s2xp*s2ym*qdtdn
      qxmym = s2xm*s2ym*qdtdn
      jym(ix ,iy ) = jym(ix ,iy ) + qxy  * vy(j,isp)
      jym(ixp,iy ) = jym(ixp,iy ) + qxpy * vy(j,isp)
      jym(ixm,iy ) = jym(ixm,iy ) + qxmy * vy(j,isp)
      jym(ix ,iyp) = jym(ix ,iyp) + qxyp * vy(j,isp)
      jym(ixp,iyp) = jym(ixp,iyp) + qxpyp* vy(j,isp)
      jym(ixm,iyp) = jym(ixm,iyp) + qxmyp* vy(j,isp)
      jym(ix ,iym) = jym(ix ,iym) + qxym * vy(j,isp)
      jym(ixp,iym) = jym(ixp,iym) + qxpym* vy(j,isp)
      jym(ixm,iym) = jym(ixm,iym) + qxmym* vy(j,isp)

      !   ix1 =x(j,isp)+1
      !   iy1 =y(j,isp)+1
      !   drx =x(j,isp)+1.d0-ix1
      !   duy =y(j,isp)+1.d0-iy1
      !   dlx =1.d0-drx
      !   dly =1.d0-duy
      !   ixh =x(j,isp)+0.5d0
      !   iyh =y(j,isp)+0.5d0
      !   drxh=x(j,isp)+0.5d0-ixh
      !   duyh=y(j,isp)+0.5d0-iyh
      !   dlxh=1.d0-drxh
      !   dlyh=1.d0-duyh
      !   ix1p=ix1+1
      !   iy1p=iy1+1
      !   ixhp=ixh+1
      !   iyhp=iyh+1
      !!
      !   qvx11 = dlxh*dly *vx(j,isp)*qdtdn
      !   qvx12 = dlxh*duy *vx(j,isp)*qdtdn
      !   qvx21 = drxh*dly *vx(j,isp)*qdtdn
      !   qvx22 = drxh*duy *vx(j,isp)*qdtdn
      !   qvy11 = dlx *dlyh*vy(j,isp)*qdtdn
      !   qvy12 = dlx *duyh*vy(j,isp)*qdtdn
      !   qvy21 = drx *dlyh*vy(j,isp)*qdtdn
      !   qvy22 = drx *duyh*vy(j,isp)*qdtdn
      !   qvz11 = dlxh*dlyh*vz(j,isp)*qdtdn
      !   qvz12 = dlxh*duyh*vz(j,isp)*qdtdn
      !   qvz21 = drxh*dlyh*vz(j,isp)*qdtdn
      !   qvz22 = drxh*duyh*vz(j,isp)*qdtdn
      !!
      !   jxm(ixh ,iy1 ) = jxm(ixh ,iy1 ) + qvx11
      !   jxm(ixh ,iy1p) = jxm(ixh ,iy1p) + qvx12
      !   jxm(ixhp,iy1 ) = jxm(ixhp,iy1 ) + qvx21
      !   jxm(ixhp,iy1p) = jxm(ixhp,iy1p) + qvx22
      !   jym(ix1 ,iyh ) = jym(ix1 ,iyh ) + qvy11
      !   jym(ix1 ,iyhp) = jym(ix1 ,iyhp) + qvy12
      !   jym(ix1p,iyh ) = jym(ix1p,iyh ) + qvy21
      !   jym(ix1p,iyhp) = jym(ix1p,iyhp) + qvy22
      !   jzm(ixh ,iyh ) = jzm(ixh ,iyh ) + qvz11
      !   jzm(ixh ,iyhp) = jzm(ixh ,iyhp) + qvz12
      !   jzm(ixhp,iyh ) = jzm(ixhp,iyh ) + qvz21
      !   jzm(ixhp,iyhp) = jzm(ixhp,iyhp) + qvz22
      !
      ! compute new particle position quantities.
      ! (x, y are in cells; one wrap is enough because |v| < c and
      ! c*dt < dx)
      !
      x(j,isp)=x(j,isp)+vx(j,isp)*dt/dx
      y(j,isp)=y(j,isp)+vy(j,isp)*dt/dy
      !   if(x(j,isp)  <  0.d0) then
      !     x(j,isp)=-x(j,isp)
      !     vx(j,isp)=-vx(j,isp)
      !   endif
      !   if(x(j,isp)  >=   xn) then
      !     x(j,isp)=2*xn-x(j,isp)
      !     vx(j,isp)=-vx(j,isp)
      !   endif
      if(x(j,isp)  <  0.d0) x(j,isp)=x(j,isp)+xn
      if(x(j,isp)  >=   xn) x(j,isp)=x(j,isp)-xn
      if(y(j,isp)  <  0.d0) y(j,isp)=y(j,isp)+yn
      if(y(j,isp)  >=   yn) y(j,isp)=y(j,isp)-yn

      ! charge density at the new position (integer nodes)
      ix = x(j,isp) + 1.5d0;    iy = y(j,isp) + 1.5d0
      ixp= ix + 1;              iyp= iy + 1
      ixm= ix - 1;              iym= iy - 1
      !
      s2x = 0.75d0 - (x(j,isp) + 1.d0 - ix)**2
      s2xp= 0.5d0 * ( 1.5d0 + x(j,isp) - ix)**2
      s2xm= 0.5d0 * (-0.5d0 - x(j,isp) + ix)**2
      s2y = 0.75d0 - (y(j,isp) + 1.d0 - iy)**2
      s2yp= 0.5d0 * ( 1.5d0 + y(j,isp) - iy)**2
      s2ym= 0.5d0 * (-0.5d0 - y(j,isp) + iy)**2
      !
      qxy   = s2x *s2y *qdxdn
      qxpy  = s2xp*s2y *qdxdn
      qxmy  = s2xm*s2y *qdxdn
      qxyp  = s2x *s2yp*qdxdn
      qxpyp = s2xp*s2yp*qdxdn
      qxmyp = s2xm*s2yp*qdxdn
      qxym  = s2x *s2ym*qdxdn
      qxpym = s2xp*s2ym*qdxdn
      qxmym = s2xm*s2ym*qdxdn
      rho(ix ,iy ) = rho(ix ,iy ) + qxy
      rho(ixp,iy ) = rho(ixp,iy ) + qxpy
      rho(ixm,iy ) = rho(ixm,iy ) + qxmy
      rho(ix ,iyp) = rho(ix ,iyp) + qxyp
      rho(ixp,iyp) = rho(ixp,iyp) + qxpyp
      rho(ixm,iyp) = rho(ixm,iyp) + qxmyp
      rho(ix ,iym) = rho(ix ,iym) + qxym
      rho(ixp,iym) = rho(ixp,iym) + qxpym
      rho(ixm,iym) = rho(ixm,iym) + qxmym
      !
      ix = x(j,isp) + 1.0d0
      ixp= ix + 1
      ixm= ix - 1
      s2x = 0.75d0 - (x(j,isp) + 0.5d0 - ix)**2
      s2xp= 0.5d0 * ( 1.0d0 + x(j,isp) - ix)**2
      s2xm= 0.5d0 * (       - x(j,isp) + ix)**2
      !
      qxy   = s2x *s2y*qdtdn
      qxpy  = s2xp*s2y*qdtdn
      qxmy  = s2xm*s2y*qdtdn
      qxyp  = s2x *s2yp*qdtdn
      qxpyp = s2xp*s2yp*qdtdn
      qxmyp = s2xm*s2yp*qdtdn
      qxym  = s2x *s2ym*qdtdn
      qxpym = s2xp*s2ym*qdtdn
      qxmym = s2xm*s2ym*qdtdn
      jxp(ix ,iy ) = jxp(ix ,iy ) + qxy  * vx(j,isp)
      jxp(ixp,iy ) = jxp(ixp,iy ) + qxpy * vx(j,isp)
      jxp(ixm,iy ) = jxp(ixm,iy ) + qxmy * vx(j,isp)
      jxp(ix ,iyp) = jxp(ix ,iyp) + qxyp * vx(j,isp)
      jxp(ixp,iyp) = jxp(ixp,iyp) + qxpyp* vx(j,isp)
      jxp(ixm,iyp) = jxp(ixm,iyp) + qxmyp* vx(j,isp)
      jxp(ix ,iym) = jxp(ix ,iym) + qxym * vx(j,isp)
      jxp(ixp,iym) = jxp(ixp,iym) + qxpym* vx(j,isp)
      jxp(ixm,iym) = jxp(ixm,iym) + qxmym* vx(j,isp)
      !
      iy = y(j,isp) + 1.0d0
      iyp= iy + 1
      iym= iy - 1
      s2y = 0.75d0 - (y(j,isp) + 0.5d0 - iy)**2
      s2yp= 0.5d0 * ( 1.0d0 + y(j,isp) - iy)**2
      s2ym= 0.5d0 * (       - y(j,isp) + iy)**2
      !
      qxy   = s2x *s2y*qdtdn
      qxpy  = s2xp*s2y*qdtdn
      qxmy  = s2xm*s2y*qdtdn
      qxyp  = s2x *s2yp*qdtdn
      qxpyp = s2xp*s2yp*qdtdn
      qxmyp = s2xm*s2yp*qdtdn
      qxym  = s2x *s2ym*qdtdn
      qxpym = s2xp*s2ym*qdtdn
      qxmym = s2xm*s2ym*qdtdn
      jzp(ix ,iy ) = jzp(ix ,iy ) + qxy  * vz(j,isp)
      jzp(ixp,iy ) = jzp(ixp,iy ) + qxpy * vz(j,isp)
      jzp(ixm,iy ) = jzp(ixm,iy ) + qxmy * vz(j,isp)
      jzp(ix ,iyp) = jzp(ix ,iyp) + qxyp * vz(j,isp)
      jzp(ixp,iyp) = jzp(ixp,iyp) + qxpyp* vz(j,isp)
      jzp(ixm,iyp) = jzp(ixm,iyp) + qxmyp* vz(j,isp)
      jzp(ix ,iym) = jzp(ix ,iym) + qxym * vz(j,isp)
      jzp(ixp,iym) = jzp(ixp,iym) + qxpym* vz(j,isp)
      jzp(ixm,iym) = jzp(ixm,iym) + qxmym* vz(j,isp)
      !
      ix = x(j,isp) + 1.5d0
      ixp= ix + 1
      ixm= ix - 1
      s2x = 0.75d0 - (x(j,isp) + 1.0d0 - ix)**2
      s2xp= 0.5d0 * ( 1.5d0 + x(j,isp) - ix)**2
      s2xm= 0.5d0 * (-0.5d0 - x(j,isp) + ix)**2
      !
      qxy   = s2x *s2y*qdtdn
      qxpy  = s2xp*s2y*qdtdn
      qxmy  = s2xm*s2y*qdtdn
      qxyp  = s2x *s2yp*qdtdn
      qxpyp = s2xp*s2yp*qdtdn
      qxmyp = s2xm*s2yp*qdtdn
      qxym  = s2x *s2ym*qdtdn
      qxpym = s2xp*s2ym*qdtdn
      qxmym = s2xm*s2ym*qdtdn
      jyp(ix ,iy ) = jyp(ix ,iy ) + qxy  * vy(j,isp)
      jyp(ixp,iy ) = jyp(ixp,iy ) + qxpy * vy(j,isp)
      jyp(ixm,iy ) = jyp(ixm,iy ) + qxmy * vy(j,isp)
      jyp(ix ,iyp) = jyp(ix ,iyp) + qxyp * vy(j,isp)
      jyp(ixp,iyp) = jyp(ixp,iyp) + qxpyp* vy(j,isp)
      jyp(ixm,iyp) = jyp(ixm,iyp) + qxmyp* vy(j,isp)
      jyp(ix ,iym) = jyp(ix ,iym) + qxym * vy(j,isp)
      jyp(ixp,iym) = jyp(ixp,iym) + qxpym* vy(j,isp)
      jyp(ixm,iym) = jyp(ixm,iym) + qxmym* vy(j,isp)

      !   ix1 =x(j,isp)+1
      !   iy1 =y(j,isp)+1
      !   drx =x(j,isp)+1.d0-ix1
      !   duy =y(j,isp)+1.d0-iy1
      !   dlx =1.d0-drx
      !   dly =1.d0-duy
      !   ixh =x(j,isp)+0.5d0
      !   iyh =y(j,isp)+0.5d0
      !   drxh=x(j,isp)+0.5d0-ixh
      !   duyh=y(j,isp)+0.5d0-iyh
      !   dlxh=1.d0-drxh
      !   dlyh=1.d0-duyh
      !   ix1p=ix1+1
      !   iy1p=iy1+1
      !   ixhp=ixh+1
      !   iyhp=iyh+1
      !!
      !   q11   = dlx*dly*qdxdn
      !   q12   = dlx*duy*qdxdn
      !   q21   = drx*dly*qdxdn
      !   q22   = drx*duy*qdxdn
      !   qvx11 = dlxh*dly *vx(j,isp)*qdtdn
      !   qvx12 = dlxh*duy *vx(j,isp)*qdtdn
      !   qvx21 = drxh*dly *vx(j,isp)*qdtdn
      !   qvx22 = drxh*duy *vx(j,isp)*qdtdn
      !   qvy11 = dlx *dlyh*vy(j,isp)*qdtdn
      !   qvy12 = dlx *duyh*vy(j,isp)*qdtdn
      !   qvy21 = drx *dlyh*vy(j,isp)*qdtdn
      !   qvy22 = drx *duyh*vy(j,isp)*qdtdn
      !   qvz11 = dlxh*dlyh*vz(j,isp)*qdtdn
      !   qvz12 = dlxh*duyh*vz(j,isp)*qdtdn
      !   qvz21 = drxh*dlyh*vz(j,isp)*qdtdn
      !   qvz22 = drxh*duyh*vz(j,isp)*qdtdn
      !!
      !   rho(ix1 ,iy1 ) = rho(ix1 ,iy1 ) + q11
      !   rho(ix1 ,iy1p) = rho(ix1 ,iy1p) + q12
      !   rho(ix1p,iy1 ) = rho(ix1p,iy1 ) + q21
      !   rho(ix1p,iy1p) = rho(ix1p,iy1p) + q22
      !   jxp(ixh ,iy1 ) = jxp(ixh ,iy1 ) + qvx11
      !   jxp(ixh ,iy1p) = jxp(ixh ,iy1p) + qvx12
      !   jxp(ixhp,iy1 ) = jxp(ixhp,iy1 ) + qvx21
      !   jxp(ixhp,iy1p) = jxp(ixhp,iy1p) + qvx22
      !   jyp(ix1 ,iyh ) = jyp(ix1 ,iyh ) + qvy11
      !   jyp(ix1 ,iyhp) = jyp(ix1 ,iyhp) + qvy12
      !   jyp(ix1p,iyh ) = jyp(ix1p,iyh ) + qvy21
      !   jyp(ix1p,iyhp) = jyp(ix1p,iyhp) + qvy22
      !   jzp(ixh ,iyh ) = jzp(ixh ,iyh ) + qvz11
      !   jzp(ixh ,iyhp) = jzp(ixh ,iyhp) + qvz12
      !   jzp(ixhp,iyh ) = jzp(ixhp,iyh ) + qvz21
      !   jzp(ixhp,iyhp) = jzp(ixhp,iyhp) + qvz22
    enddo
    !!$OMP ENDDO
    !!$OMP END PARALLEL
    !
  enddo
  !
  ! ----- MPI -----
  ! sum rho and J over all ranks (every rank then has the full arrays)
  call MPI_ALLREDUCE(rho,rho0,(nx1+2)*(ny1+2),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(jxm,jxm0,(nx1+2)*(ny1+2),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(jym,jym0,(nx1+2)*(ny1+2),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(jzm,jzm0,(nx1+2)*(ny1+2),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(jxp,jxp0,(nx1+2)*(ny1+2),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(jyp,jyp0,(nx1+2)*(ny1+2),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(jzp,jzp0,(nx1+2)*(ny1+2),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
  rho=rho0
  jxm=jxm0
  jym=jym0
  jzm=jzm0
  jxp=jxp0
  jyp=jyp0
  jzp=jzp0
  ! ---------------
  !
end subroutine curnt
