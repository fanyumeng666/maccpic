!#################
subroutine fields
!#################
  !
  !  subroutine fields solves for phi and ex, computes field energy,
  !  and advance the wave radiation fields.
  !  Order: fold the guard-cell deposits of rho and J, B by dt/2, E by dt,
  !  B by dt/2 (periodic in x and y), filter, Poisson correction of E
  !  (div E = rho), field energies, then reset rho and J.
  !
  use global, only : nx,ny,nx1,ny1,nx2,ny2,dx,dy,dt,tem,ctx,cty,ctxh,ctyh,nthreads,&
                     bx,by,bz,ex,ey,ez,bxm,bym,bzm,jxm,jxp,jym,jyp,jzm,jzp,rho,&
                     exe,eye,eze,bxe,bye,bze,nt,it,time,rnmec2i,anorml1,isav
  !$ use omp_lib
  !
  implicit none
  integer(4) ix,iy,ixm,ixp,iym,iyp
  real(8) exes,eyes,ezes,bxes,byes,bzes
  !
  !     store the previous B
  !
  !$ call omp_set_dynamic(.false.)
  !$ call omp_set_num_threads(nthreads)
  !$OMP PARALLEL PRIVATE(iy) &
  !$OMP SHARED(rho,jxm,jym,jzm,jxp,jyp,jzp)
  !$OMP DO
  do ix=0,nx2
    do iy=0,ny2
      bxm(ix,iy)=bx(ix,iy)
      bym(ix,iy)=by(ix,iy)
      bzm(ix,iy)=bz(ix,iy)
    enddo
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL
  !
  !     transform charge density.
  !     (add the deposits in the guard cells 0, nx+1, nx+2 to the periodic
  !     images, then copy back; the same is done in y below)
  !
  !$OMP PARALLEL PRIVATE(iy) &
  !$OMP SHARED(rho,jxm,jym,jzm,jxp,jyp,jzp)
  !$OMP DO
  do iy=0,ny2
    rho(1,iy)=rho(1,iy)+rho(nx1,iy)
    jxm(1,iy)=jxm(1,iy)+jxm(nx1,iy)
    jxp(1,iy)=jxp(1,iy)+jxp(nx1,iy)
    jym(1,iy)=jym(1,iy)+jym(nx1,iy)
    jyp(1,iy)=jyp(1,iy)+jyp(nx1,iy)
    jzm(1,iy)=jzm(1,iy)+jzm(nx1,iy)
    jzp(1,iy)=jzp(1,iy)+jzp(nx1,iy)
    rho(2,iy)=rho(2,iy)+rho(nx2,iy)
    jxm(2,iy)=jxm(2,iy)+jxm(nx2,iy)
    jxp(2,iy)=jxp(2,iy)+jxp(nx2,iy)
    jym(2,iy)=jym(2,iy)+jym(nx2,iy)
    jyp(2,iy)=jyp(2,iy)+jyp(nx2,iy)
    jzm(2,iy)=jzm(2,iy)+jzm(nx2,iy)
    jzp(2,iy)=jzp(2,iy)+jzp(nx2,iy)
    rho(nx,iy)=rho(nx,iy)+rho(0,iy)
    jxm(nx,iy)=jxm(nx,iy)+jxm(0,iy)
    jxp(nx,iy)=jxp(nx,iy)+jxp(0,iy)
    jym(nx,iy)=jym(nx,iy)+jym(0,iy)
    jyp(nx,iy)=jyp(nx,iy)+jyp(0,iy)
    jzm(nx,iy)=jzm(nx,iy)+jzm(0,iy)
    jzp(nx,iy)=jzp(nx,iy)+jzp(0,iy)
    rho(nx1,iy)=rho(1,iy)
    jxm(nx1,iy)=jxm(1,iy)
    jxp(nx1,iy)=jxp(1,iy)
    jym(nx1,iy)=jym(1,iy)
    jyp(nx1,iy)=jyp(1,iy)
    jzm(nx1,iy)=jzm(1,iy)
    jzp(nx1,iy)=jzp(1,iy)
    rho(nx2,iy)=rho(2,iy)
    jxm(nx2,iy)=jxm(2,iy)
    jxp(nx2,iy)=jxp(2,iy)
    jym(nx2,iy)=jym(2,iy)
    jyp(nx2,iy)=jyp(2,iy)
    jzm(nx2,iy)=jzm(2,iy)
    jzp(nx2,iy)=jzp(2,iy)
    rho(0,iy)=rho(nx,iy)
    jxm(0,iy)=jxm(nx,iy)
    jxp(0,iy)=jxp(nx,iy)
    jym(0,iy)=jym(nx,iy)
    jyp(0,iy)=jyp(nx,iy)
    jzm(0,iy)=jzm(nx,iy)
    jzp(0,iy)=jzp(nx,iy)
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL

  !$OMP PARALLEL PRIVATE(ix) &
  !$OMP SHARED(rho,jxm,jym,jzm,jxp,jyp,jzp)
  !$OMP DO
  do ix=0,nx2
    rho(ix,1)=rho(ix,1)+rho(ix,ny1)
    jxm(ix,1)=jxm(ix,1)+jxm(ix,ny1)
    jxp(ix,1)=jxp(ix,1)+jxp(ix,ny1)
    jym(ix,1)=jym(ix,1)+jym(ix,ny1)
    jyp(ix,1)=jyp(ix,1)+jyp(ix,ny1)
    jzm(ix,1)=jzm(ix,1)+jzm(ix,ny1)
    jzp(ix,1)=jzp(ix,1)+jzp(ix,ny1)
    rho(ix,2)=rho(ix,2)+rho(ix,ny2)
    jxm(ix,2)=jxm(ix,2)+jxm(ix,ny2)
    jxp(ix,2)=jxp(ix,2)+jxp(ix,ny2)
    jym(ix,2)=jym(ix,2)+jym(ix,ny2)
    jyp(ix,2)=jyp(ix,2)+jyp(ix,ny2)
    jzm(ix,2)=jzm(ix,2)+jzm(ix,ny2)
    jzp(ix,2)=jzp(ix,2)+jzp(ix,ny2)
    rho(ix,ny)=rho(ix,ny)+rho(ix,0)
    jxm(ix,ny)=jxm(ix,ny)+jxm(ix,0)
    jxp(ix,ny)=jxp(ix,ny)+jxp(ix,0)
    jym(ix,ny)=jym(ix,ny)+jym(ix,0)
    jyp(ix,ny)=jyp(ix,ny)+jyp(ix,0)
    jzm(ix,ny)=jzm(ix,ny)+jzm(ix,0)
    jzp(ix,ny)=jzp(ix,ny)+jzp(ix,0)
    rho(ix,ny1)=rho(ix,1)
    jxm(ix,ny1)=jxm(ix,1)
    jxp(ix,ny1)=jxp(ix,1)
    jym(ix,ny1)=jym(ix,1)
    jyp(ix,ny1)=jyp(ix,1)
    jzm(ix,ny1)=jzm(ix,1)
    jzp(ix,ny1)=jzp(ix,1)
    rho(ix,ny2)=rho(ix,2)
    jxm(ix,ny2)=jxm(ix,2)
    jxp(ix,ny2)=jxp(ix,2)
    jym(ix,ny2)=jym(ix,2)
    jyp(ix,ny2)=jyp(ix,2)
    jzm(ix,ny2)=jzm(ix,2)
    jzp(ix,ny2)=jzp(ix,2)
    rho(ix,0)=rho(ix,ny)
    jxm(ix,0)=jxm(ix,ny)
    jxp(ix,0)=jxp(ix,ny)
    jym(ix,0)=jym(ix,ny)
    jyp(ix,0)=jyp(ix,ny)
    jzm(ix,0)=jzm(ix,ny)
    jzp(ix,0)=jzp(ix,ny)
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL
  !do i=1,np(1)
  !write(999,'(4e13.5)') y(MYRANK*maxnp+i,1),vz(MYRANK*maxnp+i,1),y(MYRANK*maxnp+i,2),vz(MYRANK*maxnp+i,2)
  !enddo
  ! debug output: jzm at step 1 to fort.999 (written by every rank)
  if(it==1) then
    do iy=1,ny
      write(999,'(256e13.5)') (jzm(ix,iy),ix=1,nx)
    enddo
  endif

  !   rho(1,1:ny1)=rho(1,1:ny1)+rho(nx1,1:ny1)
  !   jxm(1,1:ny1)=jxm(1,1:ny1)+jxm(nx1,1:ny1)
  !   jxp(1,1:ny1)=jxp(1,1:ny1)+jxp(nx1,1:ny1)
  !   jym(1,1:ny1)=jym(1,1:ny1)+jym(nx1,1:ny1)
  !   jyp(1,1:ny1)=jyp(1,1:ny1)+jyp(nx1,1:ny1)
  !   jzm(1,1:ny1)=jzm(1,1:ny1)+jzm(nx1,1:ny1)
  !   jzp(1,1:ny1)=jzp(1,1:ny1)+jzp(nx1,1:ny1)
  !   rho(1:nx1,1)=rho(1:nx1,1)+rho(1:nx1,ny1)
  !   jxm(1:nx1,1)=jxm(1:nx1,1)+jxm(1:nx1,ny1)
  !   jxp(1:nx1,1)=jxp(1:nx1,1)+jxp(1:nx1,ny1)
  !   jym(1:nx1,1)=jym(1:nx1,1)+jym(1:nx1,ny1)
  !   jyp(1:nx1,1)=jyp(1:nx1,1)+jyp(1:nx1,ny1)
  !   jzm(1:nx1,1)=jzm(1:nx1,1)+jzm(1:nx1,ny1)
  !   jzp(1:nx1,1)=jzp(1:nx1,1)+jzp(1:nx1,ny1)
  !   rho(nx1,1:ny1)=rho(1,1:ny1)
  !   jxm(nx1,1:ny1)=jxm(1,1:ny1)
  !   jxp(nx1,1:ny1)=jxp(1,1:ny1)
  !   jym(nx1,1:ny1)=jym(1,1:ny1)
  !   jyp(nx1,1:ny1)=jyp(1,1:ny1)
  !   jzm(nx1,1:ny1)=jzm(1,1:ny1)
  !   jzp(nx1,1:ny1)=jzp(1,1:ny1)
  !   rho(1:nx1,ny1)=rho(1:nx1,1)
  !   jxm(1:nx1,ny1)=jxm(1:nx1,1)
  !   jxp(1:nx1,ny1)=jxp(1:nx1,1)
  !   jym(1:nx1,ny1)=jym(1:nx1,1)
  !   jyp(1:nx1,ny1)=jyp(1:nx1,1)
  !   jzm(1:nx1,ny1)=jzm(1:nx1,1)
  !   jzp(1:nx1,ny1)=jzp(1:nx1,1)
  !
  ! first dt/2 step for B
  !
  !$OMP PARALLEL PRIVATE(iy,iyp,iym,ix,ixp,ixm) &
  !$OMP SHARED(ex,ey,ez,bx,by,bz)
  !$OMP DO
  do iy=1,ny1
    iyp=iy+1;  iym=iy-1
    do ix=1,nx1
      ixp=ix+1;  ixm=ix-1
      bx(ix,iy)=bx(ix,iy)-ctyh*(ez(ix,iy)-ez(ix,iym))
      by(ix,iy)=by(ix,iy)+ctxh*(ez(ix,iy)-ez(ixm,iy))
      bz(ix,iy)=bz(ix,iy)-ctxh*(ey(ixp,iy)-ey(ix,iy))+ctyh*(ex(ix,iyp)-ex(ix,iy))
    enddo
    bx(0,iy)=bx(nx,iy); bx(nx2,iy)=bx(2,iy)
    by(0,iy)=by(nx,iy); by(nx2,iy)=by(2,iy)
    bz(0,iy)=bz(nx,iy); bz(nx2,iy)=bz(2,iy)
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL
  bx(0:nx2,0)=bx(0:nx2,ny); bx(0:nx2,ny2)=bx(0:nx2,2)
  by(0:nx2,0)=by(0:nx2,ny); by(0:nx2,ny2)=by(0:nx2,2)
  bz(0:nx2,0)=bz(0:nx2,ny); bz(0:nx2,ny2)=bz(0:nx2,2)
  !
  ! dt step for E
  !
  !$OMP PARALLEL PRIVATE(iy,iyp,iym,ix,ixp,ixm) &
  !$OMP SHARED(ex,ey,ez,bx,by,bz,jxm,jym,jzm,jxp,jyp,jzp)
  !$OMP DO
  do iy=1,ny1
    iyp=iy+1;  iym=iy-1
    do ix=1,nx1
      ixp=ix+1;  ixm=ix-1
      ex(ix,iy)=ex(ix,iy)+cty*(bz(ix,iy)-bz(ix,iym))-tem*(jxm(ix,iy)+jxp(ix,iy))
      ey(ix,iy)=ey(ix,iy)-ctx*(bz(ix,iy)-bz(ixm,iy))-tem*(jym(ix,iy)+jyp(ix,iy))
      ez(ix,iy)=ez(ix,iy)+ctx*(by(ixp,iy)-by(ix,iy))-&
                          cty*(bx(ix,iyp)-bx(ix,iy))-tem*(jzm(ix,iy)+jzp(ix,iy))
    enddo
    ex(0,iy)=ex(nx,iy); ex(nx2,iy)=ex(2,iy)
    ey(0,iy)=ey(nx,iy); ey(nx2,iy)=ey(2,iy)
    ez(0,iy)=ez(nx,iy); ez(nx2,iy)=ez(2,iy)
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL
  ex(0:nx2,0)=ex(0:nx2,ny); ex(0:nx2,ny2)=ex(0:nx2,2)
  ey(0:nx2,0)=ey(0:nx2,ny); ey(0:nx2,ny2)=ey(0:nx2,2)
  ez(0:nx2,0)=ez(0:nx2,ny); ez(0:nx2,ny2)=ez(0:nx2,2)
  !
  ! second dt/2 step for B
  !
  !$OMP PARALLEL PRIVATE(iy,iyp,iym,ix,ixp,ixm) &
  !$OMP SHARED(ex,ey,ez,bx,by,bz)
  !$OMP DO
  do iy=1,ny1
    iyp=iy+1;  iym=iy-1
    do ix=1,nx1
      ixp=ix+1;  ixm=ix-1
      bx(ix,iy)=bx(ix,iy)-ctyh*(ez(ix,iy)-ez(ix,iym))
      by(ix,iy)=by(ix,iy)+ctxh*(ez(ix,iy)-ez(ixm,iy))
      bz(ix,iy)=bz(ix,iy)-ctxh*(ey(ixp,iy)-ey(ix,iy))+ctyh*(ex(ix,iyp)-ex(ix,iy))
    enddo
    bx(0,iy)=bx(nx,iy); bx(nx2,iy)=bx(2,iy)
    by(0,iy)=by(nx,iy); by(nx2,iy)=by(2,iy)
    bz(0,iy)=bz(nx,iy); bz(nx2,iy)=bz(2,iy)
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL
  bx(0:nx2,0)=bx(0:nx2,ny); bx(0:nx2,ny2)=bx(0:nx2,2)
  by(0:nx2,0)=by(0:nx2,ny); by(0:nx2,ny2)=by(0:nx2,2)
  bz(0:nx2,0)=bz(0:nx2,ny); bz(0:nx2,ny2)=bz(0:nx2,2)
  !
  ! filtering rho, E and B
  ! (rho every step, E and B only every isav steps)
  !
  call smooth2(rho)
  if(mod(it,isav) == 0) then
    call smooth2(ex)
    call smooth2(ey)
    call smooth2(ez)
    call smooth2(bx)
    call smooth2(by)
    call smooth2(bz)
  endif
  !
  call Poisson
  !
  ! calculate radiation field energies.
  ! (sum over rows iy = 88..168 only, scaled by ny/80)
  !
  if(mod(it,isav) == 0) then
    exes =0.d0
    eyes =0.d0
    ezes =0.d0
    bxes =0.d0
    byes =0.d0
    bzes =0.d0
    !
    !$OMP PARALLEL PRIVATE(ix,iy) &
    !$OMP SHARED(ex,exes)
    !$OMP DO REDUCTION(+:exes)
    do ix=1,nx
      do iy=88,168
        exes =exes+ex(ix,iy)**2
      enddo
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL

    !$OMP PARALLEL PRIVATE(ix,iy) &
    !$OMP SHARED(ey,eyes)
    !$OMP DO REDUCTION(+:eyes)
    do ix=1,nx
      do iy=88,168
        eyes =eyes+ey(ix,iy)**2
      enddo
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL

    !$OMP PARALLEL PRIVATE(ix,iy) &
    !$OMP SHARED(ez,ezes)
    !$OMP DO REDUCTION(+:ezes)
    do ix=1,nx
      do iy=88,168
        ezes =ezes+ez(ix,iy)**2
      enddo
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL

    !$OMP PARALLEL PRIVATE(ix,iy) &
    !$OMP SHARED(bx,bxes)
    !$OMP DO REDUCTION(+:bxes)
    do ix=1,nx
      do iy=88,168
        bxes =bxes+bx(ix,iy)**2
      enddo
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL

    !$OMP PARALLEL PRIVATE(ix,iy) &
    !$OMP SHARED(by,byes)
    !$OMP DO REDUCTION(+:byes)
    do ix=1,nx
      do iy=88,168
        byes =byes+by(ix,iy)**2
      enddo
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL

    !$OMP PARALLEL PRIVATE(ix,iy) &
    !$OMP SHARED(bz,bzes)
    !$OMP DO REDUCTION(+:bzes)
    do ix=1,nx
      do iy=88,168
        bzes =bzes+bz(ix,iy)**2
      enddo
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL
    !
    exe =exes *anorml1 *ny/80.d0
    eye =eyes *anorml1*ny/80.d0
    eze =ezes *anorml1*ny/80.d0
    bxe =bxes *anorml1*ny/80.d0
    bye =byes *anorml1*ny/80.d0
    bze =bzes *anorml1*ny/80.d0
  endif
  !
  ! reset rho and J for the next step
  !$OMP PARALLEL PRIVATE(ix,iy) &
  !$OMP SHARED(rho,jxm,jym,jzm,jxp,jyp,jzp)
  !$OMP DO
  do ix=0,nx2
    do iy=0,ny2
      rho(ix,iy)=0.d0
      jxm(ix,iy)=0.d0
      jxp(ix,iy)=0.d0
      jym(ix,iy)=0.d0
      jyp(ix,iy)=0.d0
      jzm(ix,iy)=0.d0
      jzp(ix,iy)=0.d0
    enddo
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL
  !
end subroutine fields
!
!#################
subroutine Poisson
!#################
  !
  !  subroutine Poisson corrects E so that div E = rho: it solves
  !  lap(ph) = div E - rho (periodic 5-point Laplacian, dx = dy assumed)
  !  with cgm and subtracts grad(ph) from Ex and Ey.
  !
  use global, only : nx,ny,nx1,ny1,nx2,ny2,dx,dy,ex,ey,rho,nthreads
  !$ use omp_lib
  !
  implicit none
  integer(4) ix,ixm,ixp,iy,iym,iyp
  real(8)    dd,dx2,dy2,f(nx1,ny1),ph(nx1,ny1),php(nx1,ny1)
  !
  dx2=1.d0/dx/dx; dy2=1.d0/dy/dy
  dd=2.d0*dx2+2.d0*dy2
  !
  !$ call omp_set_dynamic(.false.)
  !$ call omp_set_num_threads(nthreads)
  !$OMP PARALLEL DEFAULT(PRIVATE) &
  !$OMP SHARED(f,ex,ey,rho)
  !$OMP DO
  do iy=1,ny
    iym=iy-1; iyp=iy+1
    do ix=2,nx
      ixm=ix-1; ixp=ix+1
      f(ix,iy) = (ex(ix,iy)-ex(ixm,iy))/dx+(ey(ix,iy)-ey(ix,iym))/dy-rho(ix,iy)
    enddo
    f(1,iy) = (ex(1,iy)-ex(nx,iy))/dx+(ey(1,iy)-ey(1,iym))/dy-rho(1,iy)
    f(nx1,iy)=f(1,iy)
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL
  f(:,ny1)=f(:,1)
  !
  f=f*dx*dx
  ph(:,:) = 0.d-4
  ! call sor(f,ph)
  call cgm(f,ph)
  !
  !$OMP PARALLEL DEFAULT(PRIVATE) &
  !$OMP SHARED(php,ex,ph)
  !$OMP DO
  do iy=1,ny
    do ix=1,nx
      ixp=ix+1
      php(ix,iy) = ex(ix,iy) - (ph(ixp,iy)-ph(ix,iy))/dx
    enddo
    php(nx1,iy) = php(1,iy)
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL
  php(:,ny1) = php(:,1)
  ex(1:nx1,1:ny1) = php(1:nx1,1:ny1)
  ex(0,1:ny1)=ex(nx,1:ny1); ex(nx2,1:ny1)=ex(2,1:ny1)
  ex(0:nx2,0)=ex(0:nx2,ny); ex(0:nx2,ny2)=ex(0:nx2,2)
  !
  !$OMP PARALLEL DEFAULT(PRIVATE) &
  !$OMP SHARED(php,ey,ph)
  !$OMP DO
  do ix=1,nx
    do iy=1,ny
      iyp=iy+1
      php(ix,iy) = ey(ix,iy) - (ph(ix,iyp)-ph(ix,iy))/dy
    enddo
    php(ix,ny1) = php(ix,1)
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL
  php(nx1,:) = php(1,:)
  ey(1:nx1,1:ny1) = php(1:nx1,1:ny1)
  ey(0,1:ny1)=ey(nx,1:ny1); ey(nx2,1:ny1)=ey(2,1:ny1)
  ey(0:nx2,0)=ey(0:nx2,ny); ey(0:nx2,ny2)=ey(0:nx2,2)
  !
end subroutine Poisson
!
!#############
SUBROUTINE sor(f,u)
!#############
  ! SOR solver for the same equation as cgm (not used; see Poisson)
  use global, only : nx,ny,nx1,ny1,nthreads
  !$ use omp_lib
  !
  implicit none
  REAL(8) rjac,f(nx1,ny1),u(nx1,ny1)
  INTEGER, PARAMETER :: MAXITS=10000
  REAL(8), PARAMETER :: EPS=1.d-5
  INTEGER ipass,j,l,n,jm
  REAL(8) anorm,anormf,omega,resid
  rjac = (cos(6.283185307179586d0/nx1)+cos(6.283185307179586d0/ny1))/2.d0
  anormf=0.d0
  !$ call omp_set_dynamic(.false.)
  !$ call omp_set_num_threads(nthreads)
  !$OMP PARALLEL DEFAULT(PRIVATE) &
  !$OMP SHARED(f,anormf)
  !$OMP DO REDUCTION(+:anormf)
  do j=1,nx
    do l=1,ny
      anormf=anormf+abs(f(j,l))
    enddo
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL
  omega=1.d0
  do n=1,MAXITS
    anorm=0.d0
    do ipass=1,2
      !$OMP PARALLEL PRIVATE(j,jm,l) &
      !$OMP SHARED(resid,anorm,omega,u,f)
      !$OMP DO REDUCTION(+:anorm)
      do j=1,nx
        jm=j-1
        if(j==1) jm=nx
        do l=2,ny
          resid=u(j+1,l)+u(jm,l)+u(j,l+1)+u(j,l-1)-4.d0*u(j,l)-f(j,l)
          anorm=anorm+abs(resid)
          u(j,l)=u(j,l)+omega*resid/4.d0
        enddo
        resid=u(j+1,1)+u(jm,1)+u(j,2)+u(j,ny)-4.d0*u(j,1)-f(j,1)
        anorm=anorm+abs(resid)
        u(j,1)=u(j,1)+omega*resid/4.d0
      enddo
      !$OMP ENDDO
      !$OMP END PARALLEL
      u(nx1,:)=u(1,:)
      u(:,ny1)=u(:,1)
      if(n==1.and.ipass==1) then
        omega=1.d0/(1.d0-.5d0*rjac**2)
      else
        omega=1.d0/(1.d0-.25d0*rjac**2*omega)
      endif
    enddo
    if(anorm<EPS*anormf)return
  enddo
  write(*,*) 'MAXITS exceeded in sor'
END SUBROUTINE sor

!#############
subroutine cgm(f,u)
!#############
  ! Conjugate-gradient solver of lap(u) = f on the periodic grid
  ! (5-point Laplacian), used by Poisson; stops when mean(r**2) < EPS
  use global, only : nx,ny,nx1,ny1,nthreads
  !$ use omp_lib
  !
  implicit none
  REAL(8)  ak,bk,r_norm,rr,pap,rr0
  REAL(8), DIMENSION(nx1,ny1) :: f,u,r,r1,p,ap
  INTEGER, PARAMETER :: MAXITS=10000
  REAL(8), PARAMETER :: EPS=1.d-9
  INTEGER j,l,n,jm,lm
  REAL(8) anorm,anormf,omega
  !
  ! initial residual r = f - lap(u) and first search direction p = r
  !$ call omp_set_dynamic(.false.)
  !$ call omp_set_num_threads(nthreads)
  !$OMP PARALLEL DEFAULT(PRIVATE) &
  !$OMP SHARED(r,f,u)
  !$OMP DO
  do j=1,nx
    jm=j-1
    if(j==1) jm=nx
    do l=1,ny
      lm=l-1
      if(l==1) lm=ny
      r(j,l)=f(j,l)-u(j+1,l)-u(jm,l)-u(j,l+1)-u(j,lm)+4.d0*u(j,l)
    enddo
    r(j,ny1)=r(j,1)
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL
  r(nx1,1:ny1)=r(1,1:ny1)
  p(:,:)=r(:,:)
  !
  ! conjugate-gradient iterations
  do n=1,MAXITS
    rr=0.d0
    !$OMP PARALLEL DEFAULT(PRIVATE) &
    !$OMP SHARED(r,rr)
    !$OMP DO REDUCTION(+:rr)
    do j=1,nx
      do l=1,ny
        rr=rr+r(j,l)**2
      enddo
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL
    if(n==1) rr0=rr

    !$OMP PARALLEL DEFAULT(PRIVATE) &
    !$OMP SHARED(p,ap)
    !$OMP DO
    do j=1,nx
      jm=j-1
      if(j==1) jm=nx
      do l=1,ny
        lm=l-1
        if(l==1) lm=ny
        ap(j,l)=p(j+1,l)+p(jm,l)+p(j,l+1)+p(j,lm)-4.d0*p(j,l)
      enddo
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL

    pap=0.d0
    !$OMP PARALLEL DEFAULT(PRIVATE) &
    !$OMP SHARED(p,ap,pap)
    !$OMP DO REDUCTION(+:pap)
    do j=1,nx
      do l=1,ny
        pap=pap+p(j,l)*ap(j,l)
      enddo
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL

    ak=rr/pap

    !$OMP PARALLEL DEFAULT(PRIVATE) &
    !$OMP SHARED(r,p,ap,u,r1,ak)
    !$OMP DO
    do j=1,nx
      do l=1,ny
        u(j,l)=u(j,l)+ak*p(j,l)
        r1(j,l)=r(j,l)-ak*ap(j,l)
      enddo
      u(j,ny1)=u(j,1)
      r1(j,ny1)=r1(j,1)
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL
    u(nx1,1:ny1)=u(1,1:ny1)
    r1(nx1,1:ny1)=r1(1,1:ny1)

    r_norm=0.d0
    !$OMP PARALLEL DEFAULT(PRIVATE) &
    !$OMP SHARED(r1,r_norm)
    !$OMP DO REDUCTION(+:r_norm)
    do j=1,nx
      do l=1,ny
        r_norm=r_norm+r1(j,l)**2
      enddo
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL

    if(r_norm/nx/ny < EPS) then
      write(*,*) 'n=',n,'rr0=',rr0,' rr1=',r_norm
      return
    endif

    ! bk=0.d0
    !!$OMP PARALLEL DEFAULT(PRIVATE) &
    !!$OMP SHARED(r,bk)
    !!$OMP DO REDUCTION(+:bk)
    ! do j=1,nx
    ! do l=1,ny
    !   bk=bk+r(j,l)**2
    ! enddo
    ! enddo
    !!$OMP ENDDO
    !!$OMP END PARALLEL
    bk=r_norm/rr

    !$OMP PARALLEL DEFAULT(PRIVATE) &
    !$OMP SHARED(p,r1,bk)
    !$OMP DO
    do j=1,nx
      do l=1,ny
        p(j,l)=r1(j,l)+bk*p(j,l)
      enddo
      p(j,ny1)=p(j,1)
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL
    p(nx1,1:ny1)=p(1,1:ny1)

    r(:,:)=r1(:,:)
  enddo
  write(*,*) 'MAXITS exceeded in cgm'
  !
end subroutine cgm
