!###############
subroutine init
!###############
  !
  ! Allocate this rank's particle arrays (slots MYRANK*maxnp+1 ...),
  ! zero all fields, set the initial magnetic field (two Harris sheets
  ! plus an X-point perturbation) and the species constants.
  !
  use global, only : com_den,nx1,ny1,nx2,ny2,nsp,maxnp,np,nproc,bx,by,bz,ex,ey,ez,&
                     jxm,jxp,jym,jyp,jzm,jzp,rho,c,dt,wce0,rd,gcsw,ny,dy,pi,csden,nx,dx,&
                     x,y,vx,vy,vz,gamma,qme,qmi,am,tx,ty,anme,IERR,MYRANK,nthreads,rmass
  !$ use omp_lib
  !
  implicit none
  ! ----- MPI -----
  include 'mpif.h'
  ! ---------------
  integer(8) j,i1,i2
  integer(4) isp,ix,iy
  real(8), dimension(nsp) :: dn,q,qm
  real(8) dbx,dby
  character(len=3) :: fnbx,fnby
  !
  write(*,*) 'number of particles per grid=',maxnp*nproc
  i1=MYRANK*maxnp+1; i2=(MYRANK+1)*maxnp
  allocate(x(i1:i2,nsp),y(i1:i2,nsp),vx(i1:i2,nsp),vy(i1:i2,nsp),vz(i1:i2,nsp),gamma(i1:i2,nsp),STAT=IERR)
  !
  do isp=1,nsp
    !$ call omp_set_dynamic(.false.)
    !$ call omp_set_num_threads(nthreads)
    !$OMP PARALLEL PRIVATE(j) &
    !$OMP SHARED(isp,i1,i2,x,y,vx,vy,vz,gamma)
    !$OMP DO
    do j=i1,i2
      x(j,isp)=0.d0
      y(j,isp)=0.d0
      vx(j,isp)=0.d0
      vy(j,isp)=0.d0
      vz(j,isp)=0.d0
      gamma(j,isp)=1.d0
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL
  enddo
  !
  !$OMP PARALLEL PRIVATE(ix,iy) &
  !$OMP SHARED(bx,by,bz,ex,ey,ez,rho,jxm,jxp,jym,jyp,jzm,jzp)
  !$OMP DO
  do ix=0,nx2
    do iy=0,ny2
      bx(ix,iy) =0.d0
      by(ix,iy) =0.d0
      bz(ix,iy) =0.d0
      ex(ix,iy) =0.d0
      ey(ix,iy) =0.d0
      ez(ix,iy) =0.d0
      rho(ix,iy)=0.d0
      jxm(ix,iy)=0.d0
      jxp(ix,iy)=0.d0
      jym(ix,iy)=0.d0
      jyp(ix,iy)=0.d0
      jzm(ix,iy)=0.d0
      jzp(ix,iy)=0.d0
    enddo
    ! Bx = -/+ B0 tanh(y/csw) for the sheets at y = ny/4 and 3ny/4, plus a
    ! perturbation of amplitude dbx, dby (one wavelength in x and in y).
    ! Each component is evaluated at its own position on the staggered grid
    ! (in cells): Bx(ix,iy) at (ix-1/2, iy-1), By(ix,iy) at (ix-1, iy-1/2),
    ! so that Bx = 0 at y = ny/4 and 3ny/4, where creator centres the
    ! particle sheets. The original lines (commented) used (ix, iy) as the
    ! position, which put the field reversal one cell below the sheets.
    dbx=0.1*abs(wce0)/sqrt(csden)*sqrt(rmass)*2*pi/(nx*dx)
    dby=0.1*abs(wce0)/sqrt(csden)*sqrt(rmass)*2*pi/(ny*dy)
    do iy=0,ny/2
      !      bx(ix,iy+ny/2+1)=abs(wce0)*tanh((iy+ny/2+1-3*ny/4)/gcsw)+&
      !         dby*cos(2*pi*(ix-nx*0.5d0)/nx+pi)*sin(2*pi*(iy+ny/2+1-3*ny/4)/ny)
      !      bx(ix,iy)=-abs(wce0)*tanh((iy-ny/4)/gcsw)-&
      !         dby*cos(2*pi*(ix-nx*0.5d0)/nx)*sin(2*pi*(iy-0.25d0*ny)/ny)
      !      by(ix,iy+ny/2+1)=-dbx*sin(2*pi*(ix-nx*0.5d0)/nx+pi)*cos(2*pi*(iy+ny/2+1-3*ny/4)/ny)
      !      by(ix,iy       ) =dbx*sin(2*pi*(ix-nx*0.5d0)/nx)*cos(2*pi*(iy       -0.25d0*ny)/ny)
      bx(ix,iy+ny/2+1)=abs(wce0)*tanh((iy+ny/2-3*ny/4)/gcsw)+&
         dby*cos(2*pi*(ix-0.5d0-nx*0.5d0)/nx+pi)*sin(2*pi*(iy+ny/2-3*ny/4)/ny)

      bx(ix,iy)=-abs(wce0)*tanh((iy-1-ny/4)/gcsw)-&
         dby*cos(2*pi*(ix-0.5d0-nx*0.5d0)/nx)*sin(2*pi*(iy-1-0.25d0*ny)/ny)

      by(ix,iy+ny/2+1)=-dbx*sin(2*pi*(ix-1-nx*0.5d0)/nx+pi)*cos(2*pi*(iy+ny/2+0.5d0-3*ny/4)/ny)
      by(ix,iy       ) =dbx*sin(2*pi*(ix-1-nx*0.5d0)/nx)*cos(2*pi*(iy-0.5d0-0.25d0*ny)/ny)
    enddo
    bx(ix,ny2)=bx(ix,0)
    by(ix,ny2)=by(ix,0)
  enddo
  !initial b-f
  ! write the initial Bx and By to the files bx0 and by0
  write(fnbx,'(a3)') 'bx0'
  open(34,file=fnbx,form='formatted',status='unknown')
  do iy=0,ny+1
    write(34,'(1024e13.5)') (bx(ix,iy),ix=0,nx+1)
  enddo
  close(34)

  write(fnby,'(a3)') 'by0'
  open(35,file=fnby,form='formatted',status='unknown')
  do iy=0,ny+1
    write(35,'(1024e13.5)') (by(ix,iy),ix=0,nx+1)
  enddo
  close(35)
  !
  !$OMP ENDDO
  !$OMP END PARALLEL
  !
  call com_den(dn,q,qm)
  !
  ! every rank uses all of its slots; am = macro-particle mass;
  ! tx, ty (background field rotation) are switched off by the *0.d0
  np(1:nsp)=maxnp
  am(1)=q(1)/qme;   tx(1)=wce0*cos(rd)*dt/2.d0/c*0.d0;      ty(1)=wce0*sin(rd)*dt/2.d0/c*0.d0
  am(2)=q(2)/qmi;   tx(2)=-am(1)*tx(1)/am(2)*0.d0;          ty(2)=-am(1)*ty(1)/am(2)*0.d0
  anme=np(1)*am(1)*nproc
  if(nsp > 2) then
    do isp=3,nsp
      if(q(isp) < 0.d0) then
        am(isp)=am(1);  tx(isp)=tx(1);  ty(isp)=ty(1)
      else
        am(isp)=am(2);  tx(isp)=tx(2);  ty(isp)=ty(2)
      endif
    enddo
  endif
  !
end subroutine init
