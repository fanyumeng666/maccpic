!#################
subroutine creator
!#################
  !
  ! subroutine creator loads particles in phase space (positions and
  ! momenta); the fields are set in init.
  !
  use input, only : rmass,rtemp
  use global, only : com_vel,nsp,npcel,maxnp,alx,aly,xn,yn,twopi,rd,parentwave,&
                     x,y,vx,vy,vz,gamma,np,csden,csw,gcsw,wce0,NPROCS,MYRANK,nthreads
  !!$ use omp_lib
  !
  implicit none
  ! ----- MPI -----
  include 'mpif.h'
  ! ---------------
  integer(8) i,j,ip,iq,nph,npback
  integer(4) isp,idum,nnx
  real(8) pxmax,pymax,dpx,dpy,ppx,ppy,fvx,fvy,dfx,dfy,ddx,ddy,va0,ve0,gamma0,gamma0a,aaa,dran2,vdp
  real(8) gam1,vvx,vvy,vvz,va0z,gamma0z,q
  !
  real(8), dimension(nsp) :: temp,va,ve,aist
  !
  call com_vel(va,ve,temp,aist)
  !
  do isp=1,nsp
    npback=nint(xn*yn)*npcel/NPROCS
    nph=np(isp)/2
    !
    ! momentum grid of nph points over +/- 5 thermal widths
    ! (aist = temperature anisotropy in x)
    pxmax=5.d0*sqrt(temp(isp)*(1.0d0+aist(isp)))
    dpx=2.d0*pxmax/dble(nph)
    pymax=5.d0*sqrt(temp(isp))
    dpy=2.d0*pymax/dble(nph)
    x(MYRANK*maxnp+1,isp)=0.d0
    y(MYRANK*maxnp+1,isp)=0.d0
    !
    ! save the indefinite integral of the cdf in the x array.
    !
    !!$ call omp_set_dynamic(.false.)
    !!$ call omp_set_num_threads(nthreads)
    !!$OMP PARALLEL DEFAULT(PRIVATE) &
    !!$OMP SHARED(isp,np,x,y,dp,pmax,temp)
    !!$OMP DO
    do i=1,nph
      ppx=(i-1.5d0)*dpx-pxmax
      fvx=exp(-(sqrt(1.d0+ppx*ppx)-1.d0)/temp(isp)/(1.0d0+aist(isp)))
      ppy=(i-1.5d0)*dpy-pymax
      fvy=exp(-(sqrt(1.d0+ppy*ppy)-1.d0)/temp(isp))
      x(MYRANK*maxnp+i+1,isp)=x(MYRANK*maxnp+i,isp)+max(fvx,0.d0)
      y(MYRANK*maxnp+i+1,isp)=y(MYRANK*maxnp+i,isp)+max(fvy,0.d0)
    enddo
    !!$OMP ENDDO
    !!$OMP END PARALLEL
    !
    ! load particles in momentum space according to the cdf.
    ! each particle is loaded with px,py=pz.
    !
    dfx=x(MYRANK*maxnp+nph+1,isp)/dble(nph)
    dfy=y(MYRANK*maxnp+nph+1,isp)/dble(nph)

    !!$OMP PARALLEL DEFAULT(PRIVATE) &
    !!$OMP SHARED(isp,np,df,x,vx,vy,vz,pmax,dp)
    !!$OMP DO
    ip=1; iq=1
    do i=1,nph
      fvx=(i-0.5d0)*dfx
      do j=1,nph
        if(fvx >= x(MYRANK*maxnp+iq+1,isp)) then
          iq=iq+1
          if(iq > nph) stop
        else
          exit
        endif
      enddo
      ppx=dpx*(iq-1+(fvx-x(MYRANK*maxnp+iq,isp))/(x(MYRANK*maxnp+iq+1,isp)-x(MYRANK*maxnp+iq,isp)))-pxmax
      vx(MYRANK*maxnp+ip,isp)=ppx

      ip=ip+1
    enddo

    ip=1; iq=1
    do i=1,nph
      fvy=(i-0.5d0)*dfy
      do j=1,nph
        if(fvy >= y(MYRANK*maxnp+iq+1,isp)) then
          iq=iq+1
          if(iq > nph) stop
        else
          exit
        endif
      enddo
      ppy=dpy*(iq-1+(fvy-y(MYRANK*maxnp+iq,isp))/(y(MYRANK*maxnp+iq+1,isp)-y(MYRANK*maxnp+iq,isp)))-pymax
      vy(MYRANK*maxnp+ip,isp)=ppy
      vz(MYRANK*maxnp+ip,isp)=ppy
      ip=ip+1
    enddo

    !!$OMP ENDDO
    !!$OMP END PARALLEL
    !
    ! load particles in configuration space uniformly,
    ! with perturbation if need be.
    ! Lower half (slots 1..nph/2): background on a uniform lattice (slots
    ! 1..npback/4), then the current sheet at y = yn/4 with density
    ! csden*sech^2((y-yn/4)/csw), loaded as rows of nnx particles spaced
    ! ddx/n(y) apart, first downwards and then upwards from the centre.
    !
    ddx=1.d0/sqrt(dble(npcel/2)); ddy=ddx
    nnx=int(xn/ddx/NPROCS)
    !$OMP PARALLEL DEFAULT(PRIVATE) &
    !$OMP SHARED(isp,nph,nnx,ddx,ddy,x,y)
    !$OMP DO
    do i=1,npback/4
      x(MYRANK*maxnp+i,isp)=mod((i-0.5d0)*ddx,xn/NPROCS) + xn*MYRANK/NPROCS
      y(MYRANK*maxnp+i,isp)=ddy*(i/nnx+0.5d0)
    enddo
    ! sheet particles below the centre (rows going down)
    ddy=0.25d0*yn+ddx*0.001
    !   do i=npback/4+1,nph/2
    do i=npback/4+1,nph/4+npback/8
      x(MYRANK*maxnp+i,isp)=mod((i-0.5d0)*ddx,xn/NPROCS) + xn*MYRANK/NPROCS
      y(MYRANK*maxnp+i,isp)=ddy
      if(mod(i-npback/4,nnx)==0) &
        ddy=y(MYRANK*maxnp+i-nnx+1,isp)-&
            ddx/csden/4*((exp((y(MYRANK*maxnp+i-nnx+1,isp)-0.25d0*yn)/gcsw)+&
                          exp((-y(MYRANK*maxnp+i-nnx+1,isp)+0.25d0*yn)/gcsw))**2)
      if(y(MYRANK*maxnp+i,isp)<0.d0) y(MYRANK*maxnp+i,isp)=y(MYRANK*maxnp+i,isp)+yn/2
    enddo
    ! sheet particles above the centre (rows going up)
    ddy=0.25d0*yn+ddx/csden
    do i=nph/4+npback/8+1,nph/2
      x(MYRANK*maxnp+i,isp)=mod((i-0.5d0)*ddx,xn/NPROCS) + xn*MYRANK/NPROCS
      y(MYRANK*maxnp+i,isp)=ddy
      if(mod(i-nph/4-npback/8,nnx)==0) &
        ddy=y(MYRANK*maxnp+i-nnx+1,isp)+&
            ddx/csden/4*((exp((y(MYRANK*maxnp+i-nnx+1,isp)-0.25d0*yn)/gcsw)+&
                          exp((-y(MYRANK*maxnp+i-nnx+1,isp)+0.25d0*yn)/gcsw))**2)
      if(y(MYRANK*maxnp+i,isp)>yn/2) y(MYRANK*maxnp+i,isp)=y(MYRANK*maxnp+i,isp)-yn/2
    enddo

    ! (alternative loading with the background mixed in; not used)
    !   do i=1,nph/2
    !     x(MYRANK*maxnp+i,isp)=mod((i-0.5d0)*ddx,xn/NPROCS) + xn*MYRANK/NPROCS
    !!     y(MYRANK*maxnp+i,isp)=((i-1)/nnx+0.5d0)*ddy
    !!     x(MYRANK*maxnp+i,isp)=(i-0.5d0)*xn/nph
    !!     y(MYRANK*maxnp+i,isp)=(i-0.5d0)*yn/nph
    !     y(MYRANK*maxnp+i,isp)=ddy
    !     if(y(MYRANK*maxnp+i,isp)==ddx) y(MYRANK*maxnp+i,isp)=y(MYRANK*maxnp+i,isp)*0.5d0
    !     if(mod(MYRANK*maxnp+i,nnx)==0) then
    !       ddy=y(MYRANK*maxnp+i-nnx+1,isp)+&
    !           ddx/(csden*(2/(exp((y(MYRANK*maxnp+i-nnx+1,isp)-0.25d0*yn)/gcsw)+&
    !                          exp(-(y(MYRANK*maxnp+i-nnx+1,isp)-0.25d0*yn)/gcsw)))**2+csden*0.2d0)
    !     end if
    !   enddo
    !!$OMP ENDDO
    !!$OMP END PARALLEL
    ! upper half (slots nph/2+1..nph): the lower half shifted by yn/2,
    ! i.e. the second sheet at y = 3yn/4
    !!$OMP PARALLEL DEFAULT(PRIVATE) &
    !!$OMP SHARED(isp,nph,x,y)
    !!$OMP DO
    do i=nph/2+1,nph
      x(MYRANK*maxnp+i,isp)=x(MYRANK*maxnp+i-nph/2,isp)
      y(MYRANK*maxnp+i,isp)=y(MYRANK*maxnp+i-nph/2,isp)+yn*0.5
    end do
    !!$OMP ENDDO
    !!$OMP END PARALLEL
    !
  enddo
  !
  ! randomize the particle momenta and positions so as
  ! to decorrelate phase space.
  ! (once for all species; smear shuffles only the velocities of slots
  ! 1..np/2)
  !
  call smear(10)

  ! second pass: mirror the momenta, then add the current-sheet drift
  do isp=1,nsp
    !
    ! second half of the slots: same positions, opposite momenta
    ! (zero net momentum)
    !$OMP PARALLEL DEFAULT(PRIVATE) &
    !$OMP SHARED(isp,MYRANK,nph,x,y,vx,vy,vz)
    !$OMP DO
    do i=1,nph
      x(MYRANK*maxnp+i+nph,isp)=  x(MYRANK*maxnp+i,isp)
      y(MYRANK*maxnp+i+nph,isp)=  y(MYRANK*maxnp+i,isp)
      vx(MYRANK*maxnp+i+nph,isp)=-vx(MYRANK*maxnp+i,isp)
      vy(MYRANK*maxnp+i+nph,isp)=-vy(MYRANK*maxnp+i,isp)
      vz(MYRANK*maxnp+i+nph,isp)=-vz(MYRANK*maxnp+i,isp)
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL

    ! drift along z: none for the background slots; sheet particles get
    ! va0 = |wce0|/(csw*csden) (the Harris current), with opposite signs for
    ! the two sheets and the two species, shared as Te:Ti (electrons
    ! va0/(1+rtemp), ions va0*rtemp/(1+rtemp)). Then momentum -> velocity
    ! with relativistic addition of vdp.
    do i=1,np(isp)
      if(i<=npback/4) then
        va0=0.d0
      elseif((i>nph/2).and.(i<=nph/2+npback/4)) then
        va0=0.d0
      elseif((i>nph).and.(i<=nph+npback/4)) then
        va0=0.d0
      elseif((i>3*nph/2).and.(i<=3*nph/2+npback/4)) then
        va0=0.d0
      elseif(y(MYRANK*maxnp+i,isp)>yn/2) then
        !    if((y(MYRANK*maxnp+i,isp)>yn/2).and.(mod(y(MYRANK*maxnp+i,isp)-0.5*ddx,ddx)<ddx*1.d-3)) then
        !       q=(2/(exp((y(MYRANK*maxnp+i,isp)-0.75d0*yn)/gcsw)+&
        !             exp(-(y(MYRANK*maxnp+i,isp)-0.75d0*yn)/gcsw)))**2
        !!        va0=(-1)**(isp+1)*0.5d0*abs(wce0)/csw*q/(csden*q+1)
        !       va0=(-1)**(isp+1)*0.5d0*abs(wce0)/csw/csden
        va0=(-1)**(isp+1)*abs(wce0)/csw/csden
      elseif(y(MYRANK*maxnp+i,isp)<yn/2) then
        !    elseif((y(MYRANK*maxnp+i,isp)<yn/2).and.(mod(y(MYRANK*maxnp+i,isp)-0.5*ddx,ddx)<ddx*1.d-3)) then
        !       q=(2/(exp((y(MYRANK*maxnp+i,isp)-0.25d0*yn)/gcsw)+&
        !             exp(-(y(MYRANK*maxnp+i,isp)-0.25d0*yn)/gcsw)))**2
        !!        va0=(-1)**(isp)*0.5d0*abs(wce0)/csw*q/(csden*q+1)
        !       va0=(-1)**(isp)*0.5d0*abs(wce0)/csw/csden
        va0=(-1)**(isp)*abs(wce0)/csw/csden
      endif
      vdp=(rtemp)**(isp-1)/(rtemp+1)*va0
      gamma0a=sqrt(1.d0-vdp*vdp)
      gam1=sqrt(1.d0+vx(MYRANK*maxnp+i,isp)**2+vy(MYRANK*maxnp+i,isp)**2+vz(MYRANK*maxnp+i,isp)**2)
      vvx=vx(MYRANK*maxnp+i,isp)/gam1
      vvy=vy(MYRANK*maxnp+i,isp)/gam1
      vvz=vz(MYRANK*maxnp+i,isp)/gam1
      vz(MYRANK*maxnp+i,isp)=(vdp+vvz)/(1.d0+vdp*vvz)
      vy(MYRANK*maxnp+i,isp)=vvy*gamma0a/(1.d0+vdp*vvz)
      vx(MYRANK*maxnp+i,isp)=vvx*gamma0a/(1.d0+vdp*vvz)
    enddo

    !   do i=1,np(isp)
    !     x(MYRANK*maxnp+i,isp) = x(MYRANK*maxnp+i,isp)+0.1d0*ddx*(dran2(idum)-0.5d0)
    !     y(MYRANK*maxnp+i,isp) = y(MYRANK*maxnp+i,isp)+0.1d0*ddx*(dran2(idum)-0.5d0)
    !   enddo
    !
  enddo
  !
end subroutine creator
