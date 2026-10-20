!################
subroutine smooth(sm)
!################
  !
  ! 1-2-1 filter in x, then in y, for an (nx1,ny1) periodic array
  ! (not called at present; fields uses smooth2)
  !
  use global, only : nx,ny,nx1,ny1
  implicit none
  integer(4) ix,iy
  !
  real(8), dimension(nx1,ny1) :: sm
  real(8), dimension(nx1)     :: dumx
  real(8), dimension(ny1)     :: dumy
  !
  ! -----  3 point filter in x -----
  !
  do iy=1,ny
    dumx(1:nx1)=sm(1:nx1,iy)
    do ix=2,nx
      sm(ix,iy)=(dumx(ix)+0.5d0*(dumx(ix-1)+dumx(ix+1)))*0.5d0
    enddo
    sm(  1,iy)=(dumx(1)+0.5d0*(dumx(nx)+dumx(2)))*0.5d0
    sm(nx1,iy)=sm(1,iy)
  enddo
  sm(1:nx1,ny1)=sm(1:nx1,1)
  !
  ! -----  3 point filter in y -----
  !
  do ix=1,nx
    dumy(1:ny1)=sm(ix,1:ny1)
    do iy=2,ny
      sm(ix,iy)=(dumy(iy)+0.5d0*(dumy(iy-1)+dumy(iy+1)))*0.5d0
    enddo
    sm(ix,  1)=(dumy(1)+0.5d0*(dumy(ny)+dumy(2)))*0.5d0
    sm(ix,ny1)=sm(ix,1)
  enddo
  sm(nx1,1:ny1)=sm(1,1:ny1)
  !!
  !! ----- 5 point filter -----
  !!
  ! do j=1,ng
  !   jm2=j-2; jm1=j-1
  !   jp2=j+2; jp1=j+1
  !   if(jm2 <  1) jm2=jm2+ng
  !   if(jm1 <  1) jm1=jm1+ng
  !   if(jp2 > ng) jp2=jp2-ng
  !   sm(j)=(dum(jm2)+dum(jp2)+4*(dum(jm1)+dum(jp1))+6*dum(j))/16
  ! enddo
  ! sm(ng1)=sm(1)
  !!
  !! ----- 7 point filter -----
  !!
  ! do j=1,ng
  !   jm3=j-3; jm2=j-2; jm1=j-1
  !   jp3=j+3; jp2=j+2; jp1=j+1
  !   if(jm3 <  1) jm3=jm3+ng
  !   if(jm2 <  1) jm2=jm2+ng
  !   if(jm1 <  1) jm1=jm1+ng
  !   if(jp3 > ng) jp3=jp3-ng
  !   if(jp2 > ng) jp2=jp2-ng
  !   sm(j)=(dum(jm3)+dum(jp3)+6*(dum(jm2)+dum(jp2))+15*(dum(jm1)+dum(jp1))+20*dum(j))/64
  ! enddo
  ! sm(ng1)=sm(1)
  !!
  !! ----- 9 point filter -----
  !!
  ! do j=1,ng
  !   jm4=j-4; jm3=j-3; jm2=j-2; jm1=j-1
  !   jp4=j+4; jp3=j+3; jp2=j+2; jp1=j+1
  !   if(jm4 <  1) jm4=jm4+ng
  !   if(jm3 <  1) jm3=jm3+ng
  !   if(jm2 <  1) jm2=jm2+ng
  !   if(jm1 <  1) jm1=jm1+ng
  !   if(jp4 > ng) jp4=jp4-ng
  !   if(jp3 > ng) jp3=jp3-ng
  !   if(jp2 > ng) jp2=jp2-ng
  !   sm(j)=(dum(jm4)+dum(jp4)+8*(dum(jm3)+dum(jp3))+28*(dum(jm2)+dum(jp2))+&
  !          56*(dum(jm1)+dum(jp1))+70*dum(j))/256
  ! enddo
  ! sm(ng1)=sm(1)
  !!
  !! ----- compensation -----
  !!
  ! dum(1:ng1)=sm(1:ng1)
  ! do j=2,ng
  !   sm(j)=1.5d0*dum(j)-(dum(j-1)+dum(j+1))/4.d0
  ! enddo
  ! sm(1)=1.5d0*dum(1)-(dum(ng)+dum(2))/4.d0
  ! sm(ng1)=sm(1)
  !
end subroutine smooth
!
!################
subroutine smooth2(sm)
!################
  !
  ! 9-point (1-2-1 x 1-2-1) binomial filter followed by a compensation
  ! step, for a (0:nx2,0:ny2) periodic array; used in fields
  !
  use global, only : nx,ny,nx1,ny1,nx2,ny2,nthreads
  !$ use omp_lib
  implicit none
  integer(4) ix,iy
  !
  real(8), dimension(0:nx2,0:ny2) :: sm,dummy
  ! real(8), dimension(0:nx2)     :: dumx
  ! real(8), dimension(0:ny2)     :: dumy
  !!
  !! -----  3 point filter in x -----
  !!
  ! do iy=1,ny1
  !   dumx(0:nx2)=sm(0:nx2,iy)
  !   do ix=1,nx1
  !     sm(ix,iy)=(dumx(ix)+0.5d0*(dumx(ix-1)+dumx(ix+1)))*0.5d0
  !   enddo
  !   sm(  0,iy)=sm(nx,iy)
  !   sm(nx2,iy)=sm( 2,iy)
  ! enddo
  ! sm(0:nx2,0)=sm(0:nx2,ny)
  ! sm(0:nx2,ny2)=sm(0:nx2,2)
  !!
  !! -----  3 point filter in y -----
  !!
  ! do ix=1,nx1
  !   dumy(0:ny2)=sm(ix,0:ny2)
  !   do iy=1,ny1
  !     sm(ix,iy)=(dumy(iy)+0.5d0*(dumy(iy-1)+dumy(iy+1)))*0.5d0
  !   enddo
  !   sm(ix,  0)=sm(ix,ny)
  !   sm(ix,ny2)=sm(ix, 2)
  ! enddo
  ! sm(0,0:ny2)=sm(nx,0:ny2)
  ! sm(nx2,0:ny2)=sm(2,0:ny2)
  !
  ! -----  9 point filter  -----
  !
  !$ call omp_set_dynamic(.false.)
  !$ call omp_set_num_threads(nthreads)
  !$OMP PARALLEL DEFAULT(PRIVATE) &
  !$OMP SHARED(dummy,sm)
  !$OMP DO
  do iy=1,ny1
    do ix=1,nx1
      dummy(ix,iy)=(4*sm(ix,iy)+2*(sm(ix-1,iy)+sm(ix+1,iy)+sm(ix,iy-1)+sm(ix,iy+1))+&
                    sm(ix-1,iy-1)+sm(ix-1,iy+1)+sm(ix+1,iy-1)+sm(ix+1,iy+1))/16.d0
    enddo
    dummy(0,iy)=dummy(nx,iy); dummy(nx2,iy)=dummy(2,iy)
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL
  ! dummy(0,1:ny1)=dummy(nx,1:ny1)
  ! dummy(nx2,1:ny1)=dummy(2,1:ny1)
  dummy(0:nx2,0)=dummy(0:nx2,ny)
  dummy(0:nx2,ny2)=dummy(0:nx2,2)
  !
  ! -----  compensation  -----
  !
  ! dummy=sm
  !$OMP PARALLEL DEFAULT(PRIVATE) &
  !$OMP SHARED(dummy,sm)
  !$OMP DO
  do iy=1,ny1
    do ix=1,nx1
      sm(ix,iy)=(20*dummy(ix,iy)-dummy(ix-1,iy)-dummy(ix+1,iy)-dummy(ix,iy-1)-dummy(ix,iy+1)-&
                    dummy(ix-1,iy-1)-dummy(ix-1,iy+1)-dummy(ix+1,iy-1)-dummy(ix+1,iy+1))/12.d0
      !   sm(ix,iy)=(36*dummy(ix,iy)-6*(dummy(ix-1,iy)+dummy(ix+1,iy)+dummy(ix,iy-1)+dummy(ix,iy+1))+&
      !                 dummy(ix-1,iy-1)+dummy(ix-1,iy+1)+dummy(ix+1,iy-1)+dummy(ix+1,iy+1))/16.d0
    enddo
    sm(0,iy)=sm(nx,iy); sm(nx2,iy)=sm(2,iy)
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL
  ! sm(0,1:ny1)=sm(nx,1:ny1)
  ! sm(nx2,1:ny1)=sm(2,1:ny1)
  sm(0:nx2,0)=sm(0:nx2,ny)
  sm(0:nx2,ny2)=sm(0:nx2,2)
  !
end subroutine smooth2
