!#################
subroutine setrho
!#################
  !
  !  subroutine setrho converts position to computer normalization and
  !  accumulates charge density at t=0.
  !  (linear weighting; not called at present: the call in main is
  !  commented out, and curnt deposits rho every step)
  !
  use global, only : com_den,nx,ny,nx1,ny1,nx2,ny2,nsp,maxnp,dx,dy,&
                     jxm,jxp,jym,jyp,jzm,jzp,rho,x,y,time,np,nthreads,IERR,MYRANK
  !$ use omp_lib
  implicit none
  ! ----- MPI -----
  include 'mpif.h'
  real(8), dimension(0:nx2,0:ny2) :: rho0
  ! ---------------
  integer(8) i
  integer(4) isp,ij,jj,j
  real(8) qdxdy,xn,yn,delrs,dells,delus,delds,qdxdn
  !
  real(8), dimension(nsp) :: dn,q,qm
  !
  call com_den(dn,q,qm)
  !
  do isp=1,nsp
    qdxdy=q(isp)/(dx*dy)*dn(isp)
    !
    !  linear weighting.
    !
    !$ call omp_set_dynamic(.false.)
    !$ call omp_set_num_threads(nthreads)
    !$OMP PARALLEL DEFAULT(PRIVATE) &
    !$OMP SHARED(isp,np,xn,yn,qdxdn,x,y,rho,MYRANK)
    !$OMP DO
    do i=MYRANK*maxnp+1,MYRANK*maxnp+np(isp)
      ij=x(i,isp); jj=y(i,isp)
      delrs=x(i,isp)-ij; dells=1.d0-delrs
      delus=y(i,isp)-jj; delds=1.d0-delus
      rho(ij+1,jj+1)=rho(ij+1,jj+1)+dells*delds*qdxdn
      rho(ij+1,jj+2)=rho(ij+1,jj+2)+dells*delus*qdxdn
      rho(ij+2,jj+1)=rho(ij+2,jj+1)+delrs*delds*qdxdn
      rho(ij+2,jj+2)=rho(ij+2,jj+2)+delrs*delus*qdxdn
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL
  enddo
  !
  rho(1,1:ny1)=rho(1,1:ny1)+rho(nx1,1:ny1)
  rho(nx1,1:ny1)=rho(1,1:ny1)
  rho(1:nx1,1)=rho(1:nx1,1)+rho(1:nx1,ny1)
  rho(1:nx1,ny1)=rho(1:nx1,1)
  !
  ! ----- MPI -----
  call MPI_ALLREDUCE(rho,rho0,(nx2+1)*(ny2+1),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
  rho=rho0
  ! ---------------
end subroutine setrho
