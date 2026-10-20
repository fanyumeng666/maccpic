!###############
subroutine cont
!###############
  !
  ! Restart: read the fields (fld.fin), rho and J (rhoj.fin) and this
  ! rank's particles (*_rankNNN.fin) written by output at the last step
  ! of the previous run, then advance the fields once.
  !
  use global, only : nx,ny,nx1,ny1,nsp,maxnp,dt,bx,by,bz,ex,ey,ez,time,it0,&
                     jxm,jxp,jym,jyp,jzm,jzp,rho,x,y,vx,vy,vz,gamma,np,nt,it,MYRANK
  !
  implicit none
  ! ----- MPI -----
  include 'mpif.h'
  ! ---------------
  integer(4) isp
  integer(8) i
  character(len=13) fnxfin,fnyfin,fngfin
  character(len=14) fnvxfin,fnvyfin,fnvzfin
  !
  open(80,file='fld.fin',form='unformatted')
  read(80) bx
  read(80) by
  read(80) bz
  read(80) ex
  read(80) ey
  read(80) ez
  close(80)
  !
  open(81,file='rhoj.fin',form='unformatted')
  read(81) rho
  read(81) jxm
  read(81) jxp
  read(81) jym
  read(81) jyp
  read(81) jzm
  read(81) jzp
  close(81)
  !
  write(fnxfin,'(a6,i3.3,a4)') 'x_rank',MYRANK,'.fin'
  open(82,file=fnxfin,form='unformatted',status='unknown')
  do isp=1,nsp
    do i=MYRANK*maxnp+1,(MYRANK+1)*maxnp
      read(82) x(i,isp)
    enddo
  enddo
  close(82)
  !
  write(fnyfin,'(a6,i3.3,a4)') 'y_rank',MYRANK,'.fin'
  open(83,file=fnyfin,form='unformatted',status='unknown')
  do isp=1,nsp
    do i=MYRANK*maxnp+1,(MYRANK+1)*maxnp
      read(83) y(i,isp)
    enddo
  enddo
  close(83)
  !
  write(fnvxfin,'(a7,i3.3,a4)') 'vx_rank',MYRANK,'.fin'
  open(84,file=fnvxfin,form='unformatted',status='unknown')
  do isp=1,nsp
    do i=MYRANK*maxnp+1,(MYRANK+1)*maxnp
      read(84) vx(i,isp)
    enddo
  enddo
  close(84)
  !
  write(fnvyfin,'(a7,i3.3,a4)') 'vy_rank',MYRANK,'.fin'
  open(85,file=fnvyfin,form='unformatted',status='unknown')
  do isp=1,nsp
    do i=MYRANK*maxnp+1,(MYRANK+1)*maxnp
      read(85) vy(i,isp)
    enddo
  enddo
  close(85)
  !
  write(fnvzfin,'(a7,i3.3,a4)') 'vz_rank',MYRANK,'.fin'
  open(86,file=fnvzfin,form='unformatted',status='unknown')
  do isp=1,nsp
    do i=MYRANK*maxnp+1,(MYRANK+1)*maxnp
      read(86) vz(i,isp)
    enddo
  enddo
  close(86)
  !
  write(fngfin,'(a6,i3.3,a4)') 'g_rank',MYRANK,'.fin'
  open(87,file=fngfin,form='unformatted',status='unknown')
  do isp=1,nsp
    do i=MYRANK*maxnp+1,(MYRANK+1)*maxnp
      read(87) gamma(i,isp)
    enddo
  enddo
  close(87)
  !
  it=it0+1
  time=it*dt
  !===========
  call fields
  !===========
  !
end subroutine cont
