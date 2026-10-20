!#################
subroutine output
!#################
  !
  ! Diagnostics and restart files (called from main every isav and iphs
  ! steps, at step 1 and at the last step). Files written:
  !   ene###.hst               field and kinetic energies vs time (rank 0)
  !   xvelc###.hst, xvion###.hst
  !                            x, y, vx, vy, vz of the first 100 electrons
  !                            on rank 1 / the first 100 ions on rank 2
  !   ex_TTTTT.dat ... bz_TTTTT.dat
  !                            2D field snapshots every iphs steps
  !                            (TTTTT = time in 1/wpe)
  !   xp<isp>_TTTTT_rank###.datap
  !                            particle snapshots every 4000 steps
  !   ex###.hst ... bz###.hst  fields along y at x = nx/2, every isav steps
  !   cf-p1_###.bin            Bz every 100 steps (rank 0), frames_per_file
  !                            frames per file; frame it/100 goes to file
  !                            it/100/frames_per_file+1, also in a
  !                            continuation run
  !   fld.fin, rhoj.fin, *_rank###.fin
  !                            restart files at the last step (read by cont)
  !
  use global, only : trace,nx,ny,nx1,ny1,nsp,maxnp,isav,iphs,idata,dx,dy,ake,&
                     bx,by,bz,ex,ey,ez,jxm,jxp,jym,jyp,jzm,jzp,rho,np,nt,it,it0,&
                     exe,eye,eze,bxe,bye,bze,anme,time,x,y,vx,vy,vz,gamma,NPROCS,MYRANK
  implicit none
  ! ----- MPI -----
  include 'mpif.h'
  ! integer(4) ISTA,IEND,irnk
  ! ---------------
  integer(8) i
  integer(4) isp,iii,itime,ij,ix,iy,ifile,irk,irk1,irk2,irk3,irk4,irk5
  real(8) ax,avx,avy,avz,agam,wem,total
  real(8) aex,aey,aez,abx,aby,abz
  !
  character(len=10) fcn1,fcn2,fcn3,fcn4
  character(len=10) fnene,fneip,fniip
  character(len=12) fnelc,fnion
  character(len=17) fnex,fney,fnez,fnbx,fnby,fnbz
  character(len=23) fnptcl
  character(len=12) fssex,fssey,fssez,fssbx,fssby,fssbz
  character(len=13) fnxfin,fnyfin,fngfin
  character(len=14) fnvxfin,fnvyfin,fnvzfin
  ! cf-p1 files: open_index = number of the file now open on unit 133 (0: none)
  integer(4), save :: open_index = 0
  integer(4), parameter :: frames_per_file = 30
  integer(4) iframe,findex
  integer(8) ipos
  character(len=13) :: fcn5
  !
  ! open files
  !
  irk=0
  if(((it0 == 0).and.(it == 0)).or.((it0 /= 0).and.(it == it0+isav))) then
    ifile=(it0+isav)/(nt-it0)+1
    !   irk=irk+1
    if(MYRANK==mod(0,NPROCS)) then
      write(fnene,'(a3,i3.3,a4)') 'ene',ifile,'.hst'
      open(10,file=fnene,status='unknown')
    elseif(MYRANK==mod(1,NPROCS)) then
      write(fnelc,'(a5,i3.3,a4)') 'xvelc',ifile,'.hst'
      open(70,file=fnelc,status='unknown',form='formatted')
    elseif(MYRANK==mod(2,NPROCS)) then
      write(fnion,'(a5,i3.3,a4)') 'xvion',ifile,'.hst'
      open(71,file=fnion,status='unknown',form='formatted')
    endif
  endif
  !
  if(mod(it,iphs) == 0) then
    !
    ! snap shot

    itime=time

    ! particle snapshot every 4000 steps (each rank writes its own slots)
    if (mod(it,4000) == 0) then
      do isp=1,nsp
        write(fnptcl,'(a2,i1.1,a1,i5.5,a5,i3.3,a6)') 'xp',isp,'_',itime,'_rank',MYRANK,'.datap'
        open(30,file=fnptcl,form='formatted',status='unknown')
        do i=maxnp*MYRANK+1,maxnp*MYRANK+np(isp)
          write(30,'(5e13.5)') x(i,isp)*dx,y(i,isp)*dy,vx(i,isp),vy(i,isp),vz(i,isp)
        enddo
        close(30)
      enddo
    end if

    ! 2D field snapshots, one component per rank (ranks 3..8 mod NPROCS)
    irk=3
    if(MYRANK==mod(irk,NPROCS)) then
      write(fssex,'(a3,i5.5,a4)') 'ex_',itime,'.dat'
      open(31,file=fssex,form='formatted',status='unknown')
      do iy=1,ny
        write(31,'(1024e13.5)') (ex(ix,iy),ix=1,nx)
      enddo
      close(31)
    endif
    irk=irk+1
    if(MYRANK==mod(irk,NPROCS)) then
      write(fssey,'(a3,i5.5,a4)') 'ey_',itime,'.dat'
      open(32,file=fssey,form='formatted',status='unknown')
      do iy=1,ny
        write(32,'(1024e13.5)') (ey(ix,iy),ix=1,nx)
      enddo
      close(32)
    endif
    irk=irk+1
    if(MYRANK==mod(irk,NPROCS)) then
      write(fssez,'(a3,i5.5,a4)') 'ez_',itime,'.dat'
      open(33,file=fssez,form='formatted',status='unknown')
      do iy=1,ny
        write(33,'(1024e13.5)') (ez(ix,iy),ix=1,nx)
      enddo
      close(33)
    endif
    irk=irk+1
    if(MYRANK==mod(irk,NPROCS)) then
      write(fssbx,'(a3,i5.5,a4)') 'bx_',itime,'.dat'
      open(34,file=fssbx,form='formatted',status='unknown')
      do iy=1,ny
        write(34,'(1024e13.5)') (bx(ix,iy),ix=1,nx)
      enddo
      close(34)
    endif
    irk=irk+1
    if(MYRANK==mod(irk,NPROCS)) then
      write(fssby,'(a3,i5.5,a4)') 'by_',itime,'.dat'
      open(35,file=fssby,form='formatted',status='unknown')
      do iy=1,ny
        write(35,'(1024e13.5)') (by(ix,iy),ix=1,nx)
      enddo
      close(35)
    endif
    irk=irk+1
    if(MYRANK==mod(irk,NPROCS)) then
      write(fssbz,'(a3,i5.5,a4)') 'bz_',itime,'.dat'
      open(36,file=fssbz,form='formatted',status='unknown')
      do iy=1,ny
        write(36,'(1024e13.5)') (bz(ix,iy),ix=1,nx)
      enddo
      close(36)
    endif
    !itime
  endif
  !
  ! p1
  ! Bz frame every 100 steps (rank 0). The file number and the position in
  ! the file follow from the step number (frame = it/100), so a continuation
  ! run writes on after the frames of the previous run instead of
  ! overwriting them.
  if ((mod(it,100) == 0) .and. (MYRANK == 0)) then
    iframe = it/100
    findex = iframe/frames_per_file + 1
    if (findex /= open_index) then
      if (open_index > 0) close(133)
      write(fcn5,'("cf-p1_",i3.3,".bin")') findex
      if (mod(iframe,frames_per_file) == 0) then
        ! first frame of this file
        open(133, file=fcn5, form='unformatted', access='stream', status='replace')
      else
        ! file begun by the previous run: keep its frames
        open(133, file=fcn5, form='unformatted', access='stream', status='unknown')
      endif
      open_index = findex
      write(*,*) 'Opened file: ', trim(fcn5)
    endif
    ipos = int(mod(iframe,frames_per_file),8)*int(size(bz),8)*8_8 + 1_8
    write(133, pos=ipos) bz
  endif

  if ((it == nt) .and. (MYRANK == 0) .and. (open_index > 0)) then
    close(133)
    open_index = 0
    write(*,*) 'Closed final output file.'
  endif
  !
  if((it == 1).or.((it /= 0).and.(mod(it,isav) == 0))) then
    !
    ! trajectory of one particle.
    ! (first 100 slots: electrons on rank 1 -> xvelc, ions on rank 2 -> xvion;
    ! 501 blank-separated columns: t, then x, y, vx, vy, vz of each particle.
    ! The format was (f10.3,500d12.5), whose columns run together.)
    !
    if(MYRANK==mod(1,NPROCS)) write(70,'(f10.3,500e13.5)') time,(x(i,1)*dx,y(i,1)*dy,&
                           vx(i,1),vy(i,1),vz(i,1),i=MYRANK*maxnp+1,MYRANK*maxnp+100)
    if(MYRANK==mod(2,NPROCS)) write(71,'(f10.3,500e13.5)') time,(x(i,2)*dx,y(i,2)*dy,&
                           vx(i,2),vy(i,2),vz(i,2),i=MYRANK*maxnp+1,MYRANK*maxnp+100)
    !
    ! time history of field, kinetic and total energies.
    !
    if(MYRANK==0) then
      wem=exe+eye+eze+bxe+bye+bze
      total=wem
      do isp=1,nsp
        ake(isp)=ake(isp)/anme
        total=total+ake(isp)
      enddo
      write(10,'(11e13.5)') time,exe,eye,eze,bxe,bye,bze,(ake(isp),isp=1,nsp),total
    endif

    ! wave patterns of each field component.
    ! (Ex..Bz along y at x = nx/2; a new file every idata records)
    !
    if(it0 == 0) ij=it/isav
    if(it0 /= 0) ij=(it-it0)/isav-1
    if(ij >= 0) then
      irk=irk+1
      irk1=irk+1
      irk2=irk1+1
      irk3=irk2+1
      irk4=irk3+1
      irk5=irk4+1
      if(mod(ij+idata,idata) == 0) then
        if(it0 == 0) iii=(ij+idata)/idata
        if(it0 /= 0) iii=it/isav/idata+1
        !       if(MYRANK==mod(irk,NPROCS)) then
        write(fnex,'(a2,i3.3,a4)') 'ex',iii,'.hst'
        open(40,file=fnex,form='formatted',status='unknown')
        !       endif
        !       if(MYRANK==mod(irk1,NPROCS)) then
        write(fney,'(a2,i3.3,a4)') 'ey',iii,'.hst'
        open(41,file=fney,form='formatted',status='unknown')
        !       endif
        !       if(MYRANK==mod(irk2,NPROCS)) then
        write(fnez,'(a2,i3.3,a4)') 'ez',iii,'.hst'
        open(42,file=fnez,form='formatted',status='unknown')
        !       endif
        !       if(MYRANK==mod(irk3,NPROCS)) then
        write(fnbx,'(a2,i3.3,a4)') 'bx',iii,'.hst'
        open(43,file=fnbx,form='formatted',status='unknown')
        !       endif
        !       if(MYRANK==mod(irk4,NPROCS)) then
        write(fnby,'(a2,i3.3,a4)') 'by',iii,'.hst'
        open(44,file=fnby,form='formatted',status='unknown')
        !       endif
        !       if(MYRANK==mod(irk5,NPROCS)) then
        write(fnbz,'(a2,i3.3,a4)') 'bz',iii,'.hst'
        open(45,file=fnbz,form='formatted',status='unknown')
        !       endif
      endif

      !     if(MYRANK==mod(irk ,NPROCS)) write(40,'(1024e13.5)') (ex(ix,1),ix=1,nx)
      !     if(MYRANK==mod(irk1,NPROCS)) write(41,'(1024e13.5)') (ey(ix,1),ix=1,nx)
      !     if(MYRANK==mod(irk2,NPROCS)) write(42,'(1024e13.5)') (ez(ix,1),ix=1,nx)
      !     if(MYRANK==mod(irk3,NPROCS)) write(43,'(1024e13.5)') (bx(ix,1),ix=1,nx)
      !     if(MYRANK==mod(irk4,NPROCS)) write(44,'(1024e13.5)') (by(ix,1),ix=1,nx)
      !     if(MYRANK==mod(irk5,NPROCS)) write(45,'(1024e13.5)') (bz(ix,1),ix=1,nx)
      write(40,'(1024e13.5)') (ex(nx/2,iy),iy=1,ny)
      write(41,'(1024e13.5)') (ey(nx/2,iy),iy=1,ny)
      write(42,'(1024e13.5)') (ez(nx/2,iy),iy=1,ny)
      write(43,'(1024e13.5)') (bx(nx/2,iy),iy=1,ny)
      write(44,'(1024e13.5)') (by(nx/2,iy),iy=1,ny)
      write(45,'(1024e13.5)') (bz(nx/2,iy),iy=1,ny)
      !     do ix=1,nx
      !       write(40) ex(ix,1)
      !       write(41) ey(ix,1)
      !       write(42) ez(ix,1)
      !       write(43) bx(ix,1)
      !       write(44) by(ix,1)
      !       write(45) bz(ix,1)
      !     enddo
      if((mod(ij+1,idata) == 0).or.(it == nt)) then
        !       if(MYRANK==mod(irk ,NPROCS)) close(40)
        !       if(MYRANK==mod(irk1,NPROCS)) close(41)
        !       if(MYRANK==mod(irk2,NPROCS)) close(42)
        !       if(MYRANK==mod(irk3,NPROCS)) close(43)
        !       if(MYRANK==mod(irk4,NPROCS)) close(44)
        !       if(MYRANK==mod(irk5,NPROCS)) close(45)
        close(40)
        close(41)
        close(42)
        close(43)
        close(44)
        close(45)
      endif
    endif
  endif
  !
  if(it == nt) then
    !
    ! last step: close the history files and write the restart files
    ! (one record per particle; read back by cont)
    if(MYRANK==0) close(10)
    if(MYRANK==mod(1,NPROCS)) close(70)
    if(MYRANK==mod(2,NPROCS)) close(71)
    !
    irk=irk+1
    if(MYRANK==mod(irk,NPROCS)) then
      open(80,file='rhoj.fin',form='unformatted')
      write(80) rho
      write(80) jxm
      write(80) jxp
      write(80) jym
      write(80) jyp
      write(80) jzm
      write(80) jzp
      close(80)
    endif
    !
    irk=irk+1
    if(MYRANK==mod(irk,NPROCS)) then
      open(81,file='fld.fin',form='unformatted')
      write(81) bx
      write(81) by
      write(81) bz
      write(81) ex
      write(81) ey
      write(81) ez
      close(81)
    endif
    !
    write(fnxfin,'(a6,i3.3,a4)') 'x_rank',MYRANK,'.fin'
    open(82,file=fnxfin,form='unformatted',status='unknown')
    do isp=1,nsp
      do i=MYRANK*maxnp+1,(MYRANK+1)*maxnp
        write(82) x(i,isp)
      enddo
    enddo
    close(82)
    !
    write(fnyfin,'(a6,i3.3,a4)') 'y_rank',MYRANK,'.fin'
    open(83,file=fnyfin,form='unformatted',status='unknown')
    do isp=1,nsp
      do i=MYRANK*maxnp+1,(MYRANK+1)*maxnp
        write(83) y(i,isp)
      enddo
    enddo
    close(83)
    !
    write(fnvxfin,'(a7,i3.3,a4)') 'vx_rank',MYRANK,'.fin'
    open(84,file=fnvxfin,form='unformatted',status='unknown')
    do isp=1,nsp
      do i=MYRANK*maxnp+1,(MYRANK+1)*maxnp
        write(84) vx(i,isp)
      enddo
    enddo
    close(84)
    !
    write(fnvyfin,'(a7,i3.3,a4)') 'vy_rank',MYRANK,'.fin'
    open(85,file=fnvyfin,form='unformatted',status='unknown')
    do isp=1,nsp
      do i=MYRANK*maxnp+1,(MYRANK+1)*maxnp
        write(85) vy(i,isp)
      enddo
    enddo
    close(85)
    !
    write(fnvzfin,'(a7,i3.3,a4)') 'vz_rank',MYRANK,'.fin'
    open(86,file=fnvzfin,form='unformatted',status='unknown')
    do isp=1,nsp
      do i=MYRANK*maxnp+1,(MYRANK+1)*maxnp
        write(86) vz(i,isp)
      enddo
    enddo
    close(86)
    !
    write(fngfin,'(a6,i3.3,a4)') 'g_rank',MYRANK,'.fin'
    open(87,file=fngfin,form='unformatted',status='unknown')
    do isp=1,nsp
      do i=MYRANK*maxnp+1,(MYRANK+1)*maxnp
        write(87) gamma(i,isp)
      enddo
    enddo
    close(87)
    !
  endif
end subroutine output
