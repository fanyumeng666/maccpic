!
!*********************************************************************
!
!  em2 -- a two-dimensional electromagnetic plasma simulation code
!         two spatial dimension and three velocity components with
!         relativistic electrons and ions
!
!  written by S. Matsukiyo, 19 Apr. 2012.
!
!*********************************************************************
!
!  Program flow
!    fresh run (irun = 1): init -> creator -> setv
!    restart   (irun > 1): init -> cont (reads the *.fin files of the
!                          previous run)
!    every step          : move (push momenta) -> curnt (push positions,
!                          deposit rho and J) -> output -> fields (advance
!                          E and B, Poisson correction of E)
!  MPI: the particles are divided among the ranks; every rank keeps the
!       whole field grid, and rho and J are summed with MPI_ALLREDUCE.
!
program main
  use global, only : nproc,irun,nt0,it0,isav,iphs,dt,parentwave,nt,it,time,&
                    IERR,NPROCS,MYRANK
  implicit none
  ! ----- MPI -----
  include 'mpif.h'
  ! ---------------
  !
  ! ===== parameters to be set in 'input.f90' =====
  !
  ! nx,ny      : number of grid cells
  ! npcel      : number of particles for each specie in a cell
  ! nsp        : number of species
  !
  ! irun       : index of the run
  ! nt0        : number of time iteration
  ! isav       : sampling interval of field time history data
  ! iphs       : sampling interval of particle snapshot data
  ! idata      : number of time sampling in a field file
  !
  ! alx,aly    : system length (in unit of c/wpe)
  ! wpe        : electron plasma freq. (fixed to 1.)
  ! wce0       : electron gyro freq. (|wce0| = B0)
  ! theta      : angle between B0 and the x-axis (B0 is in the x-z plane)
  ! rmass      : ion-to-electron mass ratio
  ! dx,dy      : spatial grid cell spacing (dx=alx/nx = dy=aly/ny)
  ! dt         : time step (=dx*0.7=dy*0.7)
  ! c          : speed of light (= 1.)
  ! qme        : charge-to-mass ratio for electrons (fixed to -1.)
  ! qmi        : charge-to-mass ratio for ions
  !
  ! parentwave : introduction of a monochromatic parent Alfven wave if .true.
  ! mno1       : mode number of the parent wave
  ! b01        : amplitude of the parent wave normalized to B0
  ! v01        : electron velocity determined by Walen relation
  ! v02        : ion velocity determined by Walen relation
  ! freq1      : frequency of the parent wave
  !
  ! va(:)      : drift velocity parallel to B0
  ! ve(:)      : ring velocity perp. to B0
  ! temp(:)    : temperature of each specie (T_j/m_j c^2)
  ! dn(:)      : relative density of each specie
  !
  ! ===== other parameters used in the code =====
  !
  !            --- field variables ---
  !
  ! jxm,jxp,jym,jyp,jzm,jzp : jx-(x), jx+(x), jy-(x), jy+(x), jz-(x), jz+(x) field values
  ! bx, by, bz : bx(x), by(x), bz(x) field values
  ! ex, ey, ez : ex(x), ey(x), ez(x) field values
  ! phi        : phix(x) potential values
  ! rho        : rhox(x) charge density values
  ! exe, eye, eze, bxe, bye, bze : energy density of each field component
  !
  !            --- particle variables ---
  !
  ! x, y       : particle axial position
  ! vx, vy, vz : particle velocity components
  ! gamma      : particle Lorentz factor
  ! np(:)      : number of superparticles for each specie (=ng*npcel)
  ! ake(:)     : kinetic energy for each specie
  !
  ! control (variables which control the flow of the program)
  !
  integer(4) isp
  !
  integer(4) t1,t2,tr,tmax,diff
  real(8) cput1,cput2
  call system_clock(t1)
  call cpu_time(cput1)
  !
  ! ----- MPI -----
  !
  call MPI_INIT(IERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERR)
  if(NPROCS /= nproc) then
    write(*,*) 'NPROCS /= nproc'
    stop
  endif
  !
  ! it0: step at which this run starts (0 for a fresh run); nt: last step
  it0    = nt0*(irun-1)-isav
  if(it0 <= 0) it0=0
  if(it0 == 0) nt=nt0-isav
  !
  ! initializing variables
  !
  !  =========
  call init
  !  =========
  !
  ! initial values of parameters for second run.
  !
  if(it0 /= 0) then
    !  ===================
    call cont
    !  ===================
    nt=nt0+it0
  else
    !
    ! initial particle distribution
    !
    !  ============
    call creator
    !   call setrho
    !  ============
    !
    ! parent wave input.
    !
    !                 ============
    !   if(parentwave) call waveint
    !                 ============
    !
    ! shift the velocities back to t = -dt/2 (leapfrog start)
    !  ================================
    !   call fields
    call setv
    !  ================================
    !
  endif
  !
  ! ----------------------------------
  ! ----- begin time step loop. ------
  ! ----------------------------------
  !
  do
    !
    if(mod(it,10) == 0) write(*,653) 'it =',it
    653 format(a8,i8)
    !
    !  ==========
    call move
    call curnt
    !  ==========
    !
    if((mod(it,isav) == 0).or.(it == 1).or.&
       (mod(it,iphs) == 0).or.(it == nt)) then
      !         ======
      call output
      !         ======
    endif
    if(it == nt) exit
    it=it+1
    time=it*dt
    !
    !  ===========
    call fields
    !  ===========
    !
  enddo
  !
  ! ----- MPI -----
  call MPI_FINALIZE(IERR)
  ! ---------------
  !
  call cpu_time(cput2)
  call system_clock(t2,tr,tmax)
  if(t2 < t1) then
    diff=tmax-t1+t2
  else
    diff=t2-t1
  endif
  write(*,*) "cpu time =",cput2-cput1,"sec"
  write(*,*) "time it took was:",diff/dble(tr)
  !
  ! end of run.
  !
  stop
end program main
