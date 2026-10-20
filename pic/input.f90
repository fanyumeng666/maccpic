!=====================================================================
!  input: run parameters (compile-time constants).
!  Every other file uses these values, so rebuild everything after
!  editing this file (make clean; make).
!  Units: length c/wpe, time 1/wpe, velocity c (wpe = c = 1).
!=====================================================================
module input
  ! grid cells, particles per cell for each species, number of species
  integer(4), parameter :: nx = 1024, ny = 512, npcel = 50, nsp = 2
  ! MPI processes (must equal the number given to mpirun), OpenMP threads
  integer(4), parameter :: nproc = 4, nthreads = 4
  ! isav=20
  ! run control (see the list at the top of main.f90)
  integer(4), parameter :: irun  =      1
  integer(4), parameter :: nt0   =   20480
  integer(4), parameter :: isav  =     20
  integer(4), parameter :: iphs  =   1000
  integer(4), parameter :: idata =   1024
  !
  ! grid spacing (dx = dy is assumed by the Poisson solver), time step,
  ! electron gyrofrequency of the outer field, B0 angle, mass ratio
  real(8),    parameter :: dx    = 0.1118D+00, dy=dx
  real(8),    parameter :: dt    = dx*0.7d0
  real(8),    parameter :: wce0  =-1.498147003616D-00
  real(8),    parameter :: theta = 0.D+00
  real(8),    parameter :: rmass = 100.D+00
  logical,    parameter :: trace = .false.
  !
  logical,    parameter :: parentwave = .false.
  integer(4), parameter :: mno1  = 512
  real(8),    parameter :: b01   = 1.D-00
  real(8),    parameter :: v01   = 0.676581209040969D-00
  real(8),    parameter :: v02   = 0.619917185536019D-01
  real(8),    parameter :: freq1 = 0.198130655636275D-00
  ! rtemp = Ti/Te; creator uses it to share the current-sheet drift
  ! between electrons and ions (keep it equal to temp(2)*rmass/temp(1))
  real(8),    parameter :: rtemp = 4.D-00

  ! Harris current sheets: csw = half-thickness (c/wpe),
  ! csden = peak sheet density / background density
  real(8),    parameter :: csw = 2.23606797749979d0, csden=5.D+00
  !.
  real(8),    parameter :: wpe=1.d0, c=1.d0, qme=-1.d0, qmi=-qme/rmass
  !
contains
  ! drift velocity va, ring velocity ve, temperature T/(m c^2) (Ti/Te = 4
  ! here) and temperature anisotropy aist of each species:
  ! T_x = temp*(1+aist), T_y = T_z = temp (aist = 0: isotropic)
  subroutine com_vel(va,ve,temp,aist)
    implicit none
    real(8), dimension(nsp) :: va,ve,temp,aist
    !
    va(1) = 0.D-00;  ve(1) = 0.00D-00;  temp(1) = 4.4888888888D-02;        aist(1) = 0.D0
    va(2) = 0.D-00;  ve(2) = 0.00D-00;  temp(2) = 1.7955555555D-01/rmass;  aist(2) = 0.D0
    !
  end subroutine com_vel
  !=
  ! mass ratio rm, relative density dn and charge-to-mass ratio qm of each
  ! species; q is the charge of one macro-particle, chosen so that npcel
  ! particles per cell give wpe = 1
  subroutine com_den(dn,q,qm)
    implicit none
    real(8), dimension(nsp) :: dn,q,qm,rm
    !
    rm(1)= 1.0D-00; rm(2)=rmass!;         rm(3)= 1.0D+00!;       rm(4)=rmass
    dn(1)= 1.0D-00; dn(2)=1.D-00!;        dn(3)= 0.5D-00!;       dn(4)=0.5D-00
    qm(1)=-1.0D-00; qm(2)=1.0D-00/rm(2)!; qm(3)=-1.0D-00/rm(3)!; qm(4)=1.0D-00/rm(4)
    !
    q(1:nsp)=dx*dy*wpe*wpe/(npcel*qm(1:nsp))/rm(1:nsp)
    ! if(nsp > 3) then
    !   write(*,*) 'nsp > 3'
    !   stop
    ! else
    !   dens=0.d0
    !   do isp=2,nsp
    !     dens=dens+dn(isp)
    !   enddo
    !   if(dens /= 1.d0) then
    !     write(*,*) 'total ion density /= 1.0'
    !     stop
    !   endif
    ! endif
    !
  end subroutine com_den
  !
end module input
