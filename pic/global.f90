!=====================================================================
!  global: variables shared by all routines (derived constants, fields,
!  particles). Grid arrays run from 0 to nx+2 / ny+2: interior points
!  are 1..nx, 1..ny, the rest are periodic guard cells.
!=====================================================================
module global
  use input
  implicit none
  !
  integer(4)  IERR,MYRANK,NPROCS    !(NPROCS = nproc)
  integer(4), parameter :: nx1 = nx+1, ny1 = ny+1, nx2 = nx+2, ny2 = ny+2
  ! particle slots per rank (maxnp) and in total (maxnpp): background
  ! plasma plus the two Harris current sheets
  ! integer(8), parameter :: maxnp = nx*ny/nproc*npcel
  integer(8), parameter :: maxnp = 2*nx*npcel*int(ny/2+2*csw/dx*csden*tanh(ny*dy/2/csw))/nproc
  integer(8), parameter :: maxnpp = 2*nx*npcel*int(ny/2+2*csw/dx*csden*tanh(ny*dy/2/csw))

  ! system size (c/wpe) and in grid units; ctx, cty = c*dt/dx, c*dt/dy;
  ! rnmec2i, anorml1 normalize the field energies; gcsw = csw in cells
  real(8),    parameter :: alx=dx*nx, aly=dy*ny
  real(8),    parameter :: xn=nx,     yn=ny
  real(8),    parameter :: pi=3.141592653589793d0, twopi=2*pi
  real(8),    parameter :: rd=theta*pi/180.d0, ak0=twopi/alx
  real(8),    parameter :: tem=0.5d0*dt, ctx=dt*c/dx, cty=dt*c/dy
  real(8),    parameter :: ctxh=ctx*0.5d0, ctyh=cty*0.5d0
  real(8),    parameter :: rnmec2i=(qme/wpe/c)**2/(alx*aly)
  real(8),    parameter :: anorml1=0.5d0*dx*dy*rnmec2i
  real(8),    parameter :: gcsw=csw/dx
  !
  ! step counters, time and energies (exe...bze: field, anme: normalization)
  integer(4), save :: nt,it,it0
  real(8),    save :: time,exe,eye,eze,bxe,bye,bze,anme
  !
  ! np: particles per rank; am: macro-particle mass; ake: kinetic energy
  integer(8), dimension(nsp),              save :: np
  real(8),    dimension(nsp),              save :: am,tx,ty,ake
  ! fields; bxm,bym,bzm keep B of the previous step (set in fields, no
  ! longer used by the particle push). Staggering
  ! (index i,j -> position in cells): rho (i-1, j-1); Ex,Bx,Jx (i-1/2, j-1);
  ! Ey,By,Jy (i-1, j-1/2); Ez,Bz,Jz (i-1/2, j-1/2)
  real(8),    dimension(0:nx2,0:ny2),      save :: bx,by,bz,ex,ey,ez,bxm,bym,bzm
  ! currents at the old (jxm...) and new (jxp...) particle positions
  real(8),    dimension(0:nx2,0:ny2),      save :: jxm,jxp,jym,jyp,jzm,jzp,rho
  ! particles (slot, species): x, y in cells, velocity v/c, Lorentz factor
  real(8),    dimension(:,:), allocatable, save :: x,y,vx,vy,vz,gamma
  !
end module global
