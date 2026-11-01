! Tail symmetry of the Box-Muller normals in the loading pattern of init_fresh
! (7 draws per particle: x, then 3 pairs for vx, vy, vz): counts beyond
! +/-2, 3, 3.5, 4, 4.5 sigma per velocity component, mean and skewness, for the
! module generator (ran2) and, as a reference, gfortran's intrinsic random_number.
program test_tails
  use random
  implicit none
  integer, parameter :: npart = 20000000, nth = 5
  real(8), parameter :: thr(nth) = (/ 2d0, 3d0, 3.5d0, 4d0, 4.5d0 /), twopi = 6.283185307179586d0
  integer(8) :: cpos(nth,3), cneg(nth,3)
  real(8) :: s1(3), s2(3), s3(3), g, u1, u2, dummy, expc, sd
  integer :: ip, k, it, gen, seed
  character(len=12) :: gname(2) = (/ 'ran2        ', 'intrinsic   ' /)
  do gen = 1, 2
    do seed = 12345, 12346
      call rng_init(seed)
      call random_seed()
      cpos = 0; cneg = 0; s1 = 0d0; s2 = 0d0; s3 = 0d0
      do ip = 1, npart
        dummy = draw()                          ! the position draw
        do k = 1, 3
          u1 = draw(); u2 = draw()
          g = sqrt(-2d0*log(u1))*cos(twopi*u2)
          s1(k) = s1(k) + g; s2(k) = s2(k) + g*g; s3(k) = s3(k) + g**3
          do it = 1, nth
            if (g > thr(it)) cpos(it,k) = cpos(it,k) + 1
            if (g < -thr(it)) cneg(it,k) = cneg(it,k) + 1
          end do
        end do
      end do
      write(*,'(a,a,i0,a,i0,a)') '==== ', trim(gname(gen)), seed, ', ', npart, ' particles'
      do k = 1, 3
        write(*,'(a,i1,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a)') '  component ', k, ': mean ', s1(k)/npart, ' (se ', 1d0/sqrt(dble(npart)), &
             ')  skew ', s3(k)/npart/(s2(k)/npart)**1.5d0, ' (se ', sqrt(6d0/npart), ')'
        do it = 1, nth
          expc = npart*0.5d0*erfc(thr(it)/sqrt(2d0))
          sd = sqrt(expc)
          write(*,'(a,f4.1,a,i9,a,f6.1,a,i9,a,f6.1,a,f12.1)') '     |g| > ', thr(it), ':  +side ', cpos(it,k), ' (', (cpos(it,k) - expc)/sd, &
               ' sd)   -side ', cneg(it,k), ' (', (cneg(it,k) - expc)/sd, ' sd)   expected', expc
        end do
      end do
    end do
  end do
contains
  function draw() result(u)
    real(8) :: u
    if (gen == 1) then
      u = rand_u()
    else
      call random_number(u)
      if (u <= 0d0) u = 1d-300
    end if
  end function draw
end program test_tails
