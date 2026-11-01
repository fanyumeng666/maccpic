! Over many seeds, reproduce exactly the loading of one rank of Case 1
! (32768 electrons then 32768 ions, 7 draws each: x, vx, vy, vz pairs) and
! collect the z-scores of six sample statistics per velocity component:
! mean, skewness, counts beyond +2, -2, +3, -3 sigma.  For independent
! N(0,1) samples every z-score is N(0,1): rms 1, 4.55 % beyond 2, 0.27 %
! beyond 3.  ran2 (module random) is compared with the intrinsic generator.
program test_seeds
  use random
  implicit none
  integer, parameter :: npart = 32768, nseed = 1000, nstat = 6
  real(8), parameter :: twopi = 6.283185307179586d0
  real(8) :: g, u1, u2, dummy, s1(3,2), s2(3,2), s3(3,2), z(nstat), e2, e3, sd2, sd3, se, m1, m2c, m3c
  integer :: cp2(3,2), cm2(3,2), cp3(3,2), cm3(3,2)
  integer :: ip, k, isp, seed, gen, ist, n, i
  real(8) :: zsum(nstat,2), zsq(nstat,2)
  integer :: n2(nstat,2), n3(nstat,2), cnt(2), plusheavier(2)
  integer, allocatable :: sdarr(:)
  character(len=10) :: sname(nstat) = (/ 'mean      ', 'skewness  ', 'N(>+2s)   ', 'N(<-2s)   ', 'N(>+3s)   ', 'N(<-3s)   ' /)
  character(len=10) :: gname(2) = (/ 'ran2      ', 'intrinsic ' /)
  e2 = npart*0.5d0*erfc(2d0/sqrt(2d0)); sd2 = sqrt(e2*(1d0 - e2/npart))
  e3 = npart*0.5d0*erfc(3d0/sqrt(2d0)); sd3 = sqrt(e3*(1d0 - e3/npart))
  se = sqrt(6d0/npart)
  call random_seed(size=n); allocate(sdarr(n))
  zsum = 0d0; zsq = 0d0; n2 = 0; n3 = 0; cnt = 0; plusheavier = 0
  do gen = 1, 2
    do seed = 1, nseed
      if (gen == 1) then
        call rng_init(seed)
      else
        do i = 1, n
          sdarr(i) = 1000003*seed + 7919*i
        end do
        call random_seed(put=sdarr)
      end if
      s1 = 0d0; s2 = 0d0; s3 = 0d0; cp2 = 0; cm2 = 0; cp3 = 0; cm3 = 0
      do isp = 1, 2
        do ip = 1, npart
          dummy = draw()
          do k = 1, 3
            u1 = draw(); u2 = draw()
            g = sqrt(-2d0*log(u1))*cos(twopi*u2)
            s1(k,isp) = s1(k,isp) + g; s2(k,isp) = s2(k,isp) + g*g; s3(k,isp) = s3(k,isp) + g**3
            if (g > 2d0) cp2(k,isp) = cp2(k,isp) + 1
            if (g < -2d0) cm2(k,isp) = cm2(k,isp) + 1
            if (g > 3d0) cp3(k,isp) = cp3(k,isp) + 1
            if (g < -3d0) cm3(k,isp) = cm3(k,isp) + 1
          end do
        end do
      end do
      do isp = 1, 2
        do k = 1, 3
          z(1) = s1(k,isp)/sqrt(dble(npart))
          ! central skewness: var = 6/N for a normal sample (the raw third moment has var 15/N)
          m1 = s1(k,isp)/npart
          m2c = s2(k,isp)/npart - m1*m1
          m3c = s3(k,isp)/npart - 3d0*m1*s2(k,isp)/npart + 2d0*m1**3
          z(2) = m3c/m2c**1.5d0/se
          z(3) = (cp2(k,isp) - e2)/sd2
          z(4) = (cm2(k,isp) - e2)/sd2
          z(5) = (cp3(k,isp) - e3)/sd3
          z(6) = (cm3(k,isp) - e3)/sd3
          cnt(gen) = cnt(gen) + 1
          if (cp3(k,isp) > cm3(k,isp)) plusheavier(gen) = plusheavier(gen) + 1
          do ist = 1, nstat
            zsum(ist,gen) = zsum(ist,gen) + z(ist)
            zsq(ist,gen) = zsq(ist,gen) + z(ist)**2
            if (abs(z(ist)) > 2d0) n2(ist,gen) = n2(ist,gen) + 1
            if (abs(z(ist)) > 3d0) n3(ist,gen) = n3(ist,gen) + 1
          end do
        end do
      end do
    end do
    write(*,'(a,a,i0,a,i0,a)') '==== ', trim(gname(gen)), nseed, ' seeds x 2 species x 3 components = ', cnt(gen), ' samples of 32768'
    write(*,'(a)') '  statistic    mean z    rms z   %|z|>2 (4.55)   %|z|>3 (0.27)'
    do ist = 1, nstat
      write(*,'(a,a,f8.3,f9.3,f12.2,f14.2)') '  ', sname(ist), zsum(ist,gen)/cnt(gen), sqrt(zsq(ist,gen)/cnt(gen)), &
           100d0*n2(ist,gen)/cnt(gen), 100d0*n3(ist,gen)/cnt(gen)
    end do
    write(*,'(a,f6.1,a)') '  +3 sigma tail heavier than -3 sigma tail in ', 100d0*plusheavier(gen)/cnt(gen), ' % of the samples (~47 % with ties)'
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
end program test_seeds
