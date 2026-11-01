!*********************************************************************
!
!  ccp1d -- one-dimensional bounded electrostatic PIC/MCC code for a
!           radio-frequency capacitive discharge (1d3v, electrons and
!           He+, Monte Carlo collisions with a fixed helium background)
!
!  Program flow
!    read_input -> allocate -> poisson_init -> (xsec_load) ->
!    fresh run: init_fresh;  restart: restart_read
!    every step: push (both species) -> mcc -> deposit_rho ->
!                poisson_solve with V(t_{n+1}) -> diag_step
!    every period: diag_period (history line), restart files
!    end: diag_final (profiles, distributions, x-t maps, info.txt)
!  MPI: particles are divided among the ranks, every rank holds the
!       whole grid; rho is summed with MPI_ALLREDUCE every step.
!  Units: SI throughout; energies in eV only in tables and output.
!
!*********************************************************************
program main
  use global
  use random
  use poisson
  use deposit
  use push
  use diag
  use init
  implicit none
  integer :: k, stop_flag, nper_stop
  real(8) :: t0
  !
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  !
  call read_input(nprocs)
  call rng_init(seed + myrank)
  call allocate_arrays()
  call poisson_init(nx, dx)
  call push_init()
  call diag_init()
  mcc_on = .false.
  if (n_gas > 0d0 .and. myrank == 0) then
    write(*,'(a)') 'ccp1d WARNING: n_gas > 0 but collisions are not implemented yet (M2); running collisionless'
  end if
  !
  if (restart) then
    call restart_read()
    call history_open(.true.)
  else
    call init_fresh()
    call write_info_header()
    call history_open(.false.)
    call diag_period()            ! period 0: initial state
  end if
  if (myrank == 0) then
    write(*,'(a,i0,a,i0,a,i0,a)') 'ccp1d: nx = ', nx, ', steps/period = ', steps_per_period, &
         ', periods ', iperiod, ' -> ', n_periods
    write(*,'(a,i0,a,i0)') 'ccp1d: particles per species per rank = ', np(1), ', maxnp = ', maxnp
    flush(6)
  end if
  t0 = MPI_WTIME()
  nper_stop = n_periods
  if (periods_this_run > 0) nper_stop = min(n_periods, iperiod + periods_this_run)
  !
  ! ---------------- time loop, one period per outer iteration ----------------
  do while (iperiod < nper_stop)
    call diag_window_check()
    do k = 1, steps_per_period
      call push_species(1)
      call push_species(2)
      ! if (mcc_on) call mcc_step()           ! M2
      call deposit_rho()
      ex_old = ex
      it = it + 1
      time = it*dt
      call poisson_solve(rho, v_drive(time), 0d0, phi, ex)
      call diag_step()
    end do
    iperiod = iperiod + 1
    call diag_period()
    ! stop cleanly when the wall-clock budget is used up (decided by rank 0)
    stop_flag = 0
    if (myrank == 0 .and. max_wall_minutes > 0d0) then
      if (MPI_WTIME() - t0 > 60d0*max_wall_minutes) stop_flag = 1
    end if
    call MPI_BCAST(stop_flag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (mod(iperiod, restart_every) == 0 .or. iperiod == nper_stop .or. stop_flag == 1) call restart_write()
    if (stop_flag == 1) then
      if (myrank == 0) write(*,'(a,i0,a)') 'ccp1d: wall-clock limit reached after period ', iperiod, &
           '; restart files written, rerun with restart = .true.'
      exit
    end if
  end do
  !
  if (iperiod == n_periods) then
    call diag_final()
  else if (myrank == 0) then
    write(*,'(a,i0,a,i0,a)') 'ccp1d: stopped after period ', iperiod, ' of ', n_periods, &
         '; rerun with restart = .true. to continue'
  end if
  if (myrank == 0) then
    write(*,'(a,f10.1,a)') 'ccp1d: done, wall time ', MPI_WTIME() - t0, ' s'
    flush(6)
  end if
  call MPI_FINALIZE(ierr)
end program main
