!=====================================================================
!  input: run parameters read from a namelist file (default ccp1d.in),
!  parameter checks and the derived quantities that do not need the
!  cross sections (dx, dt, weight, maxnp).  All quantities in SI.
!  Tests set the namelist variables directly and call derive_parameters.
!=====================================================================
module input
  implicit none
  !
  ! ---- namelist variables with Case 1 (Turner et al. 2013) defaults ----
  real(8) :: L = 0.067d0                 ! electrode gap [m]
  integer :: nx = 128                    ! number of cells
  real(8) :: freq = 13.56d6              ! drive frequency [Hz]
  real(8) :: V0 = 450d0                  ! drive amplitude [V], phi(0,t) = V0 sin(2 pi f t)
  real(8) :: freq_lf = 0d0, V0_lf = 0d0  ! second (low) frequency, added to phi(0,t)
  real(8) :: n_gas = 9.64d20             ! neutral density [m^-3] (0: collisionless)
  real(8) :: T_gas = 300d0               ! gas temperature [K]
  real(8) :: n0 = 2.56d14                ! initial plasma density [m^-3]
  real(8) :: Te0 = 30000d0, Ti0 = 300d0  ! initial temperatures [K]
  integer :: nppc = 512                  ! initial macro-particles per cell per species
  integer :: steps_per_period = 400
  integer :: n_periods = 1280            ! periods to reach (cumulative on restart)
  integer :: periods_this_run = 0        ! > 0: stop after this many periods in this invocation (restart files written)
  integer :: n_avg = 32                  ! measurement window: last n_avg periods
  integer :: n_phase = 100               ! phase bins per period for the x-t maps
  integer :: seed = 12345
  character(len=256) :: xsec_dir = '../../data/he'
  logical :: restart = .false.
  integer :: restart_every = 160         ! periods between restart files
  real(8) :: iedf_emax = 600d0, iedf_de = 1d0     ! wall energy histograms [eV]
  real(8) :: eepf_emax = 100d0, eepf_de = 0.5d0   ! bulk electron energy histogram [eV]
  real(8) :: numax_emax_e = 1000d0, numax_emax_i = 1000d0  ! energy cap of the nu_max search [eV]
  integer :: maxnp_factor = 4            ! maxnp = maxnp_factor * initial particles per rank
  real(8) :: max_wall_minutes = 0d0      ! > 0: stop cleanly (restart file) after this wall time
  logical :: dump_initial = .false.      ! write the t = 0 particles and grid (init_*.dat)
  !
  namelist /ccp1d/ L, nx, freq, V0, freq_lf, V0_lf, n_gas, T_gas, n0, Te0, Ti0, &
       nppc, steps_per_period, n_periods, periods_this_run, n_avg, n_phase, seed, xsec_dir, &
       restart, restart_every, iedf_emax, iedf_de, eepf_emax, eepf_de, &
       numax_emax_e, numax_emax_i, maxnp_factor, max_wall_minutes, dump_initial
  !
  ! ---- derived ----
  real(8)    :: dx, dt, period, weight     ! weight: real particles per m^2 per macro-particle
  integer(8) :: np_init_total              ! initial macro-particles per species (all ranks)
  integer(8) :: maxnp                      ! particle slots per species per rank
  integer    :: nbin_iedf, nbin_eepf
  character(len=256) :: input_file = 'ccp1d.in'
  !
contains
  !
  subroutine read_input(nranks)
    integer, intent(in) :: nranks
    integer :: u, ios
    if (command_argument_count() >= 1) call get_command_argument(1, input_file)
    open(newunit=u, file=trim(input_file), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ccp1d: cannot open input file ', trim(input_file)
      stop 1
    end if
    read(u, nml=ccp1d, iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ccp1d: error reading namelist &ccp1d from ', trim(input_file), ' iostat =', ios
      stop 1
    end if
    close(u)
    call derive_parameters(nranks)
  end subroutine read_input
  !
  subroutine derive_parameters(nranks)
    integer, intent(in) :: nranks
    character(len=200) :: msg
    msg = ''
    if (nx < 4)                 msg = 'nx must be >= 4'
    if (L <= 0d0)               msg = 'L must be > 0'
    if (freq <= 0d0)            msg = 'freq must be > 0'
    if (steps_per_period < 1)   msg = 'steps_per_period must be >= 1'
    if (nppc < 1)               msg = 'nppc must be >= 1'
    if (n_avg < 1 .or. n_avg > n_periods) msg = 'need 1 <= n_avg <= n_periods'
    if (n_phase < 1)            msg = 'n_phase must be >= 1'
    if (mod(steps_per_period, n_phase) /= 0) msg = 'steps_per_period must be a multiple of n_phase'
    if (n_gas < 0d0)            msg = 'n_gas must be >= 0'
    if (maxnp_factor < 1)       msg = 'maxnp_factor must be >= 1'
    if (iedf_de <= 0d0 .or. eepf_de <= 0d0) msg = 'histogram bin widths must be > 0'
    if (len_trim(msg) > 0) then
      write(*,*) 'ccp1d input error: ', trim(msg)
      stop 1
    end if
    dx = L/nx
    period = 1d0/freq
    dt = period/steps_per_period
    np_init_total = int(nx,8)*int(nppc,8)
    weight = n0*L/dble(np_init_total)
    maxnp = int(maxnp_factor,8)*((np_init_total + nranks - 1)/nranks)
    nbin_iedf = int(iedf_emax/iedf_de + 0.5d0)
    nbin_eepf = int(eepf_emax/eepf_de + 0.5d0)
  end subroutine derive_parameters
  !
  ! applied voltage at the driven electrode x = 0
  pure function v_drive(t) result(v)
    real(8), intent(in) :: t
    real(8) :: v
    real(8), parameter :: twopi = 6.283185307179586d0
    v = V0*sin(twopi*freq*t)
    if (V0_lf /= 0d0) v = v + V0_lf*sin(twopi*freq_lf*t)
  end function v_drive
  !
  subroutine print_input(u)
    integer, intent(in) :: u
    write(u,'(a)') '# ---- input parameters (namelist &ccp1d) ----'
    write(u,'(a,es14.6,a)') 'L = ', L, '  # m'
    write(u,'(a,i0)')       'nx = ', nx
    write(u,'(a,es14.6,a)') 'freq = ', freq, '  # Hz'
    write(u,'(a,es14.6,a)') 'V0 = ', V0, '  # V'
    write(u,'(a,es14.6,a)') 'freq_lf = ', freq_lf, '  # Hz'
    write(u,'(a,es14.6,a)') 'V0_lf = ', V0_lf, '  # V'
    write(u,'(a,es14.6,a)') 'n_gas = ', n_gas, '  # m^-3'
    write(u,'(a,es14.6,a)') 'T_gas = ', T_gas, '  # K'
    write(u,'(a,es14.6,a)') 'n0 = ', n0, '  # m^-3'
    write(u,'(a,es14.6,a)') 'Te0 = ', Te0, '  # K'
    write(u,'(a,es14.6,a)') 'Ti0 = ', Ti0, '  # K'
    write(u,'(a,i0)')       'nppc = ', nppc
    write(u,'(a,i0)')       'steps_per_period = ', steps_per_period
    write(u,'(a,i0)')       'n_periods = ', n_periods
    write(u,'(a,i0)')       'periods_this_run = ', periods_this_run
    write(u,'(a,i0)')       'n_avg = ', n_avg
    write(u,'(a,i0)')       'n_phase = ', n_phase
    write(u,'(a,i0)')       'seed = ', seed
    write(u,'(a,a)')        'xsec_dir = ', trim(xsec_dir)
    write(u,'(a,l1)')       'restart = ', restart
    write(u,'(a,i0)')       'restart_every = ', restart_every
    write(u,'(a,es14.6,a)') 'iedf_emax = ', iedf_emax, '  # eV'
    write(u,'(a,es14.6,a)') 'iedf_de = ', iedf_de, '  # eV'
    write(u,'(a,es14.6,a)') 'eepf_emax = ', eepf_emax, '  # eV'
    write(u,'(a,es14.6,a)') 'eepf_de = ', eepf_de, '  # eV'
    write(u,'(a,es14.6,a)') 'numax_emax_e = ', numax_emax_e, '  # eV'
    write(u,'(a,es14.6,a)') 'numax_emax_i = ', numax_emax_i, '  # eV'
    write(u,'(a,i0)')       'maxnp_factor = ', maxnp_factor
    write(u,'(a,es14.6,a)') 'max_wall_minutes = ', max_wall_minutes, '  # min'
    write(u,'(a,l1)')       'dump_initial = ', dump_initial
    write(u,'(a)') '# ---- derived ----'
    write(u,'(a,es14.6,a)') 'dx = ', dx, '  # m'
    write(u,'(a,es14.6,a)') 'dt = ', dt, '  # s'
    write(u,'(a,es14.6,a)') 'period = ', period, '  # s'
    write(u,'(a,es14.6,a)') 'weight = ', weight, '  # m^-2 (real particles per macro-particle per unit area)'
    write(u,'(a,i0)')       'np_init_total = ', np_init_total
    write(u,'(a,i0)')       'maxnp = ', maxnp
    write(u,'(a,i0)')       'nsteps_total = ', int(n_periods,8)*int(steps_per_period,8)
  end subroutine print_input
  !
end module input
