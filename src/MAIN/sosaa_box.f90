PROGRAM sosaa_box

USE second_Parameters  ! NSPEC, ind_XXX, dp, NFIX
USE second_Monitor, ONLY: SPC_NAMES, EQN_NAMES
USE second_Global, ONLY: NPHOT
USE Chemistry_Mod  ! chemistry, kpp_setup, NPHOT

IMPLICIT NONE

!------------------!
! Constants
!------------------!

! INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,307)

! General constants
REAL(dp), PARAMETER :: pi = 2*ASIN(1.0d0)
REAL(dp), PARAMETER :: ppm = 1.0d-6               ! [mol mol-1]
REAL(dp), PARAMETER :: ppb = 1.0d-9               ! [mol mol-1]
REAL(dp), PARAMETER :: ppt = 1.0d-12              ! [mol mol-1]
REAL(dp), PARAMETER :: T00 = 273.15d0             ! [K]
REAL(dp), PARAMETER :: P00 = 1.01325d5            ! [Pa]
REAL(dp), PARAMETER :: Rgas = 8.31446261815324d0  ! [J mol-1 K-1], gas constant
REAL(dp), PARAMETER :: NA = 6.02214076d23         ! [molec mol-1], Avogadro constant

! Time constants
REAL(dp), PARAMETER :: DAY_SEC = 86400.0d0      ! one day in seconds
REAL(dp), PARAMETER :: HALFDAY_SEC = 43200.0d0  ! half day in seconds
REAL(dp), PARAMETER :: HOUR_SEC = 3600.0d0      ! one hour in seconds
REAL(dp), PARAMETER :: MINUTE_SEC = 60.0d0      ! one minute in seconds

! Constants for input data
integer, parameter :: MAXLEN_NAME = 50
integer, parameter :: MAXLEN_FILE_NAME = 100
integer, parameter :: MAXNUM_INPUT_DATA = 500
integer, parameter :: MAXLEN_MESSAGE= 500

! Other parameters
integer, parameter :: LINFO = 0
integer, parameter :: LWARNING = 1
integer, parameter :: LDEBUG = 2
integer, parameter :: LERROR = 3

! Other variables
REAL(dp) :: Nair = 2.25d19  ! [molec cm-3]
character(MAXLEN_MESSAGE) :: message

!------------------!
! Station or chamber info
!------------------!

! Station, chamber or campaign name, e.g., 'hyytiala', used to name the station folder
CHARACTER(MAXLEN_NAME) :: station

! location
REAL(dp) :: latitude_deg, latitude_rad
REAL(dp) :: longitude_deg, longitude_rad

! For chamber
real(dp) :: chamber_volume  ! [m3], chamber volume  
real(dp) :: flow_rate_in, flow_rate_out  ! [m3 s-1], inflow and outflow rates
real(dp) :: flow_rate_in_Lm, flow_rate_out_Lm  ! [L min-1], inflow and outflow rates in the unit of L min-1

!------------------!
! Environment
!------------------!

REAL(dp) :: temp, RH, pres, ws
REAL(dp) :: zenith_deg  ! [deg], solar zenith angle in degree
REAL(dp) :: zenith_rad  ! [rad], solar zenith angle in radian
REAL(dp) :: cos_zenith  ! [deg], cos(zenith_deg)
REAL(dp) :: sun_par  ! [-], sunlit fraction, for box model 
real(dp) :: albedo_ground  ! [-], ground albedo used for calculating total actinic flux
REAL(dp) :: cs(NSPEC)  ! [s-1], condensation sinks
REAL(dp) :: cs_H2SO4, cs_HNO3  ! [s-1], condensation sinks for H2SO4 and HNO3
REAL(dp) :: rhov  ! [kg m-3], absolute humidity, water vapor density
REAL(dp) :: J_values(NPHOT)  ! [s-1], input photolysis reaction rates from measurements
REAL(dp) :: K_values(NKVALUES)  ! reaction rates
REAL(dp) :: Rglob  ! [W m-2], income global solar radiation
REAL(dp), DIMENSION(84) :: tmp_wl  ! [nm], temporary wavelength array
REAL(dp), DIMENSION(NWL, 2) :: lamp_spectral  ! [-], spectral intensity distribution of lamps
REAL(dp), DIMENSION(84) :: swr_distribution  ! [-], shortwave spectral distribution
REAL(dp), DIMENSION(84) :: Rglob_spec  ! [W m-2], shortwave spectral
CHARACTER(MAXLEN_FILE_NAME) :: input_dir  ! input dir, usually the path of 'sosa_in'
CHARACTER(MAXLEN_FILE_NAME) :: output_dir  ! output dir, usually './output'
REAL(dp) :: actinic_flux(NWL)  ! [photon cm-2 s-1 nm-1], calculated actinic flux

!------------------!
! Chemistry
!------------------!

real(dp) :: conc(NSPEC)
real(dp) :: conc_ovocs

real(dp) :: emis(NSPEC)  ! [molec cm-3 s-1], emission rate

real(dp) :: inflow(NSPEC)  ! [molec cm-3 s-1], inflow source rate
real(dp) :: outflow(NSPEC)  ! [molec cm-3 s-1], outflow source rate
real(dp) :: conc_inflow(NSPEC)  ! [molec cm-3], concentration in the inflow gas
real(dp) :: conc_inflow_air  ! [molec cm-3], air concentration in the inflow gas

real(dp) :: k
real(dp) :: O2, N2, MO2N2, H2O, H2
real(dp) :: RO2
real(dp) :: I_ACF_in(NWL), I_ACF_out(NWL)

real(dp) :: vd(NSPEC)
real(dp) :: zref

real(dp) :: wall_loss(NSPEC)

logical :: l_has_fix  ! if chemistry scheme has fixed concentration species

! Aqueous

real(dp) :: H_ion_aq_wall, OH_ion_aq_wall
real(dp) :: cs_aq_wall(50), cw_aq_wall

!------------------!
! Time
!------------------!

real(dp) :: time, time_start, time_end
real(dp) :: year, month, day, hour, minute, second
real(dp) :: year_start, month_start, day_start
real(dp) :: hour_start, minute_start, second_start
real(dp) :: time_of_day
integer  :: julian_day
real(dp) :: time1, time2
real(dp) :: year1, year2 
real(dp) :: month1, month2
real(dp) :: day1, day2 
real(dp) :: hour1, hour2 
real(dp) :: minute1, minute2 
real(dp) :: second1, second2 
real(dp) :: dt                          ! time step, in seconds
real(dp) :: dt_display, dt_output
integer  :: time_in_ms                  ! time in milliseconds

!------------------!
! Input data
!------------------!

character(MAXLEN_FILE_NAME ) :: input_fname
integer :: input_fid
integer :: input_ncol  ! number of data columns
integer :: input_nrow  ! number of data rows
integer :: input_nline  ! number of file rows

character(MAXLEN_NAME) :: input_vname(MAXNUM_INPUT_DATA)
character(MAXLEN_NAME) :: input_vdesc(MAXNUM_INPUT_DATA)
character(MAXLEN_NAME) :: input_vunit(MAXNUM_INPUT_DATA)

real(dp), allocatable :: input_data(:, :)

INTEGER :: col_temp, col_RH, col_pres, col_ws, col_rglob
INTEGER :: col_J_1, col_J_3, col_J_4, col_J_5, col_J_6, col_J_7, col_J_12
INTEGER :: col_fix(NFIX)  ! column numbers for fixed compounds
INTEGER :: col_O3, col_NO, col_NO2, col_SO2, col_CO  ! included in col_fix, but explicitly declared for convenience
INTEGER :: col_HONO, col_HCHO, col_CH3COCH3
INTEGER :: col_pm25, col_pm10

! Current row of input data
INTEGER :: input_line


!------------------!
! Flags of different processes
!------------------!
logical :: l_dry_deposition
logical :: l_wall_process
logical :: l_inflow
logical :: l_outflow
logical :: l_emission

!------------------!
! Others
!------------------!

INTEGER :: i, j
integer :: stat
logical :: l_debug  ! if in debug mode, show more message if true
integer :: tmp_fid

!------------------!
! Name list
!------------------!

NAMELIST /NML_MAIN/ input_dir, output_dir, &
  l_debug, l_dry_deposition, l_wall_process, l_inflow, l_outflow, l_emission
NAMELIST /NML_STATION/ station, latitude_deg, longitude_deg, &
  chamber_volume, flow_rate_in_Lm, flow_rate_out_Lm
NAMELIST /NML_TIME/ time_start, time_end, dt, dt_display, dt_output, &
                    year_start, month_start, day_start, &
                    hour_start, minute_start, second_start
NAMELIST /NML_INPUT/ input_fname, input_ncol


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! Start program
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

! Set default values
call set_default_values()

! Read input from name list in the init file
call read_input()

! Display the run information
call display_general_information()

!------------------!
! Initial state
!------------------!

! First line is the initial time
! Suppose time is time1 and time2 in line1 and line2, the current time should be:
! time1 <= time <= time2
! Find line1, line2 should be line1 + 1.
! Here input_line is actually line1.
if (col('time') >= 1) then
  input_line = MINLOC( ABS(time-input_data(:, col('time'))), DIM=1, &
    MASK=( time-input_data(:, col('time')) >= 0.0d0 ) )
else
  call print_message('Input file does not contain a time column.', LERROR)
  stop
end if

! Time
CALL update_time()

! Environment conditions
call init_environment()

! Update meteo by reading from input
call update_environment()

!----- Chemistry -----!

! Photolysis reaction rates and other reaction rates
call init_rconst()

CALL update_rconst()

! Column indices of compounds
! col_O3  = find_string_location(input_vname, 'O3')
! col_NO  = find_string_location(input_vname, 'NO')
! col_NO2 = find_string_location(input_vname, 'NO2')
! col_SO2 = find_string_location(input_vname, 'SO2')
! col_CO  = find_string_location(input_vname, 'CO')
! col_HONO     = find_string_location(input_vname, 'HONO')
! col_HCHO     = find_string_location(input_vname, 'HCHO')
! col_CH3COCH3 = find_string_location(input_vname, 'CH3COCH3')

! Column indices of compounds
l_has_fix = (NVAR < NSPEC)
IF (l_has_fix) THEN
  DO i = 1, NFIX
    col_fix(i) = col( TRIM(ADJUSTL(SPC_NAMES(NVAR+i))) )
  END DO
END IF

! Initial concentration
conc = 0.0d0
CALL update_input_concentration()

! Peroxy radical
RO2 = 0.0d0

! Condensation sink
cs = 0.0d0
CALL update_condensation_sink(cs)

! Gas dry deposition velocity
if (l_dry_deposition) then
  ! Lu2013, ACP: 0.012 m s-1 for a typical BLH 1000 m
  ! Zhu2020, ACP: 0.01 m s-1 for BLH 400 m to 1400 m
  zref = 1000.0  ! [m], reference height for deposition, use a typical BLH here
  vd = 0.0d0
  CALL update_dry_deposition(vd)
end if

! Wall process
if (l_wall_process) then
  call update_wall_process(wall_loss)
end if

! Aerosol
! col_pm25 = find_string_location(input_vname, 'PM2.5')
! col_pm10 = find_string_location(input_vname, 'PM10')

! Initialize KPP
CALL kpp_setup()

! Calculate some variables before output

! conc_ovocs
!  CALL update_gas_group(conc, &
!    (/ind_ACR, ind_C2H5CHO, ind_CH3COCH3, ind_MTBE, ind_MACR, &
!    ind_C3H7CHO, ind_MVK, ind_MEK, ind_MPRK, ind_C4H9CHO, &
!    ind_DIEK, ind_C5H11CHO/), &
!    conc_ovocs)

WRITE(message,'(A7,F15.3,A9)'), 'time = ', time, '  seconds'
call print_message(trim(message))

! Open a data file and write down the initial values

! Output initial values to files
call open_output_file()
call write_output_file()


!==================!
! Start time loop
!==================!

call print_message('Starting main loop ...')

DO WHILE (time < time_end)

  !----- Chemistry integration -----!
  call print_message('Calculate chemistry ...', LDEBUG)

  CALL CHEMISTRY(conc, time, time+dt, &
    temp, pres, &
    O2, N2, MO2N2, H2O, RO2, &
    cs_H2SO4, cs_HNO3, &
    J_values)

  time = time + dt

  IF (time >= time2) THEN
    input_line = input_line + 1
    CALL update_time()
  END IF

  !----- Emission -----!
  if (l_emission) then
    call print_message('Update and apply emission ...', LDEBUG)
    call update_emission(emis)
    call apply_emission(conc, emis, dt)
  end if

  !----- Flow in for chamber -----!
  if (l_inflow) then
    call print_message('Update and apply inflow ...', LDEBUG)
    call update_inflow(inflow, conc_inflow, flow_rate_in, chamber_volume)
    call apply_inflow(conc, inflow, dt)
  end if

  !----- Consider a general condensation sink for OVOCs -----!
  call print_message('Update and apply condensation sink ...', LDEBUG)
  ! CALL condensation_ovocs(conc, &
  !   (/ind_ACR, ind_C2H5CHO, ind_MACR, ind_C3H7CHO, ind_MVK, &
  !   ind_MEK, ind_MPRK, ind_C4H9CHO, ind_DIEK, ind_C5H11CHO/), &
  !   (/1.0d-6, 2.0d-5, 1.0d-4, 1.0d-5, 1.0d-4, &
  !   1.0d-6, 1.0d-4, 0.0d0, 2.0d-5, 0.0d0/), &
  !   dt)
  CALL update_condensation_sink(cs)
  CALL apply_condensation_sink(conc, cs, dt)

  !----- Dry deposition of all variable concentrations: 1:NVAR -----!
  if (l_dry_deposition) then
    call print_message('Update and apply dry deposition ...', LDEBUG)
    call update_dry_deposition(vd)
    CALL apply_dry_deposition(conc, vd, zref, dt)
  end if

  !----- Wall process -----!
  if (l_wall_process) then
    call print_message('Update and apply wall processes ...', LDEBUG)
    call update_wall_process(wall_loss)
    call apply_wall_process(conc, wall_loss, dt)
  end if

  !----- Dilution -----!

  !----- Flow out for chamber -----!
  if (l_outflow) then
    call print_message('Update and apply outflow ...', LDEBUG)
    ! Notice that the outflow concentration is the same as the current concentration
    ! in the chamber.
    call update_outflow(outflow, conc, flow_rate_out, chamber_volume)
    call apply_outflow(conc, outflow, dt)
  end if

  !----- Update environment -----!
  call print_message('Update environment ...', LDEBUG)
  call update_environment()

  !----- Update J_values -----!
  call print_message('Update J values ...', LDEBUG)
  CALL update_rconst()

  !----- Gas concentrations -----!
  call print_message('Update input concentrations ...', LDEBUG)
  CALL update_input_concentration()

  !----- Calculate some variables before output -----!

  ! conc_ovocs
  ! CALL update_gas_group(conc, &
  !   (/ind_ACR, ind_C2H5CHO, ind_CH3COCH3, ind_MTBE, ind_MACR, &
  !   ind_C3H7CHO, ind_MVK, ind_MEK, ind_MPRK, ind_C4H9CHO, &
  !   ind_DIEK, ind_C5H11CHO/), &
  !   conc_ovocs)

  ! time_in_ms = FLOOR(1000*time)

  if ( MODULO(FLOOR(time*1000), FLOOR(dt_display*1000)) == 0 ) then
    write(message,'(A7,F15.3,A9)'), 'time = ', time, '  seconds'
    call print_message(trim(message))
  end if

  !----- Write output every dt_output -----!
  IF ( MODULO(FLOOR(time*1000), FLOOR(dt_output*1000)) == 0 ) THEN
    ! Write data to output files
    call write_output_file()

    ! Print message of saving data
    write(message,'(A7,F15.3,A9)'), 'time = ', time, '  seconds, data saved'
    call print_message(trim(message))
  END IF

END DO  ! DO WHILE (time < time_end)

call close_output_file()


CONTAINS


subroutine set_default_values()
  year_start = 2000
  month_start = 1
  day_start = 1
  hour_start = 0
  minute_start = 0
  second_start = 0

  l_dry_deposition = .true.
  l_wall_process = .true.
  l_inflow = .true.
  l_outflow = .true.
  l_emission = .true.
end subroutine set_default_values


subroutine read_input()
  character(MAXLEN_FILE_NAME) :: fname_init

  integer :: i, nline

  real(dp), allocatable :: file_buffer(:, :)

  integer :: row1, row2

  integer :: header_line_count

  logical :: exists

  ! Get the first argument of the executable command
  CALL GETARG(1, fname_init)
  call print_message('Reading namelists from '//trim(adjustl(fname_init)))

  OPEN(UNIT=99, FILE=TRIM(ADJUSTL(fname_init)), STATUS='OLD')
  READ(UNIT=99, NML=NML_MAIN   , IOSTAT=stat); REWIND(99)  ! directories and station
  READ(UNIT=99, NML=NML_STATION, IOSTAT=stat); REWIND(99)  ! input data
  READ(UNIT=99, NML=NML_INPUT  , IOSTAT=stat); REWIND(99)  ! input header list
  READ(UNIT=99, NML=NML_TIME   , IOSTAT=stat); REWIND(99)  ! start and end time, time steps
  CLOSE(99)

  !----- Check if everything is OK -----!

  ! time_end should be larger than time_start
  if (time_start >= time_end) then
    write(*,*) '[-ERROR  -] time_end is smaller than time_start.'
    stop
  end if

  ! Check if the output directory path exists, if not, create folders needed.
  INQUIRE(FILE=TRIM(ADJUSTL(output_dir)) // '/.', EXIST=exists)
  IF (.NOT. exists) THEN
    CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(output_dir)))
  END IF

  !-----------------------------!
  ! Calculate some parameters
  !-----------------------------!

  time = time_start
  year = year_start
  month = month_start
  day = day_start
  hour = hour_start
  minute = minute_start
  second = second_start

  latitude_rad = latitude_deg * pi/180.0d0
  longitude_rad = longitude_deg * pi/180.0d0

  !-----------------------------!
  ! Read time and data from input files
  !-----------------------------!
  header_line_count = 2

  ! Get input file name
  input_fname = trim(adjustl(input_dir)) // &
      '/station/' // trim(adjustl(station)) // &
      '/' // trim(adjustl(input_fname))

  ! Get file line count
  call get_file_line_count(input_fname, input_nline)
  input_nrow = input_nline - header_line_count  ! do not count name line and unit line
  allocate(input_data(input_nrow, input_ncol))
  write(message, '(a, i10, i5)') 'nrow, ncol: ', input_nrow, input_ncol
  call print_message(trim(message))

  ! Open the input file and read the data
  OPEN(UNIT=newunit(input_fid), FILE=TRIM(ADJUSTL(input_fname)), STATUS='old')
  REWIND(input_fid)
  select case (header_line_count)
    case (3)
      READ(input_fid, *) input_vname(1:input_ncol)
      READ(input_fid, *) input_vdesc(1:input_ncol)
      READ(input_fid, *) input_vunit(1:input_ncol)
    case (2)
      READ(input_fid, *) input_vname(1:input_ncol)
      READ(input_fid, *) input_vunit(1:input_ncol)
    case default
      call print_message('Wrong header lines.', LWARNING)
  end select

  ! Check the input names and their units
  do i = 1, input_ncol
    write(message, '(i5, a)') i, &
      ' '//trim(adjustl(input_vname(i)))//' | '//trim(adjustl(input_vunit(i)))
    call print_message(trim(message))
  end do

  ! Read the input data
  READ(input_fid, *) ( (input_data(i, j), j=1, input_ncol), i=1, input_nrow )

  ! Close the input file
  CLOSE(UNIT=input_fid)

  !-----------------------------!
  ! Set variables related to flow in and out
  !-----------------------------!
  flow_rate_in = flow_rate_in_Lm * 1.0d-3/60.0d0  ! [L min] --> [m3 s-1]
  flow_rate_out = flow_rate_out_Lm * 1.0d-3/60.0d0  ! [L min] --> [m3 s-1]

  conc_inflow_air = 2.4d19  ! [molec cm-3]
  conc_inflow = 0.0d0
  conc_inflow(ind_DMS) = 1.0d3/15.0d0 * ppm * conc_inflow_air * 5.0d-3/250.0d0  ! [molec cm-3]
end subroutine read_input


! Count the non-empty lines
subroutine get_file_line_count(fname, n)
  character(len=*), intent(in) :: fname
  integer, intent(out) :: n

  character(len=1000) :: rawline
  integer :: lun
  integer :: stat

  ! initial value
  n = 0

  ! open the file
  open(unit=newunit(lun), file=trim(adjustl(fname)))
  rewind(lun)

  ! count the line number
  do
    read(lun, '(a)', iostat=stat) rawline  ! read every line including empty lines
    if (stat/=0) exit                   ! exit if eof
    if (len(trim(rawline)) > 0) n=n+1          ! add one line count if the line is not empty
  end do
  
  ! close the file
  close(lun)
end subroutine get_file_line_count


! This is a simple function to search for an available unit.
! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
! The UNIT value is returned by the function, and also by the optional
! argument. This allows the function to be used directly in an OPEN
! statement, and optionally save the result in a local variable.
! If no units are available, -1 is returned.
!
! use as:
! integer lun  ! the new unit will be saved in lun
! open(unit=newunit(lun),file='test')
!
! reference: http://fortranwiki.org/fortran/show/newunit
integer function newunit(unit)
  integer, intent(out), optional :: unit
! local
  integer, parameter :: LUN_MIN=10, LUN_MAX=1000
  logical :: opened
  integer :: lun
! begin
  newunit=-1
  do lun=LUN_MIN,LUN_MAX
    inquire(unit=lun,opened=opened)
    if (.not. opened) then
      newunit=lun
      exit
    end if
  end do
  if (present(unit)) unit=newunit
end function newunit


subroutine display_general_information()
  call print_message('Input folder: '//trim(input_dir))
  call print_message('Output folder: '//trim(output_dir))
  call print_message('Station: '//trim(station))
  call print_message('')
  call print_message('Model simulation info: ')

  write(message, '(a, f15.3, a, f15.3)') 'Simulation time: ', time_start, ' --> ', time_end
  call print_message(trim(message))

  write(message, '(a, f15.3)') 'Time step: ', dt
  call print_message(trim(message))

  write(message, '(a, f15.3)') 'Time step of display: ', dt_display
  call print_message(trim(message))

  write(message, '(a, f15.3)') 'Time step of output: ', dt_output
  call print_message(trim(message))
end subroutine display_general_information


! SUBROUTINE read_input(input_fname, input_vname, input_vdesc, input_vunit, input_data)
!   CHARACTER(MAXLEN_FILE_NAME ) :: input_fname
!   CHARACTER(MAXLEN_NAME) :: input_vname(200)
!   CHARACTER(MAXLEN_NAME) :: input_vdesc(200)
!   CHARACTER(MAXLEN_NAME) :: input_vunit(200)
!   REAL(dp) :: input_data(1000, 200)
! 
!   INTEGER :: nrow_data, ncol_data
! 
!   nrow_data = 768
!   ncol_data = 157
!   OPEN(UNIT=1000, FILE=TRIM(ADJUSTL(input_fname)), STATUS='old')
!   READ(1000, *) input_vname(1:ncol_data)
!   READ(1000, *) input_vdesc(1:ncol_data)
!   READ(1000, *) input_vunit(1:ncol_data)
!   READ(1000, *) ( (input_data(i, j), j=1, ncol_data), i=1, nrow_data )
!   CLOSE(UNIT=1000)
! END SUBROUTINE read_input


! Find the first occuracne of str in the string array str_arr
FUNCTION find_string_location(str_arr, str) RESULT(loc)
  CHARACTER(MAXLEN_NAME), INTENT(IN) :: str_arr(200)
  CHARACTER(*), INTENT(IN) :: str

  INTEGER :: loc

  INTEGER :: i, nstr

  ! Default negative value represents not found
  loc = -1

  ! Size of the string array
  nstr = SIZE(str_arr)

  DO i = 1, nstr
    IF ( TRIM(ADJUSTL(str_arr(i))) == TRIM(ADJUSTL(str)) ) THEN
      loc = i
      EXIT
    END IF
  END DO
END FUNCTION find_string_location


! Interpolate between two time points
SUBROUTINE interpolate_from_input(input_data, input_line, data_name, time_now, data_now)
  REAL(dp), INTENT(IN) :: input_data(:, :)
  INTEGER, INTENT(IN) :: input_line
  character(*), intent(in) :: data_name
  REAL(dp), INTENT(IN) :: time_now
  REAL(dp), INTENT(OUT) :: data_now

  integer :: data_col
  REAL(dp) :: time1, time2, data1, data2

  data_col = col(data_name)
  if (data_col < 0) then
    call print_message(data_name//' can not be found, interpolation cancelled.', LWARNING)
    return
  end if

  ! The first column is the seconds since
  time1 = input_data(input_line  , col('time'))
  time2 = input_data(input_line+1, col('time'))

  ! Get the data at two time points
  data1 = input_data(input_line  , data_col)
  data2 = input_data(input_line+1, data_col)

  ! Interpolate
  data_now = data1 + (data2 - data1)/(time2 - time1) * (time_now - time1)
end subroutine interpolate_from_input


subroutine update_time()
  time1 = input_data(input_line  , col('time'))
  time2 = input_data(input_line+1, col('time'))

  ! Date info
  if (col('year') > 0) then
    year1 = input_data(input_line  , col('year'))
    year2 = input_data(input_line+1, col('year'))
    year = year1
  end if
  if (col('month') > 0) then
    month1 = input_data(input_line  , col('month'))
    month2 = input_data(input_line+1, col('month'))
    month = month1
  end if
  if (col('day') > 0) then
    day1 = input_data(input_line  , col('day'))
    day2 = input_data(input_line+1, col('day'))
    day = day1
  end if
  if (col('hour') > 0) then
    hour1 = input_data(input_line  , col('hour'))
    hour2 = input_data(input_line+1, col('hour'))
    hour = hour1
  end if
  if (col('minute') > 0) then
    minute1 = input_data(input_line  , col('minute'))
    minute2 = input_data(input_line+1, col('minute'))
    minute = minute1
  end if
  if (col('second') > 0) then
    second1 = input_data(input_line  , col('second'))
    second2 = input_data(input_line+1, col('second'))
    second = second1
  end if
end subroutine update_time


subroutine init_environment()
  pres = P00
  RH = 0.7d0
  Rglob = 400.0d0
  albedo_ground = 0.2d0

  zenith_deg = 45.0d0
  sun_par = 1.0d0
  cs_H2SO4 = 1.0d-3
  cs_HNO3 = 1.0d-3
  rhov = 10.0d-3
  Rglob = 400.0d0
end subroutine init_environment


SUBROUTINE update_environment()
  ! wind speed [m s-1]
  ! CALL interpolate_from_input(input_data, input_line, 'ws', time, ws)

  ! temp [K]
  CALL interpolate_from_input(input_data, input_line, 'T', time, temp)
  if ( trim(adjustl(input_vunit(col('T')))) == 'degC' .or. &
    trim(adjustl(input_vunit(col('T')))) == 'C' ) temp = temp + T00  ! [degC] -> [K]

  ! RH [-]
  CALL interpolate_from_input(input_data, input_line, 'RH', time, RH)
  if ( trim(adjustl(input_vunit(col('RH')))) == '%' )  RH = RH / 100.0d0  ! [%] -> [-]

  ! air pressure [Pa]
  ! CALL interpolate_from_input(input_data, input_line, 'pres', time, pres)
  ! if ( trim(adjustl(input_vunit(col('pres')))) == 'hPa' )  pres = pres * 100.0d0  ! [hPa] -> [Pa]

  ! Global radiation [W m-2]
  ! CALL interpolate_from_input(input_data, input_line, 'Rglob', time, Rglob)
  ! Rglob = MAX(Rglob, 0.0d0)  ! remove negative values

  ! Short wave radiation spectral distribution
  OPEN(unit=newunit(tmp_fid), &
  FILE=TRIM(ADJUSTL(input_dir))//'/station/'//TRIM(ADJUSTL(station))// &
    '/info/swr_distribution.txt', STATUS='old')
  READ(tmp_fid, *) (swr_distribution(I), I=1,84)
  CLOSE(tmp_fid)

  OPEN(unit=newunit(tmp_fid), &
  FILE=TRIM(ADJUSTL(input_dir))//'/station/'//TRIM(ADJUSTL(station))// &
    '/info/UVH_spectral.dat', STATUS='old')
  do i = 1, NWL
    READ(tmp_fid, *) tmp_wl(i), lamp_spectral(i, 1)
  end do
  CLOSE(tmp_fid)

  OPEN(unit=newunit(tmp_fid), &
  FILE=TRIM(ADJUSTL(input_dir))//'/station/'//TRIM(ADJUSTL(station))// &
    '/info/UVX_spectral.dat', STATUS='old')
  do i = 1, NWL
    READ(tmp_fid, *) tmp_wl(i), lamp_spectral(i, 2)
  end do
  CLOSE(tmp_fid)

  Rglob_spec = Rglob * swr_distribution

  ! cos(zenith)
  julian_day = JulianDay(INT(year), INT(month), INT(day))
  time_of_day = hour1*3600.0d0 + minute*60.0d0 + second + (time - time1)
  cos_zenith = cos_solar_zenith_angle(julian_day, time_of_day, latitude_rad)
  zenith_rad = ACOS(cos_zenith)
  zenith_deg = zenith_rad*180.0d0/pi

  ! Concentration of air, H2O, etc.
  Nair = pres*NA / (Rgas*temp) * 1e-6_dp ! [molec cm-3], number concentration of air
  O2 = 0.21d0*Nair  ! [molec cm-3]
  N2 = 0.78d0*Nair  ! [molec cm-3]
  MO2N2 = N2 + O2  ! [molec cm-3]
  H2O = svp(temp-T00)*RH*NA/(Rgas*temp)*1.0d-6  ! [molec cm-3]
END SUBROUTINE update_environment


subroutine init_rconst()
  J_values = 0.0d0
end subroutine init_rconst


SUBROUTINE update_rconst()
  ! Calculate J values 
  ! call update_rconst_J(input_dir, station, zenith_deg, actinic_flux, temp, J_values)
  ! Rglob = 400.0d0
  ! call update_rconst_J(input_dir, station, zenith_deg, Rglob_spec(1:NWL), temp, J_values)
  call update_rconst_J(input_dir, station, zenith_deg, lamp_spectral(1:NWL, 2)*3.0d0/2.386d0, temp, J_values)
  ! write(*,*) 'J(1:4): ', J_values(1:4)
  ! stop

  ! J_1: O3 -> O2 + O1D
  call interpolate_from_input( input_data, input_line, 'J_1', &
    time, J_values( 1) )

  ! J_2: O3 -> O
  ! CALL interpolate_from_input( input_data, input_line, 'J_2', &
  !   time, J_values( 2) )

  ! J_3: H2O2 -> 2 OH
  ! CALL interpolate_from_input( input_data, input_line, 'J_3', &
  !   time, J_values( 3) )

  ! J_4: NO2 -> NO + O3P
  call interpolate_from_input( input_data, input_line, 'J_4', &
    time, J_values( 4) )

  ! J_5: NO3 -> NO + O2
  ! CALL interpolate_from_input( input_data, input_line, 'J_5', &
  !   time, J_values( 5) )

  ! J_6: NO3 -> NO2 + O3P
  ! CALL interpolate_from_input( input_data, input_line, 'J_6', &
  !   time, J_values( 6) )

  ! J_7: HONO -> OH + NO
  ! CALL interpolate_from_input( input_data, input_line, 'J_7', &
  !   time, J_values( 7) )

  ! J_12: HCHO → H2 + CO
  ! CALL interpolate_from_input( input_data, input_line, 'J_12', &
  !   time, J_values(12) )

  call update_chemistry_properties(temp, pres, rh, cs_aq_wall, cw_aq_wall)
  write(*,*) 'cs_aq_wall, MSA, DMS: ', cs_aq_wall(50), cs_aq_wall(33)
  call update_rconst_aq(cs_aq_wall, cw_aq_wall, temp)
END SUBROUTINE update_rconst


SUBROUTINE update_input_concentration()
  INTEGER :: i

  ! Update fix concentration species
  ! if (l_has_fix) then
  !   do i = 1, NFIX
  !     call interpolate_from_input( input_data, input_line, SPC_NAMES(NVAR+i), time, conc(NVAR+i) )
  !     conc(NVAR+i) = conc(NVAR+i) * Nair * ppb
  !   end do
  ! end if

  ! Other input species

  ! HONO
  ! conc(ind_HONO) = 5.0d0*Nair*ppb
!  CALL interpolate_from_input( input_data, input_line, 'HONO', time, conc(ind_HONO) )
  ! write(*,*) 'HONO', conc(ind_HONO)
  ! conc(ind_HONO) = conc(ind_HONO) * Nair * ppb

  ! HCHO
  ! conc(ind_HCHO) = 10.0d0*Nair*ppb
!  CALL interpolate_from_input( input_data, input_line, 'HCHO', time, conc(ind_HCHO) )
  ! write(*,*) 'HCHO', conc(ind_HCHO)
  ! conc(ind_HCHO) = conc(ind_HCHO) * Nair * ppb

  ! CH3COCH3
!  CALL interpolate_from_input( input_data, input_line, 'CH3COCH3', time, conc(ind_CH3COCH3) )
  ! write(*,*) 'CH3COCH3', conc(ind_CH3COCH3)
  ! conc(ind_CH3COCH3) = conc(ind_CH3COCH3) * Nair * ppb

  ! CH4
  ! conc(ind_CH4) = 2240.0 * Nair * ppb

  ! O3
  call interpolate_from_input( input_data, input_line, 'O3', time, conc(ind_O3) )
  if ( trim(adjustl(input_vunit(col('O3')))) == 'ppb' ) then
    conc(ind_O3) = conc(ind_O3) * Nair * ppb
  end if

  ! CO
  call interpolate_from_input( input_data, input_line, 'CO', time, conc(ind_CO) )
  if ( trim(adjustl(input_vunit(col('CO')))) == 'ppb' ) then
    conc(ind_CO) = conc(ind_CO) * Nair * ppb
  end if

  ! NO2
  call interpolate_from_input( input_data, input_line, 'NO2', time, conc(ind_NO2) )
  if ( trim(adjustl(input_vunit(col('NO2')))) == 'ppb' ) then
    conc(ind_NO2) = conc(ind_NO2) * Nair * ppb
  end if

  ! DMS
  ! call interpolate_from_input( input_data, input_line, 'DMS', time, conc(ind_DMS) )
  ! if ( trim(adjustl(input_vunit(col('DMS')))) == 'ppb' ) then
  !   conc(ind_DMS) = conc(ind_DMS) * Nair * ppb
  ! end if

  !===== Special part =====!
  ! CO in Tan2017 is about half of ours
  ! conc(ind_CO) = 0.5d0*conc(ind_CO)

  ! NO in Tan2017 is very low during afternoon
  ! conc(ind_NO) = 0.2d0*conc(ind_NO)
END SUBROUTINE update_input_concentration


SUBROUTINE update_gas_group(conc, ind_list, conc_group)
  REAL(dp), INTENT(IN) :: conc(NSPEC)
  INTEGER, INTENT(IN) :: ind_list(:)
  REAL(dp), INTENT(OUT) :: conc_group

  INTEGER :: nind

  nind = SIZE(ind_list)
  conc_group = SUM(conc(ind_list))
END SUBROUTINE update_gas_group


SUBROUTINE condensation_ovocs(conc, ind_list, cs_list, dt)
  REAL(dp), INTENT(INOUT) :: conc(NSPEC)
  INTEGER , INTENT(IN   ) :: ind_list(:)
  REAL(dp), INTENT(IN   ) :: cs_list(:)
  REAL(dp), INTENT(IN   ) :: dt

  INTEGER :: i, nind

  nind = SIZE(ind_list)
  DO i = 1, nind
    conc(ind_list(i)) = conc(ind_list(i)) * (1.0d0 - cs_list(i)*dt)
  END DO
END SUBROUTINE condensation_ovocs


subroutine update_condensation_sink(cs)
  real(dp), intent(inout) :: cs(NSPEC)

  cs = 0.0d0
end subroutine update_condensation_sink


SUBROUTINE apply_condensation_sink(conc, cs, dt)
  REAL(dp), INTENT(INOUT) :: conc(:)
  REAL(dp), INTENT(IN   ) :: cs(:)
  REAL(dp), INTENT(IN   ) :: dt

  conc(:) = conc(:) * (1.0d0 - cs(:)*dt)
END SUBROUTINE apply_condensation_sink


SUBROUTINE update_dry_deposition(vd)
  REAL(dp), INTENT(INOUT) :: vd(:)

  vd(1:NVAR) = 0.0d0  ! [s-1]
END SUBROUTINE update_dry_deposition


SUBROUTINE apply_dry_deposition(conc, vd, zref, dt)
  REAL(dp), INTENT(INOUT) :: conc(:)
  REAL(dp), INTENT(IN   ) :: vd(:)
  REAL(dp), INTENT(IN   ) :: zref
  REAL(dp), INTENT(IN   ) :: dt

  ! conc(:) = conc(:) * (1.0d0 - vd(:)*dt/zref)
  conc(:) = conc(:) * exp(- vd(:)*dt)
END SUBROUTINE apply_dry_deposition


SUBROUTINE update_wall_process(wall_loss)
  REAL(dp), INTENT(INOUT) :: wall_loss(:)

!! vd(1:NVAR) = 0.01d0  ! [s-1]
  ! wall_loss(1:NVAR) = 1.20d-4  ! [s-1]
  wall_loss(ind_MSA) = 2.28d-3  ! [s-1]
  wall_loss(ind_H2SO4) =2.13d-3 
  wall_loss(ind_OH) =2.13d-3 
  wall_loss(ind_O1D) =2.13d-3 
  wall_loss(ind_SO3) =2.13d-3 
  wall_loss(ind_NO3) =2.13d-3 
  wall_loss(ind_HSO3) =2.13d-3 
  wall_loss(ind_N2O5) =2.13d-3 
  wall_loss(ind_HODMSO2) =2.13d-3 
  wall_loss(ind_CH3SOO2) =2.13d-3 
  wall_loss(ind_HO2NO2) =2.13d-3 
  wall_loss(ind_CH3SCH2O) =2.13d-3 
  wall_loss(ind_CH3SOO) =2.13d-3 
  wall_loss(ind_CH3SO3) =2.13d-3 
  wall_loss(ind_HNO3) =2.13d-3 
  wall_loss(ind_CH3SCH2O2) =2.13d-3 
  wall_loss(ind_CH3SO2) =2.13d-3 
  wall_loss(ind_CH3SO) =2.13d-3 
  wall_loss(ind_CH3S) =2.13d-3 
  wall_loss(ind_CH3O2) =2.13d-3 
  wall_loss(ind_HO2) =2.13d-3 
  wall_loss(ind_CH3SO2O2) =2.13d-3 
  wall_loss(ind_DMSO2O2) =2.13d-3 ! for the new_dms_NO these two compounds are CH3SO2CH2O2
  wall_loss(ind_DMSO2O) =2.13d-3 ! for the new_dms_NO these two compounds are CH3SO2CH2O
  wall_loss(ind_CH3SO3) =2.13d-3 
 ! NEW_DMS_NO
 ! wall_loss(ind_CH3SO2CH2O2)=2.13d-3
 ! wall_loss(ind_CH3SO2CH2O)=2.13d-3
 ! dms
  wall_loss(ind_CH3SOHCH3) =2.13d-3 
  ! new
  wall_loss(ind_HOOCH2SCH2O2) =2.13d-3 
  wall_loss(ind_HOOCH2SCO) =2.13d-3 
  wall_loss(ind_HOOCH2S) =2.13d-3 

  ! wall_loss = 0.0d0
END SUBROUTINE update_wall_process


SUBROUTINE apply_wall_process(conc, wall_loss, dt)
  REAL(dp), INTENT(INOUT) :: conc(:)
  REAL(dp), INTENT(IN   ) :: wall_loss(:)
  REAL(dp), INTENT(IN   ) :: dt

  ! conc(:) = conc(:) * (1.0d0 - vd(:)*dt/zref)
  conc(:) = conc(:) * exp(- wall_loss(:)*dt)
END SUBROUTINE apply_wall_process


subroutine update_emission(emis)
  real(dp), intent(out) :: emis(NSPEC)

  emis = 0.0d0  ! [molec cm-3 s-1]
end subroutine update_emission


subroutine apply_emission(conc, emis, dt)
  real(dp), intent(inout) :: conc(NSPEC)
  real(dp), intent(in   ) :: emis(NSPEC)
  real(dp), intent(in   ) :: dt

  conc = conc + emis * dt
end subroutine apply_emission


subroutine update_inflow(inflow, conc_inflow, flow_rate_in, chamber_volume)
  real(dp), intent(  out) :: inflow(NSPEC)  ! [molec cm-3 s-1], inflow source

  real(dp), intent(in   ) :: conc_inflow(NSPEC)  ! [molec cm-3 s-1], inflow source concentration
  real(dp), intent(in   ) :: flow_rate_in  ! [m3 s-1], inflow and outflow rates
  real(dp), intent(in   ) :: chamber_volume  ! [m3], chamber volume  

  inflow = conc_inflow * flow_rate_in / chamber_volume  ! [molec cm-3 s-1]
end subroutine update_inflow


subroutine apply_inflow(conc, inflow, dt)
  real(dp), intent(  out) :: conc(NSPEC)
  real(dp), intent(in   ) :: inflow(NSPEC)
  real(dp), intent(in   ) :: dt

  conc = conc + inflow * dt
end subroutine apply_inflow


subroutine update_outflow(outflow, conc_outflow, flow_rate_out, chamber_volume)
  real(dp), intent(  out) :: outflow(NSPEC)  ! [molec cm-3 s-1], inflow source

  real(dp), intent(in   ) :: conc_outflow(NSPEC)  ! [molec cm-3 s-1], inflow source concentration
  real(dp), intent(in   ) :: flow_rate_out  ! [m3 s-1], inflow and outflow rates
  real(dp), intent(in   ) :: chamber_volume  ! [m3], chamber volume  

  outflow = conc_outflow * flow_rate_out / chamber_volume  ! [molec cm-3 s-1]
  ! write(*,*) 'flow_rate_out/chamber_volume = ', flow_rate_out / chamber_volume
end subroutine update_outflow


subroutine apply_outflow(conc, outflow, dt)
  real(dp), intent(  out) :: conc(NSPEC)
  real(dp), intent(in   ) :: outflow(NSPEC)
  real(dp), intent(in   ) :: dt

  conc = conc - outflow*dt
end subroutine apply_outflow



FUNCTION svp(TC)
  REAL(dp) :: TC  ! [degC]
  REAL(dp) :: svp  ! [Pa]

  REAL(dp), PARAMETER :: a0 = 6.107799961,     &
                         a1 = 4.436518524D-1,  &
                         a2 = 1.428945805D-2,  &
                         a3 = 2.650648471D-4,  &
                         a4 = 3.031240396D-6,  &
                         a5 = 2.034080948D-8,  &
                         a6 = 6.136820929D-11
  svp = (a0 + a1*TC + a2*TC**2 + a3*TC**3 + a4*TC**4 + a5*TC**5 + a6*TC**6) * 1.0d2
END FUNCTION svp


FUNCTION cos_solar_zenith_angle(julian_day, time_of_day, latitude) RESULT(cos_zenith)
  ! INPUT
  REAL(dp), INTENT(IN) :: time_of_day  ! [s], time of day
  REAL(dp), INTENT(IN) :: latitude     ! [rad], Geographic latitude
  INTEGER , INTENT(IN) :: julian_day   ! [day], Julian day (day of year), 1 -> 365

  ! OUTPUT
  REAL(dp) :: declination  ! solar declination
  REAL(dp) :: cos_zenith  ! cos(zenith angle)
  
  ! LOCAL
  REAL(dp) :: solar_hour                 ! (=zeit) solar hour angle

  ! solar declination, after Tjernstrom (1989)
  ! julian_day - 79.8 should be equal to calendar day (J), ! with 1 at March 23 (according to Sogachev et al, 2002)
  ! (Equation 4)
  declination = 0.409d0 * SIN(2.0d0*pi * (julian_day - 79.8d0)/365.24d0)
  ! wsun=0.409*sin(pi2*(day-79.8)/365.24)  ! Henderson-Sellers

  ! solar hour angle, don't know why calculated like this
  ! (Equation 5)
  solar_hour = pi*(time_of_day/HALFDAY_SEC - 1.0d0)

  ! sin(sun height angle) = cos(zenith angle) = ...  (Equation below eq. A1 in Sogachev et al.2002)
  ! (Equation 2)
  cos_zenith = SIN(latitude) * SIN(declination) + COS(latitude) * COS(declination) * COS(solar_hour)
  cos_zenith = max(1.d-2, cos_zenith)
END FUNCTION cos_solar_zenith_angle


INTEGER FUNCTION JulianDay(y, m, d)
  INTEGER, INTENT(IN) :: y, m, d
  INTEGER, PARAMETER :: md(12) = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
  INTEGER :: i

  JulianDay = 0
  !** first month
  IF (m==1) THEN
    JulianDay = d
  !** other months
  ELSE
    DO i = 1, m-1
      JulianDay = JulianDay + md(i)
    END DO
    JulianDay = JulianDay + d
    !** add 1 day for leap years
    IF (m>2 .AND. IsLeapYear(y)) THEN
      JulianDay = JulianDay + 1
    END IF
  END IF
END FUNCTION JulianDay


LOGICAL FUNCTION IsLeapYear(y)
  INTEGER, INTENT(IN) :: y

  IF (MOD(y, 100) /= 0 .AND. MOD(y, 4) == 0) THEN 
    IsLeapYear = .TRUE.
  ELSEIF (MOD(y, 400) == 0) THEN
    IsLeapYear = .TRUE.
  ELSE 
    IsLeapYear = .FALSE.
  ENDIF
END FUNCTION IsLeapYear


! A shortcut of the function find_string_location with the str_arr set to input_vname
function col(str) result(loc)
  character(*), intent(in) :: str

  integer :: loc

  loc = find_string_location(input_vname, str)

  ! if (loc < 0) then
  !   write(*,*) '[-WARNING-] ', str, ' is not found in the input names.'
  ! end if

end function col


subroutine print_message(message, level)
  character(*), intent(in)      :: message
  integer, intent(in), optional :: level

  integer :: level_

  level_ = LINFO
  if (present(level)) level_ = level

  select case (level_)
    case (LINFO)
      write(*,*) '[-INFO   -] ', message
    case (LWARNING)
      write(*,*) '[-WARNING-] ', message
    case (LDEBUG)
      if (l_debug) then
        write(*,*) '[-DEBUG  -] ', message
      end if
    case (LERROR)
      write(*,*) '[-ERROR  -] ', message
    case default
      write(*,*) '[-???????-] ', message
  end select
end subroutine print_message


subroutine open_output_file()
  OPEN(11, file=TRIM(ADJUSTL(output_dir))//"/conc_ovocs.dat", &
    status='replace', action='write')

  OPEN(12, file=TRIM(ADJUSTL(output_dir))//"/conc_inorg.dat", &
    status='replace', action='write')

  OPEN(16, file=TRIM(ADJUSTL(output_dir))//"/conc_others.dat", &
    status='replace', action='write')

  OPEN(17, file=TRIM(ADJUSTL(output_dir))//"/meteo.dat", &
    status='replace', action='write')

  OPEN(18, file=TRIM(ADJUSTL(output_dir))//"/J_values.dat", &
    status='replace', action='write')
end subroutine open_output_file


subroutine write_output_file()
  ! Concentration of OVOCs
  ! conc=conc/(Nair*ppb)
  ! conc_ovocs=conc_ovocs/(Nair*ppb)
!  WRITE(11,*) time, conc(ind_ACR), conc(ind_C2H5CHO), conc(ind_CH3COCH3), conc(ind_MTBE), conc(ind_MACR), &
!    conc(ind_C3H7CHO), conc(ind_MVK), conc(ind_MEK), conc(ind_MPRK), conc(ind_C4H9CHO), &
!    conc(ind_DIEK), conc(ind_C5H11CHO), RO2, conc_ovocs

  ! Concentration of inorganic species
  WRITE(12,*) time, &
    conc(ind_OH), conc(ind_NO3), conc(ind_HONO), conc(ind_HCHO), conc(ind_HO2), &
    conc(ind_O3), conc(ind_NO), conc(ind_NO2), conc(ind_SO2), conc(ind_CO), &
    conc(ind_DMS), conc(ind_H2SO4), conc(ind_MSA)

  ! Concentration of others
  ! WRITE(16,*) time, conc(ind_CH4), conc(ind_CL), conc(ind_BENZENE), conc(ind_CH3CHO)

  ! Environment variables
  WRITE(17,*) time, Nair, MO2N2, N2, O2, H2O, Rglob, temp, pres, RH

  ! J values
  WRITE(18,*) time, J_values
end subroutine write_output_file


subroutine close_output_file()
  CLOSE(11)
  CLOSE(12)
  CLOSE(13)
  CLOSE(14)
  CLOSE(15)
  CLOSE(16)
  CLOSE(17)
  CLOSE(18)
end subroutine close_output_file

end program sosaa_box
