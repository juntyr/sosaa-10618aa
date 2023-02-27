PROGRAM chemistry_box
  USE second_Parameters  ! NSPEC, ind_XXX, dp, NFIX
  USE second_Monitor, ONLY: SPC_NAMES, EQN_NAMES
  USE second_Global, ONLY: NPHOT
  USE Chemistry_Mod  ! chemistry, kpp_setup, NPHOT

  IMPLICIT NONE

  ! INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,307)
  REAL(dp), PARAMETER :: pi = 2*ASIN(1.0d0)
  REAL(dp), PARAMETER :: DAY_SEC = 86400.0d0
  REAL(dp), PARAMETER :: HALFDAY_SEC = 43200.0d0
  REAL(dp), PARAMETER :: HOUR_SEC = 3600.0d0  ! one hour in seconds
  REAL(dp), PARAMETER :: MINUTE_SEC = 60.0d0  ! one minute in seconds
  INTEGER, PARAMETER :: MAX_STRING_LENGTH_COMMON = 50
  INTEGER, PARAMETER :: MAX_STRING_LENGTH_FNAME  = 100

  ! Station info
  CHARACTER(MAX_STRING_LENGTH_COMMON) :: station  ! station name, e.g., 'hyytiala'
  REAL(dp) :: latitude_deg
  REAL(dp) :: latitude_rad

  ! Meteorology
  REAL(dp) :: temp, RH, pres, ws
  REAL(dp) :: zenith_deg  ! [deg], solar zenith angle in degree
  REAL(dp) :: zenith_rad  ! [rad], solar zenith angle in radian
  REAL(dp) :: cos_zenith  ! [deg], cos(zenith_deg)
  REAL(dp) :: sun_par  ! [-], sunlit fraction, for box model 
  REAL(dp) :: cs_H2SO4, cs_HNO3  ! [s-1], condensation sinks for H2SO4 and HNO3
  REAL(dp) :: rhov  ! [kg m-3], absolute humidity, water vapor density
  REAL(dp) :: J_values(NPHOT)  ! photolysis reaction rates
  REAL(dp) :: K_values(NKVALUES)  ! reaction rates
  REAL(dp) :: Rglob  ! [W m-2], income global solar radiation
  CHARACTER(MAX_STRING_LENGTH_FNAME) :: input_dir  ! input dir, usually the path of 'sosa_in'
  CHARACTER(MAX_STRING_LENGTH_FNAME) :: output_dir  ! output dir, usually './output'
  REAL(dp) :: albedo  ! [-], surface albedo
  REAL(dp) :: actinic_flux(NWL)  ! [photon cm-2 s-1 nm-1], calculated actinic flux

  ! Chemistry
  REAL(dp) :: conc(NSPEC)
  REAL(dp) :: conc_ovocs

  REAL(dp) :: k
  REAL(dp) :: O2, N2, M, H2O, H2
  REAL(dp) :: RO2
  REAL(dp) :: I_ACF_in(NWL), I_ACF_out(NWL)

  !***** time variables *****!
  REAL(dp) :: time, time_start, time_end
  REAL(dp) :: year, month, day, hour, minute, second
  REAL(dp) :: time_of_day
  INTEGER :: julian_day
  REAL(dp) :: time1, time2
  REAL(dp) :: year1, year2 
  REAL(dp) :: month1, month2
  REAL(dp) :: day1, day2 
  REAL(dp) :: hour1, hour2 
  REAL(dp) :: minute1, minute2 
  REAL(dp) :: second1, second2 
  REAL(dp) :: dt                          ! time step, in seconds
  REAL(dp) :: display_dt, output_dt
  INTEGER  :: time_in_ms                  ! time in milliseconds

  ! Other parameters
  REAL(dp), PARAMETER :: ppm = 1.0d-6   ! [mol mol-1]
  REAL(dp), PARAMETER :: ppb = 1.0d-9   ! [mol mol-1]
  REAL(dp), PARAMETER :: ppt = 1.0d-12  ! [mol mol-1]
  REAL(dp), PARAMETER :: T00 = 273.15d0  ! [K]
  REAL(dp), PARAMETER :: P00 = 1.01325d5  ! [Pa]
  REAL(dp), PARAMETER :: Rgas = 8.314d0  ! [J mol-1 K-1]
  REAL(dp), PARAMETER :: NA = 6.022d23  ! [molec mol-1]

  ! Other variables
  REAL(dp) :: Nair = 2.25d19  ! [molec cm-3]

  INTEGER :: i, j

  ! About input data
  CHARACTER(MAX_STRING_LENGTH_FNAME ) :: raw_fname
  CHARACTER(MAX_STRING_LENGTH_COMMON) :: raw_name(200)
  CHARACTER(MAX_STRING_LENGTH_COMMON) :: raw_desc(200)
  CHARACTER(MAX_STRING_LENGTH_COMMON) :: raw_unit(200)
  REAL(dp) :: raw_data(1000, 200)

  INTEGER :: col_temp, col_RH, col_pres, col_ws, col_rglob
  INTEGER :: col_J_1, col_J_3, col_J_4, col_J_5, col_J_6, col_J_7, col_J_12
  INTEGER :: col_fix(NFIX)  ! column numbers for fixed compounds
  INTEGER :: col_O3, col_NO, col_NO2, col_SO2, col_CO  ! included in col_fix, but explicitly declared for convenience
  INTEGER :: col_HONO, col_HCHO, col_CH3COCH3
  INTEGER :: col_pm25, col_pm10

  ! Current row of input data
  INTEGER :: raw_line
  REAL(dp) :: raw_sec1, raw_sec2


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  ! Start program
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

  !------------------!
  ! Read input
  !------------------!
  input_dir = '/home/local/pzzhou/Scripts/Fortran/SOSAA/sosa_in'
  output_dir = './output'
  station = 'beijing_yang'

  raw_fname = TRIM(ADJUSTL(input_dir)) // '/station/' // TRIM(ADJUSTL(station)) // '/input_interp_data.txt'
  CALL read_input(raw_fname, raw_name, raw_desc, raw_unit, raw_data)


  !------------------!
  ! Set parameters
  !------------------!
  time_start = 0.0d0
  time_end = 60.0d0
  time = time_start
  dt = 10.0d0
  display_dt = 60.0d0
  output_dt = 600.0d0

!  station = 'beijing_yang'
  latitude_deg = 39.9042d0
  latitude_rad = latitude_deg * pi/180.0d0


  !------------------!
  ! Initial state
  !------------------!

  ! First line is the initial time
  raw_line = MINLOC( ABS(time-raw_data(:,1)), DIM=1, MASK=( time-raw_data(:,1) >= 0.0d0 ) )

  !!!!! Time
  CALL update_time()


  !!!!! Meteorology conditions
  col_ws   = find_string_location(raw_name, 'ws')
  col_temp = find_string_location(raw_name, 'T')
  col_RH   = find_string_location(raw_name, 'RH')
  col_pres = find_string_location(raw_name, 'Ps')
  col_rglob = find_string_location(raw_name, 'Rglob')

  sun_par = 1.0d0
  cs_H2SO4 = 1.0d-3
  cs_HNO3 = 1.0d-3
  rhov = 10.0d-3
  albedo = 0.3d0

  ! Update meteo by reading from input
  CALL update_meteorology()


  !!!!! Chemistry

  J_values = 0.0d0

  ! Column indices of compounds
  col_O3  = find_string_location(raw_name, 'O3')
  col_NO  = find_string_location(raw_name, 'NO')
  col_NO2 = find_string_location(raw_name, 'NO2')
  col_SO2 = find_string_location(raw_name, 'SO2')
  col_CO  = find_string_location(raw_name, 'CO')
  col_HONO     = find_string_location(raw_name, 'HONO')
  col_HCHO     = find_string_location(raw_name, 'HCHO')
  col_CH3COCH3 = find_string_location(raw_name, 'CH3COCH3')

  DO i = 1, NFIX
    col_fix(i) = find_string_location( raw_name, TRIM(ADJUSTL(SPC_NAMES(NVAR+i))) )
  END DO

  ! atmospheric O2, N2, H2O
  CALL update_basic_concentration()

  ! Initial concentration
  conc = 0.0d0
  CALL update_input_concentration()




  !!!!! Initialize KPP
  CALL kpp_setup()


  !!!!! Calculate some variables before output
  
  ! conc_ovocs
  CALL update_gas_group(conc, &
    (/ind_ACR, ind_C2H5CHO, ind_CH3COCH3, ind_MTBE, ind_MACR, &
    ind_C3H7CHO, ind_MVK, ind_MEK, ind_MPRK, ind_C4H9CHO, &
    ind_DIEK, ind_C5H11CHO/), &
    conc_ovocs)

  !!!!! Open a data file and write down the initial values
  WRITE(*,'(A8,F10.3,A9)'), 'time = ', time, '  seconds'

  OPEN(11,file=TRIM(ADJUSTL(output_dir))//"/conc_ovocs.dat",status='replace',action='write')
  WRITE(11,*) time, conc(ind_ACR), conc(ind_C2H5CHO), conc(ind_CH3COCH3), conc(ind_MTBE), conc(ind_MACR), &
    conc(ind_C3H7CHO), conc(ind_MVK), conc(ind_MEK), conc(ind_MPRK), conc(ind_C4H9CHO), &
    conc(ind_DIEK), conc(ind_C5H11CHO), conc_ovocs

  OPEN(12,file=TRIM(ADJUSTL(output_dir))//"/conc_inorg.dat",status='replace',action='write')
  WRITE(12,*) time, conc(ind_OH), conc(ind_NO3), conc(ind_HONO), conc(ind_HCHO), conc(ind_HO2), conc(ind_O3), conc(ind_NO), conc(ind_NO2), conc(ind_SO2), &
    conc(ind_CO)

  OPEN(16,file=TRIM(ADJUSTL(output_dir))//"/conc_BENZENE.dat",status='replace',action='write')
  WRITE(16,*) time, conc(ind_BENZENE)

  OPEN(17,file=TRIM(ADJUSTL(output_dir))//"/meteo.dat",status='replace',action='write')
  WRITE(17,*) time, Nair, M, N2, O2, H2O, Rglob, temp, pres, RH

  OPEN(18,file=TRIM(ADJUSTL(output_dir))//"/J_values.dat",status='replace',action='write')
  WRITE(18,*) time, J_values


  !------------------!
  ! Start time loop
  !------------------!
  DO WHILE (time < time_end)

    ! chemistry integration
    CALL CHEMISTRY(conc, time, time+dt, &
      temp, pres, zenith_deg, sun_par, cs_H2SO4, cs_HNO3, rhov, &
      H2O, RO2, &
      J_values, K_values, Rglob, &
      input_dir, station, albedo, actinic_flux)

    ! Consider a general condensation sink for OVOCs 
    CALL condensation_ovocs(conc, &
      (/ind_ACR, ind_C2H5CHO, ind_CH3COCH3, ind_MTBE, ind_MACR, &
      ind_C3H7CHO, ind_MVK, ind_MEK, ind_MPRK, ind_C4H9CHO, &
      ind_DIEK, ind_C5H11CHO/), &
      0.0d-4, dt)

    time = time + dt

    IF (time >= time2) THEN
      raw_line = raw_line + 1
      CALL update_time()
    END IF

    !
    ! Update meteorology
    !
    CALL update_meteorology()

    !
    ! Update chemistry
    !

    ! air number concentration
    CALL update_basic_concentration()

    ! gas concentrations
    CALL update_input_concentration()

    !
    ! Calculate some variables before output
    !

    ! conc_ovocs
    CALL update_gas_group(conc, &
      (/ind_ACR, ind_C2H5CHO, ind_CH3COCH3, ind_MTBE, ind_MACR, &
      ind_C3H7CHO, ind_MVK, ind_MEK, ind_MPRK, ind_C4H9CHO, &
      ind_DIEK, ind_C5H11CHO/), &
      conc_ovocs)

    !
    ! Write output every output_dt
    !

    ! time_in_ms = FLOOR(1000*time)
    IF ( MODULO(FLOOR(time*1000), FLOOR(display_dt*1000)) == 0 ) THEN
      WRITE(*,'(A8,F10.3,A9)', ADVANCE='NO') 'time = ', time, '  seconds'
      IF ( MODULO(FLOOR(time*1000), FLOOR(output_dt*1000)) == 0 ) THEN
        WRITE(*,*) ', data saved'
      ELSE
        WRITE(*,*)
      END IF
    END IF


    IF ( MODULO(FLOOR(time*1000), FLOOR(output_dt*1000)) == 0 ) THEN
      WRITE(11,*) time, conc(ind_ACR), conc(ind_C2H5CHO), conc(ind_CH3COCH3), conc(ind_MTBE), conc(ind_MACR), &
        conc(ind_C3H7CHO), conc(ind_MVK), conc(ind_MEK), conc(ind_MPRK), conc(ind_C4H9CHO), &
        conc(ind_DIEK), conc(ind_C5H11CHO), conc_ovocs
      WRITE(12,*) time, conc(ind_OH), conc(ind_NO3), conc(ind_HONO), conc(ind_HCHO), conc(ind_HO2), conc(ind_O3), conc(ind_NO), conc(ind_NO2), conc(ind_SO2), &
        conc(ind_CO)
      WRITE(16,*) time, conc(ind_BENZENE)
      WRITE(17,*) time, Nair, M, N2, O2, H2O, Rglob, temp, pres, RH
      WRITE(18,*) time, J_values
    END IF

  END DO  ! DO WHILE (time < time_end)

  CLOSE(11)
  CLOSE(12)
  CLOSE(13)
  CLOSE(14)
  CLOSE(15)
  CLOSE(16)
  CLOSE(17)
  CLOSE(18)


CONTAINS


SUBROUTINE read_input(raw_fname, raw_name, raw_desc, raw_unit, raw_data)
  CHARACTER(MAX_STRING_LENGTH_FNAME ) :: raw_fname
  CHARACTER(MAX_STRING_LENGTH_COMMON) :: raw_name(200)
  CHARACTER(MAX_STRING_LENGTH_COMMON) :: raw_desc(200)
  CHARACTER(MAX_STRING_LENGTH_COMMON) :: raw_unit(200)
  REAL(dp) :: raw_data(1000, 200)

  INTEGER :: nrow_data, ncol_data

  nrow_data = 768
  ncol_data = 157
  OPEN(UNIT=1000, FILE=TRIM(ADJUSTL(raw_fname)), STATUS='old')
  READ(1000, *) raw_name(1:ncol_data)
  READ(1000, *) raw_desc(1:ncol_data)
  READ(1000, *) raw_unit(1:ncol_data)
  READ(1000, *) ( (raw_data(i, j), j=1, ncol_data), i=1, nrow_data )
  CLOSE(UNIT=1000)

END SUBROUTINE read_input


! Find the first occuracne of str in the string array str_arr
FUNCTION find_string_location(str_arr, str) RESULT(loc)
  CHARACTER(MAX_STRING_LENGTH_COMMON), INTENT(IN) :: str_arr(200)
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
SUBROUTINE interpolate_from_raw(raw_data, raw_line, col, time_now, data_now)
  REAL(dp), INTENT(IN) :: raw_data(:, :)
  INTEGER, INTENT(IN) :: raw_line, col
  REAL(dp), INTENT(IN) :: time_now
  REAL(dp), INTENT(OUT) :: data_now

  REAL(dp) :: time1, time2, data1, data2

  
  ! The first column is the seconds since
  time1 = raw_data(raw_line  , 1)
  time2 = raw_data(raw_line+1, 1)

  ! Get the data at two time points
  data1 = raw_data(raw_line  , col)
  data2 = raw_data(raw_line+1, col)

  ! Interpolate
  data_now = data1 + (data2 - data1)/(time2 - time1) * (time_now - time1)
END SUBROUTINE interpolate_from_raw


SUBROUTINE update_time()
  ! time window
  time1 = raw_data(raw_line  , 1)
  time2 = raw_data(raw_line+1, 1)

  ! Date info
  year1 = raw_data(raw_line  , 2)
  year2 = raw_data(raw_line+1, 2)
  month1 = raw_data(raw_line  , 3)
  month2 = raw_data(raw_line+1, 3)
  day1 = raw_data(raw_line  , 4)
  day2 = raw_data(raw_line+1, 4)
  hour1 = raw_data(raw_line  , 5)
  hour2 = raw_data(raw_line+1, 5)
  minute1 = raw_data(raw_line  , 6)
  minute2 = raw_data(raw_line+1, 6)
  second1 = raw_data(raw_line  , 7)
  second2 = raw_data(raw_line+1, 7)

  ! temporary solution, need to be improved later
  year = year1
  month = month1
  day = day1
END SUBROUTINE update_time


SUBROUTINE update_meteorology()
  ! wind speed [m s-1]
  CALL interpolate_from_raw(raw_data, raw_line, col_ws, time, ws)

  ! temp [K]
  CALL interpolate_from_raw(raw_data, raw_line, col_temp, time, temp)
  temp = temp + T00  ! [degC] -> [K]

  ! RH [-]
  CALL interpolate_from_raw(raw_data, raw_line, col_RH, time, RH)
  RH = RH / 100.0d0  ! [%] -> [-]

  ! air pressure [Pa]
  CALL interpolate_from_raw(raw_data, raw_line, col_pres, time, pres)
  pres = pres * 100.0d0  ! [hPa] -> [Pa]

  ! Global radiation [W m-2]
  CALL interpolate_from_raw(raw_data, raw_line, col_rglob, time, Rglob)
  Rglob = MAX(Rglob, 0.0d0)  ! remove negative values

  ! cos(zenith)
  julian_day = JulianDay(INT(year), INT(month), INT(day))
  time_of_day = hour1*3600.0d0 + minute*60.0d0 + second + (time - time1)
  cos_zenith = cos_solar_zenith_angle(julian_day, time_of_day, latitude_rad)
  zenith_rad = ACOS(cos_zenith)
  zenith_deg = zenith_rad*180.0d0/pi
END SUBROUTINE update_meteorology


SUBROUTINE update_basic_concentration()
  Nair = pres*NA / (Rgas*temp) * 1e-6_dp ! [molec cm-3], number concentration of air
  O2 = 0.21d0*Nair  ! [molec cm-3]
  N2 = 0.78d0*Nair  ! [molec cm-3]
  M = N2 + O2  ! [molec cm-3]
  H2O = svp(temp-T00)*RH*NA/(Rgas*temp)*1.0d-6  ! [molec cm-3]
END SUBROUTINE update_basic_concentration


SUBROUTINE update_input_concentration()
  INTEGER :: i

  DO i = 1, NFIX
    CALL interpolate_from_raw( raw_data, raw_line, col_fix(i), time, conc(NVAR+i) )
    conc(NVAR+i) = conc(NVAR+i) * Nair * ppb
  END DO

  ! Other input species

  ! HONO
  CALL interpolate_from_raw( raw_data, raw_line, col_HONO, time, conc(ind_HONO) )
  conc(ind_HONO) = conc(ind_HONO) * Nair * ppb

  ! HCHO
  CALL interpolate_from_raw( raw_data, raw_line, col_HCHO, time, conc(ind_HCHO) )
  conc(ind_HCHO) = conc(ind_HCHO) * Nair * ppb

  ! CH3COCH3
  CALL interpolate_from_raw( raw_data, raw_line, col_CH3COCH3, time, conc(ind_CH3COCH3) )
  conc(ind_CH3COCH3) = conc(ind_CH3COCH3) * Nair * ppb
END SUBROUTINE update_input_concentration


SUBROUTINE update_gas_group(conc, ind_list, conc_group)
  REAL(dp), INTENT(IN) :: conc(NSPEC)
  INTEGER, INTENT(IN) :: ind_list(:)
  REAL(dp), INTENT(OUT) :: conc_group

  INTEGER :: nind

  nind = SIZE(ind_list)
  conc_group = SUM(conc(ind_list))
END SUBROUTINE update_gas_group


SUBROUTINE condensation_ovocs(conc, ind_list, cs, dt)
  REAL(dp), INTENT(INOUT) :: conc(NSPEC)
  INTEGER, INTENT(IN) :: ind_list(:)
  REAL(dp), INTENT(IN) :: cs, dt

  INTEGER :: i, nind

  nind = SIZE(ind_list)
  DO i = 1, nind
    conc(ind_list(i)) = conc(ind_list(i)) * (1.0d0 - cs*dt)
  END DO
END SUBROUTINE condensation_ovocs


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


END PROGRAM main_box
