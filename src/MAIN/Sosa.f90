!==============================================================================!
!
! Comment rules
!
! level 1:
!   !===...===! until 80
!   ! xxx
!   !===...===! until 80
!
! level 2:
!   !---...---! until 60
!   ! xxx
!   !---...---! until 60
!
! level 3: !===== xxx =====!
!
! level 4: !----- xxx
!
! level 5: ! xxx
!
! level 6: !>>, !>>>>, !>>>>>>, ...
!
! level 7: code  ! xxx
!
!==============================================================================!

program Sosa

  USE SOSA_DATA
  USE SOSA_IO
  USE MEGAN_version_1
  USE MEGAN_version_2
  USE SIMBIM
  USE LAI_Hyy_month
  USE Chemistry_Mod

  USE second_Monitor, ONLY: SPC_NAMES  ! a CHARACTER(LEN=12) array of chemical names
  USE psd_interface, psd_pi=>pi, psd_dp=>dp, psd_COLLISION_RATE=>COLLISION_RATE, psd_use_dmps=>use_dmps
  USE MT_MainMet
  USE Scadis_Initial
  USE gdd_function_mod
  use traj_proc_mod

  IMPLICIT NONE
  real(dp) :: cpu2,cpu1

  !============================================================================!
  ! Set up MPI
  !============================================================================!

!#ifdef PARALLEL

  ! Initialize MPI
  CALL MPI_Init(mpi_rc)

  IF (mpi_rc /= MPI_SUCCESS) THEN
    WRITE(*,*) 'MPI initialization failed'
    STOP
  END IF

  ! Get number of tasks, i.e. how many processes are running
  CALL MPI_Comm_size(MPI_COMM_WORLD, mpi_ntasks, mpi_rc)

  ! One task is master for coordinating, output and meteorology, and others are slaves for
  ! chemistry and aerosol calculations
  mpi_nslaves = mpi_ntasks - 1

  CALL MPI_Comm_rank(MPI_COMM_WORLD, my_id, mpi_rc)  ! find out our id number

  IF (my_id == master_id) THEN
    WRITE(*,*) 'MPI: number of processes: ', mpi_ntasks
    WRITE(*,*) 'The master ID is ', my_id
  END IF

  ! Terminte the program if the number of slaves is larger than kz
  IF (mpi_nslaves > kz) THEN
    WRITE (*,*) 'Error: too many slave processes, ', mpi_nslaves, &
      'for only kz = ', kz, ' atmospheric layers.'
    CALL MPI_Finalize(mpi_rc)
    stop
  END IF

  !============================================================================!
  ! Read input data for all cores
  !============================================================================!

  WRITE(*,*) 'The core ', my_id, ' is reading input data ...'

  CALL ACDC_INIT(SPC_NAMES, stage=1)
  CALL READ_INPUTS(stage=1)
  CALL ACDC_INIT(SPC_NAMES, stage=2)

  !============================================================================!
  ! Initialize for all cores
  !============================================================================!

  if (flag_chem == 1) then
    ! CALL chemistry_init(TRIM(CHEM_DIR)//'/INPUT', STATION)
    CALL PHOTOLYSIS_ARCA_INIT(TRIM(CHEM_DIR)//'/INPUT', STATION)

    WRITE(*,*) 'The core ', my_id, ' is initiating vapours ...'

    call VAPOUR_INIT(SPC_NAMES,CHEM_DIR) ! ARCA

  end if

  ! commented because outpt seems to save PSD structure even if flag_aero == 0. This isn't too expensive to mem.
  ! IF (flag_aero == 1) THEN

      dt_aero = dt_uhma ! P.C. remove later when UHMA is gone
      WRITE(*,*) 'The core ', my_id, ' is initiating aerosol ...'
      call AERO_INIT(NSPEC, kz) ! ARCA

  ! ELSE
      ! [Future] need to consider how to process it when no uhma is used
      ! NOTE XXX replace this BS with VAPOUR INIT
      ! CALL vapor_no_uhma(trim(adjustl(CHEM_DIR)) // '/INPUT/condensable_vapors', &
      ! VAPOUR_PROP%vapour_names,VAPOUR_PROP%psat_a,VAPOUR_PROP%psat_b)
  ! END IF
  CALL READ_INPUTS(stage=2)


  !============================================================================!
  ! The sub processors go to work
  !============================================================================!
  IF (my_id /= master_id) THEN

    CALL slave_task  ! if we are not the master branch, go to slave mode and stay there
    ! From here on the code, the loop etc. is only executed by the master code (until the end and finalization)
    CALL MPI_Finalize(mpi_rc)
    stop
  END IF

  ! OPEN(UNIT=1199, FILE=TRIM(OUTPUT_DIR)//'/parts_02.txt', status='replace',action='write')
!#endif

call cpu_time(cpu1)

  !============================================================================!
  ! Initiate time and grid
  !============================================================================!

  !--------------------------------------------------------!
  ! Time
  !--------------------------------------------------------!

  first_time_loop = .true.

  nxodrad = 1

  ! Time to start the simulation.
  time = 0.0d0  ! [s], time in model, starting from 0
  time_in_day = 3600.0d0*start_date(4)+60.0d0*start_date(5)  ! [s], current time_in_day
  time_in_month = time_in_day + (start_date(3)-1)*SECONDS_IN_ONE_DAY  ! [s], current month time
  aero_start = (aero_start_date(3)-1)*SECONDS_IN_ONE_DAY + SECONDS_IN_ONE_HOUR*aero_start_date(4) + SECONDS_IN_ONE_MINUTE*aero_start_date(5) + aero_start_date(6)
  ! Check if the start time is a new day and new month
  if ( (start_date(4)==0) .and. (start_date(5)==0) .and. (start_date(6)==0) ) then
    is_newday = .true.
  else
    is_newday = .false.
  end if

  if ( (start_date(3)==1) .and. (start_date(4)==0) &
    .and. (start_date(5)==0) .and. (start_date(6) == 0) ) then
    is_newmonth = .true.
  else
    is_newmonth = .false.
  end if

  ! The clock for the finishing time
  time_end = (end_date(3)-1)*SECONDS_IN_ONE_DAY + SECONDS_IN_ONE_HOUR*end_date(4) + SECONDS_IN_ONE_MINUTE*end_date(5) + end_date(6)
  nstep = 10000000 ! Number of steps. Actually it's needed to define error in calculations. So take 1000000 for instance

  ! Date information
  julian = JulianDay(now_date(1), now_date(2), now_date(3))
  date = now_date(3)
  nclo_end = MonthDay(now_date(1), now_date(2))*48
  ! mon = start_date(2)

  ! Set pointers for input data
  if (flag_model_type==TRAJECTORY_MODE) then
    input_ptr = 1
    ! input_meteo_dt = input_meteo_2(input_meteo_pointer+1, 2) - input_meteo_2(input_meteo_pointer, 2)
    ! input_emisant_pointer = 1
    ! input_emisant_dt = input_emisant(input_emisant_pointer+1, 2) - input_emisant(input_emisant_pointer, 2)
    ! input_emisbio_pointer = 1
  end if

  ! Timers
  do k = 1, kz
    call init_timer(root_recv_aero(k))
    call init_timer(root_send_aero(k))
  end do

  !--------------------------------------------------------!
  ! Grid
  !--------------------------------------------------------!
  call generate_grid('traj', kz, hh, z, abl)
  call calculate_grid_parameters()

  ! SOME VARIABLE PARAMETERS INSIDE THE LAYER BETWEEN GROUND AND DOMAIN TOP
  if (flag_model_type==STATION_MODE) then
    CALL CONVERT_ECMWF()
  end if


  !============================================================================!
  ! Initiate canopy information
  !============================================================================!
  ! 'curved' LAI for Hyytiala pine, spruce and hardwood over 15 years: 1996 - 2010 by Janne Levola (2011 same as 2010)
  ! EM_LAI_year_Hyy = (/2.62, 2.82, 3.00, 3.12, 3.29, 3.42, 3.44, 3.47, 3.47/)
  ! This is used for the emission part.
  ! Normally the total LAI is measured.
  ! The relationships for pine and spruce are:
  ! total LAI * 0.6 = curved LAI     &     total LAI = 2.7 * projected LAI  &  projected LAI = curved LAI / 1.62
  ! The relationships for birch are:
  ! total LAI * 0.5 = EM_LAI      &      total LAI * 0.5 = projected LAI
  ! EM_LAI_year_Hyy for:
  ! spruce 15 years = 4.079
  ! spruce 30 years = 8.0
  ! spruce 50 years = 7.335
  ! pine 15 years = 1.600
  ! pine 20 years = 4.26
  ! pine 50 years = 4.160
  ! birch 10 years = 0.863
  ! birch 20 years = 4.01
  ! birch 50 years = 2.74

  EM_LAI_year_Hyy = (/2.62_dp, 2.82_dp, 3.00_dp, 3.12_dp, 3.29_dp, 3.42_dp, 3.44_dp, 3.47_dp, 3.49_dp, 3.51_dp, 3.52_dp, 3.54_dp/)  ! curved LAI
  IF (now_date(1) == 2003) EM_LAI_year = EM_LAI_year_Hyy(1)
  IF (now_date(1) == 2004) EM_LAI_year = EM_LAI_year_Hyy(2)
  IF (now_date(1) == 2005) EM_LAI_year = EM_LAI_year_Hyy(3)
  IF (now_date(1) == 2006) EM_LAI_year = EM_LAI_year_Hyy(4)
  IF (now_date(1) == 2007) EM_LAI_year = EM_LAI_year_Hyy(5)
  IF (now_date(1) == 2008) EM_LAI_year = EM_LAI_year_Hyy(6)
  IF (now_date(1) == 2009) EM_LAI_year = EM_LAI_year_Hyy(7)
  IF (now_date(1) == 2010) EM_LAI_year = EM_LAI_year_Hyy(8)
  IF (now_date(1) == 2011) EM_LAI_year = EM_LAI_year_Hyy(9)
  IF (now_date(1) == 2012) EM_LAI_year = EM_LAI_year_Hyy(10)
  IF (now_date(1) == 2013) EM_LAI_year = EM_LAI_year_Hyy(11)
  IF (now_date(1) == 2014) EM_LAI_year = EM_LAI_year_Hyy(12)
  IF (now_date(1) == 2015) EM_LAI_year = EM_LAI_year_Hyy(12)
  IF (now_date(1) >= 2016) EM_LAI_year = EM_LAI_year_Hyy(12)

  IF (hc > 0.0d0) THEN  ! if canopy exists
    !LAI_proj is ONLY used for the met part and it is the projected LAI (so the flat part of the needle).
    !So MEGAN do not get this info!
    IF ((flag_tree == 1) .OR. (flag_tree == 2)) THEN  !pine & spruce
      LAI_proj = EM_LAI_year / 1.62
      ! LAI_proj = 3.  ! Try to set projected LAI to 3 [m2 m-2], since in met part it will be divided by 2.7.
      su=1
    ELSEIF (flag_tree == 3) THEN  !birch
      LAI_proj = EM_LAI_year
      su=2
    ELSEIF (flag_tree == 4) THEN  !clearcut
      LAI_proj = 0.1
      su=1
    ENDIF
    pa  = 3 ! For the meteorology
    pa2 = 6 ! only for megan or general emission activity

    ! Find the model level at the top of the canopy (higher than the canopy, not equal to)
    k_canopy = 1
    DO k=1,kz
      if (z(k) <= hc) then
        k_canopy = k + 1
      else
        exit
      end if
    ENDDO
  else  ! no canopy exists
    k_canopy = 0
  ENDIF

  ! Initial values for simbim for sunny (S) and shadow (C) leaves
  DO k = 1,kz
    BIM_IN_S(k,1)  = 30.0   ! start value for stomatal conductance
    BIM_IN_S(k,2)  = 380.0  ! start value for CI (CO2 inside)
    BIM_IN_S(k,3)  = 1.0    ! start value for photosynthesis rate
    BIM_IN_S(k,4)  = 10.0   ! start value PGA (phosphor-glycerate-aldehyde)
    BIM_IN_S(k,5)  = 30.0   ! start value GAP (
    BIM_IN_S(k,6)  = 100.0  ! start value and constant for NADPH ()
    BIM_IN_S(k,7)  = 0.0    ! start value DXP
    BIM_IN_S(k,8)  = 0.0    ! start value MEP
    BIM_IN_S(k,9)  = 0.0    ! start value IDP
    BIM_IN_S(k,10) = 0.0    ! start value DMADP
    BIM_IN_S(k,11) = 0.0    ! start value Isoprene production rate [nmol / m2-leaf area / s]
    BIM_IN_S(k,12) = 0.0    ! start value GDP
    BIM_IN_S(k,13) = 0.0    ! start value Monoterpene production rate [nmol / m2-leaf area / s]

    BIM_IN_C(k,1)  = 30.0   ! start value for stomata conductance
    BIM_IN_C(k,2)  = 380.0  ! start value for CI (CO2 inside)
    BIM_IN_C(k,3)  = 1.0    ! start value for photosynthesis rate
    BIM_IN_C(k,4)  = 10.0   ! start value PGA (phosphor-glycerate-aldehyde)
    BIM_IN_C(k,5)  = 30.0   ! start value GAP (
    BIM_IN_C(k,6)  = 100.0  ! start value and constant for NADPH ()
    BIM_IN_C(k,7)  = 0.0    ! start value DXP
    BIM_IN_C(k,8)  = 0.0    ! start value MEP
    BIM_IN_C(k,9)  = 0.0    ! start value IDP
    BIM_IN_C(k,10) = 0.0    ! start value DMADP
    BIM_IN_C(k,11) = 0.0    ! start value Isoprene production rate [nmol / m2-leaf area / s]
    BIM_IN_C(k,12) = 0.0    ! start value GDP
    BIM_IN_C(k,13) = 0.0    ! start value Monoterpene production rate [nmol / m2-leaf area / s]
  ENDDO


  !============================================================================!
  ! Initiate meteorology
  !============================================================================!
  WRITE(*,*) 'Initiate meteorology ...'
  CALL initiate_meteorology(flag_model_type)

  ! Update latitude and longitude in traj mode
  call update_location(lat_deg, lon_deg, flag_model_type)

  ! Solar zenith angle
  cos_zenith = get_cos_zenith( time_in_day_UTC, real(julian_UTC, dp), &
    lat_deg, lon_deg)
  zenith_deg = acos(cos_zenith)*rad2deg

  ! Update sunrise and sunset time
  ! New method to calculate the sunrise and sunset time, but it is not tested
  ! yet, so use the old one upsn and downsn now in update_meteorology_stat.
  ! call get_sunrise_and_sunset_time_UTC( &
  !   real(julian_UTC, dp), lat_deg, lon_deg, &
  !   sunrise_UTC, sunset_UTC)
  ! write(*,*) 'rads: ', rads
  ! write(*,*) 'sunrise, sunset old: ', upsn, downsn
  ! write(*,*) 'sunrise, sunset new: ', &
  !   sunrise_UTC+time_zone*SECONDS_IN_ONE_HOUR, sunset_UTC+time_zone*SECONDS_IN_ONE_HOUR


  !============================================================================!
  ! Initiate chemistry
  !============================================================================!
  WRITE(*,*) 'Initiate chemistry ...'
  ! everything concerning ##!! warning!!
  !#ifndef PARALLEL
  CALL KPP_SetUp()  ! This only called once for KPP in the beginning
  !#endif

  ! Set all gas-concentrations to 0 [molec cm-3] in the beginning
  CH_CONS_ALL(:, :) = 0.0d0

  ! Set initial value for ozone
  ! write(*,*) 'time_in_month, nxodrad, dt_obs', time_in_month, nxodrad, dt_obs
  ! CH_CONS_ALL(:,ind_O3)  = linear_interp(time_in_month, nxodrad, CH_gas_hyy(2,:), dt_obs) * 2.4d19*1.0d-9

  !---------------------------------------
  ! Count how many chemical names begin with 'rOH': CH_oh_count
  ! Set up a logical array, NSPEC long, so that only for rOH-pseudochemicals we have
  ! CH_oh_flag(j) == .TRUE., j being the number of chemical. Do the same for O3 and NO3.
  !---------------------------------------
  CALL init_rOH()
  CALL init_rO3()
  CALL init_rNO3()
  CALL init_reactivity()

  !============================================================================!
  ! Initialize emissions
  !============================================================================!
  WRITE(*,*) 'Initiate emission ...'
  ! Soil emission parameters
  I = 1
  Do WHILE (z(I) < 4.7_dp)
    I = I + 1
  END DO
  temp_level = I

  ! Set all the emission rate to zero initially
  emis = 0.0d0


  !============================================================================!
  ! Initialize gas dry deposition
  !============================================================================!
  o3_weight = 0.0d0
  o3_weight(nz) = 1.0d0  ! get from measurement in layer nz

  IF (flag_gasdrydep == 1) THEN
    WRITE(*,*) 'Initiate gas dry deposition ...'

    l_drydep = .True.
    l_vpd = .True.
    l_wetskin = .TRUE.

    ! ALLOCATE(frac_veg(1:nz-2), frac_ws(1:nz-2)
    ! frac_veg(2:nz-1) = 1.0d0
    ! frac_ws(2:nz-1) = 0.0d0

    LAIl = 0.0d0
    DO k = 2, nz-1
      ! LAIl(k) = s1(k)*(dz(k)+dz(k+1))*0.5d0 * 1.62d0  ! change from meteorology LAI to emission LAI, namely projected to curved
      LAIl(k) = s1(k)*(dz(k)+dz(k+1))*0.5d0 * 2.7d0  ! change from meteorology LAI to total LAI, namely projected to total
    END DO
    LAIl(2) = LAIl(2)/2.7d0*2.0d0  ! the lowest layer (understory) is the normal-leaf shrub, not needle leaves, should be 0.5 [m2 m-2]
    DO k = 2, nz-1
      LAIl_c(k) = 0.5d0*LAIl(k) + SUM(LAIl(k+1:nz-1))
    END DO

    ! LAIl_debug = 1.0d0  ! for debug

    PARl = 0.0d0

    u_veg = 0.0d0
    rho_veg = roa
    Tc_veg = T00
    stomblock = 0.0d0

    vdep = 0.0d0
    gasdep_flux = 0.0d0



    !===== Read molar_mass, Henry's law contants, dry reactivity factors from input file H_f0_mmass.inp =====!
    OPEN(UNIT=99, FILE=TRIM(ADJUSTL(CHEM_DIR))//'/INPUT/H_f0_mmass.inp', STATUS='OLD')
    READ(99, *)
    DO
      READ(99,*,IOSTAT=dummy_io) dummy_index, dummy_MCM, dummy_SMILES, &
        dummy_HenrySE, dummy_HenryEE, dummy_HenryEG, dummy_HenryEB, dummy_HFlag, dummy_HenryA, dummy_HenryB, &
        dummy_f0Flag, dummy_f0, dummy_formula, dummy_molar_mass

      IF (dummy_io > 0) THEN  ! error
        WRITE(*,*) 'Something is wrong with H_f0_mmass.inp.'
        STOP
      ELSE IF (dummy_io < 0) THEN  ! EOF
        EXIT
      ELSE  ! everything is OK
        trname(dummy_index) = dummy_MCM
        molar_mass(dummy_index) = dummy_molar_mass
        HenrySE(dummy_index) = dummy_HenrySE
        HenryEE(dummy_index) = dummy_HenryEE
        HenryEG(dummy_index) = dummy_HenryEG
        HenryEB(dummy_index) = dummy_HenryEB
        HenryA(dummy_index) = dummy_HenryA
        HenryB(dummy_index) = dummy_HenryB
        dryreac(dummy_index) = dummy_f0
      END IF
    END DO
    CLOSE(99)
    Henry = 0.0d0
  END IF


  !============================================================================!
  ! Initiate other variables
  !============================================================================!
  WRITE(*,*) 'Initiate other variables ...'
  ! Initial values for calculating accumulative sources and sinks
  CH_CONS_0(:, :) = 0.0d0
  CH_CONS_1(:, :) = 0.0d0
  Qconc(:,:) = 0.0d0
  Qemis(:,:) = 0.0d0
  Qchem(:,:) = 0.0d0
  Qturb(:,:) = 0.0d0
  Qdepo(:,:) = 0.0d0


  !============================================================================!
  ! Open output files and write the initial values
  !============================================================================!
  WRITE(*,*) 'Open output files and write the initial values ...'
  call output_init()
  WRITE(*,*) '2 Open output files and write the initial values ...'
  call output_step()
  WRITE(*,*) '3 Open output files and write the initial values ...'
  call output_write()
  WRITE(*,*) '4 Open output files and write the initial values ...'
  ! if (flag_aero==1) call output_write_old()

  !============================================================================!
  ! Basic cycle
  !============================================================================!
  ! time

  WRITE(*,*) 'Start main loop within one month ...'
  WRITE(*,*) 'Current time', time_in_month
  WRITE(*,*) 'Ending time', time_end


  !    SSSSSSSSSSSSSSS      OOOOOOOOO        SSSSSSSSSSSSSSS              AAA                              AAA
  !  SS:::::::::::::::S   OO:::::::::OO    SS:::::::::::::::S            A:::A                            A:::A
  ! S:::::SSSSSS::::::S OO:::::::::::::OO S:::::SSSSSS::::::S           A:::::A                          A:::::A
  ! S:::::S     SSSSSSSO:::::::OOO:::::::OS:::::S     SSSSSSS          A:::::::A                        A:::::::A
  ! S:::::S            O::::::O   O::::::OS:::::S                     A:::::::::A                      A:::::::::A
  ! S:::::S            O:::::O     O:::::OS:::::S                    A:::::A:::::A                    A:::::A:::::A
  !  S::::SSSS         O:::::O     O:::::O S::::SSSS                A:::::A A:::::A                  A:::::A A:::::A
  !   SS::::::SSSSS    O:::::O     O:::::O  SS::::::SSSSS          A:::::A   A:::::A                A:::::A   A:::::A
  !     SSS::::::::SS  O:::::O     O:::::O    SSS::::::::SS       A:::::A     A:::::A              A:::::A     A:::::A
  !        SSSSSS::::S O:::::O     O:::::O       SSSSSS::::S     A:::::AAAAAAAAA:::::A            A:::::AAAAAAAAA:::::A
  !             S:::::SO:::::O     O:::::O            S:::::S   A:::::::::::::::::::::A          A:::::::::::::::::::::A
  !             S:::::SO::::::O   O::::::O            S:::::S  A:::::AAAAAAAAAAAAA:::::A        A:::::AAAAAAAAAAAAA:::::A
  ! SSSSSSS     S:::::SO:::::::OOO:::::::OSSSSSSS     S:::::S A:::::A             A:::::A      A:::::A             A:::::A
  ! S::::::SSSSSS:::::S OO:::::::::::::OO S::::::SSSSSS:::::SA:::::A               A:::::A    A:::::A               A:::::A
  ! S:::::::::::::::SS    OO:::::::::OO   S:::::::::::::::SSA:::::A                 A:::::A  A:::::A                 A:::::A
  !  SSSSSSSSSSSSSSS        OOOOOOOOO      SSSSSSSSSSSSSSS AAAAAAA                   AAAAAAAAAAAAAA                   AAAAAAA


  DO WHILE (time_end > time_in_month)

    nxodrad = int(time_in_month/dt_obs)+1

    ! Output current time which is being calculated
    IF (MOD( INT(time_in_month), INT(SECONDS_IN_HALF_HOUR) ) == 0) THEN
      WRITE(*,'("=====  julian", I5, " | ", I4.4, ".", I2.2, ".", I2.2, "-", I2.2, ":", I2.2, " | halfhour", I4, "  hrsinmth =====", I5, F10.1, F10.2)') &
        julian, now_date(1), now_date(2), now_date(3), now_date(4), now_date(5), nxodrad, nclo_end, time_in_month, time_in_month/dt_obs
    END IF

    !***** Starting gas concentrations at every time step *****!
    CH_CONS_0 = CH_CONS_ALL


    !------------------------------------------------------!
    ! Meteorology part
    !------------------------------------------------------!

    call debug_message('DEBUGING: calculating meteorology ...')

    call update_meteorology(dt_obs, flag_model_type)


    !------------------------------------------------------!
    ! Update some meteorological variables
    !------------------------------------------------------!

    ! Horizontal wind velocity
    hwind = sqrt(u1**2 + v1**2)

    ! Number concentration of air, N2 and O2
    ! N/V = n*NA/V = NA * P/(RT) [molec m-3], *1e-6 --> [molec cm-3]
    air = pres / (ta1 * Rgas) * Avog * 1e-6_dp
    N2 = air * 0.78d0
    O2 = air * 0.21d0

    ! number concentration of H2O, [molec cm-3], notice that qa1 is absolute
    ! humidity (rhov)
    H2O = qa1 / xmm_H2O * Avog * 1.0d-6
    RH = (a0 + a1 * (ta1-273.15d0)**1 + a2 * (ta1-273.15d0)**2 + a3 * (ta1-273.15d0)**3    &
             + a4 * (ta1-273.15d0)**4 + a5 * (ta1-273.15d0)**5 + a6 * (ta1-273.15d0)**6)* 100
    RH = 100d0 * H2O / (RH/pres * air)

    ! ?? Need to reorganize the input for chemistry

    select case (flag_model_type)
    case (STATION_MODE)
      print*, "Temporarly commented out for clarity by P.C."
      ! ! Getting values from input data
      !
      ! !>> downward short wave radiation at 18 m
      ! glob = linear_interp(time_in_month, nxodrad, stat_glob, dt_obs)
      !
      ! !>> downward PAR
      ! PAR = linear_interp(time_in_month, nxodrad, stat_PAR, dt_obs)
      !
      ! !>> albedo
      ! albedo = linear_interp(time_in_month, nxodrad, stat_albedo, dt_obs)
      !
      ! !>> soil moisture averaged of -26->-36 cm and -38->-61 cm
      ! soil_moist_3 = linear_interp(time_in_month, nxodrad, stat_sm_3, dt_obs)
      ! if (isnan(soil_moist_3)) soil_moist_3 = 0.3d0  ! a temporary solution
      !
      ! !>> soil moisture averaged of -2->-6 cm and -14->-25 cm
      ! soil_moist_2 = linear_interp(time_in_month, nxodrad, stat_sm_2, dt_obs)
      !
      ! !>> soil moisture in organic layer 0-5 cm
      ! soil_moist_1 = linear_interp(time_in_month, nxodrad, stat_sm_1, dt_obs)

    case (TRAJECTORY_MODE)
        call linear_interp_2p( &
          infield%var(infield%id_time)%f1d(input_ptr), &
          infield%var(infield%id_time)%f1d(input_ptr+1), &
          MAX(0d0, infield%var(infield%id_ssr)%f1d(input_ptr)), &
          MAX(0d0, infield%var(infield%id_ssr)%f1d(input_ptr+1)), &
          time, glob &
          )

        call linear_interp_2p( &
          infield%var(infield%id_time)%f1d(input_ptr), &
          infield%var(infield%id_time)%f1d(input_ptr+1), &
          MAX(0d0, infield%var(infield%id_lsm)%f1d(input_ptr)), &
          MAX(0d0, infield%var(infield%id_lsm)%f1d(input_ptr+1)), &
          time, land_sea_mask &
          )

        ! interpolate CAMS O3 concentrations from pressure levels to model height levels
        call interp_time2p_level1d( &
        infield%var(infield%id_conc_lev)%f1d*1d2, infield%var(infield%id_time)%f1d, &
        infield%var(infield%id_conc_o3)%f2d, input_ptr, &
        pres(kz:1:-1), time, & ! Because interpolation function expects monotonically increasing x-values
        CH_CONS_ALL(kz:1:-1,ind_O3) & ! So we must flip the output vector
        )
        CH_CONS_ALL(:, ind_O3) = CH_CONS_ALL(:, ind_O3) * 1d-9 * air
      ! continue

    end select


    !------------------------------------------------------!
    ! Emission part
    !------------------------------------------------------!
    CALL debug_message('DEBUGING: calculating emissions ...')

    !===== Save the concentrations before emission and chemistry =====!
    CH_CONS_temp = CH_CONS_ALL

    ! Calculate emissions every dt_emis
    IF (MOD( INT(time_in_day), INT(dt_emis) ) == 0) THEN
      ! Set the emission array to zero every time before adding up all the
      ! emissions
      emis = 0.0d0
      emi_aer = 0.0d0

      select case (flag_model_type)
      case (STATION_MODE)
      print*, "Temporarly commented out for clarity by P.C."

      !   IF (flag_tree == 1) THEN !pine
      !     ! This is to interpolate the measured monoterpene tree emission rates
      !     ! + chamber temeperature (for SEP in MEGAN)
      !
      !     !>> measured monoterpene emission rate
      !     EMI_MEAS_MONO = linear_interp(time_in_month, nxodrad, EF_meas(1,:), dt_obs)
      !
      !     !>> chamber temperature ~ leaf temperature
      !     CHAM_temp = linear_interp(time_in_month, nxodrad, EF_meas(2,:), dt_obs)
      !   ELSE
      !     EMI_MEAS_MONO = 0.0d0
      !     CHAM_temp = 0.0d0
      !   ENDIF  ! flag_tree==1
      !
      !   ! Data transfer from main-program to emission variables - make sure
      !   ! all units are correct
      !   ! Emissions using MEGAN or SIMBIM
      !   ! That is the order of emisson factors by MEGAN and SIMBIM in EM_EMI
      !   ! ISOP                  1
      !   ! MYRC                  2
      !   ! SABI                  3
      !   ! LIMO                  4
      !   ! 3CAR                  5
      !   ! OCIM                  6
      !   ! BPIN                  7
      !   ! APIN                  8
      !   ! OMTP                  9
      !   ! FARN                 10
      !   ! BCAR                 11
      !   ! OSQT                 12
      !   ! MBO                  13
      !   ! MEOH (methanol)      14
      !   ! ACTO (acetone)       15
      !   ! CH4                  16
      !   ! NO                   17
      !   ! ACTA (acetaldehyde)  18
      !   ! FORM (formaldehyde)  19
      !   ! CO                   20
      !   ! CIN (cineole, Eucalyptol) 21
      !   ! LIN (Linalool)            22
      !   IF ((flag_emis == 1) .OR. (flag_emis == 2) .OR. (flag_emis == 3)) THEN
      !     ! The input data and other parameters are only read in or set at the first
      !     IF (first_call == 1) THEN ! do only once
      !       first_call = 0
      !     ENDIF  ! first_call
      !
      !     ! First put all emissions to zero at each time step
      !     EM_EMI(:,1:22) = 0.0d0
      !
      !     ! Input
      !     EM_DATE         = ((1000*now_date(1))+julian)     ! in format YYYDDD scalar
      !     EM_Time_M2      = Timeformatter(time_in_day)    ! Function in megan module: integer output - real input
      !
      !     EM_PAR_day      = 0.0d0                             ! Incoming PAR daily average
      !     ! EM_SRAD_day     = 0.                             ! Daily average short wave radiation (W/m2)
      !     EM_TempK_day    = ta1(2)      ! Daily average temperature Kelvin, use level 2 value as a temporary solution, need to be modified in future
      !
      !     ! Monoterpene distribution data from Jaana chemotypy paper average values
      !     EM_Mon_Proc(1)  = 0.010  !  MYRC
      !     EM_Mon_Proc(2)  = 0.010  !  SABI
      !     EM_Mon_Proc(3)  = 0.023  !  LIMO
      !     EM_Mon_Proc(4)  = 0.396  !  3CAR
      !     EM_Mon_Proc(5)  = 0.010  !  OCIM
      !     EM_Mon_Proc(6)  = 0.090  !  BPIN
      !     EM_Mon_Proc(7)  = 0.437  !  APIN
      !     EM_Mon_Proc(8)  = 0.024  !  OMTP
      !
      !     select case (flag_tree)
      !     case (1)  ! pine
      !       ! LAI_Megan is for Hyytiala monthly distribution
      !       EM_LAI = EM_LAI_year * LAI_Megan(Mon)/6.5728d0     ! One sided LAI (for needle leaf trees it is the projected area of the 3D needle)
      !       IF (mon /= 1) THEN
      !         EM_LAI_past = EM_LAI_year * LAI_Megan(Mon-1)/6.5728d0
      !       ELSE
      !         EM_LAI_past = EM_LAI_year * LAI_Megan(12)/6.5728d0
      !       END IF
      !     case (2)  !spruce
      !       EM_LAI = EM_LAI_year
      !       EM_LAI_past = EM_LAI
      !     case (3)  !birch
      !       EM_LAI = EM_LAI_year * LAI_function(julian)
      !       IF (julian /= 1) THEN
      !         EM_LAI_past = EM_LAI_year * LAI_function(julian-1)
      !       ELSE
      !         EM_LAI_past = EM_LAI_year * LAI_function(1)
      !       ENDIF
      !     end select
      !
      !     ! Output values set to zero for each run
      !     IF (flag_tree == 4) THEN
      !       EM_SUN_PAR = 1.0d0
      !     ENDIF
      !
      !     EM_SUN_PAR      = 0.0d0  ! Array of sun fraction - 1 above the canopy and decreases inside the canopy
      !     EM_SunleafTK    = 0.0d0  ! Array of temparture for sun leaf
      !     EM_ShadeleafTK  = 0.0d0  ! Array of temparture for shade leaf
      !     EM_Sunfrac      = 0.0d0  ! Array of the fraction of sun leaves. i = 1 is the top canopy layer, 2 is the next layer, etc.
      !     EM_EMI          = 0.0d0  ! Emissions in molecules per cm3 per second
      !     EM_Ea1pL_M1     = 0.0d0  ! Emission factor light
      !     EM_Ea1tL_M1     = 0.0d0  ! Emission factor temperature
      !     EM_Ea1NL_M1     = 0.0d0  ! Emission factor compined
      !     EM_Ea1pL_M2     = 0.0d0  ! Emission factor light
      !     EM_Ea1tL_M2     = 0.0d0  ! Emission factor temperature
      !     EM_Ea1NL_M2     = 0.0d0  ! Emission factor compined
      !     EM_GAM_TMP      = 0.0d0  ! Gamma factor for temperature, Non isoprene
      !     EM_GAM_OTHER    = 0.0d0  ! Gamma factors for other parameters
      !
      !     ! Layer list NOTE megan used the canopy for added accuracy so "cover"
      !     ! appears in every part of the canopy, cause it just changes a few
      !     ! variables in the calculations in comparison to the other plant
      !     ! functional types, (g/s)
      !     EM_ER_HB        = 0.0d0  ! emission rate from Herbaceous
      !     EM_ER_SB        = 0.0d0  ! emission rate from Shrubland(s) , see note above
      !     EM_ER_NT        = 0.0d0  ! emission rate from Needle Trees, - || -
      !     EM_ER_BT        = 0.0d0  ! emission rate from Broadleaf Trees, - || -
      !
      !     IF (flag_emis == 1) THEN ! old Megan code, not used anymore
      !       ! CALL EMISSION_M1( &
      !       !   ta1, PAR, hwind, kz, julian, &
      !       !   EM_EW, EM_LAI, zenith_deg, time_in_day, LAD_p, &
      !       !   k_canopy, z, EM_Ea1pL_M1, EM_Ea1tL_M1, EM_Ea1NL_M1, &
      !       !   EM_Emi, EM_SUN_PAR, EM_SunleafTK, EM_ShadeleafTK, EM_Sunfrac &
      !       !   )
      !     ELSEIF (flag_emis == 2) THEN ! new Megan code
      !       Call EMISSION_M2( &
      !         julian, k_canopy, zenith_deg, lat_deg, lon_deg, &
      !         EM_Date, EM_Time_M2, PAR, EM_PAR_day, ta1, &
      !         EM_TempK_day, pres(2), EM_WVM, hwind, soil_moist_3, &
      !         LAD_p, EM_LAI, EM_LAI_past, EM_ER, EM_VAR, &
      !         EM_ER_HB, EM_ER_SB, EM_ER_NT, EM_ER_BT, EM_Ea1pL_M2, &
      !         EM_Ea1NL_M2, EM_Ea1tL_M2, EM_GAM_TMP, EM_GAM_OTHER, z, &
      !         kz, EM_EMI, EM_SUN_PAR, EM_SunleafTK, EM_ShadeleafTK, &
      !         EM_Sunfrac, STATION, EMI_MEAS_MONO, CHAM_temp, flag_tree, &
      !         now_date(1),check_rh, tsn1, tsd1, flag_syn &
      !         )
      !
      !       ! DO JJ = 1, k_canopy
      !       !   EM_Ea1pL_M1(JJ) = EM_Ea1pL_M2(2,JJ)
      !       !   EM_Ea1tL_M1(JJ) = EM_Ea1tL_M2(2,JJ)
      !       !   EM_Ea1NL_M1(JJ) = EM_Ea1NL_M2(2,JJ)
      !       ! ENDDO
      !     ELSEIF (flag_emis == 3) THEN ! SIMBIM - only for monoterpenes
      !       ! If we use SIMBIM we have to run also MEGAN to get the leaf
      !       ! temperatures and the emission rates for non-monoterpenes
      !       Call EMISSION_M2( &
      !         julian, k_canopy, zenith_deg, lat_deg, lon_deg, &
      !         EM_Date, EM_Time_M2, PAR, EM_PAR_day, ta1, &
      !         EM_TempK_day, pres(2), EM_WVM, hwind, soil_moist_3, &
      !         LAD_p, EM_LAI, EM_LAI_past, EM_ER, EM_VAR, &
      !         EM_ER_HB, EM_ER_SB, EM_ER_NT, EM_ER_BT, EM_Ea1pL_M2, &
      !         EM_Ea1NL_M2, EM_Ea1tL_M2, EM_GAM_TMP, EM_GAM_OTHER, z, &
      !         kz, EM_EMI, EM_sun_par, EM_SunleafTK, EM_ShadeleafTK, &
      !         EM_Sunfrac, STATION, EMI_MEAS_MONO, CHAM_temp, flag_tree, &
      !         now_date(1), check_rh, tsn1, tsd1, flag_syn &
      !         )
      !
      !       DO JJ = 2, k_canopy
      !         ! Initial values for simbim for sunny (S) and shadow (C) leaves
      !         ! only at time equals zero
      !         IF (time_in_day == 0.0d0) THEN
      !           EM_BIM_S(JJ,1)  = 30.0   ! start value for stomatal conductance
      !           EM_BIM_S(JJ,2)  = 380.0  ! start value for CI (CO2 inside)
      !           EM_BIM_S(JJ,3)  = 1.0    ! start value for photosynthesis rate
      !           EM_BIM_S(JJ,4)  = 10.0   ! start value PGA (phosphor-glycerate-aldehyde)
      !           EM_BIM_S(JJ,5)  = 30.0   ! start value GAP (
      !           EM_BIM_S(JJ,6)  = 100.0  ! start value and constant for NADPH ()
      !           EM_BIM_S(JJ,7)  = 0.0    ! start value DXP
      !           EM_BIM_S(JJ,8)  = 0.0    ! start value MEP
      !           EM_BIM_S(JJ,9)  = 0.0    ! start value IDP
      !           EM_BIM_S(JJ,10) = 0.0    ! start value DMADP
      !           EM_BIM_S(JJ,11) = 0.0    ! start value Isoprene production rate [nmol / m2-leaf area / s]
      !           EM_BIM_S(JJ,12) = 0.0    ! start value GDP
      !           EM_BIM_S(JJ,13) = 0.0    ! start value Monoterpene production rate [nmol / m2-leaf area / s]
      !           EM_BIM_S(JJ,14) = 0.0    ! output [nmol / m2-leaf area / s]
      !
      !           EM_BIM_C(JJ,1)  = 30.0   ! start value for stomata conductance
      !           EM_BIM_C(JJ,2)  = 380.0  ! start value for CI (CO2 inside)
      !           EM_BIM_C(JJ,3)  = 1.0    ! start value for photosynthesis rate
      !           EM_BIM_C(JJ,4)  = 10.0   ! start value PGA (phosphor-glycerate-aldehyde)
      !           EM_BIM_C(JJ,5)  = 30.0   ! start value GAP (
      !           EM_BIM_C(JJ,6)  = 100.0  ! start value and constant for NADPH ()
      !           EM_BIM_C(JJ,7)  = 0.0    ! start value DXP
      !           EM_BIM_C(JJ,8)  = 0.0    ! start value MEP
      !           EM_BIM_C(JJ,9)  = 0.0    ! start value IDP
      !           EM_BIM_C(JJ,10) = 0.0    ! start value DMADP
      !           EM_BIM_C(JJ,11) = 0.0    ! start value Isoprene production rate [nmol / m2-leaf area / s]
      !           EM_BIM_C(JJ,12) = 0.0    ! start value GDP
      !           EM_BIM_C(JJ,13) = 0.0    ! start value Monoterpene production rate [nmol / m2-leaf area / s]
      !           EM_BIM_C(JJ,14) = 0.0    ! output [nmol / m2-leaf area / s]
      !         ENDIF
      !
      !         ! Calculation of VPD (vapour pressure deficit) in kPA
      !         EM_SunleafTC(JJ)   = EM_SunleafTK(JJ) - T00
      !         EM_ShadeleafTC(JJ) = EM_ShadeleafTK(JJ) - T00
      !
      !         EM_ES_S(JJ) = (a0 + a1 * EM_SunleafTC(JJ)**1 + &
      !           a2 * EM_SunleafTC(JJ)**2 + a3 * EM_SunleafTC(JJ)**3 &
      !           + a4 * EM_SunleafTC(JJ)**4 + a5 * EM_SunleafTC(JJ)**5 + &
      !           a6 * EM_SunleafTC(JJ)**6) * 100
      !
      !         EM_ES_C(JJ) = (a0 + a1 * EM_ShadeleafTC(JJ)**1 + &
      !           a2 * EM_ShadeleafTC(JJ)**2 + a3 * EM_ShadeleafTC(JJ)**3 &
      !           + a4 * EM_ShadeleafTC(JJ)**4 + a5 * EM_ShadeleafTC(JJ)**5 + &
      !           a6 * EM_ShadeleafTC(JJ)**6) * 100
      !
      !         EM_VPD_S(JJ) = (EM_ES_S(JJ) - EM_EW(JJ)) / 1000
      !         EM_VPD_C(JJ) = (EM_ES_C(JJ) - EM_EW(JJ)) / 1000
      !
      !         IF (EM_VPD_S(JJ) .GT. 3) THEN
      !           EM_VPD_S(JJ) = 3
      !         ENDIF
      !         IF (EM_VPD_C(JJ) .GT. 3) THEN
      !           EM_VPD_C(JJ) = 3
      !         ENDIF
      !
      !         ! start value for integral of monoterpenes [nmol / m2-leaf area
      !         ! - during the time step] has to be zero before call
      !         EM_BIM_S(JJ,14) = 0.0
      !         EM_BIM_C(JJ,14) = 0.0
      !
      !         ! Call of simbim 2 times for sunny and shadow leave temperature
      !         ! at each height level inside the canopy
      !         EM_Time_in    = time_in_day/60.
      !         EM_Time_out   = EM_Time_in + EM_DT/60.
      !
      !         CALL EMISSION_SB(EM_BIM_S(JJ,:), EM_Time_in, EM_Time_out, PAR, EM_VPD_S(JJ))
      !
      !         CALL EMISSION_SB(EM_BIM_C(JJ,:), EM_Time_in, EM_Time_out, PAR, EM_VPD_C(JJ))
      !
      !         ! Emission rates from SIMBIM for all monoterpenes
      !         DO II = 2,9
      !           EM_Emi(JJ,II)  = (EM_BIM_S(JJ,14) * EM_Sunfrac(k_canopy+1-JJ) + &
      !             EM_BIM_C(JJ,14) * (1-EM_Sunfrac(k_canopy+1-JJ)))     &
      !             * 1d-9 * Avog * 1d-4 * 6.7d0 * EM_Mon_Proc(II-1) / (z(JJ)-z(JJ-1)) / 100 * LAD_p(JJ) / 1000.
      !         ENDDO
      !
      !       END DO  ! JJ = 2, k_canopy
      !
      !       ! SIMBIM does not calculate emissions for the first layer so we use the second layer twice
      !       DO II = 2, 9
      !         EM_Emi(1,II)  = EM_Emi(2,II)
      !       ENDDO
      !
      !     ENDIF  ! IF flag_emis == 1
      !   ENDIF  ! IF flag_emis is 1 or 2 or 3
      !
      ! IF (flag_tree == 4) THEN
      !   EM_SUN_PAR(1:kz) = 1.0d0
      ! ENDIF
      !
      ! IF (flag_emis_soil == 1) THEN  ! This is to read in measured VOC soil emission
      !   ! Currently we only have high resolution soil emission data for
      !   ! HUMPPA-COPEC - so parts of july and August
      !   ! This is not completely true, cause Qingyang has more data.
      !   ! I will implement this (soon).
      !   IF (((mon .EQ. 7) .OR. (mon .EQ. 8)) .AND. (k .EQ. 2)) THEN
      !     ! read in soil emission data. Unit: molecules/cm^3/s
      !     open(unit=709, &
      !       file = TRIM(ADJUSTL(input_dir_station_data))//'/soil_aaltonen.txt', &
      !       status='old')
      !     read(709,*)((EM_soil_voc(i20,j20),i20=1,13),j20=1,1488)
      !     close(unit=709)
      !
      !     DO I=1,12 !we start at time 00.00
      !       EM_soil_emi(I)  = (EM_soil_voc(I+1,nxodrad) + (time_in_month-(nxodrad-1)*dt_obs) *  &
      !         (EM_soil_voc(I+1,nxodrad+1)-EM_soil_voc(I+1,nxodrad))/dt_obs)
      !       !Some emission values are negative, and since we already take deposition into account, we here exclude the negative values.
      !       IF (EM_soil_emi(I) .LT. 0.) THEN
      !         EM_soil_emi(I) = 0.
      !       ENDIF
      !     ENDDO
      !
      !     ! Set soil emissions for chemistry
      !     CALL set_soil_emissions_for_chemistry_1(CH_CONS_ALL(k, :), EM_soil_emi(:), &
      !       dt_emis, z(k)-z(k-1))
      !   ENDIF
      ! ENDIF  ! flag_emis_soil==1
      !
      ! IF (flag_emis_soil == 2) THEN !This is in order to include Hermanni's parameterisation.
      !   !This should only be used for M05-M11
      !   IF (((mon .EQ. 5) .OR. (mon .EQ. 6) .OR. (mon .EQ. 7) .OR. (mon .EQ. 8) .OR. (mon .EQ. 9) .OR.  &
      !     (mon .EQ. 10) .OR. (mon .EQ. 11)) .AND. (k .EQ. 2)) THEN
      !     IF (ta1(temp_level) .GT. 278.15) THEN
      !       M137_soil(k) = (((0.3 * (((ta1(temp_level)-T00) * soil_moist_2)**2)) - &
      !         (0.6 * ((ta1(temp_level)-T00) * soil_moist_2)) + 2) * conversion * 1/136) / ((z(k) - z(k-1))*100)
      !       M33_soil(k) = (((0.2 * (((ta1(temp_level)-T00) * soil_moist_1)**2)) - &
      !         (0.35 * ((ta1(temp_level)-T00) * soil_moist_1)) + 0.02) * conversion * 1/32) / ((z(k) - z(k-1))*100)
      !       M69_soil(k) = (((0.005 * (((ta1(temp_level)-T00) * soil_moist_1)**2)) + &
      !         (0.015 * ((ta1(temp_level)-T00) * soil_moist_1))) * conversion * 1/68) / ((z(k) - z(k-1))*100)
      !     ELSE
      !       M137_soil(k) = 0.
      !       M33_soil(k) = 0.
      !       M69_soil(k) = 0.
      !     ENDIF
      !
      !     CALL set_soil_emissions_for_chemistry_2(CH_CONS_ALL(k, :), &
      !       M137_soil(k), M33_soil(k), M69_soil(k), dt_emis)
      !   ENDIF
      ! ENDIF  ! flag_emis_soil==2
      !
      ! IF (flag_tree == 4) THEN !CLEARCUT
      !
      !   !Added night time effect:
      !   IF ((time_in_day .GE. 0) .AND. (time_in_day .LT. 5.0*60.0*60.0)) THEN
      !     sep_factor = 1.0/3.0
      !   ELSEIF ((time_in_day .GE. 5.0*60.0*60.0) .AND. (time_in_day .LE. 15.0*60.0*60.0)) THEN
      !     sep_factor = 2.0/(3.0*36000.0)*time_in_day
      !   ELSEIF ((time_in_day .GT. 15.0*60.0*60.0)  .AND. (time_in_day .LT. 24.0*60.0*60.0)) THEN
      !     sep_factor = (-2.0/97200.0 * time_in_day) + 2.1111
      !   ENDIF
      !
      !   SEP_clear = SAMI(julian,3) * (1.0/1.0D6) * (1.0/1.0D4) * (1.0/3600.0) * sep_factor !* 1.1 !2050
      !
      !   !Top canopy layer, assuming 30% of total emission:
      !   EM_Emi(12,8) = (0.3 * 0.437 * SEP_clear * Avog) / (136 * ((z(12)-z(11))* 100)) !A-pinene
      !   EM_Emi(12,7) = (0.3 * 0.090 * SEP_clear * Avog) / (136 * ((z(12)-z(11))* 100)) !b-pinene
      !   EM_Emi(12,5) = (0.3 * 0.396 * SEP_clear * Avog) / (136 * ((z(12)-z(11))* 100)) !carene
      !   EM_Emi(12,4) = (0.3 * 0.023 * SEP_clear * Avog) / (136 * ((z(12)-z(11))* 100)) !limonene
      !   EM_Emi(12,21) = (0.3 * 0.001 * SEP_clear * Avog) / (136 * ((z(12)-z(11))* 100)) !cineole
      !   EM_Emi(12,9) = (0.3 * 0.053 * SEP_clear * Avog) / (136 * ((z(12)-z(11))* 100)) !OMT
      !   CALL set_clearcut_emissions_for_chemistry(CH_CONS_ALL(12, :), EM_Emi(12, :), dt_emis)
      !   !Second top canopy layer, assuming 40% of total emission:
      !   EM_Emi(11,8) = (0.4 * 0.437 * SEP_clear * Avog) / (136 * ((z(11)-z(10))* 100)) !A-pinene
      !   EM_Emi(11,7) = (0.4 * 0.090 * SEP_clear * Avog) / (136 * ((z(11)-z(10))* 100)) !b-pinene
      !   EM_Emi(11,5) = (0.4 * 0.396 * SEP_clear * Avog) / (136 * ((z(11)-z(10))* 100)) !carene
      !   EM_Emi(11,4) = (0.4 * 0.023 * SEP_clear * Avog) / (136 * ((z(11)-z(10))* 100)) !limonene
      !   EM_Emi(11,21) = (0.4 * 0.001 * SEP_clear * Avog) / (136 * ((z(11)-z(10))* 100)) !cineole
      !   EM_Emi(11,9) = (0.4 * 0.053 * SEP_clear * Avog) / (136 * ((z(11)-z(10))* 100)) !OMT
      !   CALL set_clearcut_emissions_for_chemistry(CH_CONS_ALL(11, :), EM_Emi(11, :), dt_emis)
      !   !Third top canopy layer, assuming 30% of total emission:
      !   EM_Emi(10,8) = (0.3 * 0.437 * SEP_clear * Avog) / (136 * ((z(10)-z(9))* 100)) !A-pinene
      !   EM_Emi(10,7) = (0.3 * 0.090 * SEP_clear * Avog) / (136 * ((z(10)-z(9))* 100)) !b-pinene
      !   EM_Emi(10,5) = (0.3 * 0.396 * SEP_clear * Avog) / (136 * ((z(10)-z(9))* 100)) !carene
      !   EM_Emi(10,4) = (0.3 * 0.023 * SEP_clear * Avog) / (136 * ((z(10)-z(9))* 100)) !limonene
      !   EM_Emi(10,21) = (0.3 * 0.001 * SEP_clear * Avog) / (136 * ((z(10)-z(9))* 100)) !cineole
      !   EM_Emi(10,9) = (0.3 * 0.053 * SEP_clear * Avog) / (136 * ((z(10)-z(9))* 100)) !OMT
      !   CALL set_clearcut_emissions_for_chemistry(CH_CONS_ALL(10, :), EM_Emi(10, :), dt_emis)
      !
      ! ENDIF  ! IF (flag_tree == 4) THEN  ! CLEARCUT
      !
      ! ! Emissions are transfered to the correct KPP-form
      ! ! DO k = 1, kz
      !   ! The subroutine is different for different chemistry schemes, and you can
      !   ! write it as what you wish for a specific chemistry scheme
      !   ! CALL set_emissions_for_chemistry(CH_CONS_ALL(k, :), EM_Emi(k, :), dt_emis)
      ! ! END DO
      ! call set_emissions(EM_Emi, emis)

      case (TRAJECTORY_MODE)
        call traj_emission_calc(infield, input_ptr, z, time, emis, emi_aer)

        ! sea salt (spray) emissions normalized with dlog(Dp)
        call calculate_sea_salt_emissions(sea_salt_psd, hwind(4),ta1(1)-273.15, land_sea_mask)
        sea_salt_psd = sea_salt_psd/10d0 ! distributed to 10 meters

        ! inter(or extra)polation seems to sometimes produce negative values (CHECK FOR ERRORS...)
        where (emis<0d0) emis=0d0
        where (emi_aer<0d0) emi_aer=0d0

      end select
      ! Add emissions to number concentrations of species
      CH_CONS_ALL = CH_CONS_ALL + emis * dt_emis
      ! Particle emissions are mixed in the PSD_INTERFACE
    end if  ! emission every dt_emis


    !------------------------------------------------------!
    ! Chemistry part
    !------------------------------------------------------!

    CALL debug_message('DEBUGING: calculating chemistry ...')

    ! Calculate chemistry every dt_chem
    IF (MOD( INT(time_in_day), INT(dt_chem) ) == 0) THEN

      !===== CHEMISTRY solved with KPP =====!

      IF (flag_chem == 1) THEN

        CALL debug_message('DEBUGING: preparing chemistry parameters')

        select case (flag_model_type)
        case (STATION_MODE)
          call prepare_chemistry_stat()
        case (TRAJECTORY_MODE)
          call prepare_chemistry_traj()
        end select

        CALL debug_message('DEBUGING: sending and receiving chemistry data with MPI')
        ! print*, 'Otsoniiiiiiiiiiiii alku',  CH_CONS_ALL(40:45,ind_O3)

        !#ifdef PARALLEL
        call calculate_chemistry_parallel()
        !#else
        ! call calculate_chemistry_serial()
        !#endif
        ! print*,time_in_day, 'T', ta1(1)
        ! print*,time_in_day, 'P', pres(1)
        ! print*,time_in_day, 'OH', CH_CONS_ALL(1, ind_OH)
        ! print*,time_in_day, 'NO', CH_CONS_ALL(1, ind_NO)
        ! print*,time_in_day, 'emiNO', CH_CONS_ALL(1, ind_emi17)
        ! print*,time_in_day, 'emiNO', emis(1,ind_NO)
        ! print*,time_in_day, 'NO2', CH_CONS_ALL(1, ind_NO2)
        ! print*,time_in_day, 'emiNO2', emis(1,ind_NO2)
        ! print*,time_in_day, 'emi ISOP', emis(1,ind_C5H8)
        ! print*,time_in_day, 'xx swr', glob
        ! print*, 'Otsoniiiiiiiiiiiii loppu',  CH_CONS_ALL(40:45,ind_O3)


        ! CALL CPU_TIME(wtime2)
        ! wwtime = wwtime + (wtime2-wtime1)
        ! wwcount = wwcount + 1

        ! Then process OH-reactivities
        IF (CH_oh_count > 0) THEN  ! only do this if we have rOH-pseudochemicals in the chemicals list
          CH_oh_cons3(:,:,1) = (CH_CONS_ALL(:,CH_oh_indices) - CH_oh_prev3(:,:,1)) / dt_chem
          CH_oh_prev3(:,:,1) = CH_CONS_ALL(:,CH_oh_indices)
          IF (MOD( INT(time_in_day), INT(CH_step2*dt_chem) ) == 0) THEN
            CH_oh_cons3(:,:,2) = (CH_CONS_ALL(:,CH_oh_indices) - CH_oh_prev3(:,:,2)) / (CH_step2 * dt_chem)
            CH_oh_prev3(:,:,2) = CH_CONS_ALL(:,CH_oh_indices)
          ENDIF
          IF (MOD( INT(time_in_day), INT(CH_step3*dt_chem) ) == 0) THEN
            CH_oh_cons3(:,:,3) = (CH_CONS_ALL(:,CH_oh_indices) - CH_oh_prev3(:,:,3)) / (CH_step3 * dt_chem)
            CH_oh_prev3(:,:,3) = CH_CONS_ALL(:,CH_oh_indices)
            CH_CONS_ALL(:,CH_oh_indices) = 0  ! this may be not necessary, maybe the numbers would not grow too large
            CH_oh_prev3 = 0
          ENDIF
        ENDIF

        ! Then process O3-reactivities
        IF (CH_o3_count > 0) THEN  ! only do this if we have rO3-pseudochemicals in the chemicals list
          CH_o3_cons3(:,:,1) = (CH_CONS_ALL(:,CH_o3_indices) - CH_o3_prev3(:,:,1)) / dt_chem
          CH_o3_prev3(:,:,1) = CH_CONS_ALL(:,CH_o3_indices)
          IF (MOD( INT(time_in_day), INT(CH_step2*dt_chem) ) == 0) THEN
            CH_o3_cons3(:,:,2) = (CH_CONS_ALL(:,CH_o3_indices) - CH_o3_prev3(:,:,2)) / (CH_step2 * dt_chem)
            CH_o3_prev3(:,:,2) = CH_CONS_ALL(:,CH_o3_indices)
          ENDIF
          IF (MOD( INT(time_in_day), INT(CH_step3*dt_chem) ) == 0) THEN
            CH_o3_cons3(:,:,3) = (CH_CONS_ALL(:,CH_o3_indices) - CH_o3_prev3(:,:,3)) / (CH_step3 * dt_chem)
            CH_o3_prev3(:,:,3) = CH_CONS_ALL(:,CH_o3_indices)
            CH_CONS_ALL(:,CH_o3_indices) = 0  ! this may be not necessary, maybe the numbers would not grow too large
            CH_o3_prev3 = 0               !
          ENDIF
        ENDIF

        ! Then process NO3-reactivities
        IF (CH_no3_count > 0) THEN  ! only do this if we have rNO3-pseudochemicals in the chemicals list
          CH_no3_cons3(:,:,1) = (CH_CONS_ALL(:,CH_no3_indices) - CH_no3_prev3(:,:,1)) / dt_chem
          CH_no3_prev3(:,:,1) = CH_CONS_ALL(:,CH_no3_indices)
          IF (MOD( INT(time_in_day), INT(CH_step2*dt_chem) ) == 0) THEN
            CH_no3_cons3(:,:,2) = (CH_CONS_ALL(:,CH_no3_indices) - CH_no3_prev3(:,:,2)) / (CH_step2 * dt_chem)
            CH_no3_prev3(:,:,2) = CH_CONS_ALL(:,CH_no3_indices)
          ENDIF
          IF (MOD( INT(time_in_day), INT(CH_step3*dt_chem) ) == 0) THEN
            CH_no3_cons3(:,:,3) = (CH_CONS_ALL(:,CH_no3_indices) - CH_no3_prev3(:,:,3)) / (CH_step3 * dt_chem)
            CH_no3_prev3(:,:,3) = CH_CONS_ALL(:,CH_no3_indices)
            CH_CONS_ALL(:,CH_no3_indices) = 0  ! this may be not necessary, maybe the numbers would not grow too large
            CH_no3_prev3 = 0               !
          ENDIF
        ENDIF

        ! Then process reactivity
        IF (CH_reactivity_count > 0) THEN  ! only do this if we have reactivity-pseudochemicals in the chemicals list
          CH_reactivity_cons3(:,:,1) = (CH_CONS_ALL(:,CH_reactivity_indices) - CH_reactivity_prev3(:,:,1)) / dt_chem
          CH_reactivity_prev3(:,:,1) = CH_CONS_ALL(:,CH_reactivity_indices)
          IF (MOD( INT(time_in_day), INT(CH_step2*dt_chem) ) == 0) THEN
            CH_reactivity_cons3(:,:,2) = (CH_CONS_ALL(:,CH_reactivity_indices) - CH_reactivity_prev3(:,:,2)) / (CH_step2 * dt_chem)
            CH_reactivity_prev3(:,:,2) = CH_CONS_ALL(:,CH_reactivity_indices)
          ENDIF
          IF (MOD( INT(time_in_day), INT(CH_step3*dt_chem) ) == 0) THEN
            CH_reactivity_cons3(:,:,3) = (CH_CONS_ALL(:,CH_reactivity_indices) - CH_reactivity_prev3(:,:,3)) / (CH_step3 * dt_chem)
            CH_reactivity_prev3(:,:,3) = CH_CONS_ALL(:,CH_reactivity_indices)
            CH_CONS_ALL(:,CH_reactivity_indices) = 0  ! this may be not necessary, maybe the numbers would not grow too large
            CH_reactivity_prev3 = 0               !
          ENDIF
        ENDIF

      ENDIF  ! should this be after the time_par? End of flag_chem .EQ. 1.

      !***** Cumulative emissions *****!
      ! DO I=1, outemi_count
      !   Qemis(:, outemi_cheminds(I)) = Qemis(:, outemi_cheminds(I)) + EM_Emi(:, outemi_meganinds(I)) * dt_chem
      ! END DO
      ! Qemis(:,ind_C5H8)      =  Qemis(:,ind_C5H8)      + EM_Emi(:,1)  * dt_chem  ! C5H8 (Isoprene)
      ! Qemis(:,ind_Myrcene)   =  Qemis(:,ind_Myrcene)   + EM_Emi(:,2)  * dt_chem  ! Myrcene
      ! Qemis(:,ind_Sabinene)  =  Qemis(:,ind_Sabinene)  + EM_Emi(:,3)  * dt_chem  ! Sabinene
      ! Qemis(:,ind_LIMONENE)  =  Qemis(:,ind_LIMONENE)  + EM_Emi(:,4)  * dt_chem  ! LIMONENE
      ! Qemis(:,ind_Carene)    =  Qemis(:,ind_Carene)    + EM_Emi(:,5)  * dt_chem  ! Carene
      ! Qemis(:,ind_Ocimene)   =  Qemis(:,ind_Ocimene)   + EM_Emi(:,6)  * dt_chem  ! Ocimene
      ! Qemis(:,ind_BPINENE)   =  Qemis(:,ind_BPINENE)   + EM_Emi(:,7)  * dt_chem  ! Bpinene
      ! Qemis(:,ind_APINENE)   =  Qemis(:,ind_APINENE)   + EM_Emi(:,8)  * dt_chem  ! Apinene
      ! Qemis(:,ind_OMT)       =  Qemis(:,ind_OMT)       + EM_Emi(:,9)  * dt_chem  ! Other monoterpenes
      ! Qemis(:,ind_Farnesene) =  Qemis(:,ind_Farnesene) + EM_Emi(:,10) * dt_chem  ! Farnesene
      ! Qemis(:,ind_BCARY)     =  Qemis(:,ind_BCARY)     + EM_Emi(:,11) * dt_chem  ! BCARY (Beta-Carophyllene)
      ! Qemis(:,ind_OSQ)       =  Qemis(:,ind_OSQ)       + EM_Emi(:,12) * dt_chem  ! Other sesquiterpenes
      ! Qemis(:,ind_MBO)       =  Qemis(:,ind_MBO)       + EM_Emi(:,13) * dt_chem  ! MBO (2methyl-3buten-2ol)
      ! Qemis(:,ind_CH3OH)     =  Qemis(:,ind_CH3OH)     + EM_Emi(:,14) * dt_chem  ! CH3OH (Methanol)
      ! Qemis(:,ind_CH3COCH3)  =  Qemis(:,ind_CH3COCH3)  + EM_Emi(:,15) * dt_chem  ! CH3COCH3 (Aceton)
      ! Qemis(:,ind_CH4)       =  Qemis(:,ind_CH4)       + EM_Emi(:,16) * dt_chem  ! CH4 (Methane)
      ! Qemis(:,ind_NO)        =  Qemis(:,ind_NO)        + EM_Emi(:,17) * dt_chem  ! NO
      ! Qemis(:,ind_CH3CHO)    =  Qemis(:,ind_CH3CHO)    + EM_Emi(:,18) * dt_chem  ! CH3CHO (Acetaldehyde)
      ! Qemis(:,ind_HCHO)      =  Qemis(:,ind_HCHO)      + EM_Emi(:,19) * dt_chem  ! HCHO (Formaldehyde)
      ! Qemis(:,ind_CO)        =  Qemis(:,ind_CO)        + EM_Emi(:,20) * dt_chem  ! CO
      ! Qemis(:,ind_Cineole)   =  Qemis(:,ind_Cineole)   + EM_Emi(:,21) * dt_chem  ! Cineole (Eucalyptol) (not included in the new megan code - used is the same values as Sabinene)
      ! Qemis(:,ind_Linalool)  =  Qemis(:,ind_Linalool)  + EM_Emi(:,22) * dt_chem  ! Linalool

      !***** Cumulative chemistry sources and sinks *****!
      Qchem = Qchem + CH_CONS_ALL - CH_CONS_temp  ! remember that here the chemistry also includes the emissions
      ! substract Qemis is wrong since it is cumulative quantity, not the current one
      ! if (flag_emis==0)then
      !   CH_CONS_ALL(2,ind_APINENE) = linear_interp(time_in_month, nxodrad, CH_VOC_hyy(1,:), dt_obs) * air(2)*1.0d-9*0.51
      !   CH_CONS_ALL(2,ind_BPINENE) = linear_interp(time_in_month, nxodrad, CH_VOC_hyy(1,:), dt_obs) * air(2)*1.0d-9*0.12
      !   CH_CONS_ALL(2,ind_LIMONENE) = linear_interp(time_in_month, nxodrad, CH_VOC_hyy(1,:), dt_obs) * air(2)*1.0d-9*0.09
      !   CH_CONS_ALL(2,ind_CARENE) = linear_interp(time_in_month, nxodrad, CH_VOC_hyy(1,:), dt_obs) * air(2)*1.0d-9*0.28
      !   CH_CONS_ALL(2,ind_C5H8) = linear_interp(time_in_month, nxodrad, CH_VOC_hyy(2,:), dt_obs) * air(2)*1.0d-9
      !   ! CH_CONS_ALL(2,ind_MVK) = linear_interp(time_in_month, nxodrad, CH_VOC_hyy(3,:), dt_obs) * air(2)*1.0d-9
      !   !        CH_CONS_ALL(2,ind_MEK) = linear_interp(time_in_month, nxodrad, CH_VOC_hyy(4,:), dt_obs) * air(2)*1.0d-9
      !   CH_CONS_ALL(2,ind_CH3OH) = linear_interp(time_in_month, nxodrad, CH_VOC_hyy(5,:), dt_obs) * air(2)*1.0d-9
      !   CH_CONS_ALL(2,ind_CH3CHO) = linear_interp(time_in_month, nxodrad, CH_VOC_hyy(6,:), dt_obs) * air(2)*1.0d-9
      !   !        CH_CONS_ALL(2,ind_C2H5OH) = linear_interp(time_in_month, nxodrad, CH_VOC_hyy(7,:), dt_obs) * air(2)*1.0d-9
      !   CH_CONS_ALL(2,ind_CH3COCH3) = linear_interp(time_in_month, nxodrad, CH_VOC_hyy(8,:), dt_obs) * air(2)*1.0d-9
      !   CH_CONS_ALL(2,ind_CH3CO2H) = linear_interp(time_in_month, nxodrad, CH_VOC_hyy(9,:), dt_obs) * air(2)*1.0d-9
      !   CH_CONS_ALL(2,ind_H2SO4) = linear_interp(time_in_month, nxodrad, CH_H2SO4_hyy, dt_obs)
      ! endif
    ENDIF  ! IF (MOD( INT(time_in_day), INT(dt_chem) ) == 0)


    !==========================================================================!
    !
    ! GAS DRY DEPOSITION PART
    !
    !==========================================================================!

    if (MOD( INT(time_in_day), INT(dt_depo) ) == 0) then

      IF (flag_gasdrydep == 1) THEN
        CALL debug_message('DEBUGING: calculating gas dry deposition')
        LAIl_sl(2:nz-1) = LAIl(2:nz-1)*psn(2:nz-1)
        LAIl_sh(2:nz-1) = LAIl(2:nz-1) - LAIl_sh(2:nz-1)

        qnud_dep = qnud
        qnud_dep(1:11)=qnud_dep(12)
        Tc_dep = ta1-T00
        RH_dep = qnud_dep*8.314d0/0.0180153d0*ta1 / &
          (a0+a1*Tc_dep+a2*Tc_dep**2+a3*Tc_dep**3+a4*Tc_dep**4+a5*Tc_dep**5+a6*Tc_dep**6) / 100.0d0

        !===== Method from Ullar et al. (2012, ACP) =====!
        ! taub = EXP(-0.7*LAIl_c*0.83/cos_zenith)
        ! taud = EXP(-0.7*LAIl_c*0.83)
        ! PARl_sl(2:nz-1) = rsnt+taud(2:nz-1)*rskt
        ! PARl_sh(2:nz-1) = taud(2:nz-1)*rskt
        ! PARl(2:nz-1) = PARl_sl(2:nz-1)*taub(2:nz-1) + PARl_sh(2:nz-1)*(1.0d0-taub(2:nz-1))

        !===== Only consider the top radiation =====!
        ! PARl(2:nz-1) = rsnt*phsn + rskt*phsk
        ! PARl(2:nz-1) = rsnt*0.42+ rskt*0.60

        !===== Consider the PAR attenuation =====!
        !!!!! Consider the upward scattered PAR
        PARl(2:nz-1) = (rsnt*gl(2:nz-1)/cos_zenith*phsn &
          + rskt*psk(2:nz-1)*gd(2:nz-1)/cos_zenith*phsk &
          + fphd(2:nz-1) + fphu(2:nz-1))*psn(2:nz-1) &
          + (rskt*psk(2:nz-1)*gd(2:nz-1)/cos_zenith*phsk &
          + fphd(2:nz-1) &
          + fphu(2:nz-1))*(1.0d0-psn(2:nz-1))

        !!!!! Do not consider the upward scattered PAR
        ! PARl(2:nz-1) = (rsnt*gl(2:nz-1)/cos_zenith*phsn &
        !   + rskt*psk(2:nz-1)*gd(2:nz-1)/cos_zenith*phsk &
        !   + fphd(2:nz-1) )*psn(2:nz-1) &
        !   + (rskt*psk(2:nz-1)*gd(2:nz-1)/cos_zenith*phsk &
        !   + fphd(2:nz-1))*(1.0d0-psn(2:nz-1))
        ! PARl(nz-1) = PARl(nz-2)
        ! PARl(nz-2) = PARl(nz-2) - 30.0d0
        ! PARl(2:nz-1) = PARl(2:nz-1) * (1.0d0-tphu-rphu)

        !!!!! Consider the scattered radiation from below
        PARl_sl(2:nz-1) = rsnt*gl(2:nz-1)/cos_zenith*phsn &
          + rskt*psk(2:nz-1)*gd(2:nz-1)/cos_zenith*phsk &
          + fphd(2:nz-1) + fphu(2:nz-1)

        PARl_sh(2:nz-1) = rskt*psk(2:nz-1)*gd(2:nz-1)/cos_zenith*phsk &
          + fphd(2:nz-1) + fphu(2:nz-1)

        !!!!! Do not consider the scattered radiation from below
        ! PARl_sl(2:nz-1) = rsnt*gl(2:nz-1)/cos_zenith*phsn &
        !   + rskt*psk(2:nz-1)*gd(2:nz-1)/cos_zenith*phsk + fphu(2:nz-1)
        ! PARl_sh(2:nz-1) = rskt*psk(2:nz-1)*gd(2:nz-1)/cos_zenith*phsk &
        !   + fphu(2:nz-1)
        ! write(*,*) 'PARl:', PARl(1:nz)
        u_veg(2:nz-1) = sqrt(u1(2:nz-1)**2 + v1(2:nz-1)**2)
        rho_veg(2:nz-1) = roa
        Tc_veg(2:nz-1) = tsn1(2:nz-1)*psn(2:nz-1) + tsd1(2:nz-1)*(1.0d0-psn(2:nz-1)) ! Tc_veg is the averaged leaf temperature
        ! Tc_veg(2:nz-1) = ta1(2:nz-1)  ! Tc_veg = air T
        DO k=2,nz-1
          ! Henry(:,k) = HenryA * EXP(HenryB*(1.0d0/298.15d0 - 1.0d0/Tc_veg(k)))
          Henry(:,k) = HenryA(:)
        END DO

        CALL get_gas_dry_deposition_velocity(time_in_day, &
          nz-2, NSPEC, &  ! it starts from level 2 (level 1 is the ground, level nz is just over the canopy)
          l_drydep, l_vpd, l_wetskin, &
          frac_veg(2:nz-1), frac_ws(2:nz-1), &
          z(2:nz-1), hc, LAIl(2:nz-1), LAIl_sl(2:nz-1), LAIl_sh(2:nz-1), &
          PARl(2:nz-1), PARl_sl(2:nz-1), PARl_sh(2:nz-1), psn(2:nz-1), &
          u_veg(2:nz-1), ur(2:nz-1), rho_veg(2:nz-1), pres(2:nz-1), Tc_veg(2:nz-1), RH_dep(2:nz-1), &
          wg1(1), wgmax, wgwilt, &
          stomblock, &
          dvsn(2:nz-1), dvsd(2:nz-1), psn(2:nz-1)*dvsn(2:nz-1)+(1.0d0-psn(2:nz-1))*dvsd(2:nz-1), &
          gstm_h2o_sn(2:nz-1), gstm_h2o_sd(2:nz-1), psn(2:nz-1)*gstm_h2o_sn(2:nz-1)+(1.0d0-psn(2:nz-1))*gstm_h2o_sd(2:nz-1), &
          SPC_NAMES, molar_mass, SUM(Henry(:, 2:nz-1), DIM=2)/(nz-2.0d0), dryreac, ind_O3, ind_SO2, &  ! average of Henry at all levels
          vdep(2:nz-1, :), dep_output(2:nz-1,:,:))
        DO k=2, nz-1
          gasdep_flux(k,:) = -vdep(k,:)*CH_CONS_ALL(k,:)/( 0.5d0*(dz(k)+dz(k+1)) )
        END DO

        ! Update concentrations
        CH_CONS_ALL = CH_CONS_ALL + dt_mete*gasdep_flux

        !***** Cumulative deposition flux *****!
        Qdepo = Qdepo + dt_mete*gasdep_flux
      END IF  ! IF (flag_gasdrydep == 1) THEN
    end if  ! if (MOD( INT(time_in_day), INT(dt_depo) ) == 0) then


    !==========================================================================!
    !
    ! AEROSOL PART
    !
    !==========================================================================!

    CALL debug_message('DEBUGING: calculating aerosol processes')


    ! Calculation of nucleation rates by limononic acid and sulfuric acid

    IF (flag_aero == 1 .AND. time_in_month>=aero_start) THEN
      IF (MOD( INT(time_in_day), INT(dt_aero) ) == 0) THEN

        ! In the beginning of everyday or the first time loop, open a new dmps file
        ! (since they are saved as one file per day) to read the data.
        IF ((flag_model_type/=TRAJECTORY_MODE).and.(first_time_loop .or. is_newday)) THEN
          ! dmps_file_name = TRIM(ADJUSTL(input_dir_station)) &
          !   //'/dmps/'// year_str//'/dm'//year_str2//''//month_str//''//day_str//'.sum'
          dmps_file_name = '/xxx/xxx/xxx/xxxx/dm180504.sum'

          ! write(*,*) '[Putian Debug] Reading new dmps file.', dmps_file_name
          line = 1
          ! G_options%month = mon

          ! try to find the number of sizebins (ccount) in the measurement
          open(12345, file = TRIM(ADJUSTL(dmps_file_name)))

          ! read the first line to dummy
          read(12345, '(A)') dummy_line
          close(12345)

          ! read values from dummy
          read(dummy_line, *, END = 8888) (val(ii), ii=1, maxcol)
          8888    continue

          ! Calculate number of sizebin from measurement ccount
          ii=3 ! since the first two numbers are zero.
          ccount = 0
          do while (val(ii)>0)
            ccount = ccount +1
            ii = ii +1
          end do
        END IF

        ! dmps data file time interval is 10 minutes
        IF (time_in_month/600.0 == FLOOR(time_in_month/(600.0))) THEN
          line=line+1
        ENDIF


        ! For the first 10 minutes every day, the particle concentration inside
        ! the mixing layer is set to measurement data, while the concentration
        ! above the mixing layer height is set as a portion of the measurement
        ! values.
        ! IF ( (now_date(4) < 1) .AND. (time_in_month/600.0 == FLOOR(time_in_month/600.0)) ) then
        IF ( time_in_month == aero_start) then ! Only initialize once P.C.
            ! print*, time_in_month, 'reading in initial PSD, line, # columns', line, ccount
          ! call use_dmps(TRIM(ADJUSTL(dmps_file_name)), &
          !   ccount, line,G_particles,G_ambient%density,initial_dist%mfractions)

          do k = 1, pblh_ind
            N_CONC(k,:) = current_psd%conc_fs
            Radius(k,:) = current_psd%diameter_fs * 0.5d0
            Rdry(k,:)   = current_psd%dp_dry_fs * 0.5d0
            Mass(k,:)   = current_psd%particle_mass_fs
            MASS_COMPO(k,:,:) = current_psd%composition_fs
          end do

          do k = pblh_ind+1, kz
              N_CONC(k,:) = current_psd%conc_fs * 0.05d0
              Radius(k,:) = current_psd%diameter_fs * 0.5d0
              Rdry(k,:)   = current_psd%dp_dry_fs * 0.5d0
              Mass(k,:)   = current_psd%particle_mass_fs
              MASS_COMPO(k,:,:) = current_psd%composition_fs
          end do
        end if

        call calculate_aerosol_parallel()

      endif  ! IF (MOD( INT(time_in_day), INT(dt_aero) ) == 0) THEN

      ! Calculate aerosol deposition velocity
      do k = 1, kz
        call DRY_DEPLAYER( k,ta1(k),pres(k),hwind(k),ur(k),s1(k), &
          face,Rk_vel(k,:),RADIUS(k,:))
      end do
      ! print*, 'avg Rk vel @1', sum(Rk_vel(1,:))/n_bins_par
      IF (ANY(Rk_vel<0d0)) print*, 'WTF XXXXXXXXXXXXx', Rk_vel
    end if  ! IF (flag_aero == 1 .AND. time_in_month>=aero_start) THEN

    !==========================================================================!
    ! Mixing of scalars
    !==========================================================================!

    !------------------------------------------------------!
    ! Mixing of chemical species
    !------------------------------------------------------!

    ! The chemical species will be mixed when either chemistry or emission is
    ! switched on
    if (flag_chem == 1 .or. flag_emis > 0) then

      CALL debug_message('DEBUGING: mixing chemical compounds')

      CH_CONS_temp = CH_CONS_ALL

      !===== Get the input O3 concentration from input file when the chemistry is off =====!
      ! IF (.FALSE.) THEN
      !   CALL debug_message('DEBUGGING: mixing O3')
      !   o3_new = linear_interp(time_in_month, nxodrad, CH_gas_hyy(2,:), dt_obs) * air(nz)*1.0d-9
      !
      !   DO k=2,kz-1
      !     terma(k) = (1.0d0-o3_weight(k)) * dt_mete*(2.0d0*Kht(k)-w(k)*dz(k))/da(k)
      !     termc(k) = (1.0d0-o3_weight(k)) * dt_mete*(2.0d0*Khb(k)+w(k)*dz(k+1))/dc(k)
      !     termb(k) = (1.0d0-o3_weight(k)) * dt_mete*(2.0d0*(dz(k)*Kht(k)+ dz(k+1)*Khb(k))/(dz(k+1)+dz(k)) + &
      !       w(k)*(dz(k+1)-dz(k)))/db(k) + 1.0d0
      !   ENDDO
      !
      !   DO k = 1,kz
      !     ! CH_CONS_VER(k) = CH_CONS_ALL(k,ind_O3) + (1.0d0-o3_weight(k))*dt_mete*gasdep_flux(k,ind_O3) &
      !     !                  + o3_weight(k)*(o3_new-CH_CONS_ALL(k,ind_O3))
      !     CH_CONS_VER(k) = CH_CONS_ALL(k,ind_O3) + o3_weight(k)*(o3_new-CH_CONS_ALL(k,ind_O3))
      !   ENDDO
      !
      !   !!!!! calculate once for all layers !!!!!
      !   CALL gtri(terma, termc, termb, CH_CONS_VER, CH_CONS_VER_NEW, 2, 0.d0, 3.d0, 2, 0.d0, 3.d0, kz, 2)
      !
      !   DO k=2,kz-1
      !     !!!!! kt: [m2 s-1], 100*kt: [cm m s-1], the flux unit is then [molec cm-2 s-1], upward is positive
      !     CH_CONS_FLUX(k,ind_O3) = -0.5*100.0d0*( Kht(k)*(CH_CONS_VER_NEW(k+1) - CH_CONS_VER_NEW(k))/dz(k+1) +  &
      !       Khb(k)*(CH_CONS_VER_NEW(k) - CH_CONS_VER_NEW(k-1))/dz(k) )
      !   ENDDO
      !
      !   DO k=1,kz
      !     CH_CONS_ALL(k,ind_O3) = MAX(0.0d0, CH_CONS_VER_NEW(k))
      !   ENDDO
      ! END IF  ! mixing of ozone

      !==========================================================================!
      !
      ! Mixing of chemicals
      !
      !==========================================================================!


      IF (flag_mix_chem == 1 .or.flag_aero == 1 ) THEN
        !----------------------------------------------------------------
        !  solution for solving the vertical transport of the gas species
        !----------------------------------------------------------------
        DO k=2,kz-1
          terma(k) = dt_mete*(2.0d0*Kht(k)-w(k)*dz(k))/da(k)
          termc(k) = dt_mete*(2.0d0*Khb(k)+w(k)*dz(k+1))/dc(k)
          termb(k) = dt_mete*(2.0d0*(dz(k)*Kht(k) + dz(k+1)*Khb(k)) /(dz(k+1)+dz(k)) +w(k)*(dz(k+1)-dz(k)))/db(k)+ 1.0d0
        ENDDO
      ENDIF  ! flag_mix_chem == 1 .or.flag_aero == 1

        IF (flag_mix_chem == 1) THEN
        DO J = 1, NSPEC
          ! IF (J==ind_O3) CYCLE  ! do not mix O3 since it is alreaday calculated above
          IF (CH_oh_flag(J)) CYCLE  ! do not mix rOH-chemicals, go back to start next do loop cycle right away
          IF (CH_o3_flag(J)) CYCLE  ! do not mix rO3-chemicals, go back to start next do loop cycle right away
          IF (CH_no3_flag(J)) CYCLE  ! do not mix rNO3-chemicals, go back to start next do loop cycle right away
          IF (CH_reactivity_flag(J)) CYCLE  ! do not mix reactivity-chemicals, go back to start next do loop cycle right away

          CH_CONS_VER = CH_CONS_ALL(:,J)  ! + dt_mete*gasdep_flux(:,J)

          CALL gtri(terma, termc, termb, CH_CONS_VER, CH_CONS_VER_NEW, 2, 0.0d0, 3.0d0, 2, 0.0d0, 3.0d0, kz, 2)

          DO k=2,kz-1
            !!!!! kt: [m2 s-1], 100*kt: [cm m s-1], so the unit of flux is [molec cm-2 s-1], upward is positive
            CH_CONS_FLUX(k,j) = -0.5d0*100.0d0*( Kht(k)*(CH_CONS_VER_NEW(k+1) - CH_CONS_VER_NEW(k))/dz(k+1) +  &
              Khb(k)*(CH_CONS_VER_NEW(k) - CH_CONS_VER_NEW(k-1))/dz(k) )
          ENDDO

          DO k=1,kz
            CH_CONS_ALL(k,J) = MAX(CH_CONS_VER_NEW(k), 0.0d0)
          ENDDO
        ENDDO
      ENDIF  ! flag_mix_chem == 1

      Qturb = Qturb + CH_CONS_ALL - CH_CONS_temp
      Qturb_now = CH_CONS_ALL - CH_CONS_temp  ! concentration change caused by turbulent transport at current time point.
    end if  ! if (flag_chem == 1 .or. flag_emis > 0)


    !***** Calculate the concentration change *****!
    ! CH_CONS_1 = CH_CONS_ALL
    ! Qconc = CH_CONS_1 - CH_CONS_0
    ! Qfluxdep = CH_CONS_1 - CH_CONS_temp
    ! write(*,*) 'Qconc(APINENE)', Qconc(1:20, ind_APINENE)

    !------------------------------------------------------!
    ! Mixing of aerosol particles
    !------------------------------------------------------!

    if (flag_aero == 1 .AND. time_in_month>=aero_start) then
      if (flag_mix_aero == 1) then
          ! print*, 'total column PM2 before mixing: ', sum(SUM(MASS_COMPO,3) * N_CONC)
          do k=1,kz
              do i=1,n_bins_par
                VOL_COMPO(k,i,:) = MASS_COMPO(k,i,:)/VAPOUR_PROP%density * N_CONC(k,i)
              end do
          end do
          ! print*, 'total column PV2 before mixing: ', sum(VOL_COMPO)
        CALL SCPARFLUX(kz,z,dz,dt_mete,alt,w1,kt,N_CONC,PAR_FLUX,VOL_COMPO,Rk_vel,ta1,pres, CLUSTERS)
        ! print*, 'total column PV2 after mixing: ', sum(VOL_COMPO)
        do k=1,kz
            do i=1,n_bins_par
                IF (N_CONC(k,i) > 0d0) THEN
                    MASS_COMPO(k,i,:) = VOL_COMPO(k,i,:)*VAPOUR_PROP%density / N_CONC(k,i)
                ELSE
                    MASS_COMPO(k,i,:) = 0d0
                end if
            end do



        end do
        ! print*,   'total column PM2 after mixing : ', sum(SUM(MASS_COMPO,3) * N_CONC)
      end if
    end if
    ! do k=1,kz
    !     IF (MOD(INT(time_in_day), 600)==0) WRITE(1199, *) k, N_CONC(k, :)
    !     IF (MOD(INT(time_in_day), 600)==0) FLUSH(1199)
    ! end do

    !------------------------------------------------------!
    ! Update time
    !------------------------------------------------------!

    CALL debug_message('DEBUGING: updating time')

    ! Update time_in_day
    time_in_day=time_in_day+dt_mete


    ! Update month and day
    IF (time_in_day >= SECONDS_IN_ONE_DAY) THEN  ! move to next day
      is_newday = .true.
      julian = julian + 1
      now_date(3) = now_date(3) + 1
      IF ( now_date(3) > MONTH_DAYS(now_date(2)) ) THEN  ! move to next month
        is_newmonth = .true.
        now_date(2) = now_date(2) + 1
        now_date(3) = 1
      ELSE
        is_newmonth = .false.
      END IF
      time_in_day = MOD(time_in_day, SECONDS_IN_ONE_DAY)
    ELSE
      is_newday = .false.
    END IF

    ! Update input time pointer
    if (flag_model_type==TRAJECTORY_MODE) then
      if (time >= infield%var(infield%id_time)%f1d(input_ptr+1)) then
        input_ptr = input_ptr + 1
      end if
    end if

    ! Update hour, minute and second
    now_date(4) = INT(time_in_day/3600.0d0)
    now_date(5) = INT((time_in_day-now_date(4)*3600.d0)/60.0d0)
    now_date(6) = INT(time_in_day-now_date(4)*3600.0d0-now_date(5)*60.0d0)

    ! Update model time and time_in_month
    ! Here time_in_month needs to be modified if it is withtin next month, this
    ! needs to be done in future.
    time = time + dt_mete
    time_in_month = time_in_month + dt_mete

    ! Update UTC time
    time_in_day_UTC = time_in_day - time_zone*SECONDS_IN_ONE_HOUR
    if (time_in_day_UTC < 0) then  ! UTC time is in previous day
      time_in_day_UTC = time_in_day_UTC + c864
      julian_UTC = julian - 1
    else if (time_in_day_UTC > c864) then  ! UTC time is in next day
      time_in_day_UTC = time_in_day_UTC - c864
      julian_UTC = julian + 1
    else
      julian_UTC = julian
    end if

    ! Write the output data
    IF (MOD(INT(time_in_day),1800) == 0) THEN
      CALL debug_message('DEBUGING: writing output data')
      call output_step()
      CALL output_write()
      ! if (flag_aero==1) call output_write_old()

      !***** Reset the cumulative quantities *****!
      Qemis = 0.0d0
      Qchem = 0.0d0
      Qdepo = 0.0d0
      Qturb = 0.0d0
    END IF

    ! Update time string
    ! Now the year time is not updated, need to be added in future
    WRITE(year_str,'(I4)') now_date(1)  ! '1999', '2000', '2001', ...
    WRITE(year_str2, '(I2.2)') MOD(now_date(1), 100)  ! '99', '00', '01', ...
    WRITE(month_str_Mxx, '(A1, I2.2)') 'M', now_date(2)  ! 'M01', 'M02', ..., 'M12'
    WRITE(month_str, '(I2.2)') now_date(2)  ! '01', '02', ..., '12'
    WRITE(day_str, '(I2.2)') now_date(3)  ! '01', '02', ..., '31'


    !------------------------------------------------------!
    ! Update latitude and longitude, and local variables
    !------------------------------------------------------!

    ! Update latitude and longitude in traj mode
    call update_location(lat_deg, lon_deg, flag_model_type)

    ! Solar zenith angle
    cos_zenith = get_cos_zenith( time_in_day_UTC, real(julian_UTC, dp), &
      lat_deg, lon_deg)
    zenith_deg = acos(cos_zenith)*rad2deg

    ! Update sunrise and sunset time
    ! New method to calculate the sunrise and sunset time, but it is not tested
    ! yet, so use the old one upsn and downsn now in update_meteorology_stat.
    ! call get_sunrise_and_sunset_time_UTC( &
    !   real(julian_UTC, dp), lat_deg, lon_deg, &
    !   sunrise_UTC, sunset_UTC)
    ! write(*,*) 'rads: ', rads
    ! write(*,*) 'sunrise, sunset old: ', upsn, downsn
    ! write(*,*) 'sunrise, sunset new: ', &
    !   sunrise_UTC+time_zone*SECONDS_IN_ONE_HOUR, sunset_UTC+time_zone*SECONDS_IN_ONE_HOUR


    !------------------------------------------------------!
    ! Other variables
    !------------------------------------------------------!

    first_time_loop = .false.

  END DO  ! DO WHILE (time_end > time_in_month)



  !  SSSSSSSSSSSSSSS        OOOOOOOOO      SSSSSSSSSSSSSSS AAAAAAA                   AAAAAAAAAAAAAA                   AAAAAAA
  ! S:::::::::::::::SS    OO:::::::::OO   S:::::::::::::::SSA:::::A                 A:::::A  A:::::A                 A:::::A
! S::::::SSSSSS:::::S OO:::::::::::::OO S::::::SSSSSS:::::SA:::::A               A:::::A    A:::::A               A:::::A
! SSSSSSS     S:::::SO:::::::OOO:::::::OSSSSSSS     S:::::S A:::::A             A:::::A      A:::::A             A:::::A
!             S:::::SO::::::O   O::::::O            S:::::S  A:::::AAAAAAAAAAAAA:::::A        A:::::AAAAAAAAAAAAA:::::A
!             S:::::SO:::::O     O:::::O            S:::::S   A:::::::::::::::::::::A          A:::::::::::::::::::::A
!        SSSSSS::::S O:::::O     O:::::O       SSSSSS::::S     A:::::AAAAAAAAA:::::A            A:::::AAAAAAAAA:::::A
!     SSS::::::::SS  O:::::O     O:::::O    SSS::::::::SS       A:::::A     A:::::A              A:::::A     A:::::A
!   SS::::::SSSSS    O:::::O     O:::::O  SS::::::SSSSS          A:::::A   A:::::A                A:::::A   A:::::A
!  S::::SSSS         O:::::O     O:::::O S::::SSSS                A:::::A A:::::A                  A:::::A A:::::A
! S:::::S            O:::::O     O:::::OS:::::S                    A:::::A:::::A                    A:::::A:::::A
! S:::::S            O::::::O   O::::::OS:::::S                     A:::::::::A                      A:::::::::A
! S:::::S     SSSSSSSO:::::::OOO:::::::OS:::::S     SSSSSSS          A:::::::A                        A:::::::A
! S:::::SSSSSS::::::S OO:::::::::::::OO S:::::SSSSSS::::::S           A:::::A                          A:::::A
!  SS:::::::::::::::S   OO:::::::::OO    SS:::::::::::::::S            A:::A                            A:::A
!    SSSSSSSSSSSSSSS      OOOOOOOOO        SSSSSSSSSSSSSSS              AAA                              AAA









  ! if (flag_aero==1) call output_netcdf_done
  call output_done()

  write(*,*) 'balans = ',balans

  !#ifdef PARALLEL
  IF (my_id == master_id) THEN
    ! send message to slaves to stop, so the whole program ends cleanly
    DO slave_i = 1,mpi_nslaves
      CALL MPI_SEND( mpi_end_task_code, 1, MPI_INTEGER, &
        slave_i, mpi_task_code_tag, &
        MPI_COMM_WORLD, mpi_rc&
        )
    END DO
  END IF

  CALL MPI_Finalize(mpi_rc)
  !#endif
  call cpu_time(cpu2)
  print*, 'RUNTIME:', cpu2 - cpu1


CONTAINS

!    ▄▄▄▄▄▄▄▄▄▄▄  ▄            ▄▄▄▄▄▄▄▄▄▄▄  ▄               ▄  ▄▄▄▄▄▄▄▄▄▄▄
!   ▐░░░░░░░░░░░▌▐░▌          ▐░░░░░░░░░░░▌▐░▌             ▐░▌▐░░░░░░░░░░░▌
!   ▐░█▀▀▀▀▀▀▀▀▀ ▐░▌          ▐░█▀▀▀▀▀▀▀█░▌ ▐░▌           ▐░▌ ▐░█▀▀▀▀▀▀▀▀▀
!   ▐░▌          ▐░▌          ▐░▌       ▐░▌  ▐░▌         ▐░▌  ▐░▌
!   ▐░█▄▄▄▄▄▄▄▄▄ ▐░▌          ▐░█▄▄▄▄▄▄▄█░▌   ▐░▌       ▐░▌   ▐░█▄▄▄▄▄▄▄▄▄
!   ▐░░░░░░░░░░░▌▐░▌          ▐░░░░░░░░░░░▌    ▐░▌     ▐░▌    ▐░░░░░░░░░░░▌
!    ▀▀▀▀▀▀▀▀▀█░▌▐░▌          ▐░█▀▀▀▀▀▀▀█░▌     ▐░▌   ▐░▌     ▐░█▀▀▀▀▀▀▀▀▀
!             ▐░▌▐░▌          ▐░▌       ▐░▌      ▐░▌ ▐░▌      ▐░▌
!    ▄▄▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄▄▄ ▐░▌       ▐░▌       ▐░▐░▌       ▐░█▄▄▄▄▄▄▄▄▄
!   ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░▌       ▐░▌        ▐░▌        ▐░░░░░░░░░░░▌
!    ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀         ▀          ▀          ▀▀▀▀▀▀▀▀▀▀▀

  !#ifdef PARALLEL
  SUBROUTINE slave_task
    ! This is where all other branches, except the master branch, do their work
    !
    ! And their work is: waiting for packets containing input for CHEMISTRY
    ! calling CHEMISTRY, and sending results back
    IMPLICIT NONE
    INTEGER :: task_code,iis, ii0
    REAL(kind=dp) :: time1,time2, CH_RO2, day_of_year, dt_aero
    ! RO2 here will hide the visibility of the global variable RO2(kz)
    REAL(kind=dp) :: tmcontr_local
    REAL(kind=dp) :: ta1_local, pres1_local, zenith_deg_local, sun_par_local
    real(dp) :: cs_H2SO4_local, cs_HNO3_local,U10m, lsm
    real(dp) :: air_local, N2_local, O2_local, H2O_local
    ! variables local to this subroutine, to hold the received values from the corresponding variables from the main program
    INTEGER :: year_int, month_int
    REAL(kind=dp) :: qa1

    REAL(kind=dp), DIMENSION(NPHOT)    :: J_values_local
    REAL(kind=dp), DIMENSION(NKVALUES) :: K_values_local
    REAL(kind=dp), DIMENSION(75)       :: ACF_local
    REAL(kind=dp), DIMENSION(9)        :: aero_emissions
    REAL(kind=dp)                      :: seasalt_emissions(n_bins_par)

    CALL KPP_SetUp()
    ! OPEN(UNIT=1999, FILE=TRIM(OUTPUT_DIR)//'/emi_aero.txt', status='replace',action='write')

    DO

      CALL MPI_RECV(task_code,1,MPI_INTEGER,master_id,mpi_task_code_tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_rc)

      IF (task_code == mpi_do_chemistry_code) THEN

        CALL MPI_RECV( &
          mpi_send_buffer, mpi_send_buffer_size, MPI_DOUBLE_PRECISION, &
          master_id, mpi_buffer_tag, &
          MPI_COMM_WORLD, MPI_STATUS_IGNORE,mpi_rc &
          )

        kz_tag_dp        = mpi_send_buffer(index_kz_tag)
        time1            = mpi_send_buffer(index_time1)
        time2            = mpi_send_buffer(index_time2)
        tmcontr_local    = mpi_send_buffer(index_tmcontr)
        nxodrad          = INT(mpi_send_buffer(index_nxodrad))
        year_int         = INT(mpi_send_buffer(index_year))
        month_int        = INT(mpi_send_buffer(index_month))
        ta1_local        = mpi_send_buffer(index_ta1)
        pres1_local      = mpi_send_buffer(index_pr1)
        zenith_deg_local = mpi_send_buffer(index_zenith_deg)
        sun_par_local    = 1d0 ! mpi_send_buffer(index_sun_par)
        cs_H2SO4_local   = mpi_send_buffer(index_res1)
        cs_HNO3_local    = mpi_send_buffer(index_res2)
        qa1              = mpi_send_buffer(index_qa1)
        air_local        = mpi_send_buffer(index_air)
        N2_local         = mpi_send_buffer(index_N2)
        O2_local         = mpi_send_buffer(index_O2)
        H2O_local        = mpi_send_buffer(index_H2O)
        CH_CONS          = mpi_send_buffer(index1_CONS:index2_CONS)
        albedo           = mpi_send_buffer(index_albedo)
        glob             = mpi_send_buffer(index_glob)
        ! CH_RES_org       = mpi_send_buffer(index1_res_org:index2_res_org)


        ! In CHEMISTRY, any value above -1 is replaced by MCM calculated photorates. This
        ! behaviour is to enable sending in user defined rates (will be implemented later)
        J_values_local = -999d0

        CALL CHEMISTRY(CH_CONS, time1, time2, &
          ta1_local, pres1_local, zenith_deg_local, sun_par_local*glob, albedo, &
          cs_H2SO4_local, cs_HNO3_local, &
          air_local, N2_local, O2_local, H2O_local, &
          CH_RO2, J_values_local, K_values_local, ACF_local, int(kz_tag_dp) )

        mpi_recv_buffer(index_kz_tag)                              = kz_tag_dp
        mpi_recv_buffer(index_RO2)                                 = CH_RO2
        mpi_recv_buffer(index_H2O_recv)                            = H2O_local
        mpi_recv_buffer(index1_CONS_recv:index2_CONS_recv)         = CH_CONS
        mpi_recv_buffer(index1_J_values_recv:index2_J_values_recv) = J_values_local
        mpi_recv_buffer(index1_K_values_recv:index2_K_values_recv) = K_values_local
        mpi_recv_buffer(index1_ACF_recv:index2_ACF_recv)           = ACF_local

        CALL MPI_SEND( &
          mpi_recv_buffer, mpi_recv_buffer_size, MPI_DOUBLE_PRECISION, &
          master_id, mpi_buffer_tag, &
          MPI_COMM_WORLD, mpi_rc &
          )

      ELSE IF (task_code == mpi_do_aerosol_code) THEN

        CALL MPI_RECV(new_mpi_send_aer_buf,new_mpi_send_aer_buf_size,MPI_DOUBLE_PRECISION,master_id,mpi_buffer_tag, &
          MPI_COMM_WORLD, MPI_STATUS_IGNORE,mpi_rc)

        kz_tag_dp                  = new_mpi_send_aer_buf(ind_kz_tag)
        time1                      = new_mpi_send_aer_buf(index_time1)
        time2                      = new_mpi_send_aer_buf(index_time2)
        dt_aero                    = new_mpi_send_aer_buf(ind_aer_ts)
        U10m                       = new_mpi_send_aer_buf(ind_U10m)
        lsm                        = new_mpi_send_aer_buf(ind_lsm)
        grh                        = new_mpi_send_aer_buf(ind_rh)
        gtempk                     = new_mpi_send_aer_buf(ind_tempk)
        gpres                      = new_mpi_send_aer_buf(ind_pres)
        aero_emissions             = new_mpi_send_aer_buf(ind1_emiaer:ind2_emiaer)
        seasalt_emissions          = new_mpi_send_aer_buf(ind1_seasalt:ind2_seasalt)
        CURRENT_PSD%conc_fs        = new_mpi_send_aer_buf(ind1_n_conc:ind2_n_conc)
        CH_CONS                    = new_mpi_send_aer_buf(ind1_CONS_aer:ind2_CONS_aer)
        CURRENT_PSD%composition_fs = RESHAPE(new_mpi_send_aer_buf(ind1_composition:ind2_composition), [n_bins_par,n_cond_tot])

        DO iis=1, size(G_ACDC)
            if (iis==1) THEN
                ii0 = ind1_CONS_clust
            ELSE
                ii0 = ii0 + G_ACDC(iis-1)%neq_syst
            END IF
            G_ACDC(iis)%ACDC_clusConc = new_mpi_send_aer_buf(ii0:ii0+G_ACDC(iis)%neq_syst-1)
        END DO

        ! if (int(kz_tag_dp)<=4) THEN
        !     ! print'(2(a,es8.1))', 'lsm',lsm, 'U',U10M
        !     call calculate_sea_salt_emissions(seasalt_emissions, U10m)
        !     ! print*, 'PN tot',sum(0.1d0*seasalt_emissions*(1-lsm))*bin_ratio_lg
        ! ELSE
        !     seasalt_emissions = 0d0
        ! END IF

        ! do iis=1,n_bins_par
        !     print* ,lsm, current_PSD%diameter_fs(iis)*1d6, seasalt_emissions(iis)*(1-lsm)
        ! end do
        ! stop

        call run_arca_psd(CH_CONS, dt_aero, gtempk, gpres, aero_emissions,seasalt_emissions, new_mpi_recv_aer_buf(ind1_oth_nuc_rates))

        DO iis=1, size(G_ACDC)
            if (iis==1) THEN
                ii0 = ind1_CONS_clust
            ELSE
                ii0 = ii0 + G_ACDC(iis-1)%neq_syst
            END IF
            new_mpi_recv_aer_buf(ii0:ii0+G_ACDC(iis)%neq_syst-1) = G_ACDC(iis)%ACDC_clusConc
        END DO

        ! IF (MOD(INT(time1), 3600)==0) write(1999, *) int(kz_tag_dp), sum(dconc_dep_mix/dt_aero), dconc_dep_mix/dt_aero/bin_ratio_lg
        ! IF (MOD(INT(time1), 3600)==0) FLUSH(1999)
        !
        new_mpi_recv_aer_buf(ind_kz_tag)                        = kz_tag_dp
        new_mpi_recv_aer_buf(ind1_n_conc:ind2_n_conc)           = CURRENT_PSD%conc_fs
        new_mpi_recv_aer_buf(ind1_composition:ind2_composition) = RESHAPE(CURRENT_PSD%composition_fs, [n_bins_par*n_cond_tot])
        new_mpi_recv_aer_buf(ind1_CONS_aer:ind2_CONS_aer)       = CH_CONS
        new_mpi_recv_aer_buf(ind1_sink:ind2_sink)               = VAPOUR_PROP%sink
        new_mpi_recv_aer_buf(ind1_acdc_rates:ind2_acdc_rates)   = [(G_ACDC(i)%J_OUT_CM3, i=1,size(G_ACDC))] !G_ACDC(1)%J_OUT_M3

        CALL MPI_SEND(new_mpi_recv_aer_buf,new_mpi_recv_aer_buf_size,&
          MPI_DOUBLE_PRECISION,master_id,mpi_buffer_tag, &
          MPI_COMM_WORLD,mpi_rc)

      ELSE  ! task_code must be mpi_end_task_code
          EXIT

      END IF

    END DO
  END SUBROUTINE slave_task

  !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !  Calculates the frequency of collision (1/cm3/s) from h2so4 and limononic acid based on Seinfeld page 152
  !
  !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  FUNCTION COLLISION(Temp)
    IMPLICIT NONE
    REAL(kind=dp) :: collision
    REAL(kind=dp) :: Temp, m_reduced, sigma, mass1, mass2, r_gas1, r_gas2
    REAL(kind=dp), PARAMETER :: k          = 5.67051e-8,               &  ! Boltzman constant (W/m2/K4)
      k2         = 1.3807d-23,               &  ! Stefan Boltzman constant (J/K)
      Avog       = 6.0221E23,                &  ! Avogadro-number (molecules / mol)
      molarmass1 = 98.08,                    &  ! Molar mass of H2SO4
      molarmass2 = 184.11,                   &  ! Molar mass of limononic acid
      molecvol1  = 8.804e-29,                &  ! Molecular volume of H2SO4 (m3)
      molecvol2  = 3.534e-28                    ! Molecular volume of limononic acid (m3)
    r_gas1 = (molecvol1*3./(4.*PI))**(1./3.)       ! Radius of one H2SO4 molecule (m)
    r_gas2 = (molecvol2*3./(4.*PI))**(1./3.)       ! Radius of one limononic acid molecules (m)
    ! Collision cross-section, sigma = in simple hard sphere collision theory the area of the circle with radius
    ! equal to the collision diameter. G.B. 56; see also 1996, 68, 158
    sigma = r_gas1**2 + r_gas2**2 * 0.4            ! collision cross section
    mass1 = molarmass1 / Avog                      ! mass of one sulfuric acid molecule in g
    mass2 = molarmass2 / Avog                      ! mass of one limononic acid molecule in g
    m_reduced = mass1 * mass2 / (mass1 + mass2)    ! reduced mass
    COLLISION =  (8 * k2 * 1E7 * Temp / PI / m_reduced )**0.5 * sigma*1E4 * PI
  END FUNCTION COLLISION


  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   Calculates the collision rate (m^3/s) of vapour molecules
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  FUNCTION COLLISION_RATE(TS7, PS7)
    IMPLICIT NONE

    !INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
    REAL(kind=dp) :: collision_rate
    ! INTEGER :: ii
    REAL(kind=dp) :: TS7, air_free_path, dif_vap1, dif_vap2, r_vap1, r_vap2, viscosity
    REAL(kind=dp) :: PS7
    !REAL(real_x), DIMENSION(sections) :: knudsen, corr, mean_path

    REAL(kind=dp), PARAMETER :: &
      molecvol1  = 8.804e-29,   &   ! Molecular volume of H2SO4
      molecvol2  = 3.534e-28,   &   ! Molecular volume of limononic acid
      molarmass1 = 98.08,       &   ! Molar mass of H2SO4
      molarmass2 = 184.11,      &   ! Molar mass of limononic acid
      diff_vol1  = 51.96,       &   ! diffusion volume of H2SO4
      diff_vol2  = 60.00            ! diffusion volume of limononic acid - assumption

    ! Air mean free path (m) and viscosity (kg s^-1 m^-1)
    air_free_path = (6.73d-8 * TS7 * (1. + 110.4 / TS7)) / (296. * PS7 * 100 / Patm * 1.373)

    viscosity     = (1.832d-5*406.4*TS7**1.5)/(5093*(TS7+110.4))

    r_vap1 = (molecvol1*3./(4.*PI))**(1./3.)  ! radius of condensable molecules (m)
    r_vap2 = (molecvol2*3./(4.*PI))**(1./3.)  ! radius of condensable molecules (m)

    ! molecular diffusion coefficient (m^2/s)
    dif_vap1 = DIFFUSION(TS7, PS7*100, molarmass1, diff_vol1)
    dif_vap2 = DIFFUSION(TS7, PS7*100, molarmass2, diff_vol2)

    ! mean free path (m)
    ! mean_path = 3*(dif_vap + dif_part)/sqrt(8*Sigma2*TS7/PI*(1./molecmass(jj) + 1./mass))

    !knudsen = mean_path/(radius + r_vap)
    !corr = knudsen + 1.      ! transitional correction factor (Fuchs and Sutugin, 1971)
    !corr = corr/(0.377*knudsen + 1 + 4./(3.*alpha(jj))*(knudsen+knudsen**2))

    collision_rate = 4*PI*(r_vap1 + r_vap1)*(dif_vap1 + dif_vap2)     ! (m^3/s)
  END FUNCTION COLLISION_RATE

  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   Evaluates molecular diffusion coefficient (m^2/s) for condensing vapour (see Reid et al. Properties of Gases and Liquids (1987))
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  REAL FUNCTION diffusion(TS7, PS7, molar_gas, d_gas)
    IMPLICIT NONE

    !INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
    REAL(kind=dp) :: PS7
    REAL(kind=dp) :: TS7, molar_gas, d_gas
    REAL(kind=dp) :: d_air

    REAL(kind=dp), PARAMETER :: MAir = 28.97  ! Molar mass of Air (Seinfeld page1292)

    d_air = 19.70  ! diffusion volume of air (???)

    diffusion = 1.e-7*TS7**1.75*SQRT(1./MAir + 1./molar_gas)
    diffusion = diffusion/((PS7/Patm)*(d_air**(1./3.)+d_gas**(1./3.))**2)
  END FUNCTION diffusion

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  SUBROUTINE DRY_DEPLAYER(k,TS7,PS7,US7,USTR2,s1S7,face,Rk_veloc, Radius)

    IMPLICIT NONE
    !argments
    REAL(kind=dp), DIMENSION(n_bins_par), INTENT(OUT) :: Rk_veloc
    REAL(kind=dp),INTENT(IN) :: TS7,PS7,US7,USTR2,s1S7,face
    INTEGER,                            INTENT(IN) :: k
    ! local
    INTEGER :: ii, zxc=0
    REAL(kind=dp) :: Air_mean_path,Air_viscosity,Air_density,dneedle,Air_kinvis,Re_needle,CB,nB
    REAL(kind=dp), DIMENSION(n_bins_par) :: Corr_coeff,Diff_part,Schmidt,sherwood, IB, &
      taurel, betas, kx, stm, interc, impact, Brownian, Settling_veloc, Radius
    REAL(kind=dp), PARAMETER :: boltzmann = 1.3807d-23
    dneedle=2.d-3  ! diameter of a needle, 2 mm


    ! Ambient air properties
    Air_mean_path = (6.73d-8 * TS7 * (1. + 110.4 / TS7)) / (296. * PS7 / Patm * 1.373)

    Air_viscosity = (1.832d-5 * 406.4 * TS7**1.5) / (5093 * (TS7 + 110.4))

    Air_density   = 1.2929 * T00 * PS7 / (TS7 * Patm)  ! [kg/m^3]

    Air_kinvis=Air_viscosity/Air_density ! [m2 s-1]

    ! Cunningham correction factor (based on Allen and Raabe (1982))
    Corr_coeff = 1. + Air_mean_path / (2. * Radius) * &
      (2.514 + 0.8 * EXP(-0.55 * (2. * Radius / Air_mean_path)))

    ! Collection efficiency of particles due to Brownian diffusion
    Diff_part = boltzmann * TS7 * Corr_coeff / (6. * PI * Air_viscosity * Radius)  ! Diffusivity of particles [m^2/s]

    Schmidt   = Air_kinvis/ Diff_part         ! _Particle_ Schmidt number

    taurel=Corr_coeff * Mass(k,:) / (6. * PI * Air_viscosity * Radius) ! relaxation time
    taurel=Corr_coeff * (4.*pi/3.*Radius**(3.)*1000.) / (6. * PI * Air_viscosity * Radius) ! relaxation time

    if (s1s7==0.) then ! deposition to ground
      Brownian=3.*sqrt(3.)/(29*PI)*Schmidt**(-2./3.)*ustr2  ! no hay like hesston

      Settling_veloc = Corr_coeff * Mass(k,:) * grav/ (6. * PI * Air_viscosity * Radius)

    else
      Re_needle=us7*dneedle/Air_kinvis  ! Reynolds number for a needle

      if (Re_needle<4.d3) then
        CB=0.467
        nB=0.5
        IB=0.99  ! table 3, normal distr
      elseif ((Re_needle .lt. 4.d4) .and.( Re_needle .ge. 4.d3)) then
        CB=0.203
        nB=0.6
        IB=0.99  ! table 3, normal distr
      else
        CB=0.025
        nB=0.8
        IB=0.995  ! table 3, normal distr
      endif

      sherwood=CB*Schmidt**(1./3.)*Re_needle**nB

      Brownian=sherwood*Diff_part/dneedle

      Brownian=s1S7*face*IB*Brownian
      Settling_veloc = Corr_coeff * Mass(k,:) * grav/ (6. * PI * Air_viscosity * Radius)

      kx=0.27
      betas=0.6
      Interc=s1S7*face*2.*kx*2.*Radius/dneedle*us7

      Stm=us7*taurel/dneedle

      Impact=s1S7*face*us7*kx*Stm**2./(2.*betas**2.)*(1./(1.+2.*betas/stm)+log(1.+2.*betas/stm)-1.)

      Brownian=Brownian+Interc+Impact
    end if
    Rk_veloc = Settling_veloc+ Brownian


    ! NOTICE Rk_veloc(k=1) is deposition velocity (m/s).

    if (k==1) then !A2
      Rk_veloc=Rk_veloc/z(2)*2.  ! from deposition velocity (m/s) -> 1/s
    end if

  END SUBROUTINE DRY_DEPLAYER

  !===============================================================================
  ! SUBROUTINE SCPARFLUX(znmax,z,dz,tstep,alt,wup,kt,N_conc_S7,PAR_FLUX,Vol_conc_S7,Rk_vel,ta1,p)
  !
  !   IMPLICIT NONE
  !     REAL(kind=dp), DIMENSION(kz,n_bins_par,n_cond_tot), INTENT(INOUT) :: Vol_conc_S7
  !     REAL(kind=dp), DIMENSION(kz,n_bins_par),            INTENT(INOUT) :: N_conc_S7
  !     REAL(kind=dp), DIMENSION(kz,n_bins_par),              INTENT(OUT) :: PAR_FLUX
  !
  !     REAL(kind=dp), DIMENSION(kz,n_bins_par), INTENT(IN) :: Rk_vel    ! m/s
  !     INTEGER,INTENT(IN) :: znmax
  !     REAL(kind=dp),DIMENSION(kz), INTENT(IN) :: z,dz,alt,kt,wup,ta1,p
  !     REAL(kind=dp), INTENT(IN) :: tstep
  !
  !   ! local
  !   INTEGER :: J,k,I
  !   REAL(kind=dp),DIMENSION(kz):: a,b,c,d_N,d_V, da,db,dc,df
  !   REAL(kind=dp), DIMENSION(kz) :: PAR_VER,PAR_VER_NEW,VOL_VER,VOL_VER_NEW
  !   REAL(kind=dp) :: dkz1,daz,dkz,fktt,fktd
  !
  !   !SIIRRETTY YLEMPaa, CANOPY-laskuista tahan
  !   df(1)=0.
  !
  !   do k=2,znmax-1
  !     da(k)=dz(k+1)*(dz(k)+dz(k+1))
  !     db(k)=dz(k)*dz(k+1)
  !     dc(k)=dz(k)*(dz(k)+dz(k+1))
  !     df(k)=0.5
  !   enddo
  !   df(1)=df(2)
  !
  !   do J = 1, n_bins_par
  !     PAR_VER(1)= N_conc_S7(2,J) - tstep * Rk_vel(1,J)*N_conc_S7(2,j)  ! deposition
  !     PAR_VER(znmax)= N_conc_S7(znmax,J)
  !     if (PAR_VER(1)<0.) then
  !       PAR_VER(1)=0.
  !     end if
  !     do k=2,znmax-1
  !       PAR_VER(k) = N_conc_S7(k,j)   &
  !         - tstep * Rk_vel(k,J)*N_conc_S7(k,J) !deposition
  !
  !
  !       fktt   = 2.*(df(k)*alt(k+1)*kt(k+1) +(1.-df(k))*alt(k)*kt(k))
  !       fktd   = 2.*((1.-df(k-1))*alt(k-1)*kt(k-1) +df(k-1)*alt(k)*kt(k))
  !
  !       ! fktt   = df(k)*alt(k+1)*kt(k+1) +(1.-df(k))*alt(k)*kt(k)
  !       ! fktd   = (1.-df(k-1))*alt(k-1)*kt(k-1) +df(k-1)*alt(k)*kt(k)
  !
  !       a(k)   = tstep*(fktt-wup(k)*dz(k))/da(k)
  !       c(k)   = tstep*(fktd+wup(k)*dz(k+1))/dc(k)
  !       b(k)   = tstep*((dz(k)*fktt+ dz(k+1)*fktd) /(dz(k+1)+dz(k)) &
  !         +wup(k)*(dz(k+1)-dz(k)))/db(k)+ 1.0d0  &
  !         - tstep * Rk_vel(k,J) !deposition
  !     enddo
  !
  !     call gtri(a, c, b, PAR_VER, PAR_VER_NEW, 1, PAR_VER(1), 3.d0, 2, 0.d0, 3.d0, znmax, 2)
  !
  !     ! No flux due to deposition at the lowest level. Flux=0.
  !     !       call gtri(a, c, b, PAR_VER, PAR_VER_NEW, 2,0.d0, 3.d0, 2, 0.d0, 3.d0, znmax, 2)
  !
  !     do k=1,znmax
  !       if (PAR_VER_NEW(k) < 0.) PAR_VER_NEW(k) = 0. !JL
  !       N_conc_S7(k,J) = PAR_VER_NEW(k)
  !     enddo
  !
  !     PAR_FLUX(1,J)=0.
  !     PAR_FLUX(znmax,J)=0.
  !     do  k=2,znmax-1
  !       fktt   = (df(k)*alt(k+1)*kt(k+1) +(1.-df(k))*alt(k)*kt(k))
  !       fktd   = ((1.-df(k-1))*alt(k-1)*kt(k-1) +df(k-1)*alt(k)*kt(k))
  !       !fktt   = 100*100*(df(k)*alt(k+1)*kt(k+1) +(1.-df(k))*alt(k)*kt(k))
  !       !fktd   = 100*100*((1.-df(k-1))*alt(k-1)*kt(k-1) +df(k-1)*alt(k)*kt(k))
  !       !JL mean of two levels, added minus, unit #/cm²/s
  !       PAR_FLUX(k,J)=-0.5*((fktt*(PAR_VER_NEW(k+1) - PAR_VER_NEW(k))   /dz(k+1)) +  &
  !         (fktd*(PAR_VER_NEW(k)   - PAR_VER_NEW(k-1)) /dz(k)))
  !     enddo
  !
  !     do I =1, n_cond_tot
  !       VOL_VER(1)= Vol_conc_S7(2,J,I) &
  !         - tstep * Rk_vel(1,J)*Vol_conc_S7(2,J,I)  ! deposition
  !       do k=2,znmax
  !         VOL_VER(k) = Vol_conc_S7(k,J,I) &
  !           - tstep * Rk_vel(k,J)*Vol_conc_S7(k,J,I) !deposition
  !       enddo
  !       call gtri(a, c, b, VOL_VER, VOL_VER_NEW, 1, VOL_VER(1), 3.d0, 2, 0.d0, 3.d0, znmax, 2)
  !       do k=1,znmax
  !         Vol_conc_S7(k,J,I) = VOL_VER_NEW(k) !_NEW(k) !JL
  !         if (Vol_conc_S7(k,J,I) < 0.) Vol_conc_S7(k,J,I) = 0. !JL
  !       enddo
  !     enddo
  !
  !   enddo
  !
  ! END SUBROUTINE SCPARFLUX

  ! ===============================================================================
SUBROUTINE SCPARFLUX(znmax,z,dz,tstep,alt,wup,kt,N_conc_S7,PAR_FLUX,Vol_conc_S7,Rk_vel,ta1,p,N_conc_Clust)

    IMPLICIT NONE
    REAL(kind=dp), DIMENSION(kz,n_bins_par,n_cond_tot), INTENT(INOUT) :: Vol_conc_S7
    REAL(kind=dp), DIMENSION(kz,n_bins_par),            INTENT(INOUT) :: N_conc_S7
    REAL(kind=dp), DIMENSION(:,:),                      INTENT(INOUT) :: N_conc_Clust
    REAL(kind=dp), DIMENSION(kz,n_bins_par),              INTENT(OUT) :: PAR_FLUX

    REAL(kind=dp), DIMENSION(kz,n_bins_par), INTENT(IN) :: Rk_vel    ! m/s
    INTEGER,INTENT(IN) :: znmax
    REAL(kind=dp),DIMENSION(kz), INTENT(IN) :: z,dz,alt,kt,wup,ta1,p
    REAL(kind=dp), INTENT(IN) :: tstep

    ! local
    INTEGER :: J,k,I
    REAL(kind=dp), DIMENSION(kz) :: PAR_VER,PAR_VER_NEW,VOL_VER,VOL_VER_NEW

    do J = 1, n_bins_par
        ! Boundaries
        PAR_VER(1)      = MAX(0d0, N_conc_S7(2,J) - tstep * Rk_vel(1,J)*N_conc_S7(2,j))  ! settling depos
        PAR_VER(znmax)  = N_conc_S7(znmax,J)

        ! middle layers
        do k=2,znmax-1 !or not?!
            PAR_VER(k)  = N_conc_S7(k,j) - tstep * Rk_vel(k,J)*N_conc_S7(k,J) ! settling depos
        enddo

        call gtri(terma, termc, termb, PAR_VER, PAR_VER_NEW, 1, PAR_VER(1), 3.d0, 2, 0.d0, 3.d0, znmax, 2)
        N_conc_S7(1:znmax,J) = MAX(PAR_VER_NEW(1:znmax), 0d0)

        PAR_FLUX(1,J)=0.
        PAR_FLUX(znmax,J)=0.
        do  k=2,znmax-1
            ! Same form as in chemistry mixing, except not multiplying by 100 since C is m⁻³
            PAR_FLUX(k,j) = -0.5d0 * ( Kht(k)*(PAR_VER_NEW(k+1) &
                            - PAR_VER_NEW(k))/dz(k+1) + Khb(k)*(PAR_VER_NEW(k) &
                            - PAR_VER_NEW(k-1))/dz(k) )
        enddo

        do I =1, n_cond_tot
            VOL_VER(1) = Vol_conc_S7(2,J,I) - tstep * Rk_vel(1,J)*Vol_conc_S7(2,J,I)  ! settling depos
            VOL_VER(znmax) = Vol_conc_S7(znmax,J,I)  ! settling depos

            do k=2,znmax-1
                VOL_VER(k) = Vol_conc_S7(k,J,I) - tstep * Rk_vel(k,J)*Vol_conc_S7(k,J,I) ! settling depos
            enddo

            call gtri(terma, termc, termb, VOL_VER, VOL_VER_NEW, 1, VOL_VER(1), 3.d0, 2, 0.d0, 3.d0, znmax, 2)

            Vol_conc_S7(1:znmax,J,I) = MAX(VOL_VER_NEW(1:znmax), 0d0) !_NEW(k) !JL
        enddo

    enddo

    do J = 1, size(N_conc_Clust, 2)
        ! Boundaries
        PAR_VER(1)      = MAX(0d0, N_conc_Clust(2,J) - tstep * Rk_vel(1,1)*N_conc_Clust(2,j))  ! settling depos, using smallest particle for depo rate
        PAR_VER(znmax)  = N_conc_Clust(znmax,J)

        ! middle layers
        do k=2,znmax-1 !or not?!
            PAR_VER(k)  = N_conc_Clust(k,j) - tstep * Rk_vel(k,1)*N_conc_Clust(k,J) ! settling depos, using smallest particle for depo rate
        enddo

        call gtri(terma, termc, termb, PAR_VER, PAR_VER_NEW, 1, PAR_VER(1), 3.d0, 2, 0.d0, 3.d0, znmax, 2)
        N_conc_Clust(1:znmax,J) = MAX(PAR_VER_NEW(1:znmax), 0d0)

        ! PAR_FLUX(1,J)=0.
        ! PAR_FLUX(znmax,J)=0.
        !
        ! do  k=2,znmax-1
        !     ! Same form as in chemistry mixing, except not multiplying by 100 since C is m⁻³
        !     PAR_FLUX(k,j) = -0.5d0 * ( Kht(k)*(PAR_VER_NEW(k+1) &
        !                     - PAR_VER_NEW(k))/dz(k+1) + Khb(k)*(PAR_VER_NEW(k) &
        !                     - PAR_VER_NEW(k-1))/dz(k) )
        ! enddo


    enddo

END SUBROUTINE SCPARFLUX


  ! subroutine for defining the vertical grid
  subroutine generate_grid(method, kz, hh, z, abl)
    !----- Input -----!

    character(*), intent(in)  :: method  ! method to generate the grid blocks
    integer     , intent(in)  :: kz      ! number of levels
    real(dp)    , intent(in)  :: hh      ! [m], model domain height

    !----- Output -----!

    real(dp), intent(out) :: z(kz)  ! [m], vertical coordinate

    ! Grid type, will be removed in future
    ! 1: linear grid
    ! 2: logarithmic grid
    integer, intent(out)  :: abl

    !----- Local -----!

    real(dp) :: dz  ! distance of the model levels in linear case
    integer  :: k   ! layer index

    !----- Main -----!

    select case ( trim(adjustl(method)) )
    case ('LIN')  ! linear grids
!      hh=23.5d0
      dz = hh/(kz-1)
      z(1) = 0.0_dp
      do k=2,kz
        z(k)=(k-1)*dz
      end do
      z(kz)=hh
    case ('LOG')  ! logarithmic grids
!      hh = 3000.0_dp
      z(1) = 0.0_dp
      do k=2,kz
        z(k)=z(1)+exp(log(hh)*(k-1.0d0)/(kz-1.0d0))-1.0d0
      enddo
      z(kz)=hh
    case ('MAN_ASAM')  ! manually set grids, ASAM case
      z(1) = 0.0_dp
      z(2:kz) = (/2.0d0, 4.0d0, 6.0d0, 8.0d0, 10.0d0, &
        12.0d0, 14.0d0, 16.0d0, 18.0d0, 20.0d0, &
        23.595d0, 27.866d0, 32.879d0, 38.762d0, 45.667d0, &
        53.772d0, 63.284d0, 74.447d0, 87.549d0, 102.92d0, &
        120.97d0, 142.15d0, 167.01d0, 196.19d0, 230.44d0, &
        270.63d0, 317.80d0, 373.17d0, 438.14d0, 514.41d0, &
        603.91d0, 708.96d0, 832.26d0, 976.96d0, 1146.8d0, &
        1346.1d0, 1580.0d0, 1854.6d0, 2176.8d0, 2555.1d0, &
        3000.0d0 /)
    case ('MAN_1')  ! unknown grid settings
      z(1) = 0.0_dp
      z(2:kz) = (/2.5     , 7.5    , 12.5   , 17.5    , 22.75   , 28.525 , 34.88  , &
        41.87   , 49.555 , 58.01  , 67.31   , 77.54   , 88.795 , 101.175, &
        114.79  , 129.765, 146.24 , 164.365 , 184.305 , 206.235, 230.36 , &
        256.9   , 286.085, 318.19 , 353.51  , 392.365 , 435.105, 482.115, &
        533.825 , 590.705, 653.275, 722.1   , 797.81  , 881.095, 972.705, &
        1073.475, 1184.32, 1306.25, 1440.375, 1587.915, 1750.21, 1928.73 /)
    case ('traj')  ! for trajectory case
      z(1:5) = (/0.0d0, 2.0d0, 5.0d0, 10.0d0, 15.0d0/)

      ! 20, 30, 40, 50, 60, 70, 80, 90
      do k=6, 13
        z(k) = 20.0d0 + 10.0d0*(k-6)
      end do

      ! 100 to hh with logarithm grids, here kz should be <= 50 to make the level
      ! just above 100 m higher than 110 m, which means dz is always increasing.
      do k=14, kz
        z(k) = exp( log(100.0d0) + (k-14)*log(hh/100.0d0)/(kz-14) )
      end do
    end select

    abl = 2  ! set it to 2 currently, need to reorganize it in future
  end subroutine generate_grid


  subroutine calculate_grid_parameters()
    ! Grid intervals
    do k=2,kz
      dz(k)=z(k)-z(k-1)
    enddo

    do k=2,kz-1
       da(k)=dz(k+1)*(dz(k)+dz(k+1))
       db(k)=dz(k)*dz(k+1)
       dc(k)=dz(k)*(dz(k)+dz(k+1))
       ! df(k)=(dz(k+1)+dz(k))/dz(k+1)/4.0d0
       df(k)=0.5d0
    enddo
    df(1)=df(2)
  end subroutine calculate_grid_parameters


  subroutine init_rOH()
    CH_oh_flag = .FALSE.
    CH_oh_count = 0

    DO CH_oh_i = 1, NSPEC
      IF (SPC_NAMES(CH_oh_i)(1:3) == 'rOH') THEN
        CH_oh_flag(CH_oh_i) = .TRUE.
        CH_oh_count = CH_oh_count + 1
      ENDIF
    ENDDO

    IF (CH_oh_count > 0) THEN
      ! collect indices of the OH-reactivity pseudochemicals
      ALLOCATE(CH_oh_indices(CH_oh_count))
      CH_oh_count = 0
      DO CH_oh_i = 1, NSPEC
        IF (SPC_NAMES(CH_oh_i)(1:3) == 'rOH') THEN
          CH_oh_count = CH_oh_count + 1
          CH_oh_indices(CH_oh_count) = CH_oh_i
        ENDIF
      ENDDO
      ALLOCATE(CH_oh_cons3(kz,CH_oh_count,3))
      ALLOCATE(CH_oh_prev3(kz,CH_oh_count,3))
      CH_oh_cons3 = 0
      CH_oh_prev3 = 0
    ENDIF
  end subroutine init_rOH


  subroutine init_rO3()
    CH_o3_flag = .FALSE.
    CH_o3_count = 0

    DO CH_o3_i = 1, NSPEC
      IF (SPC_NAMES(CH_o3_i)(1:3) == 'rO3') THEN
        CH_o3_flag(CH_o3_i) = .TRUE.
        CH_o3_count = CH_o3_count + 1
      ENDIF
    ENDDO

    IF (CH_o3_count > 0) THEN
      ! collect indices of the O3-reactivity pseudochemicals
      ALLOCATE(CH_o3_indices(CH_o3_count))
      CH_o3_count = 0
      DO CH_o3_i = 1, NSPEC
        IF (SPC_NAMES(CH_o3_i)(1:3) == 'rO3') THEN
          CH_o3_count = CH_o3_count + 1
          CH_o3_indices(CH_o3_count) = CH_o3_i
        ENDIF
      ENDDO
      ALLOCATE(CH_o3_cons3(kz,CH_o3_count,3))
      ALLOCATE(CH_o3_prev3(kz,CH_o3_count,3))
      CH_o3_cons3 = 0
      CH_o3_prev3 = 0
    ENDIF
  end subroutine init_rO3

  subroutine init_rNO3()
    CH_no3_flag = .FALSE.
    CH_no3_count = 0

    DO CH_no3_i = 1, NSPEC
      IF (SPC_NAMES(CH_no3_i)(1:4) == 'rNO3') THEN
        CH_no3_flag(CH_no3_i) = .TRUE.
        CH_no3_count = CH_no3_count + 1
      ENDIF
    ENDDO

    IF (CH_no3_count > 0) THEN
      ! collect indices of the NO3-reactivity pseudochemicals
      ALLOCATE(CH_no3_indices(CH_no3_count))
      CH_no3_count = 0
      DO CH_no3_i = 1, NSPEC
        IF (SPC_NAMES(CH_no3_i)(1:4) == 'rNO3') THEN
          CH_no3_count = CH_no3_count + 1
          CH_no3_indices(CH_no3_count) = CH_no3_i
        ENDIF
      ENDDO
      ALLOCATE(CH_no3_cons3(kz,CH_no3_count,3))
      ALLOCATE(CH_no3_prev3(kz,CH_no3_count,3))
      CH_no3_cons3 = 0
      CH_no3_prev3 = 0
    ENDIF
  end subroutine init_rNO3

  subroutine init_reactivity()
    CH_reactivity_flag = .FALSE.
    CH_reactivity_count = 0

      DO CH_reactivity_i = 1, NSPEC
        IF (SPC_NAMES(CH_reactivity_i)(1:2) == 'r_') THEN
          CH_reactivity_flag(CH_reactivity_i) = .TRUE.
          CH_reactivity_count = CH_reactivity_count + 1
        ENDIF
      ENDDO

      ALLOCATE(CH_reactivity_indices(CH_reactivity_count))
      CH_reactivity_count = 0
      ! collect indices of the reactivity pseudochemicals
      DO CH_reactivity_i = 1, NSPEC
        IF (SPC_NAMES(CH_reactivity_i)(1:2) == 'r_') THEN
          CH_reactivity_count = CH_reactivity_count + 1
          CH_reactivity_indices(CH_reactivity_count) = CH_reactivity_i
        ENDIF
      ENDDO
      ALLOCATE(CH_reactivity_cons3(kz,CH_reactivity_count,3))
      ALLOCATE(CH_reactivity_prev3(kz,CH_reactivity_count,3))
      CH_reactivity_cons3 = 0
      CH_reactivity_prev3 = 0

  end subroutine init_reactivity



  subroutine prepare_chemistry_stat()
    ! select case (flag_model_type)
    ! case (1)
    ! Get the number concentrations of specific species
    CH_CONS_ALL(:,ind_O3)  = linear_interp(time_in_month, nxodrad, &
      CH_gas_hyy(2,:), dt_obs) * air(:)*1.0d-9  ! only use observation at level nz

    CH_CONS_ALL(:,ind_SO2) = linear_interp(time_in_month, nxodrad, &
      CH_gas_hyy(3,:), dt_obs) * air(:)*1.0d-9

    CH_CONS_ALL(:,ind_NO) = linear_interp(time_in_month, nxodrad, &
      CH_gas_hyy(4,:), dt_obs) * air(:)*1.0d-9

    CH_CONS_ALL(:,ind_NO2) = linear_interp(time_in_month, nxodrad, &
      CH_gas_hyy(5,:), dt_obs) * air(:)*1.0d-9  ! NO2=NOx-NO, already done in preprocessing
    CH_CONS_ALL(:,ind_CO) = linear_interp(time_in_month, nxodrad, &
      CH_gas_hyy(6,:), dt_obs) * air(:)*1.0d-9

    if (flag_NH3==1)then
      NH3_hyy_ALL = linear_interp(time_in_month, nxodrad, &
        NH3_hyy, dt_obs) * air(:)*1.0d-9
    else if (flag_NH3==2)then
      ! if number concentration is used
      NH3_hyy_ALL = linear_interp(time_in_month, nxodrad, NH3_hyy, dt_obs)
    endif

    ! if (flag_emis==0) then
    !   CH_CONS_ALL(2,ind_APINENE) = linear_interp(time_in_month, nxodrad, &
    !     CH_VOC_hyy(1,:), dt_obs) * air(2)*1.0d-9*0.51
    !   CH_CONS_ALL(2,ind_BPINENE) = linear_interp(time_in_month, nxodrad, &
    !     CH_VOC_hyy(1,:), dt_obs) * air(2)*1.0d-9*0.12
    !   CH_CONS_ALL(2,ind_LIMONENE) = linear_interp(time_in_month, nxodrad, &
    !     CH_VOC_hyy(1,:), dt_obs) * air(2)*1.0d-9*0.09
    !   CH_CONS_ALL(2,ind_CARENE) = linear_interp(time_in_month, nxodrad, &
    !     CH_VOC_hyy(1,:), dt_obs) * air(2)*1.0d-9*0.28
    !   CH_CONS_ALL(2,ind_C5H8) = linear_interp(time_in_month, nxodrad, &
    !     CH_VOC_hyy(2,:), dt_obs) * air(2)*1.0d-9
    !   ! CH_CONS_ALL(2,ind_MVK) = linear_interp(time_in_month, nxodrad, &
    !   !   CH_VOC_hyy(3,:), dt_obs) * air(2)*1.0d-9
    !   ! CH_CONS_ALL(2,ind_MEK) = linear_interp(time_in_month, nxodrad, &
    !   !   CH_VOC_hyy(4,:), dt_obs) * air(2)*1.0d-9
    !   CH_CONS_ALL(2,ind_CH3OH) = linear_interp(time_in_month, nxodrad, &
    !     CH_VOC_hyy(5,:), dt_obs) * air(2)*1.0d-9
    !   CH_CONS_ALL(2,ind_CH3CHO) = linear_interp(time_in_month, nxodrad, &
    !     CH_VOC_hyy(6,:), dt_obs) * air(2)*1.0d-9
    !   ! CH_CONS_ALL(2,ind_C2H5OH) = linear_interp(time_in_month, nxodrad, &
    !   !   CH_VOC_hyy(7,:), dt_obs) * air(2)*1.0d-9
    !   CH_CONS_ALL(2,ind_CH3COCH3) = linear_interp(time_in_month, nxodrad, &
    !     CH_VOC_hyy(8,:), dt_obs) * air(2)*1.0d-9
    !   CH_CONS_ALL(2,ind_CH3CO2H) = linear_interp(time_in_month, nxodrad, &
    !     CH_VOC_hyy(9,:), dt_obs) * air(2)*1.0d-9
    !   CH_CONS_ALL(2,ind_H2SO4) = linear_interp(time_in_month, nxodrad, &
    !     CH_H2SO4_hyy, dt_obs)
    ! endif
    ! Read CH4 concentration

    IF ( TRIM(ADJUSTL(STATION)) == 'hyytiala' )  THEN
      ! The CH4 concentration was measured with Picarro G1301 (Picarro Inc.,
      ! USA) and it has been corrected for dilution and spectroscopic
      ! effects by the instrument itself. Data processed by Olli Peltola.
      ! There exists measurements from 16.8, 67.2 and 125 m, but since they
      ! are mostly all most identical, for now I only use the values from
      ! 16.8 m. Only the measurements from year 2014 are good.
      ! For Aug I used the same input as for July, since no data is avail.
      ! According to IPPC then the atmospheric growth rate of methane has
      ! been 6 ppb yr-1 the past 10 years or so. unit of input: ppm.

      ! Read CH4 concentration in 2014 and add trend to obtain data in other years
      CH_CONS_all(:, ind_CH4) = linear_interp(time_in_month, nxodrad, &
        hyy_ch4 + (now_date(1)-2014) * 0.006_dp, dt_obs) * air(:) / 1.0d6
    ELSE
      CH_CONS_all(:,ind_CH4) = (1.8d0 * air(:)) / 1.0d6  ! * 1.0562 !year 2050
    ENDIF

    CH_CONS_ALL(:,ind_H2)  = (0.5d0 * air(:)) * 1.0d-6

    CALL debug_message('DEBUGING: interpolating condensation sinks')

    IF (flag_aero==0)THEN
      ! index_vapor=lookup_no_uhma(vapor_names,SPC_NAMES)
      do i = 1, n_cond_tot
        if (VAPOUR_PROP%ind_ch(i) > 0) then
          vapor(:,i)=max(CH_CONS_ALL(:,VAPOUR_PROP%ind_ch(i))*1d6, 1.0d0)
        end if
      end do
      ! vapor(:,1:n_cond_tot)=max(CH_CONS_ALL(:,index_vapor)*1e6,1.d0)
    ELSE
      CH_cs_H2SO4 = sink(:,VAPOUR_PROP%ind_H2SO4)
    END IF

    CH_cs_HNO3 = linear_interp(time_in_month, nxodrad, stat_cs_HNO3, dt_obs)

    DO k = 1, kz  ! loop over kz
      IF (flag_aero==0)THEN
        !if aero is not used
        CH_cs_H2SO4(k) = linear_interp(time_in_month, nxodrad, stat_cs_H2SO4, dt_obs)
        VAPOUR_PROP%c_sat = saturation_conc_m3(VAPOUR_PROP%psat_a, VAPOUR_PROP%psat_b, ta(k))
        ! CH_RES_org=(vapor(k,1:n_cond_tot)-VAPOUR_PROP%c_sat)/vapor(k,1:n_cond_tot)*CH_cs_H2SO4(k)*0.5d0
        ! CH_RES_org = CH_cs_H2SO4(k)* 0.5d0  ! ?? temporary use for testing
      ! else
        ! CH_RES_org=0.d0
      end if
    end do
  end subroutine prepare_chemistry_stat


  subroutine prepare_chemistry_traj()
    ! Set several species to a fixed mixing ratio for testing, this needs to be
    ! modify or read as input in future
    CH_CONS_ALL(:,ind_O3 ) = 50.0d0 * air(:) * 1.0d-9
    CH_CONS_ALL(:,ind_CH4) = 1800.0d0 * air(:) * 1.0d-9
    CH_CONS_ALL(:,ind_H2 ) = 0.5d0 * air(:) * 1.0d-6
  end subroutine prepare_chemistry_traj


  subroutine calculate_chemistry_parallel()

    mpi_sendcount = 0
    mpi_recvcount = 0
    ! CALL CPU_TIME(wtime1)

    DO k = 1, kz  ! loop over kz
      IF (mpi_sendcount >= mpi_nslaves) THEN
        ! receive data from a slave
        CALL MPI_RECV( &
          mpi_recv_buffer, mpi_recv_buffer_size, MPI_DOUBLE_PRECISION, &
          MPI_ANY_SOURCE, MPI_ANY_TAG, &
          MPI_COMM_WORLD, mpi_status, mpi_rc &
          )

        call store_data_from_slaves(1)

        mpi_source_id = mpi_status(MPI_SOURCE)
        mpi_dest_id = mpi_source_id
        mpi_recvcount = mpi_recvcount + 1
      ELSE
        mpi_dest_id = mpi_sendcount + 1
      END IF  ! mpi_sendcound >= mpi_nslaves

      ! send message "chemicals databuffer coming" to slave
      CALL MPI_SEND( &
        mpi_do_chemistry_code, 1, MPI_INTEGER, &
        mpi_dest_id, mpi_task_code_tag, &
        MPI_COMM_WORLD, mpi_rc &
        )

      ! Prepare all things for sending data

      ! CH_CONS(:) = CH_CONS_ALL(k,:)

      mpi_send_buffer(index_kz_tag                 ) = REAL(k,KIND=dp)
      mpi_send_buffer(index_time1                  ) = time_in_day
      mpi_send_buffer(index_time2                  ) = time_in_day+ dt_chem
      mpi_send_buffer(index_tmcontr                ) = time_in_month
      mpi_send_buffer(index_year                   ) = REAL(now_date(1), dp)
      mpi_send_buffer(index_month                  ) = REAL(now_date(2), dp)
      mpi_send_buffer(index_nxodrad                ) = REAL(nxodrad,KIND=dp)
      mpi_send_buffer(index_ta1                    ) = ta1(k)
      mpi_send_buffer(index_pr1                    ) = pres(k)
      mpi_send_buffer(index_zenith_deg             ) = zenith_deg
      mpi_send_buffer(index_sun_par                ) = EM_SUN_PAR(k)
      mpi_send_buffer(index_res1                   ) = SINK(k,VAPOUR_PROP%ind_H2SO4) !CH_cs_H2SO4(k)
      mpi_send_buffer(index_res2                   ) = CH_cs_HNO3(k)
      mpi_send_buffer(index_qa1                    ) = qa1(k)
      mpi_send_buffer(index_air                    ) = air(k)
      mpi_send_buffer(index_N2                     ) = N2(k)
      mpi_send_buffer(index_O2                     ) = O2(k)
      mpi_send_buffer(index_H2O                    ) = H2O(k)
      mpi_send_buffer(index1_CONS:index2_CONS      ) = CH_CONS_ALL(k,:) !CH_CONS
      mpi_send_buffer(index_albedo                 ) = albedo
      mpi_send_buffer(index_glob                   ) = glob
      ! mpi_send_buffer(index1_RES_org:index2_RES_org) = CH_RES_org

      ! send chemicals databuffer to slave
      CALL MPI_SEND(mpi_send_buffer,mpi_send_buffer_size,MPI_DOUBLE_PRECISION,mpi_dest_id,mpi_buffer_tag, &
        MPI_COMM_WORLD,mpi_rc)
      mpi_sendcount = mpi_sendcount + 1
    ENDDO  ! loop over kz

    ! reveive data from the rest of the slaves - these were not yet gathered in the previous loop
    DO slave_i = 1,mpi_nslaves

      ! receive data from a slave
      CALL MPI_RECV(mpi_recv_buffer,mpi_recv_buffer_size, &
        MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG, &
        MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_rc)

        call store_data_from_slaves(1)

    END DO

    ! To remove small fluctuations around zero
    where (CH_CONS_ALL<0d0.and.CH_CONS_ALL>-1d-3) CH_CONS_ALL = 0d0


  end subroutine calculate_chemistry_parallel


  subroutine calculate_aerosol_parallel()
    ! call start_timer(timer_total_aero)
    mpi_sendcount = 0
    mpi_recvcount = 0

    ! Condensation sink name in chemistry is CH_CS, while in aerosol it's called
    ! G_ambient%sink(1)
    ! G_ambient%sink(1)=CH_CS

    DO k = 1, kz  ! loop over kz
      ! Assign environmental data

      IF (mpi_sendcount >= mpi_nslaves) THEN

        ! receive data from a slave
        CALL MPI_RECV(new_mpi_recv_aer_buf, new_mpi_recv_aer_buf_size, &
          MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
          MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS, &
          mpi_rc)

        call store_data_from_slaves(2)

        mpi_source_id = mpi_status(MPI_SOURCE)
        mpi_recvcount = mpi_recvcount + 1
        mpi_dest_id   = mpi_source_id
      ELSE
        mpi_dest_id = mpi_sendcount + 1
      END IF

      ! send message "aerosol databuffer coming" to slave
      CALL MPI_SEND(mpi_do_aerosol_code, 1, MPI_INTEGER, &
        mpi_dest_id, mpi_task_code_tag, &
        MPI_COMM_WORLD, &
        mpi_rc)

      new_mpi_send_aer_buf(ind_kz_tag)                          = REAL(k, KIND=dp)
      new_mpi_send_aer_buf(index_time1)                         = time_in_day
      new_mpi_send_aer_buf(index_time2)                         = time_in_day + dt_aero
      new_mpi_send_aer_buf(ind_aer_ts)                          = dt_aero
      new_mpi_send_aer_buf(ind_tempk)                           = ta1(k)
      new_mpi_send_aer_buf(ind_pres)                            = pres(k)
      new_mpi_send_aer_buf(ind_rh)                              = RH(k)
      new_mpi_send_aer_buf(ind_U10m)                            = hwind(4)
      new_mpi_send_aer_buf(ind_lsm)                             = land_sea_mask
      new_mpi_send_aer_buf(ind_firstcall)                       = 0
      new_mpi_send_aer_buf(ind_IPR)                             = 3d0
      new_mpi_send_aer_buf(ind1_emiaer:ind2_emiaer)             = emi_aer(k, :)
      if (k<=4) THEN ! Sea spray emitted on first 10 meters
          new_mpi_send_aer_buf(ind1_seasalt:ind2_seasalt)       = sea_salt_psd(:)
      else
          new_mpi_send_aer_buf(ind1_seasalt:ind2_seasalt)       = 0d0
      end if
      new_mpi_send_aer_buf(ind1_n_conc:ind2_n_conc)             = N_CONC(k, :)
      new_mpi_send_aer_buf(ind1_CONS_clust:ind2_CONS_clust)     = CLUSTERS(k, :)
      new_mpi_send_aer_buf(ind1_CONS_aer:ind2_CONS_aer)         = CH_CONS_ALL(k, :)
      new_mpi_send_aer_buf(ind1_composition:ind2_composition)   = RESHAPE(MASS_COMPO(k,:,:), [n_bins_par*n_cond_tot])

      CALL MPI_SEND(new_mpi_send_aer_buf, new_mpi_send_aer_buf_size, &
        MPI_DOUBLE_PRECISION, mpi_dest_id, mpi_buffer_tag, &
        MPI_COMM_WORLD, mpi_rc)

      mpi_sendcount = mpi_sendcount + 1

    ENDDO  ! loop over kz

    ! reveive data from the cores which are still working
    DO slave_i = 1, mpi_nslaves
      ! receive data from a slave
      CALL MPI_RECV(new_mpi_recv_aer_buf,new_mpi_recv_aer_buf_size, &
        MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG, &
        MPI_COMM_WORLD, MPI_STATUS,mpi_rc)

        call store_data_from_slaves(2)

    END DO  ! slave_i = 1, mpi_nslaves

    ! if (MOD(int(time_in_day), 600) == 0) print'(2(a,es9.2),a,f4.1," ",a,9(es9.2))', &
    !     'NH3', sum(CH_CONS_ALL(1:2,ind_NH3)),' DMA', sum(CH_CONS_ALL(1:2,ind_DMA)), ' J (sum, n, +, -, ...,org) cm⁻³ @ ',z(2),'m:', Formation_rates(2,:)
    ! print*,'time', time_in_day, 'OH', CH_CONS_ALL(3,ind_OH), 'CS', SINK(3,VAPOUR_PROP%ind_H2SO4), 'swr', glob
    ! print*, '   ', 'SO2', CH_CONS_ALL(3,ind_SO2), 'H2SO4', CH_CONS_ALL(3,ind_H2SO4), 'PN tot', SUM(N_CONC)


  end subroutine calculate_aerosol_parallel

  subroutine store_data_from_slaves(process)
      IMPLICIT NONE
      integer, intent(in) :: process
      integer :: i

      select case (process)

      case (1) ! 1: Chemistry
          kz_tag = INT(mpi_recv_buffer(index_kz_tag))
          CH_RO2(kz_tag) = mpi_recv_buffer(index_RO2)
          CH_J_values_ALL(kz_tag,:) = mpi_recv_buffer(index1_J_values_recv:index2_J_values_recv)
          CH_K_values_ALL(kz_tag,:) = mpi_recv_buffer(index1_K_values_recv:index2_K_values_recv)
          CH_ACF_ALL(kz_tag,:) = mpi_recv_buffer(index1_ACF_recv:index2_ACF_recv)
          CH_CONS_ALL(kz_tag,:) = mpi_recv_buffer(index1_CONS_recv:index2_CONS_recv)

      case (2) ! 2: Aerosol

        kz_tag                 = INT(new_mpi_recv_aer_buf(ind_kz_tag))
        MASS_COMPO(kz_tag,:,:) = RESHAPE(new_mpi_recv_aer_buf(ind1_composition:ind2_composition), [n_bins_par, n_cond_tot])
        N_CONC(kz_tag,:)       = new_mpi_recv_aer_buf(ind1_n_conc:ind2_n_conc)
        CH_CONS_ALL(kz_tag, :) = new_mpi_recv_aer_buf(ind1_CONS_aer:ind2_CONS_aer)
        ! NUC_RATE(kz_tag)       = new_mpi_recv_aer_buf(ind1_acdc_rates)
        GR(kz_tag,:)           = new_mpi_recv_aer_buf(ind1_gr:ind2_gr)
        SINK(kz_tag,:)         = new_mpi_recv_aer_buf(ind1_sink:ind2_sink)
        CLUSTERS(kz_tag,:)     = new_mpi_recv_aer_buf(ind1_CONS_clust:ind2_CONS_clust)
        Formation_rates(kz_tag,1:ind2_acdc_rates-ind1_acdc_rates+1) = new_mpi_recv_aer_buf(ind1_acdc_rates:ind2_acdc_rates)
        Formation_rates(kz_tag,ind2_acdc_rates-ind1_acdc_rates+2)   = new_mpi_recv_aer_buf(ind1_oth_nuc_rates)
        NUC_RATE(kz_tag) = sum(Formation_rates(kz_tag,[(1+(i*4), i=1,size(G_ACDC))])) + Formation_rates(kz_tag, size(Formation_rates,2))
      end select

  end subroutine store_data_from_slaves


  subroutine update_location(lat, lon, flag_model_type)
    real(dp), intent(inout) :: lat, lon  ! [deg]
    integer , intent(in) :: flag_model_type

    select case (flag_model_type)
    ! The station is not moving
    case (STATION_MODE)
      continue
    ! Interpolate lat and lon
    case (TRAJECTORY_MODE)
      call linear_interp_2p( &
        infield%var(infield%id_time)%f1d(input_ptr), &
        infield%var(infield%id_time)%f1d(input_ptr+1), &
        infield%var(infield%id_lat)%f1d(input_ptr), &
        infield%var(infield%id_lat)%f1d(input_ptr+1), &
        time, lat_deg &
        )
      call linear_interp_2p( &
        infield%var(infield%id_time)%f1d(input_ptr), &
        infield%var(infield%id_time)%f1d(input_ptr+1), &
        infield%var(infield%id_lon)%f1d(input_ptr), &
        infield%var(infield%id_lon)%f1d(input_ptr+1), &
        time, lon_deg &
        )
    end select
  end subroutine update_location


! Two methods:
!   - Salter et al. doi:10.5194/acp-15-11047-2015
!   - Ovadnevaite et al. doi:10.5194/acp-14-1837-2014 This method is clever in the way that SST is encorporated in the viscosity, but
!     it needs significant wave height (Hw) as input. I converted Hw to U10, so this method is somewhat truncated. Both give
!     very similar results, but Salter has less particles in smallest size (<20 nm)
subroutine calculate_sea_salt_emissions(dFlogdD,U10,T,mask)
    implicit none

    real(dp) :: ReHw,F(5), dFlogdD(:)
    real(dp), save      :: Sigma(5)                                   ! [-] geometric standard deviation of the mode
    real(dp), save      :: CMD(5)                                     ! [µm] Count median diameter (CMD)
    real(dp), PARAMETER :: A(3) = [-5.2168d5, 0d0,0d0]                ! Fitting parameters in the 3 modes in Salter method
    real(dp), PARAMETER :: B(3) = [3.31725d7,7.374d5,1.4210d4]        ! Fitting parameters in the 3 modes in Salter method
    real(dp), PARAMETER :: C(3) = [-6.95275d8, -2.4803d7, 1.4662d7]   ! Fitting parameters in the 3 modes in Salter method
    real(dp), PARAMETER :: D(3) = [1.0684d10, 7.7373d8, 1.7075d8]     ! Fitting parameters in the 3 modes in Salter method
    real(dp), PARAMETER :: aa = 1d0/3.80595262d-07                    ! To convert U10 to ReHw. Crude, but we don't know the wave height! Obtained by
    real(dp), PARAMETER :: bb = 0.010845103482110157d0/1.11531369d+01 ! Transforming U10 -> Re using figure 1 in doi:10.5194/acp-14-1837-2014
    real(dp), PARAMETER :: cc = 2.3383131276948776d0
    real(dp)            :: U10, mask, T                               ! Wind speed at 10 m and Land-Sea mask (1=land), SST (°C)
    integer             :: ii, n_modes
    LOGICAL, PARAMETER  :: Salter_method = .true.                     ! from doi:10.5194/acp-15-11047-2015

    IF (Salter_method) THEN
        n_modes = 3
        CMD(1:n_modes) = [0.095,0.6,1.5]
        Sigma(1:n_modes) = [2.10,1.72,1.60]
        do ii=1,n_modes
            F(ii) = (2d-8 * U10**3.41) * (A(ii) * T**3 + B(ii) * T**2 + C(ii) * T + D(ii) )
        end do

    ELSE
        n_modes = 5
        ReHw  = aa * (exp(bb*U10**cc) - 1d0 )
        Sigma(1:n_modes) = [1.37, 1.5, 1.42, 1.53, 1.85]
        CMD(1:n_modes) = [0.018,0.041,0.09,0.23,0.83]

        F(1) = 104.5  * (max(ReHw - 1e5,0d0))**0.556
        F(2) = 0.0442 * (max(ReHw - 1e5,0d0))**1.08
        F(3) = 149.6  * (max(ReHw - 1e5,0d0))**0.545
        F(4) = 2.96   * (max(ReHw - 1e5,0d0))**0.784
        F(5) = 0.51   * (max(ReHw - 2e5,0d0))**0.87

    END IF

    dFlogdD = 0d0

    do ii=1,n_modes
        dFlogdD = dFlogdD + F(ii) /(SQRT(2*PI)*log10(Sigma(ii))) * EXP( -0.5 * (log10(1d6*current_psd%diameter_fs/CMD(ii))/log10(Sigma(ii)))**2 )
    end do

    dFlogdD = dFlogdD * (1-mask)
    ! if (mask > 0.9) dFlogdD = dFlogdD * 0.01 ! since then the location is not open ocean

end subroutine calculate_sea_salt_emissions

end program Sosa
