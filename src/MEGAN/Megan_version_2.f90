MODULE MEGAN_version_2

  !***********************************************************************
  !
  ! MEGAN version 2.1 as according to Guenther et al., Geosci. Model Dev., 5, 1471-1492, 2012
  !
  !***********************************************************************
  !
  !   Scientific algorithm
  !
  !             Emission = [EF] [GAMMA]
  !                        where 
  !                        [EF]    = emission factor (ug/m2h)
  !                        [GAMMA] = emission activity factor (non-dimension)
  !
  !             GAMMA    = [GAMME_CE] [GAMMA_age] [GAMMA_SM]
  !                        where 
  !                        [GAMME_CE]  = canopy correction factor
  !                        [GAMMA_A] = leaf age correction factor
  !                        [GAMMA_S]  = soil moisture correction factor
  !                        Assumption: [GAMMA_SM]  = 1  (11/27/06)
  !             
  !             GAMME_CE = [GAMMA_LAI] [GAMMA_P] [GAMMA_T]
  !                        where 
  !                        [GAMMA_LAI] = leaf area index factor
  !                        [GAMMA_P]   = PPFD emission activity factor
  !                        [GAMMA_T]   = temperature response factor
  !
  !             Emission = [EF] [GAMMA_LAI] [GAMMA_P] [GAMMA_T] [GAMMA_A]
  !                        Derivation:
  !                        Emission = [EF] [GAMMA_etc] (1-LDF) + [EF] [GAMMA_etc] [LDF] [GAMMA_P]
  !                        Emission = [EF] [GAMMA_etc] {(1-LDF) + [LDF] [GAMMA_P]}
  !                        Emission = [EF] [GAMMA_etc] {(1-LDF) + [LDF] [GAMMA_P]}
  !                        where LDF = light dependent function (non-dimension)
  !                                    (See LD_FCT.EXT)
  !
  !             Final Equation
  !             
  !             Emission = [EF] [GAMMA_LAI] [GAMMA_T] [GAMMA_age] * { (1-LDF) + [LDF] [GAMMA_P] }
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE M2_Canopy
  USE M2_GAMMA_ETC        
  USE M2_LD_FCT
  USE M2_SPC_MGN
  USE constants_mod, only : dp
  USE Sosa_data, ONLY: input_dir_station_info

  IMPLICIT NONE

  INTEGER, PARAMETER :: TSTLEN = 1        !Time step using for calculating GAM_AGE [days]
  INTEGER, PARAMETER :: NCOLS =1          ! Number of columns
  INTEGER, PARAMETER :: NROWS =1          ! Number of rows

CONTAINS 

  SUBROUTINE EMISSION_M2(Day, Layersin, Beta, LATin, LONGin, DATEin, TIMEin, PPFDin, D_PPFDin,         &
       TEMPin, DTEMPin, PRESin, HUMin, WINDin, SMOISTin, LADpin, LAIcin, LAIpin, ERout,          &
       VARout, ER_HB_out, ER_SB_out, ER_NT_out, ER_BT_out, GAM_PHO_out, GAM_CE_out, GAM_T_out,   &
       GAM_TMP_out, GAM_OTHER_out, Zin, kzin, EMI, sun_par, Sunleaftk, Shadeleaftk, sunfrac, station, &
       EMI_MEAS_MONO, CHAM_temp, Treeflag, year_int,check_rh, Sunleaftk_in, Shadeleaftk_in, &
       Synflag)

    REAL(kind=dp) :: BIOMASS !Foliage mass
    REAL(kind=dp) :: BIOMASS_MAX !Seasonal maximum foliage mass (if needle tree: BIOMASS = BIOMASS_MAX)
    real :: EMI_MEAS_MONO, CHAM_temp !Measured emission rates and temperature from chambers
    INTEGER, INTENT(IN) ::  Layersin ,kzin        ! Number of Layers with vegetation in them, presumed from ground 0 to Layersin.
    INTEGER, INTENT(IN) ::  DATEin, TIMEin        ! Date (yyyddd) and time (hhmmss)
    CHARACTER(*), INTENT(IN) ::  station      
    INTEGER, INTENT(IN) :: Treeflag, year_int, Synflag !Treeflag is which tree species
    REAL(kind=dp) :: CAN_HEI !canopy height [m]
    REAL(kind=dp) :: CAN_DEP !canopy depth [m]
    REAL(kind=dp), PARAMETER     ::  Avog = 6.0221E23      ! Avogadro-number [molecules cm3]
    REAL(kind=dp), SAVE                ::  BIOMASS_func(366)
    REAL(kind=dp), SAVE                ::  TDAY(366,3)
    REAL(kind=dp), SAVE                ::  EFin(366,22)          ! emission potentials
    REAL(kind=dp), INTENT(IN)    ::  LATin, LONGin         ! LATitutde and LONGitude
    REAL(kind=dp), INTENT(IN)    ::  DTEMPin               ! TEMPerature (K) and Daily average TEMPerature (K)
    REAL(kind=dp), INTENT(IN)    ::  PPFDin, D_PPFDin      ! PAR and Daily average of PAR in (umol/m2/s)
    REAL(kind=dp), INTENT(IN)    ::  PRESin                ! Pressure (Pa)
    REAL(kind=dp), INTENT(IN)    ::  LAIcin, LAIpin       ! LAIc is LAI for current timestep or month, LAIp for previous
    REAL(kind=dp), INTENT(IN)    ::  SMOISTin              ! Soil moisture [m^3 m^-3]
    REAL(kind=dp), INTENT(IN), DIMENSION(kzin) ::  LADpin  ! Leaf area density as [0,1] in the canopy, or in other words the 
                                                           ! fraction of leaf available in each layer with respect to the total leaf 
                                                           ! amount in the whole column
    REAL(kind=dp), INTENT(IN), DIMENSION(:)    ::  TEMPin
    REAL(kind=dp), INTENT(IN), DIMENSION(:)    ::  HUMin   ! Humidity (g/kg)
    REAL(kind=dp), INTENT(IN), DIMENSION(:)    ::  Zin     ! Layer top array
    REAL(kind=dp), DIMENSION(:), INTENT(IN) :: WINDin !Wind speed (m/s)? 
    REAL(kind=dp), INTENT(INOUT) :: ERout(N_MGN_SPC, kzin) ! Emission rates output (normally g/s, can be set with a boolean to ton/hr)
    CHARACTER*16       VARout(N_MGN_SPC)       !VOC names in the same order as ERout
    REAL(kind=dp), DIMENSION(N_MGN_SPC, kzin),    INTENT(OUT) ::  ER_HB_out, ER_SB_out, ER_NT_out, ER_BT_out  !emissions by plant functional type
    REAL(kind=dp), DIMENSION(4, kzin), INTENT(OUT)            ::  GAM_PHO_out, GAM_CE_out, GAM_T_out          ! Gamma factors for light, canopy environment and temperature
    REAL(kind=dp), DIMENSION(N_MGN_SPC, kzin),    INTENT(OUT) ::  GAM_TMP_out                                 ! Gamma factor for what?
    REAL(kind=dp), DIMENSION(N_MGN_SPC, 3),       INTENT(OUT) ::  GAM_OTHER_out                               !gamma factors soil moisture, leaf age and lai correction
    REAL(kind=dp),DIMENSION(kzin) ::  sun_par, Sunleaftk, Shadeleaftk, sunfrac, Sunleaftk_in, Shadeleaftk_in
    REAL(kind=dp), DIMENSION(Layersin) :: midpoints !Layer centerpoints for gamme_CE
    REAL(kind=dp), DIMENSION(30) :: check_rh
    REAL(kind=dp), DIMENSION(N_MGN_SPC, kzin) :: ERtemporary
    REAL(kind=dp) :: ISO_EMI, MBO_EMI, MYR_EMI, SAB_EMI, LIM_EMI, CAR_EMI, OCI_EMI, BPI_EMI, API_EMI, FAR_EMI, BCA_EMI, & 
         MET_EMI, ACT_EMI, ACA_EMI, FOR_EMI, CH4_EMI, NO_EMI,  OMT_EMI, OSQ_EMI, CO_EMI,  LIN_EMI, CIN_EMI, &
         emi_factor, Beta, TSTEP2
    REAL(kind=dp) :: SEP_mono, SEP_birch
    REAL(kind=dp), DIMENSION(kzin,22) :: EMI

    ! Program name
    CHARACTER*16  :: PROGNAME = 'MEGAN'

    !..  model parameters
    INTEGER       TS
    INTEGER       SDATE                            ! Start date YYYYDDD
    INTEGER       STIME                            ! Start time HHMMSS
    INTEGER       TSTEP                            ! time step

    !..  Internal parameters
    ! internal paramters (status and buffer)
    INTEGER       IOS                              ! i/o status
    CHARACTER*256 MESG                             ! message buffer

    INTEGER, PARAMETER :: NEMIS = N_MGN_SPC ! number of MEGAN species

    ! local variables and their descriptions:
    REAL(kind=dp)          LDF                              ! Light dependent factor

    REAL(kind=dp), ALLOCATABLE    :: ER_BT( :,: )           ! output emission buffer
    REAL(kind=dp), ALLOCATABLE    :: ER_NT( :,: )           ! output emission buffer
    REAL(kind=dp), ALLOCATABLE    :: ER_SB( :,: )           ! output emission buffer
    REAL(kind=dp), ALLOCATABLE    :: ER_HB( :,: )           ! output emission buffer  
    REAL(kind=dp), ALLOCATABLE    :: LAT( :,: )             ! input latitude of grid cell
    REAL(kind=dp), ALLOCATABLE    :: LONG( :,: )            ! input longitude of grid cell
    REAL(kind=dp), ALLOCATABLE    :: LAIp( :,: )            ! previous timestep or month LAI
    REAL(kind=dp), ALLOCATABLE    :: LAIc( :,: )            ! current LAI
    REAL(kind=dp), ALLOCATABLE    :: TEMP(:)                ! input hourly temperature (K)
    REAL(kind=dp), ALLOCATABLE    :: dummyTEMP(:, :, :)     ! for gamma_TNSIP
    REAL(kind=dp), ALLOCATABLE    :: PPFD( :,: )            ! calculated PAR (umol/m2.s)
    REAL(kind=dp), ALLOCATABLE    :: D_PPFD( :,: )          ! daily PAR (umol/m2.s)
    REAL(kind=dp), ALLOCATABLE    :: D_TEMP( :,: )          ! input daily temperature (K) - currently 0K
    REAL(kind=dp), ALLOCATABLE    :: GAM_PHO(:, :,: )       ! light correction factor
    REAL(kind=dp), ALLOCATABLE    :: GAM_TMP( :,:,: )       ! temperature correction factor
    REAL(kind=dp), ALLOCATABLE    :: GAM_LHT( :,: )         ! LAI correction factor
    REAL(kind=dp), ALLOCATABLE    :: GAM_AGE( :,: )         ! leaf age correction factor
    REAL(kind=dp), ALLOCATABLE    :: GAM_SMT( :,: )         ! Soil moisture correction factor
    REAL(kind=dp), ALLOCATABLE    :: GAM_T(:, :,: )
    REAL(kind=dp), ALLOCATABLE    :: GAM_CE (:, :,: )
    REAL(kind=dp), ALLOCATABLE    :: Wind(:) 
    REAL(kind=dp), ALLOCATABLE    :: Pres( :,: ) 
    REAL(kind=dp), ALLOCATABLE    :: Humidity(:) 
    REAL(kind=dp), ALLOCATABLE    :: DI( :,: ) 
    REAL(kind=dp), ALLOCATABLE    :: SMOIS( :,: )
    REAL(kind=dp), ALLOCATABLE    :: GAMMAfactors( :,:, :)  ! gamma factors (  GAMMA_T, GAMMA_PHO, GAMMA_CE), In layers, for different plant types 
    REAL(kind=dp), ALLOCATABLE    :: LADp(:)
    INTEGER, ALLOCATABLE :: Cantype( :,: ) 
    INTEGER, ALLOCATABLE :: ISLTYP( :,: )

    INTEGER          ILINE                         ! current line
    CHARACTER(LEN=1000) LINE                       ! input line buffer
    INTEGER       CID, INX, INY                    ! Input grid x and y
    INTEGER, PARAMETER :: MXTCOL = 9               ! Columns in an input line
    CHARACTER*30     SEGMENT( MXTCOL )             ! Input line fields

    INTEGER, PARAMETER :: NVARS = MXTCOL - 3       ! Number of LAT, LONG, and
    ! PFT factor variables
    CHARACTER*16 VNAME( NVARS )                    ! Variable names
    INTEGER, PARAMETER :: NPFT = 4                 ! guess, number of plant functional types
    REAL(kind=dp), ALLOCATABLE :: PFTFBUF( :, :, : )        ! PFT factor array+Lat and LONG
    REAL(kind=dp), ALLOCATABLE :: PFTF( :, :, : )           ! PFT factor array

    CHARACTER*16   VARIABLENAMES( 30 ) 

    INTEGER, PARAMETER:: NrTyp = 4, NrCha = 16
    REAL(kind=dp),DIMENSION(NrCha,NrTyp ),SAVE:: Canopychar

    ! REAL(kind=dp), DIMENSION(:,:) :: Wind, Pres,Humidity,DI
    ! REAL(kind=dp),DIMENSION(NrCha,NrTyp) :: Canopychar 
    ! INTEGER,DIMENSION(:,:) :: Cantype

    ! loop indices
    INTEGER       T, S, I, J ,K, VAR, VAR2, lay     ! Counters
    INTEGER       NMAP            ! Index
    INTEGER       IDATE           ! Looping date YYYYDDD
    INTEGER       ITIME           ! Looping time HHMMSS

    ! times
    INTEGER       MON             ! Month from YYYYDDD
    INTEGER       DAY             ! Day from YYYYDDD

    INTEGER :: test1,test2,test3

    LOGICAL, SAVE :: first_call = .TRUE.

    ! months
    CHARACTER*3   MONTHS( 12 )
    DATA          MONTHS &
         &  / 'JAN' , 'FEB' , 'MAR' , 'APR' , 'MAY' , 'JUN' , &
         &    'JUL' , 'AUG' , 'SEP' , 'OCT' , 'NOV' , 'DEC'   /
    CHARACTER*2   MONNUM( 12 )
    DATA          MONNUM &
         &  /  '1 ' ,  '2 ' ,  '3 ' ,  '4 ' ,  '5 ' ,  '6 ' , &
         &     '7 ' ,  '8 ' ,  '9 ' ,  '10' ,  '11' ,  '12'   /


    !**********************************************************************

    !======================================================================
    !..  Begin program
    !======================================================================
    !----------------------------------------------------------------------
    !....1) File set up and assign I/O parameters
    !----------------------------------------------------------------------




    !..  Get input parameters

    SDATE= DATEin !2001187 !'Model start date (YYYYDDD)'
    STIME= TIMEin !'Model start time (HHMMSS)'
    !TSTEP= 1000 !'Model time step (HHMMSS)' ! not used

    DO S = 1, NEMIS
       VARIABLENAMES( S ) = TRIM( MGN_SPC( S ) )
    ENDDO

    VARIABLENAMES(NEMIS + 1) = 'D_TEMP'

    VARIABLENAMES(NEMIS + 2) = 'D_PPFD'

    !----------------------------------------------------------------------
    !....2) Process emission rates
    !----------------------------------------------------------------------
    !..  Allocate memory

    ALLOCATE ( ER_BT( NCOLS, NROWS ), STAT = IOS )
    ALLOCATE ( ER_NT( NCOLS, NROWS ), STAT = IOS )
    ALLOCATE ( ER_SB( NCOLS, NROWS ), STAT = IOS )
    ALLOCATE ( ER_HB( NCOLS, NROWS ), STAT = IOS )
    ALLOCATE ( LAT( NCOLS, NROWS ), STAT = IOS ) 
    ALLOCATE ( LONG( NCOLS, NROWS ), STAT = IOS ) 
    ALLOCATE ( LAIp( NCOLS, NROWS ), STAT = IOS )
    ALLOCATE ( LAIc( NCOLS, NROWS ), STAT = IOS )
    ALLOCATE ( D_PPFD( NCOLS, NROWS ), STAT = IOS )
    ALLOCATE ( D_TEMP( NCOLS, NROWS ), STAT = IOS )
    ALLOCATE ( PPFD( NCOLS, NROWS ), STAT = IOS )
    ALLOCATE ( TEMP((Layersin+1) ), STAT = IOS ) 
    ALLOCATE (  dummyTEMP(NCOLS, NROWS, Layersin) )
    ALLOCATE ( GAM_PHO(NPFT, NCOLS, NROWS ), STAT = IOS )
    ALLOCATE ( GAM_TMP( NCOLS, NROWS, Layersin ), STAT = IOS )
    ALLOCATE ( GAM_LHT( NCOLS, NROWS ), STAT = IOS )
    ALLOCATE ( GAM_AGE( NCOLS, NROWS ), STAT = IOS )
    ALLOCATE ( GAM_SMT( NCOLS, NROWS ), STAT = IOS )
    ALLOCATE ( GAM_T(NPFT,NCOLS, NROWS ), STAT = IOS )      
    ALLOCATE ( Wind((Layersin+1)), STAT = IOS )
    ALLOCATE ( Pres( NCOLS, NROWS ), STAT = IOS ) 
    ALLOCATE ( Humidity((Layersin+1) ), STAT = IOS )
    ALLOCATE ( DI( NCOLS, NROWS ), STAT = IOS )
    ALLOCATE ( Cantype ( NCOLS, NROWS ), STAT = IOS )
    ALLOCATE ( GAM_CE (NPFT, NCOLS, NROWS ), STAT = IOS )
    ALLOCATE ( ISLTYP ( NCOLS, NROWS ), STAT = IOS )
    ALLOCATE ( SMOIS ( NCOLS, NROWS ), STAT = IOS )

    ! Allocate memory for PFTF
    ALLOCATE ( PFTF( NPFT, NCOLS, NROWS ), STAT = IOS )
    ALLOCATE (GAMMAfactors(3,4, Layersin))
    ALLOCATE (LADp(Layersin))

    !...  Get input file unit

    ! 1 = LAT, 2 = LONG,
    ! 3 = PFTF_BT, 4 = PFTF_NT,
    ! 5 = PFTF_SB, 6 = PFTF_HB
!!! Order does matter !!!
    !well, not the same order anymore


    !.. Partial fractions of different plant types 
    !.. (doesn't need to sum to 100%, because parts of the grid can be barren)
    !.. You should hardcode appropriate values for your model here, or this can be added to input parameters
    if (TRIM(ADJUSTL(station)) == 'WLG')  then
       PFTF(1, 1, 1) = 43.7  ! percentage of Disturbed grassland
       PFTF(2, 1, 1) = 26.5  ! percentage of Thornveld
       PFTF(3, 1, 1) = 15.9  ! percentage of Maize field 
       PFTF(4, 1, 1) = 13.9  ! percentage of  Disturbed grassland
       
    else
       PFTF(1, 1, 1) = 0 !PFTF_BT  !percentage of Broadleaf trees
       PFTF(2, 1, 1) = 100 !PFTF_NT  !percentage of Needle leaf trees
       PFTF(3, 1, 1) = 0 !PFTF_SB  !percentage of Shurblands
       PFTF(4, 1, 1) = 0 !PFTF_HB  !percentage of Herbaceous cover
    endif

    !..  Get D_PPFD, D_TEMP, LAT and LONG (constant for each month)

    LAT    = LATin   ! 0 !latitude
    LONG   = LONGin  ! 0 !longitude
    D_TEMP = DTEMPin ! Daily average temperature K - currently 0K
    D_PPFD = D_PPFDin ! 100.0 ! Daily average PAR

    !--------INPUT OTHER VARIABLES----------------

    DI     = 0.0       ! something related to DIstomata, value 0 came with model
    SMOIS  = SMOISTin  ! Soil moisture
    ISLTYP = 6         ! Soiltype. This is only important for isoprene drought response algorithm.
                       ! For soil types, see table 2 in Chen, F. & J. Dudhia
                       ! (2001) Coupling an advanced land surface-hydrology
                       ! model with the Penn State-NCAR MM5 modeling system.
                       ! Part I: Model implementation and sensitivity. Monthly
                       ! Weather Review, 129, 569-585. 
 
    midpoints= centralizer(Zin, Layersin)

    !rescaling to [0,1]
    midpoints= midpoints/(Zin(Layersin))

    !BIOMASS. Or FOLIARMASS. DRY WEIGT. Unit is g/cm^2
    IF (Treeflag == 1) THEN !pine
       ! BIOMASS = 0.05043    !50 y, HENVI
       ! BIOMASS = 0.05161  ! 20 y
       ! BIOMASS = 0.01939  !15 y, HENVI
       BIOMASS = 0.0538   ! usually used in SOSA-this can be deleted if no one is using this value for other stations
       IF (TRIM(ADJUSTL(station)) == 'hyytiala') THEN
          !The values from 1995-2011 are from forest inventories by Janne
          !Levula. Values for later years are artificial values by Ditte
          !following trends. This biomass is the sum of the foliage mass from
          !both pine and spruce. The reason for the decrease in 2010 is due to
          !snow damage on the trees.
          IF (year_int .EQ. 2018) BIOMASS = 0.0599
          IF (year_int .EQ. 2017) BIOMASS = 0.0588
          IF (year_int .EQ. 2016) BIOMASS = 0.0577
          IF (year_int .EQ. 2015) BIOMASS = 0.0566
          IF (year_int .EQ. 2014) BIOMASS = 0.0555
          IF (year_int .EQ. 2013) BIOMASS = 0.0544
          IF (year_int .EQ. 2012) BIOMASS = 0.0533 
          IF (year_int .EQ. 2011) BIOMASS = 0.0523
          IF (year_int .EQ. 2010) BIOMASS = 0.0509
          IF (year_int .EQ. 2009) BIOMASS = 0.0519
          IF (year_int .EQ. 2008) BIOMASS = 0.0519
          IF (year_int .EQ. 2007) BIOMASS = 0.0509
          IF (year_int .EQ. 2006) BIOMASS = 0.0498
          IF (year_int .EQ. 2005) BIOMASS = 0.0489
          IF (year_int .EQ. 2004) BIOMASS = 0.0477
          IF (year_int .EQ. 2003) BIOMASS = 0.0463 
       ENDIF
    ELSEIF (Treeflag == 2) THEN !spruce 
       !BIOMASS = 0.12225 !50 y, HENVI
       ! BIOMASS = 0.13344 ! 30 y
        BIOMASS = 0.06798 ! 15 y, HENVI 
    ELSEIF (Treeflag == 3) THEN !birch - this is the maximum biomass throughout the growing season
      ! BIOMASS_MAX = 0.01957 !50 y, HENVI
       ! BIOMASS_MAX = 0.02865 ! 20 y
        BIOMASS_MAX = 0.00616 ! 10y, HENVI
    ENDIF
 
    !--------INPUT standard emission potential

    !   LOGICAL, SAVE :: first_call = .TRUE.

    IF (first_call) THEN ! do only once
       first_call = .FALSE.

       if (TRIM(ADJUSTL(station)) .eq. 'hyytiala') then
          IF (Treeflag == 1) THEN ! pine
             ! Input of the standard emissions potentials from Jaana Baeck
             OPEN(unit=5,file=TRIM(ADJUSTL(input_dir_station_info))//'/EF_day.txt', status='old')
          ELSEIF (Treeflag == 2) THEN !spruce
             OPEN(unit=5,file=TRIM(ADJUSTL(input_dir_station_info))//'/EF_spruce.txt', status='old')
          ELSEIF (Treeflag == 3) THEN !birch
             OPEN(unit=5,file=TRIM(ADJUSTL(input_dir_station_info))//'/EF_birch.txt', status='old')
          ENDIF
      elseif (station .eq. 'Manitou') then
          OPEN(unit=5,file=TRIM(ADJUSTL(input_dir_station_info))//'/EF_day.txt')
      endif
       

       DO J = 1,366
          READ(5,*) (EFin(J,I), I = 1,22)
       ENDDO
       CLOSE(5)
       If (TRIM(ADJUSTL(station)) == 'hyytiala') THEN
          ! Input of the Canopy characteristics for NrTyp with (Manitou shares the same file with Hyytiala)
          !For canopy.txt file, the different columns contain the following:
          !C1: pine, C2: mixed forest, C3: broadleaf forest, C4: mixed vegetation
          !C5: shrubs, C6: grass, C7: crop, C8: test, C9: spruce
          !For canopy.txt file, the different columns contain the following:
          !C1: Scots pine, 15 year, but I should change the other stuff as well.
          open(4,file=''//TRIM(ADJUSTL(input_dir_station_info))//'/canopy.txt')  
          Canopychar = 0.0
          DO J = 1,NrCha
             READ(4,*) (Canopychar(J,I), I = 1, NrTyp)
          ENDDO
          CLOSE(4)
          !The canopy height and depth is given for pine and spruce weighted by their basal area.
          !Values from 1995-2011 are from Janne Levula. The later years are based on
          !trends. Unit: m
          !CAN_HEI = Canopychar(4,1)
          IF (year_int .EQ. 2014) CAN_HEI = 19.32
          IF (year_int .EQ. 2013) CAN_HEI = 18.97
          IF (year_int .EQ. 2012) CAN_HEI = 18.63
          IF (year_int .EQ. 2011) CAN_HEI = 18.29
          IF (year_int .EQ. 2010) CAN_HEI = 17.98
          IF (year_int .EQ. 2009) CAN_HEI = 17.58
          IF (year_int .EQ. 2008) CAN_HEI = 17.16
          IF (year_int .EQ. 2007) CAN_HEI = 16.86
          IF (year_int .EQ. 2006) CAN_HEI = 16.57
          IF (year_int .EQ. 2005) CAN_HEI = 16.27
          IF (year_int .EQ. 2004) CAN_HEI = 15.97
          IF (year_int .EQ. 2003) CAN_HEI = 15.68
          Canopychar(4,1) = CAN_HEI
          IF (year_int .EQ. 2014) CAN_DEP = 8.63
          IF (year_int .EQ. 2013) CAN_DEP = 8.49
          IF (year_int .EQ. 2012) CAN_DEP = 8.36 
          IF (year_int .EQ. 2011) CAN_DEP = 8.22 
          IF (year_int .EQ. 2010) CAN_DEP = 8.22
          IF (year_int .EQ. 2009) CAN_DEP = 8.16
          IF (year_int .EQ. 2008) CAN_DEP = 8.03
          IF (year_int .EQ. 2007) CAN_DEP = 8.13
          IF (year_int .EQ. 2006) CAN_DEP = 8.24
          IF (year_int .EQ. 2005) CAN_DEP = 8.06
          IF (year_int .EQ. 2004) CAN_DEP = 7.88
          IF (year_int .EQ. 2003) CAN_DEP = 7.70
          Canopychar(1,1) = CAN_DEP
         ! write(*,*) 'this is the canopy(4,1) value in megan: ', Canopychar(4,1)
       ENDIF  

       IF (TRIM(ADJUSTL(station)) .eq. 'hyytiala') then
          IF (Treeflag == 3) THEN !birch
             !This is used for the tree emission of birch
             !C1: the daily average measured temperature at 4.2 m in unit C.
             !C2: if the temperature in C1 is above 5C, then T-5 is written in C2. unit is C.
             !C3: cumulative sum of daily averaged temperature above 5C.
             !open(2710,file=''//TRIM(ADJUSTL(input_dir_station_info)) // '/TDAY2010.txt')
             open(2710,file=''//TRIM(ADJUSTL(input_dir_station_info)) // '/TDAYCLIMATE.txt')
             do J = 1,366
                read(2710,*) (TDAY(J,I), I=1,3)
             enddo
             close(2710)
             open(2711,file=''//TRIM(ADJUSTL(input_dir_station_info)) // '/biomass_function.txt')
             do J = 1,366
                read(2711,*) BIOMASS_func(J)
             enddo
             close(2711)
          ENDIF
       ENDIF


    ENDIF !(we only read in once)
  
    IDATE = SDATE
    ITIME = STIME
  
    !..Initialize hourly variables
    ! TEMP(:) = TEMPin(1:(Layersin+1))!280 ! hourly temperature K

    !TEMP(1:Layersin)  = TEMpin(Layersin : 1 : -1)
    TEMP(1:Layersin)  = TEMpin(1:Layersin)
    TEMP( Layersin+1) = TEMpin(Layersin+1)

    !  dummyTEMP(1,1,:) = TEMPin(1:Layersin) !because calculation of Gamma_tmp requires a dimension(:,:,:) array

    dummyTEMP(1, 1 , 1 : Layersin) = TEMPin( Layersin : 1 : -1)
    PPFD     =  PPFDin                   ! hourly PAR
    LAIp     =  LAIpin                   ! Previous month or timestep LAI
    LAIc     =  LAIcin                   ! Current month or timestep LAI

    Wind(:)  =  WINDin( 1: (Layersin+1)) ! m/s

    Wind(1:Layersin) = Windin( Layersin : 1 : -1)
    Wind(Layersin+1) = Windin(Layersin+1)

    Pres =  PRESin !1013E3 ! Pa? probably
    Humidity(:) = HUMin(1: (Layersin+1)) ! [g/kg]

    !Humidity(1:Layersin) = HUMin( Layersin : 1 : -1)
    Humidity(1:Layersin) = HUMin(1:Layersin)
    Humidity(Layersin+1) = HUMin(Layersin+1)

    ! LADp = LADpin(1:Layersin)
    LADp= LADpin( LAyersin : 1 : -1 )

    ! leaf temperature from scadis as input 
    Sunleaftk_in(1:Layersin)  = Sunleaftk_in(Layersin : 1 : -1)
    Shadeleaftk_in(1:Layersin)  = Shadeleaftk_in(Layersin : 1 : -1)
    
    !..Go over all the chemical species
    DO S = 1, NEMIS

       !..Initialize variables
       ER_BT = 0.
       ER_NT = 0.
       ER_SB = 0. 
       ER_HB = 0.
       GAM_PHO = 0.
       GAM_TMP = 0.
       GAM_LHT = 0.
       GAM_AGE = 0.
       GAM_SMT = 0.
       GAM_T = 0.
       GAM_CE = 0.

       Gammafactors = -1.0

       CALL  GAMME_CE(IDATE,ITIME,NCOLS,NROWS,BETA,LAT, &
            LONG,PPFD,TEMP,LAIc,Wind,Pres,Humidity, &
            DI, PFTF,Canopychar,NrCha,NrTyp,NPFT, &
            GAM_T,GAM_PHO,GAM_CE, LADp, midpoints, &
            Layersin, GAMMAfactors, sun_par, &
            Sunleaftk, Shadeleaftk, sunfrac,check_rh, &
            station, Sunleaftk_in, Shadeleaftk_in)
            
          
       CALL GAMMA_TNISP(NCOLS, NROWS, VARIABLENAMES(S), TEMP, Layersin, GAM_TMP)

       CALL GAMMA_LAI(NCOLS, NROWS, LAIc, GAM_LHT)

       IF(Treeflag == 3) THEN !birch
          CALL GAMMA_A(NCOLS, NROWS, VARIABLENAMES(S), LAIp, LAIc, TSTLEN, D_TEMP, GAM_AGE) !correction factor for age of leaves [0 1] - canopy specific
       ELSE
          DO I=1,NCOLS
             DO J=1,NROWS
                GAM_AGE(I,J) = 1.0
             ENDDO
          ENDDO
       ENDIF

       CALL  GAMMA_S(NCOLS, NROWS, SMOIS, ISLTYP, VARIABLENAMES(S), GAM_SMT) !correction factor for soil moisture [0 1] - canopy specific, but only relevant for isoprene

       DO VAR = 1, NEMIS !Looping over all emitted compounds
          IF((TRIM(VARIABLENAMES(S))) .EQ. (TRIM(LDF_SPC(VAR)))) THEN
             LDF= LDF_FCT(VAR)
          END IF
       END DO

       !Same thing but separately for every layer.

       SELECT CASE ( VARIABLENAMES(S) )

       CASE ('ISOP')

          DO I=1,NCOLS
             DO J=1,NROWS
                DO lay=1, Layersin
                      ER_BT_out(S,lay) = (GAMMAfactors(3,1,lay) * GAM_AGE(I,J) * GAM_SMT(I,J))
                      ER_NT_out(S,lay) = (GAMMAfactors(3,2,lay) * GAM_AGE(I,J) * GAM_SMT(I,J))
                      ER_SB_out(S,lay) = (GAMMAfactors(3,3,lay) * GAM_AGE(I,J) * GAM_SMT(I,J))
                      ER_HB_out(S,lay) = (GAMMAfactors(3,4,lay) * GAM_AGE(I,J) * GAM_SMT(I,J))
                ENDDO
             ENDDO
          ENDDO

       ! write(*,*) 'GAM_AGE for isoprene: ', GAM_AGE(1,1)


       CASE ('MBO','MYRC','SABI','LIMO','3CAR','OCIM','BPIN', 'APIN','FARN','BCAR','MEOH','ACTO','ACTA','FORM', &
            'CH4',  'NO','OMTP','OSQT',  'CO','CINE','LINA')


          DO I=1,NCOLS
             DO J=1,NROWS
                Do lay=1, Layersin
                      ER_BT_out(S,lay) = (GAM_TMP(I,J,lay) * GAM_AGE(I,J) * GAM_LHT(I,J) * GAM_SMT(I,J)) * &
                           ((1-LDF) + (GAMMAfactors(2,1,lay)*LDF))

                      ER_NT_out(S,lay) = (GAM_TMP(I,J,lay) * GAM_AGE(I,J) * GAM_LHT(I,J) * GAM_SMT(I,J)) * &
                           ((1-LDF) + (GAMMAfactors(2,2,lay)*LDF))

                      ER_SB_out(S,lay) = (GAM_TMP(I,J,lay) * GAM_AGE(I,J) * GAM_LHT(I,J) * GAM_SMT(I,J)) * &
                           ((1-LDF) + (GAMMAfactors(2,3,lay)*LDF))

                      ER_HB_out(S,lay) = (GAM_TMP(I,J,lay) * GAM_AGE(I,J) * GAM_LHT(I,J) * GAM_SMT(I,J)) * &
                           ((1-LDF) + (GAMMAfactors(2,4,lay)*LDF))
                ENDDO
             ENDDO
          ENDDO

        !write(*,*) 'GAM_AGE for other compounds: ', GAM_AGE(1,1)

       CASE DEFAULT
          WRITE(*,*) 'Error: Chemical species, invalid variable: Luxi'
          WRITE(*,*) TRIM(VARIABLENAMES(S))
          STOP

       ENDSELECT

       !----------------------------------------------------------------------
       !....3) Write out the calculated ER and met data
       !----------------------------------------------------------------------
       ! Write emission to file

       IF (S <= N_MGN_SPC) THEN
          VARout(S) = VARIABLENAMES( S)
       END IF





       IF (S .EQ. 9) THEN
       
          DO lay = 1, 4 !looping over plant functional types
             DO J = 1,layersin
                GAM_T_out(lay, J)     = GAMMAfactors(1, lay, (layersin-J+1)) !Note GAM_T is not actually used anywhere
                GAM_PHO_out(lay, J)   = GAMMAfactors(2, lay, (layersin-J+1))
                GAM_CE_out(lay, J)    = GAMMAfactors(3, lay, (layersin-J+1))
             ENDDO
          END DO
       ENDIF

       GAM_OTHER_out(S, 1) = GAM_AGE(1,1)
       GAM_OTHER_out(S, 2) = GAM_LHT(1,1)
       GAM_OTHER_out(S, 3) = GAM_SMT(1,1)

       !GAM_TMP_out(S, 1 : Layersin) = GAM_TMP(1,1, Layersin : 1 : -1)
       GAM_TMP_out(S, 1 : Layersin) = GAM_TMP(1,1, 1:Layersin)

    ENDDO ! End loop for emission species (S)


    if (TRIM(ADJUSTL(station)) .eq. 'hyytiala') then 
       ERtemporary = ER_NT_out    
    else   
       ERtemporary = ER_HB_out + ER_NT_out + ER_SB_out + ER_BT_out
    end if
    


    !Erout= Ertemporary

    ERout( : , 1 : Layersin) = Ertemporary( :, Layersin : 1 : -1)
    !Rout = Ertemporary

    if (TRIM(ADJUSTL(station)) .eq. 'hyytiala') then
     
        IF (Treeflag == 1) THEN !pine 
           !Standard emission potential for (sum of) monoterpenes:
           SEP_mono = EMI_MEAS_MONO / EXP(0.09 * (CHAM_temp - 30.0))*5
           !SEP_mono and EMI_MEAS_MONO is in the unit ng g-1 s-1
           !Unit of CHAM_temp is celcius. CHAM_temp is actually the temperature in the chamber, 
           !but we assume that it is the measured leaf temeperature. Max difference should be ~2C.
           !0.09 is an emprical coefficient
           !30 (C) is the leaf temperature at standard conditions
           !This equation is from Tarvainen et al., ACP, 2005.

           !Standard emission potential for the monoterpenes. 
           !I have devided the monoterpene speciation based on Baeck et al., Biogeosciences, 2012.
           ! ng g(ndw)-1 s-1 --> g cm-2 s-1
           API_EMI = 0.437d0 * SEP_mono * 1.0d-9 * BIOMASS
           CAR_EMI = 0.396d0 * SEP_mono * 1.0d-9 * BIOMASS
           BPI_EMI = 0.090d0 * SEP_mono * 1.0d-9 * BIOMASS
           LIM_EMI = 0.023d0 * SEP_mono * 1.0d-9 * BIOMASS
           OMT_EMI = 0.053d0 * SEP_mono * 1.0d-9 * BIOMASS
           MYR_EMI = 0.0d0
           SAB_EMI = 0.0d0
           OCI_EMI = 0.0d0
           CIN_EMI = 0.001d0 * SEP_mono * 1.0d-9 * BIOMASS
           LIN_EMI = 0.0d0

           ! Standard emission potential [ng/g(needledryweight)/h into g/cm2/s] for Hyytiala spring from Hakola et al. Biogeosc., 3, 2006
           ! (with g(needledryweight into cm2 by CarboEuroflux data: Foliage biomass are from Janne Levula (1995-2011) and around 0.5 kg/m2 or 0.05 g/cm2)
           !The SEP for methanol, acetone and acetaldehyde are based on Juho's measurements from year 2009-2011
           ISO_EMI = EFin(Day,1)  * 1.0d-9 / 3600 * BIOMASS*4.0*1.5
           ! MYR_EMI = EFin(Day,2)  * 1E-9 / 3600 * BIOMASS
           ! SAB_EMI = EFin(Day,3)  * 1E-9 / 3600 * BIOMASS
           ! LIM_EMI = EFin(Day,4)  * 1E-9 / 3600 * BIOMASS
           ! CAR_EMI = EFin(Day,5)  * 1E-9 / 3600 * BIOMASS
           ! OCI_EMI = EFin(Day,6)  * 1E-9 / 3600 * BIOMASS
           ! BPI_EMI = EFin(Day,7)  * 1E-9 / 3600 * BIOMASS
           ! API_EMI = EFin(Day,8)  * 1E-9 / 3600 * BIOMASS
           ! OMT_EMI = EFin(Day,9)  * 1E-9 / 3600 * BIOMASS
           FAR_EMI = EFin(Day,10) * 1.0d-9 / 3600 * BIOMASS*3
           BCA_EMI = EFin(Day,11) * 1.0d-9 / 3600 * BIOMASS*3
           OSQ_EMI = EFin(Day,12) * 1.0d-9 / 3600 * BIOMASS*3
           MBO_EMI = EFin(Day,13) * 1.0d-9 / 3600 * BIOMASS
           !  MET_EMI = EFin(Day,14) * 1.0d-9 / 3600 * BIOMASS
           MET_EMI = 75.0d0*1.0d-9*1.0d-4  ! 75 [ng m-2 s-1] (Rantala et al., 2015), then --> [g cm-2 s-1]
           ACT_EMI = EFin(Day,15) * 1.0d-9 / 3600 * BIOMASS*4.5
           CH4_EMI = EFin(Day,16) * 1.0d-9 / 3600 * BIOMASS  ! 0 currently, read from input in Sosa.f90
           NO_EMI  = EFin(Day,17) * 1.0d-9 / 3600 * BIOMASS  ! 0 currently, read from input in Sosa.f90
           ACA_EMI = EFin(Day,18) * 1.0d-9 / 3600 * BIOMASS*6
           ! FOR_EMI = EFin(Day,19) * 1.0d-9 / 3600 * BIOMASS*2
           FOR_EMI = 75.0d0 * 1.0d-9 / 3600 * BIOMASS
           CO_EMI  = EFin(Day,20) * 1.0d-9 / 3600 * BIOMASS  ! 0 currently, read from input in Sosa.f90
           ! CIN_EMI = EFin(Day,21) * 1E-9 / 3600 * BIOMASS
           ! LIN_EMI = EFin(Day,22) * 1E-9 / 3600 * BIOMASS

        ELSEIF (Treeflag == 2) THEN !spruce
           ! HERE COMES THE SPRUCE SPECIFIC SEP
           ISO_EMI = EFin(Day,1)  * 1.0d-6 * BIOMASS * 1/3600
           MYR_EMI = EFin(Day,2)  * 1.0d-6 * BIOMASS * 1/3600
           SAB_EMI = EFin(Day,3)  * 1.0d-6 * BIOMASS * 1/3600
           LIM_EMI = EFin(Day,4)  * 1.0d-6 * BIOMASS * 1/3600
           CAR_EMI = EFin(Day,5)  * 1.0d-6 * BIOMASS * 1/3600
           OCI_EMI = EFin(Day,6)  * 1.0d-6 * BIOMASS * 1/3600
           BPI_EMI = EFin(Day,7)  * 1.0d-6 * BIOMASS * 1/3600
           API_EMI = EFin(Day,8)  * 1.0d-6 * BIOMASS * 1/3600
           OMT_EMI = EFin(Day,9)  * 1.0d-6 * BIOMASS * 1/3600
           FAR_EMI = EFin(Day,10) * 1.0d-6 * BIOMASS * 1/3600
           BCA_EMI = EFin(Day,11) * 1.0d-6 * BIOMASS * 1/3600
           OSQ_EMI = EFin(Day,12) * 1.0d-6 * BIOMASS * 1/3600
           MBO_EMI = EFin(Day,13) * 1.0d-6 * BIOMASS * 1/3600
           MET_EMI = EFin(Day,14) * 1.0d-6 * BIOMASS * 1/3600
           ACT_EMI = EFin(Day,15) * 1.0d-6 * BIOMASS * 1/3600
           CH4_EMI = EFin(Day,16) * 1.0d-6 * BIOMASS * 1/3600
           NO_EMI  = EFin(Day,17) * 1.0d-6 * BIOMASS * 1/3600
           ACA_EMI = EFin(Day,18) * 1.0d-6 * BIOMASS * 1/3600
           FOR_EMI = EFin(Day,19) * 1.0d-6 * BIOMASS * 1/3600
           CO_EMI  = EFin(Day,20) * 1.0d-6 * BIOMASS * 1/3600
           CIN_EMI = EFin(Day,21) * 1.0d-6 * BIOMASS * 1/3600
           LIN_EMI = EFin(Day,22) * 1.0d-6 * BIOMASS * 1/3600

        ELSEIF (Treeflag == 3) THEN !birch
          !HERE COMES THE BIRCH SPECIFIC SEP
          IF (Day .GT. 274) THEN
             SEP_birch = 0.0
          ELSE
             IF (TDAY(Day,3) == 0.0) THEN
                SEP_birch = 0.0
             ELSEIF ((TDAY(Day,3) .GT. 0.0) .AND. (TDAY(Day,3) .LT. 80.0)) THEN
                SEP_birch = 3.63 * 1D-6 * BIOMASS * 1/3600
             ELSEIF ((TDAY(Day,3) .GE. 80.0) .AND. (TDAY(Day,3) .LT. 400.0)) THEN
                SEP_birch = 0.68 * 1D-6 * BIOMASS * 1/3600
             ELSEIF (TDAY(DAY,3) .GE. 400) THEN
                SEP_birch = 7.71 * 1D-6 * BIOMASS * 1/3600
             ENDIF
          ENDIF
          ISO_EMI = EFin(Day,1)  * SEP_birch * 1/100
          MYR_EMI = EFin(Day,2)  * SEP_birch * 1/100
          SAB_EMI = EFin(Day,3)  * SEP_birch * 1/100
          LIM_EMI = EFin(Day,4)  * SEP_birch * 1/100
          CAR_EMI = EFin(Day,5)  * SEP_birch * 1/100
          OCI_EMI = EFin(Day,6)  * SEP_birch * 1/100
          BPI_EMI = EFin(Day,7)  * SEP_birch * 1/100
          API_EMI = EFin(Day,8)  * SEP_birch * 1/100
          OMT_EMI = EFin(Day,9)  * SEP_birch * 1/100
          FAR_EMI = EFin(Day,10) * SEP_birch * 1/100
          BCA_EMI = EFin(Day,11) * SEP_birch * 1/100
          OSQ_EMI = EFin(Day,12) * SEP_birch * 1/100
          MBO_EMI = EFin(Day,13) * SEP_birch * 1/100
          MET_EMI = EFin(Day,14) * SEP_birch * 1/100
          ACT_EMI = EFin(Day,15) * SEP_birch * 1/100
          CH4_EMI = EFin(Day,16) * SEP_birch * 1/100
          NO_EMI  = EFin(Day,17) * SEP_birch * 1/100
          ACA_EMI = EFin(Day,18) * SEP_birch * 1/100
          FOR_EMI = EFin(Day,19) * SEP_birch * 1/100
          CO_EMI  = EFin(Day,20) * SEP_birch * 1/100
          CIN_EMI = EFin(Day,21) * SEP_birch * 1/100
          LIN_EMI = EFin(Day,22) * SEP_birch * 1/100
          !write(*,*) 'API_EMI = ', API_EMI
          !write(*,*) 'Day = ', Day
          !write(*,*) 'SEP_birch = ', SEP_birch
          !write(*,*) 'EFin(apin) = ',EFin(Day,8)
          !write(*,*) 'TDAY = ', TDAY(Day,3)

        ENDIF
       
       elseif (TRIM(ADJUSTL(station)) .eq. 'Manitou') then 


       ISO_EMI = EFin(Day,1)  
       MYR_EMI = EFin(Day,2)  
       SAB_EMI = EFin(Day,3)  
       LIM_EMI = EFin(Day,4)  
       CAR_EMI = EFin(Day,5)  
       OCI_EMI = EFin(Day,6)  
       BPI_EMI = EFin(Day,7)  
       API_EMI = EFin(Day,8)  
       OMT_EMI = EFin(Day,9)  
       FAR_EMI = EFin(Day,10) 
       BCA_EMI = EFin(Day,11) 
       OSQ_EMI = EFin(Day,12) 
       MBO_EMI = EFin(Day,13) 
       MET_EMI = EFin(Day,14) 
       ACT_EMI = EFin(Day,15) 
       CH4_EMI = EFin(Day,16) 
       NO_EMI  = EFin(Day,17) 
       ACA_EMI = EFin(Day,18) 
       FOR_EMI = EFin(Day,19) 
       CO_EMI  = EFin(Day,20) 
       CIN_EMI = EFin(Day,21) 
       LIN_EMI = EFin(Day,22) 
    endif
    
    ! rosa add emission for africa here!

    ! Emission with EF in g/cm2/s with ER into molecules/cm3/s  
    ! (ERout is the correction factor which is calcuated based on radiation, temperature, LAI, bla bla bal.... LAI should be 
    ! projected leaf area index. )
      DO lay = 2,Layersin
      ! rosa: when you have the emissions for africa, check that these correct too!
       EMI(lay,1)  = ERout(1,lay)  * ISO_EMI /  68 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,2) = 0.
       !EMI(lay,2)  = ERout(2,lay)  * MYR_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,3) = 0.
       !EMI(lay,3)  = ERout(3,lay)  * SAB_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       !EMI(lay,4)  = LIM_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,4)  = ERout(4,lay)  * LIM_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
      ! EMI(lay,5)  = CAR_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,5)  = ERout(5,lay)  * CAR_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,6) = 0.
       !EMI(lay,6)  = ERout(6,lay)  * OCI_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       !EMI(lay,7)  = BPI_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,7)  = ERout(7,lay)  * BPI_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       !EMI(lay,8)  = API_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,8)  = ERout(8,lay)  * API_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       !EMI(lay,9)  = OMT_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,9)  = ERout(9,lay)  * OMT_EMI / 136 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,10) = ERout(10,lay) * FAR_EMI / 204 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,11) = ERout(11,lay) * BCA_EMI / 204 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,12) = ERout(12,lay) * OSQ_EMI / 204 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,13) = ERout(13,lay) * MBO_EMI /  86 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,14) = ERout(14,lay) * MET_EMI /  32 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,15) = ERout(15,lay) * ACT_EMI /  58 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,16) = ERout(16,lay) * CH4_EMI /  16 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,17) = ERout(17,lay) * NO_EMI  /  30 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,18) = ERout(18,lay) * ACA_EMI /  44 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,19) = ERout(19,lay) * FOR_EMI /  30 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,20) = ERout(20,lay) * CO_EMI  /  28 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       !EMI(lay,21) = CIN_EMI / 154 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,21) = ERout(21,lay) * CIN_EMI / 154 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
       EMI(lay,22) = 0.
       !EMI(lay,22) = ERout(22,lay) * LIN_EMI / 154 * Avog / (Zin(lay)-Zin(lay-1)) / 100 * Ladpin(lay)
      ENDDO

    IF (Zin(1) .EQ. 0.) THEN

       DO J = 1,22
          EMI(1,J) = 0.
       ENDDO

    else
! rosa: check that correct here too
       EMI(1,1)  = ERout(1,1)  * ISO_EMI /  68 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,2)  = 0.
       !EMI(1,2)  = ERout(2,1)  * MYR_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,3)  = 0.
       !EMI(1,3)  = ERout(3,1)  * SAB_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
       !EMI(1,4)  = LIM_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,4)  = ERout(4,1)  * LIM_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
       !EMI(1,5)  = CAR_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,5)  = ERout(5,1)  * CAR_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,6)  = 0.
       !EMI(1,6)  = ERout(6,1)  * OCI_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
       !EMI(1,7)  = BPI_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,7)  = ERout(7,1)  * BPI_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
       !EMI(1,8)  = API_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,8)  = ERout(8,1)  * API_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
       !EMI(1,9)  = OMT_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,9)  = ERout(9,1)  * OMT_EMI / 136 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,10) = ERout(10,1) * FAR_EMI / 204 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,11) = ERout(11,1) * BCA_EMI / 204 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,12) = ERout(12,1) * OSQ_EMI / 204 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,13) = ERout(13,1) * MBO_EMI /  86 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,14) = ERout(14,1) * MET_EMI /  32 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,15) = ERout(15,1) * ACT_EMI /  58 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,16) = ERout(16,1) * CH4_EMI /  16 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,17) = ERout(17,1) * NO_EMI  /  30 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,18) = ERout(18,1) * ACA_EMI /  44 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,19) = ERout(19,1) * FOR_EMI /  30 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,20) = ERout(20,1) * CO_EMI  /  28 * Avog / Zin(1) / 100 * Ladpin(1)
       !EMI(1,21) = CIN_EMI / 154 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,21) = ERout(21,1) * CIN_EMI / 154 * Avog / Zin(1) / 100 * Ladpin(1)
       EMI(1,22) = 0.
       !EMI(1,22) = ERout(22,1) * LIN_EMI / 154 * Avog / Zin(1) / 100 * Ladpin(1)
    ENDIF

    DO lay = Layersin+1, kzin
       sun_par(lay) = 1.
       DO J = 1,22
          EMI(lay,J) = 0.
       ENDDO
    ENDDO

    DEALLOCATE ( ER_BT      )   ! output emission buffer
    DEALLOCATE ( ER_NT      )   ! output emission buffer
    DEALLOCATE ( ER_SB      )   ! output emission buffer
    DEALLOCATE ( ER_HB      )   ! output emission buffer
    DEALLOCATE ( LAT     )   ! input latitude of grid cell
    DEALLOCATE ( LONG    )   ! input longitude of grid cell
    DEALLOCATE ( LAIp    )   ! previous monthly LAI
    DEALLOCATE ( LAIc    )   ! current monthly LAI
    DEALLOCATE ( TEMP    )   ! input hourly temperature (K)
    DEALLOCATE ( PPFD    )   ! calculated PAR (umol/m2.s)
    DEALLOCATE ( D_PPFD  )   ! daily PAR (umol/m2.s)
    DEALLOCATE ( D_TEMP  )   ! input daily temperature (K) - currently 0K
    DEALLOCATE ( GAM_PHO )   ! light correction factor
    DEALLOCATE ( GAM_TMP )   ! temperature correction factor
    DEALLOCATE ( GAM_LHT )   ! LAI correction factor
    DEALLOCATE ( GAM_AGE )   ! leaf age correction factor
    DEALLOCATE ( GAM_SMT )   ! Soil moilture correction factor
    DEALLOCATE ( GAM_T )
    DEALLOCATE ( Wind )   
    DEALLOCATE ( Pres )   
    DEALLOCATE ( Humidity )
    DEALLOCATE ( DI )   
    DEALLOCATE ( Cantype )
    DEALLOCATE ( GAM_CE )
   DEALLOCATE ( ISLTYP  )
    DEALLOCATE ( SMOIS )
    DEALLOCATE ( PFTF)
    DEALLOCATE ( GAMMAfactors)

    !======================================================================
    !..  FORMAT
    !======================================================================
1000 FORMAT( A )
1010 FORMAT( 40( A, :, I8, :, 1X ) )

    !======================================================================
    !..  End program
    !======================================================================

  END SUBROUTINE EMISSION_M2


  !======================================================================
  !..  Some added functionality, can be moved to a new module
  !======================================================================

  Pure function Timeformatter(input) result(output)

    IMPLICIT NONE

    REAL(kind=dp), INTENT(IN) :: input
    INTEGER :: output

    REAL(kind=dp) :: tmp
    INTEGER :: h, m ,s

    tmp= input
    output=0

    h   = INT(tmp/3600) ! truncation of values should work as intented
    tmp = tmp - (h*3600)
    m   = INT(tmp/60)
    tmp = tmp - (m*60)
    s   = INT(tmp)

    output = (10000*h) + (100*m) + s

  END function timeformatter


  Pure function centralizer(Inputheightarray, Inputlayers) result(output)

    IMPLICIT NONE

    REAL(kind=dp), DIMENSION(:), INTENT(IN) :: Inputheightarray
    INTEGER, INTENT(IN) :: Inputlayers
    REAL(kind=dp), DIMENSION(Inputlayers) :: output

    INTEGER :: ind
    REAL(kind=dp) :: temp


    output= Inputheightarray(1:Inputlayers)


    DO ind=1, Inputlayers
       temp=0

       IF(ind == 1) THEN
          temp= Inputheightarray(ind)
       ELSE
          temp = Inputheightarray(ind) - Inputheightarray(ind-1)
       END IF

       output(ind) = output(ind) - (temp*0.5)

    END DO

  END function centralizer

END MODULE MEGAN_version_2
