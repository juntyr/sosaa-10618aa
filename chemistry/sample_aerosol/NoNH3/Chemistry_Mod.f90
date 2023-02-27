!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   A TEMPLATE CHEMISTRY MODULE
!
!   This module includes:
!       - Conversion of solar irradiance to actinic flux
!       - Calculation of J-values
!       - Calculation of complex reaction rates
!       - Calls for KPP that calculates concentrations of the included compounds
!
!   This version is for MCM version 3.2 (J-values & complex reaction rate)
!
!   See version control for comments on changes
!
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo



MODULE Chemistry_Mod  ! chemistry_mod is the name of the chemistry module in SOSA. In MALTE-BOX the chemistry module is called CHEMISTRY.


  !Include the following second_*.f90 files that are specific for different chemistries:
  USE second_Main

  USE second_Precision, ONLY : dp      ! KPP Numerical type
  USE second_Parameters

  ! Use the variables from second_Global, and set the values for them
  ! NPHOT: Number of photochemical reactions used in KPP
  USE second_Global, ONLY : TEMP, J, NPHOT, M, N2, O2, H2O, RO2

  IMPLICIT NONE

  PRIVATE

  !Public subroutines:
  PUBLIC :: KPP_SetUp, CHEMISTRY, chemistry_init
  PUBLIC :: set_emissions, set_emissions_for_chemistry
  PUBLIC :: set_soil_emissions_for_chemistry_1, set_soil_emissions_for_chemistry_2
  PUBLIC :: set_clearcut_emissions_for_chemistry

  !Public variables:
  PUBLIC :: NKVALUES, NWL


  !Global variables:
  INTEGER, PARAMETER :: NKVALUES = 46   ! Number of rate coefficients used in KPP. Hand-copied from second_Main.f90
                                        ! NKVALUES = 46 if MCMv3.3 or MCMv3.3.1, NKVALUES = 42 if MCMv3.2

  ! File path for the in- and output values:
  CHARACTER(LEN=*), PARAMETER :: &
       filename1 = '../../sosa_in', &
       filename2 = '../../sosa_out'

  ! Number of wavelength used here
  INTEGER, PARAMETER :: NWL = 75

  ! Wavelength
  INTEGER, DIMENSION(NWL) :: WL

  ! Delta lambda (so wl band)
  REAL(dp), PARAMETER :: DL = 5.0d0

  ! ACS: absorption cross section
  ! QY: quantum yield for photolysis
  REAL(kind=dp), DIMENSION(NWL) :: &
    ACSA_o3, ACSB_o3, QY_o1d, QY_o3p, &
    ACS_h2o2, QY_h2o2, &
    ACS_no2, QY_no2, &
    ACS_no3, QY_no3_no, QY_no3_no2, &
    ACS_hono, QY_hono, &
    ACS_hno3, QY_hno3, &
    ACSA_hcho, ACSB_hcho, QY_hcho_h_hco, QY_hcho_h2_co, &
    ACS_ch3cho, QY_ch3cho, &
    ACS_c2h5cho, QY_c2h5cho, &
    ACS_nc3h7cho, QY1_nc3h7cho, QY2_nc3h7cho, &
    ACS_ic3h7cho, QY_ic3h7cho, &
    ACS_macr, QY1_macr, QY2_macr, &
    ACSA_ch3coch3, ACSB_ch3coch3, ACSC_ch3coch3, ACSD_ch3coch3, QY_ch3coch3, &
    ACS_mek, QY_mek, &
    ACS_mvk, QY1_mvk, QY2_mvk, &
    ACS_glyox, QY1_glyox, QY3_glyox, QY2_glyox, &
    ACS_mglyox, QY_mglyox, &
    ACS_biacet, QY_biacet, &
    ACS_ch3ooh, QY_ch3ooh, &
    ACSA_ch3no3, ACSB_ch3no3, QY_ch3no3, &
    ACSA_c2h5no3, ACSB_c2h5no3, QY_c2h5no3, &
    ACS_n_c3h7no3, QY_n_c3h7no3, &
    ACSA_i_c3h7no3, ACSB_i_c3h7no3, QY_i_c3h7no3, &
    ACS_t_c4h9no3, QY_t_c4h9no3, &
    ACS_noa, QY1_noa, QY2_noa, &
    ACS_ho2no2, QY1_ho2no2, QY2_ho2no2, &
    ACSA_n2o5, ACSB_n2o5, QY1_n2o5, QY2_n2o5, &
    ACS_test

  REAL(dp) :: ACS_c2h5no3, ACS_ch3coch3, ACS_ch3no3, ACS_hcho, ACS_i_c3h7no3, ACS_n2o5
  ! Constants
  REAL(kind=dp), PARAMETER :: Avog  = 6.0223d23, & ! Avogrado's constant in molecules / mol
                              RG    = 8.314d0,   & ! Gas constant in J / mol / K
                              M_Air = 28.97d-3,   & ! Average molar mass of air in g / mol
                              M_H2O = 18.0153d-3    ! Molar mass of water in g / mol
  REAL(dp), PARAMETER :: &
    h_planck   = 6.62606957D-34, &  ! [J s], Planck constant
    c_light_nm = 2.99792458D+17     ! [nm s-1], light speed

  ! Short wave radiation distribution
  REAL(dp), DIMENSION(84) :: swr_distribution

  INTEGER :: i


CONTAINS


  SUBROUTINE chemistry_init(INPUT_DIR, STATION)
    CHARACTER(*), INTENT(IN) :: INPUT_DIR
    CHARACTER(*), INTENT(IN) :: STATION
    REAL :: ignore
    CHARACTER(100) :: input_dir_photolysis

    input_dir_photolysis = TRIM(ADJUSTL(INPUT_DIR))//'/general/mcm/photolysis'

    ! Input of absorption cross spectrum and quantum yields for different molecules:

    ! Open data files
    OPEN(900,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/o3_mcm_v2.dat',          STATUS = 'OLD')
    OPEN(901,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/o1d_mcm_v2.dat',         STATUS = 'OLD')
    OPEN(902,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/o3p_mcm_v2.dat',         STATUS = 'OLD')
    OPEN(903,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/h2o2_mcm_v2.dat',        STATUS = 'OLD')
    OPEN(904,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/no2_mcm_v2.dat',         STATUS = 'OLD')
    OPEN(905,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/no2_qy_mcm_v2.dat',      STATUS = 'OLD')
    OPEN(906,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/no3_mcm_v2.dat',         STATUS = 'OLD')
    OPEN(907,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/no3_no_mcm_v2.dat',      STATUS = 'OLD')
    OPEN(908,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/no3_no2_mcm_v2.dat',     STATUS = 'OLD')
    OPEN(909,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/hono_mcm_v2.dat',        STATUS = 'OLD')
    OPEN(910,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/hno3_mcm_v2_298.dat',    STATUS = 'OLD')
    OPEN(911,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/hcho_mcm_v2.dat',        STATUS = 'OLD')
    OPEN(912,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/hcho_h_hco_mcm_v2.dat',  STATUS = 'OLD')
    OPEN(913,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/hcho_h2_co_mcm_v2.dat',  STATUS = 'OLD')
    OPEN(914,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/ch3cho_mcm_v2.dat',      STATUS = 'OLD')
    OPEN(915,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/ch3cho_qy_mcm_v2.dat',   STATUS = 'OLD')
    OPEN(916,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/c2h5cho_mcm_v2.dat',     STATUS = 'OLD')
    OPEN(917,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/c2h5cho_qy_mcm_v2.dat',  STATUS = 'OLD')
    OPEN(918,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/nc3h7cho_mcm_v2.dat',    STATUS = 'OLD')
    OPEN(919,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/ic3h7cho_mcm_v2.dat',    STATUS = 'OLD')
    OPEN(920,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/ic3h7cho_qy_mcm_v2.dat', STATUS = 'OLD')
    OPEN(921,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/macr_mcm_v2.dat',        STATUS = 'OLD')
    OPEN(922,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/ch3coch3_mcm_v2.dat',    STATUS = 'OLD')
    OPEN(923,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/ch3coch3_qy_mcm_v2.dat', STATUS = 'OLD')
    OPEN(924,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/mek_mcm_v2.dat',         STATUS = 'OLD')
    OPEN(925,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/mvk_mcm_v2.dat',         STATUS = 'OLD')
    OPEN(926,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/mvk_qy_mcm_v2.dat',      STATUS = 'OLD')
    OPEN(927,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/glyox_mcm_v2.dat',       STATUS = 'OLD')
    OPEN(928,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/glyox_qy_mcm_v2.dat',    STATUS = 'OLD')
    OPEN(929,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/mglyox_mcm_v2.dat',      STATUS = 'OLD')
    OPEN(930,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/mglyox_qy_mcm_v2.dat',   STATUS = 'OLD')
    OPEN(931,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/biacet_mcm_v2.dat',      STATUS = 'OLD')
    OPEN(932,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/ch3ooh_mcm_v2.dat',      STATUS = 'OLD')
    OPEN(933,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/ch3no3_mcm_v2.dat',      STATUS = 'OLD')
    OPEN(934,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/c2h5no3_mcm_v2.dat',     STATUS = 'OLD')
    OPEN(935,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/n_c3h7no3_mcm_v2.dat',   STATUS = 'OLD')
    OPEN(936,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/i_c3h7no3_mcm_v2.dat',   STATUS = 'OLD')
    OPEN(937,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/tc4h9no3_mcm_v2.dat',    STATUS = 'OLD')
    OPEN(938,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/noa_mcm_v2.dat',         STATUS = 'OLD')
    OPEN(939,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/ho2no2_aitkinson.dat',   STATUS = 'OLD')
    OPEN(940,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/n2o5_aitkinson.dat',     STATUS = 'OLD')
    OPEN(941,  FILE = TRIM(ADJUSTL(input_dir_photolysis))//'/test_O3.dat',            STATUS = 'OLD')

    ! Read photolysis data
    DO I=1,NWL
       READ(900,*)  WL(I), ACSA_o3(I), ACSB_o3(I)
       READ(901,*)  WL(I), QY_o1d(I)
       READ(902,*)  WL(I), QY_o3p(I)
       READ(903,*)  WL(I), ACS_h2o2(I), QY_h2o2(I)
       READ(904,*)  WL(I), ACS_no2(I)
       READ(905,*)  WL(I), QY_no2(I)
       READ(906,*)  WL(I), ACS_no3(I)
       READ(907,*)  WL(I), QY_no3_no(I)
       READ(908,*)  WL(I), QY_no3_no2(I)
       READ(909,*)  WL(I), ACS_hono(I), QY_hono(I)
       READ(910,*)  WL(I), ACS_hno3(I), QY_hno3(I)
       READ(911,*)  WL(I), ACSA_hcho(I), ACSB_hcho(I)
       READ(912,*)  WL(I), QY_hcho_h_hco(I)
       READ(913,*)  WL(I), QY_hcho_h2_co(I)
       READ(914,*)  WL(I), ACS_ch3cho(I)
       READ(915,*)  WL(I), QY_ch3cho(I)
       READ(916,*)  WL(I), ACS_c2h5cho(I)
       READ(917,*)  WL(I), QY_c2h5cho(I)
       READ(918,*)  WL(I), ACS_nc3h7cho(I), QY1_nc3h7cho(I), QY2_nc3h7cho(I)
       READ(919,*)  WL(I), ACS_ic3h7cho(I)
       READ(920,*)  WL(I), QY_ic3h7cho(I)
       READ(921,*)  WL(I), ACS_macr(I), QY1_macr(I), QY2_macr(I)
       READ(922,*)  WL(I), ACSA_ch3coch3(I), ACSB_ch3coch3(I), ACSC_ch3coch3(I), ACSD_ch3coch3(I)
       READ(923,*)  WL(I), QY_ch3coch3(I)
       READ(924,*)  WL(I), ACS_mek(I), QY_mek(I)
       READ(925,*)  WL(I), ACS_mvk(I)
       READ(926,*)  WL(I), QY1_mvk(I), QY2_mvk(I)
       READ(927,*)  WL(I), ACS_glyox(I)
       READ(928,*)  WL(I), QY1_glyox(I), QY3_glyox(I), QY2_glyox(I)
       READ(929,*)  WL(I), ACS_mglyox(I)
       READ(930,*)  WL(I), QY_mglyox(I)
       READ(931,*)  WL(I), ACS_biacet(I), QY_biacet(I)
       READ(932,*)  WL(I), ACS_ch3ooh(I), QY_ch3ooh(I)
       READ(933,*)  WL(I), ACSA_ch3no3(I), ACSB_ch3no3(I), QY_ch3no3(I)
       READ(934,*)  WL(I), ACSA_c2h5no3(I), ACSB_c2h5no3(I), QY_c2h5no3(I)
       READ(935,*)  WL(I), ACS_n_c3h7no3(I), QY_n_c3h7no3(I)
       READ(936,*)  WL(I), ACSA_i_c3h7no3(I), ACSB_i_c3h7no3(I), QY_i_c3h7no3(I)
       READ(937,*)  WL(I), ACS_t_c4h9no3(I), QY_t_c4h9no3(I)
       READ(938,*)  WL(I), ACS_noa(I), QY1_noa(I), QY2_noa(I)
       READ(939,*)  WL(I), ACS_ho2no2(I), QY1_ho2no2(I), QY2_ho2no2(I)
       READ(940,*)  WL(I), ACSA_n2o5(I), ACSB_n2o5(I), QY1_n2o5(I), QY2_n2o5(I)
       READ(941,*)  WL(I), ACS_test(I)
    END DO

    ! Close all the files
    DO I = 900, 941
      CLOSE(I)
    END DO

    ! Read short wave radiation distribution
    IF (TRIM(ADJUSTL(STATION)) == 'hyytiala') THEN
      OPEN(900, FILE=TRIM(ADJUSTL(INPUT_DIR))//'/station/hyytiala/info/swr_distribution.txt', STATUS='old')
      READ(900, *) (swr_distribution(I), I=1,84)
      CLOSE(900)

  ELSE IF (TRIM(ADJUSTL(STATION)) == 'traj') THEN
      OPEN(900, FILE=TRIM(ADJUSTL(INPUT_DIR))//'/general/mcm/photolysis/glob_swr_distr.txt', STATUS='old')
      do i=1,84
        READ(900, *) ignore, swr_distribution(I)
      end do
      CLOSE(900)
    END IF
  END SUBROUTINE chemistry_init


  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   CHEMIRTY
  !
  !   Calculates the chemistry for each layer by using KPP
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


  ! SUBROUTINE CHEMISTRY(CONS, time1_in, time2_in, temp_in, pres_in, zenith_deg_in, SUN_PAR_LAYER, &
  !      cs_H2SO4_in, cs_HNO3_in, q_H2O, H2O, RO2_out, J_values, K_values, chem_glo, INPUT_DIR, STATION, Aground, I_ACF_out)
  ! SUBROUTINE CHEMISTRY(CONS, &
  !   time1_in, time2_in, &
  !   temp_in, pres_in, zenith_deg_in, radiation_in, albedo_in, &
  !   cs_H2SO4_in, cs_HNO3_in, M_in, N2_in, O2_in, H2O_in, &
  !   RO2_out, J_values_out, K_values_out, I_ACF_out)

  SUBROUTINE CHEMISTRY(CONS, &
    time1_in, time2_in, &
    temp_in, pres_in, zenith_deg_in, radiation_in, albedo_in, &
    cs_H2SO4_in, cs_HNO3_in, &
    M_in, N2_in, O2_in, H2O_in, &
    RO2_out, J_values_out, K_values_out, I_ACF_out)

    REAL(kind=dp), INTENT(inout)                    :: CONS(NSPEC)       ! Gas concentrations

    REAL(kind=dp), INTENT(in)                       :: time1_in,      &  ! start of time step in KPP
                                                       time2_in,      &  ! end of time step in KPP
                                                       temp_in,       &  ! [K], temperature at this level of atmosphere
                                                       pres_in,       &  ! [Pa], pressure at this level of atmosphere
                                                       zenith_deg_in, &  ! [deg], solar zenith angle
                                                       radiation_in,  &  ! [-], incoming radiation already considering penetration trough the canopy
                                                       albedo_in,     &  ! ground albedo from measurements (not sure for what wl)
                                                       cs_H2SO4_in,   &  ! condensation sink of H2SO4
                                                       cs_HNO3_in,    &  ! condensation sink of HNO3
                                                       M_in,          &  ! [molec cm-3], number concentration of M (air)
                                                       N2_in,         &  ! [molec cm-3], number concentration of N2
                                                       O2_in,         &  ! [molec cm-3], number concentration of O2
                                                       H2O_in            ! [molec cm-3], number concentration of water vapor

    REAL(kind=dp)                                   :: zenith_deg,    &  ! [deg], solar zenith angle in degree
                                                       RH                ! Relative humidity
    REAL(kind=dp), INTENT(out)                      :: RO2_out           ! Concentration of peroxy radical

    REAL(kind=dp), INTENT(out), DIMENSION(NPHOT)    :: J_values_out      ! J_values, so 'rate constants' for photolysis reactions
    REAL(kind=dp), INTENT(out), DIMENSION(NKVALUES) :: K_values_out      ! K_values; complex rate constants
    REAL(kind=dp), INTENT(out), DIMENSION(NWL)      :: I_ACF_out         ! Calculated actinic flux. Unit: photon * cm^-2 * s^-1 * nm^-1

    REAL(kind=dp), DIMENSION(NWL)                   :: I_ACF             ! [W m-2], Calculated actinic flux. Different unit than I_ACF_out

!----------------------------------------------------------------------------------------------------------------------

    ! Calculation of Air, H2O, O2, N2, and in molecule / cm3
    ! M = pres_in / (temp_in * RG) / M_Air * Avog * 1.0d-6
    ! N2 = M * 0.78d0
    ! O2 = M * 0.21d0

    ! Set the values for the variables in second_Global
    M  = M_in
    N2 = N2_in
    O2 = O2_in
    H2O = H2O_in

    TEMP = temp_in

    zenith_deg = zenith_deg_in

    ! Calculated actinic flux from spectral data measured in Hyytiala
    IF (zenith_deg .eq. 0.) THEN
      zenith_deg = 90.
    ENDIF

    ! This considers the radiation reflected from below layers
    I_ACF(1:NWL) = radiation_in/(1d0-albedo_in) * swr_distribution(1:NWL) * (1 + 2 * albedo_in * COS(zenith_deg*3.14/180) + albedo_in)

    ! Calculates J-values and set J values declared in second_Global
    CALL PHOTOLYSIS_KPP(zenith_deg, I_ACF, I_ACF_out, temp_in, J_values_out)
    J = J_values_out

    ! IF (zenith_deg >= 90) THEN
    !    DO i = 1, NPHOT
    !       J_values_out(i) = 0.
    !    ENDDO
    ! ENDIF

    ! Some safety check
    DO I = 1,NPHOT
       IF (J_values_out(I) < -1d9 .OR. J_values_out(I) > 1d9 ) then
          WRITE(*,*) 'Note by Sampo: J_value(I) has bad value.'
          WRITE(*,*) 'I, J_value(I) = ', I, J_values_out(I)
          STOP
       ENDIF
    ENDDO

    CALL KPP_Proceed(CONS, time1_in, time2_in)

    ! [future] Set K values for testing, need to change this in future
    K_values_out = 0.0d0

    ! Return RO2 value to main program
    RO2_out = RO2
  END SUBROUTINE CHEMISTRY


  !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   Subroutine for the photodissociation reactions
  !
  !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  SUBROUTINE PHOTOLYSIS_KPP(Beta, SAFN, SAFN_out, T, Jpr)

    IMPLICIT NONE

    REAL(kind=dp), DIMENSION(NWL) :: SAFN, SAFN_out ! Actinic flux in two different units

    REAL(kind=dp) :: T,   &  ! [K], temperature
                     Beta    ! solar zenit angle

    ! Dimension parameters
    REAL(dp), PARAMETER :: VAA=1.0d-20, VAB=1.0d-19, VAC=1.0d-21, VAD=1.0d-24, VAAAA=1.0d-20*1.0d14

    REAL(dp), INTENT(OUT) :: Jpr(NPHOT)  ! photolysis rate


    REAL :: T1, AAAA, BBBB, CCCC, LOOO, QY_o3a, KULMA


    ! Put all J-values to 0 when new time step begins:
    Jpr = 0.0d0

    T1 = T - 230.0d0
    AAAA =  0.332d0 + 2.5650D-4 * T1 + 1.152D-5 * T1**2 + 2.3130D-8 * T1**3
    BBBB = -0.575d0 + 5.5900D-3 * T1 - 1.439D-5 * T1**2 + 3.2700D-8 * T1**3
    CCCC =  0.466d0 + 8.8830D-4 * T1 - 3.546D-5 * T1**2 - 3.5190D-7 * T1**3
    LOOO =  308.2d0 + 4.4871D-2 * T1 + 6.938D-5 * T1**2 - 2.5452D-6 * T1**3

    DO I = 1, NWL
      SAFN_out(I) = SAFN(I) * WL(I)/(h_planck*c_light_nm) / 1.0d4
    ENDDO

    DO I = 1, NWL
      IF (WL(I) .LE. 300) THEN
        QY_o3a = 0.9
      ELSE
        QY_o3a = AAAA * ATAN (BBBB * (REAL(WL(I)) - LOOO)) + CCCC
      ENDIF
      IF (QY_o3a .GT. 0.9) QY_o3a = 0.9
      IF (QY_o3a .LT. 0.)  QY_o3a = 0.
      Jpr(1) = Jpr(1) + ACS_test(I) * QY_o3a * SAFN((WL(I)-280)/5+1) * DL * VAAAA
      IF (Beta .GT. 89.9) THEN
        KULMA = 89.9
        Jpr(2) = 1.23D-3 * EXP(-0.6 / COS(KULMA * 3.14159265 / 180.))
      ELSE
        Jpr(2) = 1.23D-3 * EXP(-0.6 / COS(Beta * 3.14159265 / 180.))
      ENDIF
      IF (Beta .GE. 90.) Jpr(2) = 0.

      ! !     ACS_o3 = ACSA_o3(I) * EXP(ACSB_o3(I)/T)
      ! !     PR1 = PR1 + ACS_o3 * QY_o1d(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light_nm) * 1/10**4
      ! !     PR2 = PR2 + ACS_o3 * QY_o3p(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light_nm) * 1/10**4
      ! if (PR3 < -1.0d3 .or. PR3 > 1.0d3 .and. I>70) then
      !   write(*,*) I, PR3, ACS_h2o2(I), QY_h2o2(I), SAFN( (WL(I)-280)/5+1 ), DL, WL(I), (WL(I)-280)/5+1
      ! endif
      ACS_hcho = (ACSA_hcho(I) * VAC) + (ACSB_hcho(I) * VAD * (T-298))
      ACS_ch3coch3 = ACSA_ch3coch3(I) * (1 + (ACSB_ch3coch3(I) * T) + (ACSC_ch3coch3(I) * (T**2)) + (ACSD_ch3coch3(I) * (T**3)))
      ACS_ch3no3 = ACSA_ch3no3(I) * VAA * EXP(ACSB_ch3no3(I) * 1.0d-3 * (T-298))
      ACS_c2h5no3 = ACSA_c2h5no3(I) * EXP(ACSB_c2h5no3(I) * (T-298))
      ACS_i_c3h7no3 = ACSA_i_c3h7no3(I) * EXP(ACSB_i_c3h7no3(I) * (T-298))
      ACS_n2o5 = ACSA_n2o5(I) * 10**(1000*ACSB_n2o5(I)*((1/T)-(1.0d0/298.0d0)))

      ! In the old method this is multiplied: SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light_nm) * 1/10**4 at the end,
      ! but it is just equal to SAFN_out(I).
      ! So Putian modified to a simpler style.
      Jpr(3)  = Jpr(3) + ACS_h2o2(I) * QY_h2o2(I) * SAFN_out(I)
      Jpr(4)  = Jpr(4) + ACS_no2(I) * QY_no2(I) * SAFN_out(I) * VAA
      Jpr(5)  = Jpr(5) + ACS_no3(I) * VAB * QY_no3_no(I) *SAFN_out(I)
      Jpr(6)  = Jpr(6) + ACS_no3(I) * VAB * QY_no3_no2(I) *SAFN_out(I)
      Jpr(7)  = Jpr(7) + ACS_hono(I) * QY_hono(I) *SAFN_out(I)
      Jpr(8)  = Jpr(8) + ACS_hno3(I) * VAA * QY_hno3(I) * SAFN_out(I)

      Jpr(11) = Jpr(11) + ACS_hcho * QY_hcho_h_hco(I) * SAFN_out(I)
      Jpr(12) = Jpr(12) + ACS_hcho * QY_hcho_h2_co(I) * SAFN_out(I)
      Jpr(13) = Jpr(13) + ACS_ch3cho(I) * QY_ch3cho(I) * SAFN_out(I)
      Jpr(14) = Jpr(14) + ACS_c2h5cho(I) * QY_c2h5cho(I) * SAFN_out(I)
      Jpr(15) = Jpr(15) + ACS_nc3h7cho(I) * VAC * QY1_nc3h7cho(I) * SAFN_out(I)
      Jpr(16) = Jpr(16) + ACS_nc3h7cho(I) * VAC * QY2_nc3h7cho(I) * SAFN_out(I)
      Jpr(17) = Jpr(17) + ACS_ic3h7cho(I) * VAC * QY_ic3h7cho(I) * SAFN_out(I)
      Jpr(18) = Jpr(18) + ACS_macr(I) * QY1_macr(I) * SAFN_out(I)
      Jpr(19) = Jpr(19) + ACS_macr(I) * QY2_macr(I) * SAFN_out(I)
      Jpr(21) = Jpr(21) + ACS_ch3coch3 * QY_ch3coch3(I) * SAFN_out(I)
      Jpr(22) = Jpr(22) + ACS_mek(I) * QY_mek(I) * SAFN_out(I)
      Jpr(23) = Jpr(23) + ACS_mvk(I) * QY1_mvk(I) * SAFN_out(I)
      Jpr(24) = Jpr(24) + ACS_mvk(I) * QY2_mvk(I) * SAFN_out(I)

      Jpr(31) = Jpr(31) + ACS_glyox(I) * VAA * QY1_glyox(I) * SAFN_out(I)
      Jpr(32) = Jpr(32) + ACS_glyox(I) * VAA * QY2_glyox(I) * SAFN_out(I)
      Jpr(33) = Jpr(33) + ACS_glyox(I) * VAA * QY3_glyox(I) * SAFN_out(I)
      Jpr(34) = Jpr(34) + ACS_mglyox(I) * VAA * QY_mglyox(I) * SAFN_out(I)
      Jpr(35) = Jpr(35) + ACS_biacet(I) * VAA * QY_biacet(I) * SAFN_out(I)

      Jpr(41) = Jpr(41) + ACS_ch3ooh(I) * VAA * QY_ch3ooh(I) * SAFN_out(I)

      Jpr(51) = Jpr(51) + ACS_ch3no3 * QY_ch3no3(I) * SAFN_out(I)
      Jpr(52) = Jpr(52) + ACS_c2h5no3 * QY_c2h5no3(I) * SAFN_out(I)
      Jpr(53) = Jpr(53) + ACS_n_c3h7no3(I) * QY_n_c3h7no3(I) * SAFN_out(I)
      Jpr(54) = Jpr(54) + ACS_i_c3h7no3 * QY_i_c3h7no3(I) * SAFN_out(I)
      Jpr(55) = Jpr(55) + ACS_t_c4h9no3(I) * QY_t_c4h9no3(I) * SAFN_out(I)
      Jpr(56) = Jpr(56) + ACS_noa(I) * QY1_noa(I) * SAFN_out(I)
      Jpr(57) = Jpr(57) + ACS_noa(I) * QY2_noa(I) * SAFN_out(I)

      Jpr(67) = Jpr(67) + ACS_ho2no2(I) * QY1_ho2no2(I) * SAFN_out(I) * VAA
      Jpr(68) = Jpr(68) + ACS_ho2no2(I) * QY2_ho2no2(I) * SAFN_out(I) * VAA

      Jpr(70) = Jpr(70) + ACS_n2o5 * QY1_n2o5(I)  * SAFN_out(I) * VAA
      Jpr(71) = Jpr(71) + ACS_n2o5 * QY2_n2o5(I)  * SAFN_out(I) * VAA
    ENDDO

  END SUBROUTINE PHOTOLYSIS_KPP


  subroutine set_emissions(em_emi, emis)
    real(dp), intent(in ) :: em_emi(:, :)
    real(dp), intent(out) :: emis(:, :)

    ! write(*,*) 'em_emi d1: ', size(em_emi, dim=1)
    ! write(*,*) 'em_emi d2: ', size(em_emi, dim=2)
    ! write(*,*) 'emis d1: ', size(emis, dim=1)
    ! write(*,*) 'emis d2: ', size(emis, dim=2)

    ! Set MEGAN emissions

    emis(:,ind_C5H8     ) =  EM_Emi(:,1)   ! C5H8 (Isoprene)
    emis(:,ind_myrcene  ) =  EM_Emi(:,2)   ! Myrcene
    emis(:,ind_sabinene ) =  EM_Emi(:,3)   ! Sabinene
    emis(:,ind_LIMONENE ) =  EM_Emi(:,4)   ! LIMONENE
    emis(:,ind_carene   ) =  EM_Emi(:,5)   ! Carene
    emis(:,ind_ocimene  ) =  EM_Emi(:,6)   ! Ocimene
    emis(:,ind_BPINENE  ) =  EM_Emi(:,7)   ! Bpinene
    emis(:,ind_APINENE  ) =  EM_Emi(:,8)   ! Apinene
    emis(:,ind_OMT      ) =  EM_Emi(:,9)   ! Other monoterpenes
    emis(:,ind_farnesene) =  EM_Emi(:,10)  ! Farnesene
    emis(:,ind_BCARY    ) =  EM_Emi(:,11)  ! BCARY (Beta-Carophyllene)
    emis(:,ind_OSQ      ) =  EM_Emi(:,12)  ! Other sesquiterpenes
    emis(:,ind_MBO      ) =  EM_Emi(:,13)  ! MBO (2methyl-3buten-2ol)
    emis(:,ind_CH3OH    ) =  EM_Emi(:,14)  ! CH3OH (Methanol)
    emis(:,ind_CH3COCH3 ) =  EM_Emi(:,15)  ! CH3COCH3 (Aceton)
    emis(:,ind_CH4      ) =  EM_Emi(:,16)  ! CH4 (Methane)
    emis(:,ind_NO       ) =  EM_Emi(:,17)  ! NO
    emis(:,ind_CH3CHO   ) =  EM_Emi(:,18)  ! CH3CHO (Acetaldehyde)
    emis(:,ind_HCHO     ) =  EM_Emi(:,19)  ! HCHO (Formaldehyde)
    emis(:,ind_CO       ) =  EM_Emi(:,20)  ! CO
    emis(:,ind_cineole  ) =  EM_Emi(:,21)  ! Cineole (Eucalyptol) (not
                                           ! included in the new megan code
                                           ! - used is the same values as sabinene)
    emis(:,ind_linalool ) =  EM_Emi(:,22)  ! Linalool
  end subroutine set_emissions


SUBROUTINE set_emissions_for_chemistry(cons, em_emi, dt_chem)
  REAL(kind=dp), INTENT(INOUT) :: cons(:)    ! (NSPEC) species concentrations
  REAL(kind=dp), INTENT(IN   ) :: em_emi(:)  ! (22) emission at one layer for each compound
  REAl(kind=dp), INTENT(IN   ) :: dt_chem    ! time step for chemistry

  cons(ind_Emi1)  =  em_emi(1)  * dt_chem  ! C5H8 (Isoprene)
  cons(ind_Emi2)  =  em_emi(2)  * dt_chem  ! Myrcene
  cons(ind_Emi3)  =  em_emi(3)  * dt_chem  ! Sabinene
  cons(ind_Emi4)  =  em_emi(4)  * dt_chem  ! LIMONENE
  cons(ind_Emi5)  =  em_emi(5)  * dt_chem  ! Carene
  cons(ind_Emi6)  =  em_emi(6)  * dt_chem  ! Ocimene
  cons(ind_Emi7)  =  em_emi(7)  * dt_chem  ! Bpinene
  cons(ind_Emi8)  =  em_emi(8)  * dt_chem  ! Apinene
  cons(ind_Emi9)  =  em_emi(9)  * dt_chem  ! Other monoterpenes
  cons(ind_Emi10) =  em_emi(10) * dt_chem  ! Farnesene
  cons(ind_Emi11) =  em_emi(11) * dt_chem  ! BCARY (Beta-Carophyllene)
  cons(ind_Emi12) =  em_emi(12) * dt_chem  ! Other sesquiterpenes
  cons(ind_Emi13) =  em_emi(13) * dt_chem  ! MBO (2methyl-3buten-2ol)
  cons(ind_Emi14) =  em_emi(14) * dt_chem  ! CH3OH (Methanol)
  cons(ind_Emi15) =  em_emi(15) * dt_chem  ! CH3COCH3 (Aceton)
  cons(ind_Emi16) =  em_emi(16) * dt_chem  ! CH4 (Methane)
  cons(ind_Emi17) =  em_emi(17) * dt_chem  ! NO
  cons(ind_Emi18) =  em_emi(18) * dt_chem  ! CH3CHO (Acetaldehyde)
  cons(ind_Emi19) =  em_emi(19) * dt_chem  ! HCHO (Formaldehyde)
  cons(ind_Emi20) =  em_emi(20) * dt_chem  ! CO
  cons(ind_Emi21) =  em_emi(21) * dt_chem  ! Cineole (Eucalyptol) (not included in the new megan code - used is the same values as Sabinene)
  cons(ind_Emi22) =  em_emi(22) * dt_chem  ! Linalool
END SUBROUTINE set_emissions_for_chemistry


SUBROUTINE set_soil_emissions_for_chemistry_1(cons, em_soil_emi, dt_chem, dz)
  REAL(kind=dp), INTENT(INOUT) :: cons(:)         ! (NSPEC) species concentrations
  REAL(kind=dp), INTENT(IN   ) :: em_soil_emi(:)  ! (12) soil emission at one layer for each compound
  REAl(kind=dp), INTENT(IN   ) :: dt_chem         ! time step for chemistry
  REAl(kind=dp), INTENT(IN   ) :: dz              ! height interval

  cons(ind_Emi1)  = cons(ind_Emi1)  + (em_soil_emi(4)  * dt_chem) / dz/100  ! C5H8 (Isoprene)
  cons(ind_Emi3)  = cons(ind_Emi3)  + (em_soil_emi(7)  * dt_chem) / dz/100  ! Camphene
  cons(ind_Emi4)  = cons(ind_Emi4)  + (em_soil_emi(11) * dt_chem) / dz/100  ! Limonene
  cons(ind_Emi5)  = cons(ind_Emi5)  + (em_soil_emi(9)  * dt_chem) / dz/100  ! Carene
  cons(ind_Emi7)  = cons(ind_Emi7)  + (em_soil_emi(8)  * dt_chem) / dz/100  ! Bpinene
  cons(ind_Emi8)  = cons(ind_Emi8)  + (em_soil_emi(6)  * dt_chem) / dz/100  ! Apinene
  cons(ind_Emi9)  = cons(ind_Emi9)  + (em_soil_emi(12) * dt_chem) / dz/100  ! OMT
  cons(ind_Emi14) = cons(ind_Emi14) + (em_soil_emi(1)  * dt_chem) / dz/100  ! CH3OH (Methanol)
  cons(ind_Emi15) = cons(ind_Emi15) + (em_soil_emi(3)  * dt_chem) / dz/100  ! CH3COCH3 (Aceton)
  cons(ind_Emi18) = cons(ind_Emi18) + (em_soil_emi(2)  * dt_chem) / dz/100  ! CH3CHO (Acetaldehyde)
  cons(ind_Emi21) = cons(ind_Emi21) + (em_soil_emi(10) * dt_chem) / dz/100  ! Cineole
END SUBROUTINE set_soil_emissions_for_chemistry_1


SUBROUTINE set_soil_emissions_for_chemistry_2(cons, M137_soil, M33_soil, M69_soil, dt_chem)
  REAL(kind=dp), INTENT(INOUT) :: cons(:)         ! (NSPEC) species concentrations
  REAL(kind=dp), INTENT(IN   ) :: M137_soil, M33_soil, M69_soil  ! soil emission
  REAl(kind=dp), INTENT(IN   ) :: dt_chem         ! time step for chemistry

  cons(ind_Emi1 ) = cons(ind_Emi1 ) + (M69_soil  * dt_chem          )  ! C5H8 (Isoprene)
  cons(ind_Emi4 ) = cons(ind_Emi4 ) + (M137_soil * dt_chem * 0.023d0)  ! Limonene
  cons(ind_Emi5 ) = cons(ind_Emi5 ) + (M137_soil * dt_chem * 0.396d0)  ! Carene
  cons(ind_Emi7 ) = cons(ind_Emi7 ) + (M137_soil * dt_chem * 0.090d0)  ! Bpinene
  cons(ind_Emi8 ) = cons(ind_Emi8 ) + (M137_soil * dt_chem * 0.437d0)  ! Apinene
  cons(ind_Emi9 ) = cons(ind_Emi9 ) + (M137_soil * dt_chem * 0.053d0)  ! OMT
  cons(ind_Emi14) = cons(ind_Emi14) + (M33_soil  * dt_chem          )  ! CH3OH (Methanol)
  cons(ind_Emi21) = cons(ind_Emi21) + (M137_soil * dt_chem * 0.001d0)  ! Cineole
END SUBROUTINE set_soil_emissions_for_chemistry_2


SUBROUTINE set_clearcut_emissions_for_chemistry(cons, em_emi, dt_chem)
  REAL(kind=dp), INTENT(INOUT) :: cons(:)    ! (NSPEC) species concentrations
  REAL(kind=dp), INTENT(IN   ) :: em_emi(:)  ! (22) emission at one layer for each compound
  REAl(kind=dp), INTENT(IN   ) :: dt_chem    ! time step for chemistry

  cons(ind_Emi4)  = dt_chem * em_emi(4)
  cons(ind_Emi5)  = dt_chem * em_emi(5)
  cons(ind_Emi7)  = dt_chem * em_emi(7)
  cons(ind_Emi8)  = dt_chem * em_emi(8)
  cons(ind_Emi9)  = dt_chem * em_emi(9)
  cons(ind_Emi21) = dt_chem * em_emi(21)
END SUBROUTINE set_clearcut_emissions_for_chemistry


END MODULE Chemistry_Mod
