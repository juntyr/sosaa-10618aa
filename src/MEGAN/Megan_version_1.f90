!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!                                                                        
!   EMISSION                                                     
!                          
!   by Alex Guenther (transferred into F90 by Michael Boy) - April 2005                                             
!                                                                        
!   Calculates the emission activity for the selected canopytype 
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

MODULE Megan_version_1

  USE constants_mod, only: dp
  
  IMPLICIT NONE

  PUBLIC :: EMISSION_M1

CONTAINS

!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   SUBROUTINE EMIACTIVITY
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!   
!   Input varibles
!
!   Day						Julian day
!   Lat						Latitude
!   Hour					Hour of the day
!   Tc						Temperature [C]
!   PAR 					Incoming photosynthetic active radiation [umol/m2/s1]
!   Wind					Wind speed [m s-1]
!   Humidity				Relative humidity [%]
!   Cantypye				Defines set of canopy characteristics
!   LAI						Leaf area index [m2 per m2 ground area]
!   DI						???
!   Pres					Pressure [Pa]
!
!
!   Used variables:
!
!   PPFDfrac				Fraction of total solar radiation that is PPFD     
!   Solar					Solar radiation [W/m2]
!   Maxsolar				Maximum of solar radiation
!   Beta					Sin of solar angle above horizon
!   Sinbeta					Solar angle above horizon
!	TairK					Array of canopy air temperature [K]
!   Ws						Array of canopy wind speed [m/s]
!	HumidairPa				Array of canopy ambient humidity in [Pa]
!   StomataDI				An index indicating water status of leaves. used to modify stomatal conductance
!   Transmis				Transmission of PPFD that is diffuse
!   Difffrac				Fraction of PPFD that is diffuse
!   PPFDfrac				Fraction of solar rad that is PPFD
!   Trate					Stability of boundary ???
!   SH						Sensible heat flux ???
!	VPgausWt				Array of gaussian weighting factors 
!   VPgausDis 				Array of gaussian weighting factors
!   VPslwWT					Array of gaussian weighting factors					
!   SunFrac					Array of the fraction of sun leaves. i = 1 is the top canopy layer, 2 is the next layer, etc.
!   SunPPFD					Array of incoming (NOT absorbed) PPFD on a sun leaf [umol/m2/s]
!   ShadePPFD				Array of incoming (NOT absorbed) PPFD on a shade leaf [umol/m2/s]
!   SunQv					Array of visible radiation (in and out) fluxes on sun leaves
!   ShadeQv					Array of absorbed visible radiation (in and out) fluxes on shade leaves
!   SunQn					Array of absorbed near IR radiation (in and out) fluxes on sun leaves
!   ShadeQn					Array of absorbed near IR radiation (in and out) fluxes on shade leaves
!   SunleafTK				Array of leaf temperature for sun leaves [K]
!   SunleafSH				Array of sensible heat flux for sun leaves [W/m2]
!   SunleafLH				Array of latent heat flux for sun leaves [W/m2]
!   SunleafIR				Array of infrared flux for sun leaves [W/m2]
!   ShadeleafTK				Array of leaf temperature for shade leaves [K]
!   ShadeleafSH				Array of sensible heat flux for shade leaves [W/m2]
!   ShadeleafLH				Array of latent heat flux for shade leaves [W/m2]
!   ShadeleafIR				Array of infrared flux for shade leaves [W/m2]
!   QbAbsV, QbAbsN			Absorbed direct beam light for visible and near infra red
!   QdAbsV, QdAbsN			Array of absorbed diffuse light for visible and near infra red
!   QsAbsV, QsAbsN			Array of absorbed scattered light for visible and near infra red
!   QBeamV, QBeamN			Above canopy beam (direct) light for visible and near infra red
!   QDiffV, QDiffN			Above canopy diffuse light for visible and near infra red
!	Ea1pLayer				Array of emission activity of light per layer 
!	Ea1tLayer				Array of emission activity of temperature per layer
!	Ea1Layer				Array of companied emission activity
!	Ea1pCanopy				Total emission activity of light 
!   Ea1tCanopy				Total emission activity of temperature
!   Ea1Canopy				Total companied emission activity
!	
!   Calcbeta				Function: Calculation of solar zenith angle
!   WaterVapPres			Function: Convert water mixing ratio (g/kg) to water vapor pressure
!   Stability				Function: Temperature lapse rate
!   Ea1t99					Function: Temperature dependence activity factor for emission type 1
!   Ea1p99 					Function:
!   DIstomata				Function:
!   CalcEccentricity		Function:
!						
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  SUBROUTINE EMISSION_M1(TairK, PAR, Wind, layers, day, HumidairPa, LAI, BETA, time, s2, k_canopy, z,                &
                         Ea1pLayer_new, Ea1tLayer_new, Ea1NLayer_new, EMI, SUN_PAR_NEW, SunleafTK, ShadeleafTK, Sunfrac)

   
    IMPLICIT NONE

    !INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  

    INTEGER, PARAMETER ::        &
         NrTyp    = 7,           &   ! Number of different biosphere typs
         NrCha    = 16,          &   ! Number of canopy characteristics
         Cantype  = 1	             ! Defines set of canopy characteristics

    INTEGER :: II, JJ, KK, Day, kz, layers, k_canopy, time

    REAL(kind=dp), PARAMETER ::                             &
         Td                  = 301,                &           ! Temperature used in Ea1t99
         S0                  = 1367.,              &           ! Solar constant (W/m2)
         Avog                = 6.0221E23,          &           ! Avogadro-number [molecules cm3]
         WaterAirRatio       = 18.016 / 28.97,     &           ! Ratio between water and air molecules
         ConvertWm2toUmolm2s = 4.766                           ! Convert radiation from [W/m2] to [umol/m2/s1] - new megan

    REAL(kind=dp), SAVE, DIMENSION(NrCha,NrTyp) :: Canopychar                     ! Array of canopy characteristics for NrTyp of canopy type
    REAL(kind=dp), SAVE, DIMENSION(1106,5)      :: SEP_janne                      ! Standard emission potential monoterpenes and others for one 12 days by Janne
    REAL(kind=dp), DIMENSION(366,22)            :: EFin                           ! Standard emission potential monoterpenes and others each day

!    REAL(kind=dp), DIMENSION(NrCha,NrTyp) :: Canopychar                     ! Array of canopy characteristics for NrTyp of canopy type
!    REAL(kind=dp), DIMENSION(366,15)      :: SEP                            ! Standard emission potential monoterpenes and others each day
!    REAL(kind=dp), DIMENSION(1106,5)      :: SEP_janne                      ! Standard emission potential monoterpenes and others for one 12 days by Janne

    REAL(kind=dp) :: LAI, emi_factor, emi_par, ISO_EMI, MBO_EMI, MYR_EMI, SAB_EMI, LIM_EMI, CAR_EMI, OCI_EMI, &
         BPI_EMI, API_EMI, FAR_EMI, BCA_EMI, MET_EMI, ACT_EMI, ACA_EMI, FOR_EMI, CH4_EMI, NO_EMI, &
         OMT_EMI, OSQ_EMI, CO_EMI,  LIN_EMI, CIN_EMI

    REAL(kind=dp), DIMENSION(k_canopy) :: VPgausWt, VPgausDis, VPslwWT, QdAbsV, QsAbsV, QdAbsn, QsAbsn, SunQv, &
         ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD, TairK, HumidairPa, Ws, Sunleaftk, SunleafSH, &
         SunleafLH, SunleafIR, Shadeleaftk, ShadeleafSH, ShadeleafLH, ShadeleafIR, s2, Wind, sun_par, &
         Sunfrac

    REAL(kind=dp), DIMENSION(Layers) :: z, Ea1pLayer, Ea1tLayer, Ea1Layer, Emethanol, Ea1pLayer_new, &
         Ea1tLayer_new, Ea1Layer_new, Ea1NLayer_new, Emethanol_new, sun_par_new

    REAL(kind=dp), DIMENSION(Layers,22)   :: EMI
    REAL(kind=dp), DIMENSION(k_canopy)    :: midpoints   !Layer centerpoints for gamme_CE

    REAL(kind=dp) :: Sinbeta, PAR, Solar, Maxsolar, Transmis, Difffrac, PPFDfrac, QbAbsn, Trate, StomataDI, Qbeamv, & 
            Qdiffv, Qbeamn, Qdiffn, QbAbsV, Ea1tCanopy, Ea1pCanopy, Ea1Canopy, &
            SH, DI, LH, BETA, TSTEP

!#ifdef PARALLEL
    CHARACTER(LEN=*), PARAMETER :: &
         filename1 = 'malte_in'
!#else
 !   CHARACTER(LEN=*), PARAMETER :: &
  !       filename1 = '../malte_in'
!#endif

    !--------INPUT standard emission potential
    ! 1  =  Isoprene
    ! 2  =  MYRC
    ! 3  =  Sabinene
    ! 4  =  Limonen
    ! 5  =  3-Carene
    ! 6  =  Ocimene
    ! 7  =  Beta-pinene
    ! 8  =  Alpha-pinene
    ! 9  =  Other monoterpene
    ! 10 = FARN
    ! 11 = Betacarophylene
    ! 12 = Other sesquiterpenes
    ! 13 =  MBO (2methyl-3buten-2ol)
    ! 14 = Methanol
    ! 15 = Aceton
    ! 16 = Methan
    ! 17 = NO
    ! 18 = Acetaldehyde
    ! 19 = Formaldehyde
    ! 20 = CO
    ! 21 = Cineole
    ! 22 = Linalool

!    IF (time .EQ. TSTEP) THEN
       ! Input of the Canopy characteristics for NrTyp with:
       !open(unit=4,file=''//filename1//'/General/Hyytiala/canopy.txt')
       open(unit=4,file= 'sosa_in/General/station/Manitou/canopy.txt')
       DO JJ = 1,16
          READ(4,*) (Canopychar(JJ,II), II = 1, 4)
       ENDDO
       CLOSE(4)

       OPEN(unit=5,file='sosa_in/General/station/Manitou/EF_day.txt')
       DO JJ = 1,366
          READ(5,*) (EFin(JJ,II), II = 1,22)
       ENDDO
       CLOSE(5)
!    ENDIF


    ! Standard emission potential [ng/g(needledryweight)/h into g/cm2/s] for Hyytiala spring from Hakola et al. Biogeosc., 3, 2006
    ! (with g(needledryweight into cm2 by CarboEuroflux data: Standing leaf biomass Hyytiala 0.538 kg/m2 or 0.0538 g/cm2)

    ISO_EMI = EFin(Day,1)  * 1E-9 / 3600 * 0.0538
    MYR_EMI = EFin(Day,2)  * 1E-9 / 3600 * 0.0538
    SAB_EMI = EFin(Day,3)  * 1E-9 / 3600 * 0.0538
    LIM_EMI = EFin(Day,4)  * 1E-9 / 3600 * 0.0538
    CAR_EMI = EFin(Day,5)  * 1E-9 / 3600 * 0.0538
    OCI_EMI = EFin(Day,6)  * 1E-9 / 3600 * 0.0538
    BPI_EMI = EFin(Day,7)  * 1E-9 / 3600 * 0.0538
    API_EMI = EFin(Day,8)  * 1E-9 / 3600 * 0.0538
    OMT_EMI = EFin(Day,9)  * 1E-9 / 3600 * 0.0538
    FAR_EMI = EFin(Day,10) * 1E-9 / 3600 * 0.0538
    BCA_EMI = EFin(Day,11) * 1E-9 / 3600 * 0.0538
    OSQ_EMI = EFin(Day,12) * 1E-9 / 3600 * 0.0538
    MBO_EMI = EFin(Day,13) * 1E-9 / 3600 * 0.0538
    MET_EMI = EFin(Day,14) * 1E-9 / 3600 * 0.0538
    ACT_EMI = EFin(Day,15) * 1E-9 / 3600 * 0.0538
    CH4_EMI = EFin(Day,16) * 1E-9 / 3600 * 0.0538
    NO_EMI  = EFin(Day,17) * 1E-9 / 3600 * 0.0538
    ACA_EMI = EFin(Day,18) * 1E-9 / 3600 * 0.0538
    FOR_EMI = EFin(Day,19) * 1E-9 / 3600 * 0.0538
    CO_EMI  = EFin(Day,20) * 1E-9 / 3600 * 0.0538
    CIN_EMI = EFin(Day,21) * 1E-9 / 3600 * 0.0538
    LIN_EMI = EFin(Day,22) * 1E-9 / 3600 * 0.0538


    Sinbeta   = SIN((90-Beta) / 57.29578)
    IF (Sinbeta .EQ. 1) THEN  
       Sinbeta = 0.
    ENDIF
    StomataDI = DIstomata(DI)                                 
    Maxsolar  = Sinbeta * S0 * CalcEccentricity(Day)
    Solar     = PAR / ConvertWm2toUmolm2s * 2

    Call SolarFractions(0, Solar, Maxsolar, Transmis, Difffrac, PPFDfrac)

    !Call GaussianIntegration(VPgausWt, VPgausDis, k_canopy)              !Gaussian weighting factors
	midpoints = centralizer(z, k_canopy)

	!rescaling to [0,1]
	midpoints= midpoints/(z(k_canopy))

    DO II = 1,k_canopy 
       IF (s2(k_canopy+1-II) .LT. 0.001) THEN
          s2(k_canopy+1-II) = 0.0005
       ENDIF
       Ws(II)        = Wind(k_canopy+1-II)                                                   
       VPgausWt(II)  = s2(k_canopy+1-II)
       !VPgausWt(II)  = s2(II)
       !VPgausDis(II) = z(II) / z(k_canopy)
       !VPgausDis(II) = (II - 0.5) / k_canopy
       VPgausDis(II) = midpoints(II)
    ENDDO

    Qdiffv    = PPFDfrac * Solar * Difffrac
    Qbeamv    = PPFDfrac * Solar * (1 - Difffrac)
    Qdiffn    = (1 - PPFDfrac) * Solar * Difffrac
    Qbeamn    = (1 - PPFDfrac) * Solar * (1 - Difffrac)

    !Call WeightSLW(VPgausDis, VPgausWt, VPslwWT, LAI, k_canopy)


    Call CanopyRad(VPgausDis, Sinbeta, Qbeamv, Qdiffv, Qbeamn, Qdiffn, Cantype, Canopychar,   &
         Sunfrac, QbAbsV, QdAbsV, QsAbsV, QbAbsn, QdAbsn, QsAbsn, SunQv,         &
         ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD, NrCha, NrTyp, LAI, k_canopy)


    Trate = Stability(Canopychar, Cantype, Solar, NrCha, NrTyp)


    Call CanopyEB(Trate, VPgausDis, Canopychar, Cantype, StomataDI, TairK, HumidairPa, Ws, SunPPFD, ShadePPFD, SunQv,   &
         ShadeQv, SunQn, ShadeQn, Sunleaftk, SunleafSH, SunleafLH, SunleafIR, Shadeleaftk, ShadeleafSH,        &
         ShadeleafLH, ShadeleafIR, NrCha, NrTyp, k_canopy)

    Ea1tlayer        = 0.
    Ea1player        = 0.
    Ea1layer         = 0.
    SH                = 0.
    LH                = 0.



    DO II = 1, k_canopy    
       Ea1tLayer(II) = Ea1t99(Sunleaftk(II), Td) * Sunfrac(II) +                                                              &
            Ea1t99(Shadeleaftk(II), Td) * (1 - Sunfrac(II))

       Ea1pLayer(II) = Ea1p99(VPgausDis(II) * LAI, SunPPFD(II)) * Sunfrac(II) +                                               & 
            Ea1p99(VPgausDis(II) * LAI, ShadePPFD(II)) * (1 - Sunfrac(II))

       Ea1Layer(II)  = Ea1t99(Sunleaftk(II), Td) * Ea1p99(VPgausDis(II) * LAI, SunPPFD(II)) * Sunfrac(II) +                   &
            Ea1t99(Shadeleaftk(II), Td) * Ea1p99(VPgausDis(II) * LAI, ShadePPFD(II)) * (1 - Sunfrac(II))

       IF (SunPPFD(1) .GT. 0.) THEN
          sun_par(II) = SunPPFD(II) / SunPPFD(1)
       ELSE
          sun_par(II) = 1.
       ENDIF
       

       ! Emission for methanol from Harley (Biogeosciences 2007) for mature gray pine
       Emethanol(II) = EXP(0.116 * (Sunleaftk(II) - 273.15 - 30)) + EXP(0.116 * (Shadeleaftk(II) - 273.15 - 30))

    ENDDO


    IF (Ea1tLayer(II) .LT. 0.003) THEN
       Ea1tLayer(II) = 0.003
    ENDIF
    IF (Ea1Layer(II) .LT. 0.003) THEN
       Ea1Layer(II) = 0.003
    ENDIF


    DO II = 1,k_canopy
       Ea1tLayer_new(II) = Ea1tLayer(k_canopy+1-II)
       Ea1pLayer_new(II) = Ea1pLayer(k_canopy+1-II)
       Ea1Layer_new(II)  = Ea1Layer(k_canopy+1-II)
       Emethanol_new(II) = Emethanol(k_canopy+1-II)
       sun_par_new(II)   = sun_par(k_canopy+1-II)
       Ea1NLayer_new(II) = 0.8*Ea1tLayer_new(II) + 0.05*Ea1pLayer_new(II)
    ENDDO

    DO II = 2,k_canopy
       ! Emission with SEP in g/cm2/s with Ea1 into molecules/cm3/s
       EMI(II,1)  = Ea1Layer_new(II)  * ISO_EMI /  68 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,2)  = Ea1NLayer_new(II) * MYR_EMI / 136 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,3)  = Ea1NLayer_new(II) * SAB_EMI / 136 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,4)  = Ea1NLayer_new(II) * LIM_EMI / 136 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,5)  = Ea1NLayer_new(II) * CAR_EMI / 136 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,6)  = Ea1NLayer_new(II) * OCI_EMI / 136 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,7)  = Ea1NLayer_new(II) * BPI_EMI / 136 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,8)  = Ea1NLayer_new(II) * API_EMI / 136 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,9)  = Ea1NLayer_new(II) * OMT_EMI / 136 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,10) = Ea1NLayer_new(II) * FAR_EMI / 204 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,11) = Ea1NLayer_new(II) * BCA_EMI / 204 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,12) = Ea1NLayer_new(II) * OSQ_EMI / 204 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,13) = Ea1NLayer_new(II) * MBO_EMI /  86 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,14) = Ea1NLayer_new(II) * MET_EMI /  32 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,15) = Ea1NLayer_new(II) * ACT_EMI /  58 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,16) = Ea1NLayer_new(II) * CH4_EMI /  16 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,17) = Ea1NLayer_new(II) * NO_EMI  /  30 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,18) = Ea1NLayer_new(II) * ACA_EMI /  44 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,19) = Ea1NLayer_new(II) * FOR_EMI /  30 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,20) = Ea1NLayer_new(II) * CO_EMI  /  28 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,21) = Ea1NLayer_new(II) * CIN_EMI / 154 * Avog / (z(II)-z(II-1)) / 100 * s2(II)
       EMI(II,22) = Ea1NLayer_new(II) * LIN_EMI / 154 * Avog / (z(II)-z(II-1)) / 100 * s2(II)

    ENDDO
!!$
!!$      EMI(1,1)  = Ea1Layer_new(1)   * ISO_EMI /  68 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,2)  = Ea1NLayer_new(1)  * MYR_EMI / 136 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,3)  = Ea1NLayer_new(1)  * SAB_EMI / 136 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,4)  = Ea1NLayer_new(1)  * LIM_EMI / 136 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,5)  = Ea1NLayer_new(1)  * CAR_EMI / 136 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,6)  = Ea1NLayer_new(1)  * OCI_EMI / 136 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,7)  = Ea1NLayer_new(1)  * BPI_EMI / 136 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,8)  = Ea1NLayer_new(1)  * API_EMI / 136 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,9)  = Ea1NLayer_new(1)  * OMT_EMI / 136 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,10) = Ea1NLayer_new(1)  * FAR_EMI / 204 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,11) = Ea1NLayer_new(1)  * BCA_EMI / 204 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,12) = Ea1NLayer_new(1)  * OSQ_EMI / 204 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,13) = Ea1NLayer_new(1)  * MBO_EMI /  86 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,14) = Ea1NLayer_new(1)  * MET_EMI /  32 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,15) = Ea1NLayer_new(1)  * ACT_EMI /  58 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,16) = Ea1NLayer_new(1)  * CH4_EMI /  16 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,17) = Ea1NLayer_new(1)  * NO_EMI  /  30 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,18) = Ea1NLayer_new(1)  * ACA_EMI /  44 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,19) = Ea1NLayer_new(1)  * FOR_EMI /  30 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,20) = Ea1NLayer_new(1)  * CO_EMI  /  28 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,21) = 0.                * CIN_EMI / 154 * Avog / z(1) / 100 * s2(1)
!!$      EMI(1,22) = 0.                * LIN_EMI / 154 * Avog / z(1) / 100 * s2(1)
!!$


    DO II = k_canopy+1, layers
       sun_par_new(II) = 1.
       DO JJ = 1,22
          EMI(II,JJ) = 0.
       ENDDO
       Ea1tLayer_new(II) = 0.
       Ea1pLayer_new(II) = 0.
       Ea1Layer_new(II)  = 0.
    ENDDO

    sun_par_new(1) = sun_par_new(2)

  END SUBROUTINE EMISSION_M1
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!   Calculates the solar zenith angle
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    FUNCTION Calcbeta(Day, Lat, Hour)


    IMPLICIT NONE

    INTEGER :: Day
!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
 
    REAL(kind=dp) :: Pi, Rpi, Rpi180, Hour, Lat, SinDelta, CosDelta, A, B, Sinbeta, Calcbeta


    Pi       = 3.14159
	Rpi180   = 57.29578

    SinDelta = -SIN(0.40907) * COS(6.28 * (Day + 10) / (365))
    CosDelta = (1 - SinDelta**2.)**0.5

    A = SIN(Lat / Rpi180) * SinDelta
    B = COS(Lat / Rpi180) * Cosdelta

    Sinbeta = A + B * COS(2 * Pi * (Hour - 12) / 24)

    Calcbeta = ASIN(Sinbeta) * 57.29578

    END FUNCTION


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION DIstomata
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    FUNCTION DIstomata(DI)

    
	IMPLICIT NONE

!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  

    REAL(kind=dp) :: DIhigh, DIlow, DIstomata, DI


    ! > -.5 incipient,  mild or no drought; < -4 extreme drought
    DIhigh = -0.5
	DIlow  = -5
   
    IF (DI > DIhigh) THEN
       DIstomata = 1  ! no drought
    ELSEIF (DI > DIlow) THEN
       DIstomata = 1 - (0.9 * ((DI - DIhigh) / (DIlow - DIhigh))) ! interpolate
    ELSE
       DIstomata = 0  ! Maximum drought, maximum stomatal resistance
    ENDIF

    END FUNCTION DIstomata


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION CalcEccentricity
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    Function CalcEccentricity(Day)


	IMPLICIT NONE

!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
    INTEGER :: Day

    REAL(kind=dp) :: CalcEccentricity

   
	CalcEccentricity = 1 + 0.033 * COS(2 * 3.14 * (Day - 10) / 365)

    END FUNCTION CalcEccentricity


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   SUBROUTINE GaussianIntegration
!
!   Calculate weighting factors and distances
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    SUBROUTINE GaussianIntegration(Weightgauss, Distgauss, k_canopy)


	IMPLICIT NONE

!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  

    INTEGER :: Day, II, k_canopy 

    REAL(kind=dp), DIMENSION(k_canopy) :: Weightgauss, Distgauss


    IF (k_canopy .EQ. 1) THEN
       Weightgauss(1) = 1
	   Distgauss(1)   = 0.5
    ELSEIF (k_canopy .EQ. 3) THEN
       Weightgauss(1) = 0.277778
	   Weightgauss(2) = 0.444444
	   Weightgauss(3) = 0.277778
       Distgauss(1)   = 0.112702
	   Distgauss(2)   = 0.5
	   Distgauss(3)   = 0.887298
    ELSEIF (k_canopy .EQ. 5) THEN
       Weightgauss(1) = 0.1184635
	   Weightgauss(2) = 0.2393144
	   Weightgauss(3) = 0.284444444
	   Weightgauss(4) = 0.2393144
	   Weightgauss(5) = 0.1184635
       Distgauss(1)   = 0.0469101
	   Distgauss(2)   = 0.2307534
	   Distgauss(3)   = 0.5 
	   Distgauss(4)   = 0.7692465
	   Distgauss(5)   = 0.9530899
    ELSE
       DO II = 1, k_canopy
          Weightgauss(II) = 1. / k_canopy
          Distgauss(II)   = (II - 0.5) / k_canopy
       ENDDO
    ENDIF

    END SUBROUTINE GaussianIntegration


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   SUBROUTINE SolarFractions
!
!   Transmission, fraction of PPFD that is diffuse, fraction of solar rad that is PPFD
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    SUBROUTINE SolarFractions(Timeperiod, Solar, Maxsolar, Transmis, FracDiff, PPFDfrac)


	IMPLICIT NONE

!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
    INTEGER :: Timeperiod

    REAL(kind=dp) :: Solar, Maxsolar, Transmis, FracDiff, PPFDfrac, TransMin, TransSlope


    IF (Timeperiod .EQ. 1) THEN     ! Daily transmission
       TransMin   = 0.26
       TransSlope = 1.655
    ELSE                       ! Hourly transmission
       TransMin   = 0.26
       TransSlope = 1.655
    ENDIF
  
    IF (Maxsolar <= 0) THEN
       Transmis = 0.5
    ELSE
       Transmis = Solar / Maxsolar
    ENDIF
  
  
    ! Estimate diffuse fraction based on daily transmission (Roderick 1999, Goudriann and Van Laar 1994- P.33)
  
    IF (Transmis > 0.81) THEN  
       FracDiff = 0.05
    ELSEIF (Transmis > TransMin) THEN
       FracDiff = 0.96 - TransSlope * (Transmis - TransMin)
    ELSE
       FracDiff = 0.96
    ENDIF

    !The fraction of total solar radiation that is PPFD (43% to 55%) G. and L. 84
	PPFDfrac = 0.43 + FracDiff * 0.12

    END SUBROUTINE


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   SUBROUTINE WeightSLW
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    SUBROUTINE WeightSLW(Distgauss, Weightgauss, SLW, LAI, k_canopy)


	IMPLICIT NONE

!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  

    INTEGER :: II, k_canopy  

    REAL(kind=dp) :: SLWsum, RealII, LAI

    REAL(kind=dp), DIMENSION(k_canopy) :: Distgauss, Weightgauss, SLW


    SLWsum = 0
  
    DO II = 1, k_canopy
	   RealII = II
       SLW(II) = (0.63 + 0.37 * EXP(-(LAI * Distgauss(II)) * ((RealII - 1) / k_canopy)))
       SLWsum = SLWsum + SLW(II) * Weightgauss(II)
    ENDDO
  
    DO II = 1, k_canopy 
	   SLW(II) = SLW(II) / SLWsum
	ENDDO

    END SUBROUTINE WeightSLW


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine CanopyRad
!
!   Canopy light environment model
!   Code developed by Alex Guenther, based on Spitters et al. (1986), Goudrian and Laar (1994), Leuning (1997)
!   Initial code 8-99, modified 7-2000 and 12-2001
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
    SUBROUTINE CanopyRad(Distgauss, Sinbeta, Qbeamv, Qdiffv, Qbeamn, Qdiffn, Cantype, Canopychar,  &
                         Sunfrac, QbAbsV, QdAbsV, QsAbsV, QbAbsn, QdAbsn, QsAbsn, SunQv,           &
	  				     ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD, NrCha, NrTyp, LAI, k_canopy)



	IMPLICIT NONE

!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  

    INTEGER :: II, JJ, NrCha, NrTyp, Cantype, k_canopy 

    REAL(kind=dp), DIMENSION(k_canopy) :: Sunfrac, ShadePPFD, SunPPFD, Qdabsv, Qsabsv, Qsabsn, Shadeqv,        &
	                                   SunQn, Distgauss, QdAbsn, SunQv, ShadeQn 

    REAL(kind=dp) :: ScatV, ScatN, RefldV, RefldN, ReflbV, ReflbN, Kb, Kd, KbpV, KbpN, KdpV, KdpN,           &
	                LAIdepth, Cluster, ConvertPPFD, QdAbsVL, QsAbsVL, QdAbsNL, QsAbsNL, Qbeamv, Qdiffv,     & 
					Sinbeta, QbAbsV, Qbeamn, Qbabsn, Qdiffn, LAI 
					 
    REAL(kind=dp), DIMENSION(NrCha,NrTyp)  :: Canopychar

 
    If (((Qbeamv + Qdiffv) > 0.001) .AND. (Sinbeta > 0.00002) .AND. (LAI > 0.001)) Then    ! Daytime
       ConvertPPFD = 4.766  !4.55
  
       ! Scattering coefficients (scatV,scatN), diffuse and beam reflection coefficients (ref..) for visible or near IR
       ScatV   = Canopychar(5,Cantype)
       ScatN   = Canopychar(6,Cantype)
       RefldV  = Canopychar(7,Cantype)
       RefldN  = Canopychar(8,Cantype)
       Cluster = Canopychar(9,Cantype)

       ! Extinction coefficients for black leaves for beam (kb) or diffuse (kd)
       ! Kb = Cluster * 1 / Sinbeta                                                 ! (0.5 assumes a spherical leaf angle distribution (0.5 = cos (60 deg))
       Kb = Cluster * 0.5 / Sinbeta                                                 ! (0.5 assumes a spherical leaf angle distribution (0.5 = cos (60 deg))
       Kd = 0.8 * Cluster                                                           ! (0.8 assumes a spherical leaf angle distribution)

       Call CalcExtCoeff(Qbeamv, ScatV, Kb, Kd, ReflbV, KbpV, KdpV, QbAbsV)

       Call CalcExtCoeff(Qbeamn, ScatN, Kb, Kd, ReflbN, KbpN, KdpN, QbAbsN)

       DO II = 1, k_canopy
          LAIdepth   = LAI * Distgauss(II)     ! LAI depth at this layer
          Sunfrac(II) = EXP(-Kb * LAIdepth)    !fraction of leaves that are sunlit
          !write(*,*) 'xxx', ii, laidepth, lai, distgauss(II), Sunfrac(II)
          Call CalcRadComponents(Qdiffv, Qbeamv, kdpV, kbpV, kb, scatV, refldV, reflbV, LAIdepth, QbAbsV, QdAbsVL, QsAbsVL)
   
          Call CalcRadComponents(Qdiffn, Qbeamn, kdpN, kbpN, kb, scatN, refldN, reflbN, LAIdepth, QbAbsn, QdAbsNL, QsAbsNL)
   
          ShadePPFD(II) = (QdAbsVL + QsAbsVL) * ConvertPPFD / (1 - scatV)
          SunPPFD(II)   = ShadePPFD(II) + (QbAbsV * ConvertPPFD / (1 - scatV))
          QdAbsV(II)    = QdAbsVL 
		  QsAbsV(II)    = QsAbsVL 
		  QdAbsn(II)    = QdAbsNL
		  QsAbsn(II)    = QsAbsNL
          ShadeQv(II)   = QdAbsVL + QsAbsVL
          SunQv(II)     = ShadeQv(II) + QbAbsV
          ShadeQn(II)   = QdAbsNL + QsAbsNL
          SunQn(II)     = ShadeQn(II) + QbAbsn
       ENDDO

    Else                                                                            ! Night time
       QbAbsV = 0
	   QbAbsn = 0
       
	   DO II = 1, k_canopy
          Sunfrac(II)   = 0.2
		  SunQn(II)     = 0
		  ShadeQn(II)   = 0
		  SunQv(II)     = 0
		  ShadeQv(II)   = 0
		  SunPPFD(II)   = 0
		  ShadePPFD(II) = 0
          QdAbsV(II)    = 0
		  QsAbsV(II)    = 0
		  QdAbsn(II)    = 0
		  QsAbsn(II)    = 0
       ENDDO
    ENDIF
    
	End SUBROUTINE CanopyRad


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine CalcExtCoeff
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    SUBROUTINE CalcExtCoeff(Qbeam, scat, kb, kd, reflb, kbp, kdp, QbeamAbsorb)


    IMPLICIT NONE

!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
    REAL(kind=dp) :: Qbeam, scat, Kb, Kd, Reflb, Kbp, Kdp, QbeamAbsorb, P


    P     = (1 - scat)**0.5
    Reflb = 1 - Exp((-2 * ((1 - P) / (1 + P)) * kb) / (1 + kb))
    
	! Extinction coefficients
	Kbp   = Kb * P 
    Kdp   = Kd * P
    
	QbeamAbsorb = kb * Qbeam * (1 - scat) !calculate Absorbed beam radiation

    END SUBROUTINE CalcExtCoeff


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine CalcRadComponents
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    SUBROUTINE CalcRadComponents(Qdiff, Qbeam, kdp, kbp, kb, scat, refld, reflb, LAIdepth, QbAbs, QdAbs, QsAbs)


    IMPLICIT NONE

!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
    REAL(kind=dp) :: Qdiff, Qbeam, kdp, kbp, kb, scat, refld, reflb, LAIdepth, QbAbs, QdAbs, QsAbs
    

	QdAbs = Qdiff *   Kdp * (1 - Refld) * Exp(-Kdp * LAIdepth)
    QsAbs = Qbeam * ((Kbp * (1 - Reflb) * Exp(-Kbp * LAIdepth)) - (Kb * (1 - Scat) * Exp(-Kb * LAIdepth)))

    END SUBROUTINE CalcRadComponents


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!   Convert water mixing ratio (g/kg) to water vapor pressure (Pa or Kpa depending on units of input )
!   Mixing ratio (g/kg), temp (C), pressure (KPa)!
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
    FUNCTION WaterVapPres(Dens, Pres, WaterAirRatio) 


    IMPLICIT NONE

!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
    REAL(kind=dp) :: Dens, Dens001, Pres, WaterVapPres, WaterAirRatio


    Dens001 = Dens * 0.001

    WaterVapPres = (Dens001 / (Dens001 + WaterAirRatio)) * Pres

    END FUNCTION


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    FUNCTION Stability(Canopychar, Cantype, Solar, NrCha, NrTyp)


    IMPLICIT NONE
!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
    INTEGER :: Cantype, NrCha, NrTyp  

    REAL(kind=dp) :: Solar, Trateboundary, Stability

    REAL(kind=dp), DIMENSION(NrCha, NrTyp)  :: Canopychar


    Trateboundary = 500
 
    IF (Solar > Trateboundary) THEN
	    ! Daytime temperature lapse rate
       Stability = Canopychar(12, Cantype)
    ELSEIF (Solar > 0) THEN
       Stability = Canopychar(12, Cantype) - ((Trateboundary - Solar) / Trateboundary) * & 
	              (Canopychar(12, Cantype) - Canopychar(13, Cantype))
    ELSE
       ! Nightime temperature lapse rate
	   Stability = Canopychar(13, Cantype) 
    ENDIF

    END FUNCTION Stability


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine CanopyEB
!
!   Canopy energy balance model for estimating leaf temperature
!   Code developed by Alex Guenther, based on Goudrian and Laar (1994), Leuning (1997)
!   Initial code 8-99, modified 7-2000 and 12-2001
!
!   Note: i denotes an array containing a vertical profile through the canopy with 0 (above canopy conditions) plus 1 to number of canopy layers
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    SUBROUTINE CanopyEB(Trate, Distgauss, Canopychar, Cantype, StomataDI, TairK, HumidairPa, Ws, SunPPFD,    &
	                    ShadePPFD, SunQv, ShadeQv, SunQn, ShadeQn, Sunleaftk, SunleafSH, SunleafLH,          &
						SunleafIR, Shadeleaftk, ShadeleafSH, ShadeleafLH, ShadeleafIR, NrCha, NrTyp,         &
						k_canopy)


    IMPLICIT NONE

!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  

    INTEGER :: II, Cantype, NrCha, NrTyp, k_canopy

    REAL(kind=dp) :: 	Cdepth, Lwidth, Llength, Cheight, Eps, TranspireType, Deltah,            &
	                    IRin, Trate, IRout, StomataDI  !,                              & 
                      !  UnexposedLeafIRin, ExposedLeafIRin

    REAL(kind=dp), DIMENSION(k_canopy) :: Distgauss, TairK, HumidairPa, Ws, SunPPFD, SunQv, SunQn, SunleafIR,    &
	                                   Sunleaftk, SunleafSH, SunleafLH, ShadePPFD, ShadeQv, ShadeQn,          &
									   ShadeleafIR, Shadeleaftk, ShadeleafSH, ShadeleafLH                     

    REAL(kind=dp), DIMENSION(NrCha, NrTyp)  :: Canopychar
									 

 
	Cdepth        = Canopychar(1, Cantype)
    Lwidth        = Canopychar(2, Cantype)
	Llength       = Canopychar(3, Cantype)
    Cheight       = Canopychar(4, Cantype)
	Eps           = Canopychar(10,Cantype)
    TranspireType = Canopychar(11,Cantype)


    DO II = 1, k_canopy
  
       IRin            = UnexposedLeafIRin(TairK(k_canopy+1-II), Eps)
	   ShadeleafIR(II) = 2 * IRin

	   SunleafIR(II)   = 0.5 * ExposedLeafIRin(HumidairPa(k_canopy+1-II), TairK(k_canopy+1-II)) + 1.5 * IRin


       ! Sun
       CALL LeafEB(SunPPFD(II), SunQv(II) + SunQn(II), SunleafIR(II), Eps, TranspireType, Lwidth, Llength, TairK(k_canopy+1-II), &
            HumidairPa(k_canopy+1-II), Ws(II), Sunleaftk(II), SunleafSH(II), SunleafLH(II), IRout, StomataDI)
  
       SunleafIR(II) = SunleafIR(II) - IRout
  
       ! Shade
       CALL LeafEB(ShadePPFD(II), ShadeQv(II) + ShadeQn(II), ShadeleafIR(II), Eps, TranspireType, Lwidth, Llength, &
            TairK(k_canopy+1-II), HumidairPa(k_canopy+1-II), Ws(II), Shadeleaftk(II), ShadeleafSH(II), &
            ShadeleafLH(II), IRout, StomataDI)
 
       ShadeleafIR(II) = ShadeleafIR(II) - IRout
 
    ENDDO
    END SUBROUTINE CanopyEB


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION UnexposedLeafIRin
!
!   Calculate IR into leaf that is not exposed to the sky
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    FUNCTION UnexposedLeafIRin(Tk, Eps)


    IMPLICIT NONE

!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
    REAL(kind=dp) :: Eps, Tk, UnexposedLeafIRin

    REAL(kind=dp), PARAMETER ::  Sigma1    = 5.67051e-8        ! Boltzman constant (W/m2/K4)

   
	UnexposedLeafIRin = Eps * Sigma1 * (Tk**4.) !

    END FUNCTION UnexposedLeafIRin

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION ExposedLeafIRin
!
!   Calculate IR into leaf that is exposed to the sky
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    FUNCTION ExposedLeafIRin(HumidPa, Tk)


    IMPLICIT NONE

!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
    REAL(kind=dp) :: Tk, HumidPa, EmissAtm, ExposedLeafIRin
    REAL(kind=dp), PARAMETER ::  Sigma1    = 5.67051e-8        ! Boltzman constant (W/m2/K4)


      
	! Apparent atmospheric emissivity for clear skies: function of water vapor pressure (Pa) and ambient 
	! Temperature (K) based on Brutsaert(1975) referenced in Leuning 1997
    
	EmissAtm        = 0.642 * (HumidPa / Tk)**(1./7.)
	ExposedLeafIRin = EmissAtm * Sigma1 * (Tk**4.)  

    END FUNCTION ExposedLeafIRin


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine LeafEB
!
!   Leaf energy balance
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    SUBROUTINE LeafEB(PPFD, Q, IRin, Eps, TranspireType, Lwidth, Llength, TairK, HumidairPa, Ws,   &
	                  Tleaf, SH, LH, IRout, StomataDI)


    IMPLICIT NONE

!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
    INTEGER :: II 

    REAL(kind=dp) :: PPFD, Q, IRin, Eps, TranspireType, Lwidth, Llength, TairK, HumidairPa, Ws, Tleaf, SH, LH, IRout, StomataDI,  &
                    HumidAirKgm3, GHforced, StomRes, IRoutairT, LatHv, LHairT,    &
					Tdelt, Balance, E1, GH1, SH1, LH1, IRout1, GH, Pa, Tk


    IF (Ws <= 0) Then 
	   Ws = 0.001
    ENDIF

	! Air vapor density kg m-3
    HumidAirKgm3 = 0.002165 * HumidairPa / TairK

	! Heat convection coefficient (W m-2 K-1) for forced convection. Nobel page 366   
	GHforced = 0.0259 / (0.004 * ((Llength / Ws)**0.5)) 

	! Stomatal resistence s m-1
	StomRes  = ResSC(PPFD, stomataDI)
    
	IRoutairT = LeafIROut(tairK, eps)
    
	! Latent heat of vaporization (J Kg-1)
	LatHv = LHV(TairK)
    
	! Latent heat flux
	LHairT = LeafLE(TairK, HumidAirKgm3, LatHv, GHforced, StomRes, TranspireType)
    
	Tdelt = 1
	Balance = 10
	II = 1
    E1 = (Q + IRin - IRoutairT - LHairT)
    IF (E1 .EQ. 0.) THEN
	   E1 = -1.
	ENDIF

    DO II = 1, 10
	   IF (ABS(Balance) > 2) THEN
          GH1     = LeafBLC(GHforced, Tdelt, Llength)       ! Boundary layer conductance
          SH1     = LeafH(Tdelt, GH1)                       ! Convective heat flux
          LatHv   = LHV(TairK + Tdelt)                      ! Latent heat of vaporization (J Kg-1)
          LH      = LeafLE(TairK + Tdelt, HumidAirKgm3, LatHv, GH1, StomRes, TranspireType)
          LH1     = LH - LHairT
          IRout   = LeafIROut(TairK + Tdelt, Eps)
          IRout1  = IRout - IRoutairT
          Tdelt   = E1 / ((SH1 + LH1 + IRout1) / Tdelt)
          Balance = Q + IRin - IRout - SH1 - LH
	   ENDIF
	ENDDO

	If (Tdelt > 10)  Tdelt = 10
    If (Tdelt < -10) Tdelt = -10
    
	Tleaf = TairK + Tdelt

    GH    = LeafBLC(GHforced, Tleaf - TairK, Llength)
    SH    = LeafH(Tleaf - TairK, GH)
    LH    = LeafLE(Tleaf, HumidAirKgm3, LatHv, GH, StomRes, TranspireType)
    IRout = LeafIROut(Tleaf, Eps)

	END SUBROUTINE LeafEB


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!   Leaf stomatal cond. resistence s m-1
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    FUNCTION ResSC(Par, StomataDI) 


    IMPLICIT NONE

!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
    REAL(kind=dp) :: Par, StomataDI, SCadj, ResSC

 
    SCadj = StomataDI * ((0.0027 * 1.066 * Par) / ((1 + 0.0027 * 0.0027 * Par**2.)**0.5))

    IF (SCadj < 0.1) THEN
       ResSC = 2000
    ELSE
       ResSC = 200 / SCadj
    ENDIF
    
	END FUNCTION 


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!   IR thermal radiation energy output by leaf
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    FUNCTION LeafIROut(Tleaf, Eps) 


    IMPLICIT NONE

!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
    REAL(kind=dp) :: Tleaf, Eps, LeafIROut 
    REAL(kind=dp), PARAMETER ::  Sigma1    = 5.67051e-8        ! Boltzman constant (W/m2/K4)


 
    LeafIROut = Eps * Sigma1 * (2 * (Tleaf**4.))

    END FUNCTION 


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!   Latent Heat of vaporization(J Kg-1)from Stull p641
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    
	FUNCTION LHV(Tk)


    IMPLICIT NONE

!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
    REAL(kind=dp) :: Tk, LHV

 
    LHV = 2501000 - (2370 * (Tk - 273)) 
    
	END FUNCTION 


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!   Latent energy term in Energy balance
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    FUNCTION LeafLE(Tleaf, Ambvap, LatHv, GH, StomRes, TranspireType)


    IMPLICIT NONE

!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
    REAL(kind=dp) :: Tleaf, Ambvap, LatHv, GH, StomRes, TranspireType, LeafRes, Vapdeficit, LeafLE, LE
 
    LeafRes    = (1 / (1.075 * (GH / 1231))) + StomRes
    Vapdeficit = (SvdTk(Tleaf) - Ambvap)

	! Latent heat of vap (J Kg-1) * vap deficit(Kg m-3) / leaf resistence (s m-1)
    LE         = TranspireType * (1 / LeafRes) * LatHv * Vapdeficit
    
	IF (LE < 0) THEN 
	   LeafLE = 0 
	ELSE
	   LeafLE = LE
    ENDIF

	END FUNCTION 


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!   Boundary layer conductance
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    FUNCTION LeafBLC(GHforced, Tdelta, Llength)


    IMPLICIT NONE

!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
    REAL(kind=dp) :: GHforced, Tdelta, Llength, Ghfree, LeafBLC

 
    ! This is based on Leuning 1995 p.1198 except using molecular conductivity (.00253 W m-1 K-1 Stull p 640) instead of molecular
    ! diffusivity so that you end up with a heat convection coefficient (W m-2 K-1) instead of a conductance for free convection
    
	IF (Tdelta >= 0) THEN
       GhFree = 0.5 * 0.00253 * ((160000000 * Tdelta / (Llength**3.))**0.25) / Llength 
    ELSE
       GhFree = 0
    ENDIF

    LeafBLC = GHforced + GhFree
    
	END FUNCTION 


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!   Convective energy term in Energy balance (W m-2 heat flux from both sides of leaf)
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    FUNCTION LeafH(Tdelta, GH) 


    IMPLICIT NONE
 
!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
    REAL(kind=dp) :: Tdelta, GH, LeafH


    ! 2 sides X conductance X Temperature gradient
    LeafH = 2 * GH * Tdelta
    
	END FUNCTION 


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!   Saturation vapor density  (kg/m3)
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    Function SvdTk(Tk)


    IMPLICIT NONE
 
!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
    REAL(kind=dp) :: Tk, Svp, SvdTk
    

	Svp = 10**((-2937.4 / Tk) - (4.9283 * LOG10(Tk)) + 23.5518)  ! Saturation vapor pressure (millibars)
    
	SvdTk = 0.2165 * Svp / Tk 

    END FUNCTION


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!   Temperature dependence activity factor for emission type 1 (e.g. isoprene, MBO)
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    FUNCTION Ea1t99(Tm, Td)


    IMPLICIT NONE
 
!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
    REAL(kind=dp) :: Tm, Td, Ea1t99, Ctm1, Ctm2, Topt, X, Eopt
    
    
	IF (Tm < 260) THEN
       Ea1t99 = 0
    ELSE
	   ! Energy of activation and deactivation 
       Ctm1 = 95
	   Ctm2 = 230

	   ! Temperature at which maximum emission occurs
       Topt = 312.5 + 0.5 * (Td - 301) 
       X    = ((1 / Topt) - (1 / Tm)) / 0.00831
  
       ! Maximum emission (relative to emission at 30 C)
	   Eopt   = 1.9 * EXP(0.125 * (Td - 301)) 
       Ea1t99 = Eopt * Ctm2 * Exp(Ctm1 * X) / (Ctm2 - Ctm1 * (1 - EXP(Ctm2 * X)))
    ENDIF

    END FUNCTION


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    FUNCTION Ea1p99(LAI2, PPFD)


    IMPLICIT NONE
 
!    INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)  
    REAL(kind=dp) :: PPFD, Alpha, C1, Ea1p99, LAI2
    

    Alpha  = 0.001 + 0.00085 * LAI2
    C1     = 1.42 * EXP(-0.3 * LAI2)
    Ea1p99 = (Alpha * C1 * PPFD) / ((1 + Alpha**2. * PPFD**2.)**0.5)

    END FUNCTION


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

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

    END MODULE Megan_version_1

!!!!!! Not in use


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!Function ConvertHumidityGG2Pa(Pair, gg) !input air Pressure (Pa),humidity (g/g)
!ConvertHumidityGG2Pa = Pair * gg / 0.622
!End Function

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!Calculates water vapor Pressure in Pa
!Function CalcwaterVPpa(RH, tk) !input temperature in K
!CalcwaterVPpa = RH * 0.01 * (0.6112 * Exp((17.67 * (tk - 273.16)) / (tk - 29.66)))
!End Function

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!Calculates the saturation vapor Pressure in kPa, page 276 of Stull
!Function svpKPa(tk) !input temperature in K
!svpKPa = 0.6112 * Exp((17.67 * (tk - 273.16)) / (tk - 29.66))
!End Function

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!Function Log10(x) As Single
!   Log10 = Log(x) / Log(10#)
!End Function

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!Function ResBL(lwidth, llength, ws) As Single !boundary layer resistance m s-1
!ResBL = 200 * (lwidth ^ 0.2 * llength ^ 0.3 / ws ^ 0.5)
!End Function

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!Function DIstomata(DI)
!DIhigh = -0.5: DIlow = -5 !>-.5 incipient,  mild or no drought; <-4 extreme drought
! If DI > DIhigh Then
!  DIstomata = 1  ! no drought
! ElseIf DI > DIlow Then
!  DIstomata = 1 - (0.9 * ((DI - DIhigh) / (DIlow - DIhigh))) ! interpolate
! Else
!  DIstomata = 0  ! maximum drought: maximum stomatal resistance
! End If
!End Function

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!Function RH2Pa(RH, tc, P)
!Function RH2Pa()
!gkg = waterdens(RH, tc)
!RH2Pa = WaterVapPres(gkg, tc, P)
!End Function

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!Convert RH (%) to g/Kg,
!Function waterdens(RH, tc) !RH(%), Temp (C)
!Function waterdens() !RH(%), Temp (C)
!RH = 68: tc = 30
!waterdens = 0.01 * RH * (-0.082 + Exp(1.279 + 0.0687 * tc))
!End Function

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!Function CanopyNoEB(airt, k_canopy, Suntk, SunSH, SunLH, SunIR, Shadetk, ShadeSH, ShadeLH, ShadeIR)
!For i = 1 To k_canopy
! Suntk(i) = airt
! SunSH(i) = 0
! SunLH(i) = 0
! SunIR(i) = 0
! Shadetk(i) = airt
! ShadeSH(i) = 0
! ShadeLH(i) = 0
! ShadeIR(i) = 0
!Next i
!End Function

!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!Function Asin(x)
!Asin = Atn(x / Sqr(-x * x + 1))
!End Function
