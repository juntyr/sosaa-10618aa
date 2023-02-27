!   By Alex Guenther --April 2005
!
!   Modified by Xuemei -July 2007  
!
! Modified for different IO
!¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
!
!   Input and output files have to be selected before starting the program
!
!¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
!!
!   Input varibles
!
!   Day                                         Julian day
!   Lat                                         Latitude
!   Hour                                        Hour of the day
!   Tc                                          Temperature [C] -> values look like [K] (rosa 12.9.2014)
!   PPFD                                        Incoming photosynthetic active radiation [umol/m2/s1]
!   Wind                                        Wind speed [m s-1]
!   Humidity                                    Mixing ratio [g/kg]
!   Cantype                                    Defines set of canopy characteristics
!   LAI                                         Leaf area index [m2 per m2 ground area]
!   DI                                          ???
!   Pres                                        Pressure [Pa]
!
!
!   Used variables:
!
!   PPFDfrac                            Fraction of total solar radiation that is PPFD
!   Solar                                       Solar radiation [W/m2]
!   Maxsolar                            Maximum of solar radiation
!   Beta                                        Sin of solar angle above horizon
!   Sinbeta                                     Solar angle above horizon
!   TairK0                                      Above canopy air temperature [K]
!       TairK                                   Array of canopy air temperature [K]
!   Ws0                                         Above canopy wind speed [m/s]
!   Ws                                          Array of canopy wind speed [m/s]
!   HumidairPA0                         Above canopy ambient humidity [Pa]
!       HumidairPa                              Array of canopy ambient humidity in [Pa]
!   StomataDI                           An index indicating water status of leaves. used to modify stomatal conductance
!   Transmis                            Transmission of PPFD that is diffuse
!   Difffrac                            Fraction of PPFD that is diffuse
!   PPFDfrac                            Fraction of solar rad that is PPFD
!   Trate                                       Stability of boundary ???
!   SH                                          Sensible heat flux ???
!       VPgausWt                                Array of gaussian weighting factors
!   VPgausDis                           Array of gaussian weighting factors
!   VPslwWT                                     Array of gaussian weighting factors
!   SunFrac                                     Array of the fraction of sun leaves. i = 1 is the top canopy layer, 2 is the next layer, etc.
!   SunPPFD                                     Array of incoming (NOT absorbed) PPFD on a sun leaf [umol/m2/s]
!   ShadePPFD                           Array of incoming (NOT absorbed) PPFD on a shade leaf [umol/m2/s]
!   SunQv                                       Array of visible radiation (in and out) fluxes on sun leaves
!   ShadeQv                                     Array of absorbed visible radiation (in and out) fluxes on shade leaves
!   SunQn                                       Array of absorbed near IR radiation (in and out) fluxes on sun leaves
!   ShadeQn                                     Array of absorbed near IR radiation (in and out) fluxes on shade leaves
!   SunleafTK                           Array of leaf temperature for sun leaves [K]
!   SunleafSH                           Array of sensible heat flux for sun leaves [W/m2]
!   SunleafLH                           Array of latent heat flux for sun leaves [W/m2]
!   SunleafIR                           Array of infrared flux for sun leaves [W/m2]
!   ShadeleafTK                         Array of leaf temperature for shade leaves [K]
!   ShadeleafSH                         Array of sensible heat flux for shade leaves [W/m2]
!   ShadeleafLH                         Array of latent heat flux for shade leaves [W/m2]
!   ShadeleafIR                         Array of infrared flux for shade leaves [W/m2]
!   QbAbsV, QbAbsN                      Absorbed direct beam light for visible and near infra red
!   QdAbsV, QdAbsN                      Array of absorbed diffuse light for visible and near infra red
!   QsAbsV, QsAbsN                      Array of absorbed scattered light for visible and near infra red
!   QBeamV, QBeamN                      Above canopy beam (direct) light for visible and near infra red
!   QDiffV, QDiffN                      Above canopy diffuse light for visible and near infra red
!       Ea1pLayer                               Array of emission activity of light per layer
!       Ea1tLayer                               Array of emission activity of temperature per layer
!       Ea1Layer                                Array of companied emission activity
!       Ea1pCanopy                              Total emission activity of light
!   Ea1tCanopy                          Total emission activity of temperature
!   Ea1Canopy                           Total companied emission activity
!
!   Calcbeta                            Function: Calculation of solar zenith angle
!   WaterVapPres                        Function: Convert water mixing ratio (kg/kg) to water vapor pressure
!   Stability                           Function: Temperature lapse rate
!   Ea1t99                                      Function: Temperature dependence activity factor for emission type 1
!   Ea1p99                                      Function:
!   DIstomata                           Function:
!   CalcEccentricity            		Function:
!
!¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
MODULE M2_canopy

  USE constants_mod, ONLY : dp

  IMPLICIT NONE
 
  public

CONTAINS

  SUBROUTINE GAMME_CE(JDATE,JTIME,NCLOS,NROWS,BETA,LAT, LONG,PPFD0,Tc,LAI,Wind,Pres,Humidity,DI,  &
       PFTF,Canopychar,NrCha,NrTyp,NPFT, Ealt,Ealp, Ealpt, LAdp, midpoint,    &
       Layers, gammasout, sun_par, Sunleaftk, Shadeleaftk, Sunfrac,check_rh, station, Sunleaftk_in, Shadeleaftk_in)

    IMPLICIT NONE

    !	   INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
    
    ! added by Rosa 12.8.2014
    CHARACTER(len=3), INTENT(IN) :: station 
    INTEGER                      :: start_ind, end_ind

    INTEGER :: Layers   ! Number of layers inside the canopy

    REAL(kind=dp), PARAMETER :: ConvertWm2toUmolm2s = 4.766, &                        ! Convert radiation from [W/m2] to [umol/m2/s1]
         SolarConstant       = 1367,  &                        ! Solar constant [W/m2]
         Td                  = 301,   &                        ! Temperature used in Ea1t99
         WaterAirRatio       = 18.016 / 28.97, &               ! Ratio between water and air molecules
         Cce  = 0.57

    INTEGER :: II, JJ, KK,I,J,JDATE,JTIME,K, NCLOS,NROWS,NrCha,NrTyp,NPFT, counter

    INTEGER  Day(NCLOS,NROWS)

    REAL(kind=dp), DIMENSION(NCLOS,NROWS) :: Pres, LAI, DI,PPFD
    REAL(kind=dp), DIMENSION(NCLOS,NROWS) :: PPFD0
    REAL(kind=dp), DIMENSION(NrCha,NrTyp) :: Canopychar             ! Array of canopy characteristics for NrTyp of canopy type

    REAL(kind=dp) LONG(NCLOS,NROWS),Hour(NCLOS,NROWS)

    REAL(kind=dp) :: Beta
    REAL(kind=dp),DIMENSION(NCLOS,NROWS) :: Lat,Sinbeta
    REAL(kind=dp),DIMENSION(NCLOS,NROWS,Layers) ::  VPgausWt, VPgausDis, VPslwWT, Sunfrac, QdAbsV, QsAbsV, QdAbsn,           &  
         QsAbsn, SunQv, ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD, TairK,       &
         HumidairPa, Ws, Sunleaftk, SunleafSH, SunleafLH,SunleafIR, Shadeleaftk,  &
         ShadeleafSH,  ShadeleafLH,ShadeleafIR,Ea1pLayer, Ea1tLayer, Ea1Layer

    REAL(kind=dp), DIMENSION(NCLOS,NROWS) :: Solar, Maxsolar, Transmis,Difffrac, PPFDfrac, QbAbsn, Trate, StomataDI, Qbeamv, &
         Qdiffv, Qbeamn, Qdiffn, QbAbsV,Ea1tCanopy, Ea1pCanopy, Ea1Canopy, TairK0,       &
         HumidairPa0,Ws0, SH,PFTP

    REAL(kind=dp), DIMENSION(NPFT,NCLOS,NROWS)::  Ealpt,PFTF,Ealt,Ealp

    REAL(kind=dp)::  Calcbeta

    INTEGER, DIMENSION(NCLOS,NROWS) ::  Cantype

    REAL(kind=dp), DIMENSION(Layers)     :: LADp, midpoint, sun_par, Sunleaftk_in, Shadeleaftk_in
    REAL(kind=dp), DIMENSION(Layers)     :: check_rh
    REAL(kind=dp), DIMENSION(Layers+1)   :: Tc, Wind, Humidity !tempreature, wind, humidity
    REAL(kind=dp), DIMENSION(3,NPFT,Layers) :: gammasout ! gamma factors bundled

    TairK(1,1,:) = Tc( 1 : Layers ) ! [K]
    Ws(1,1,:)    = Wind( 1 : Layers )

    DO counter=1, Layers
       HumidairPa(1,1, counter) = WaterVapPres(Humidity(counter), Pres(1,1), WaterAirRatio)
        check_rh(counter) = WaterVapPres(Humidity(counter), Pres(1,1), WaterAirRatio)
     END DO
 

    Day(:,:)  =  MOD(JDATE,1000)

    ! Convert from XXXXXX format to XX.XX (solar hour)
    ! HOUR = 0 -> 23.xx
    DO I = 1, NCLOS
       DO J = 1, NROWS
          Hour(I,J) = JTIME/10000 + LONG(I,J)/15   ! Solar hour
          IF ( Hour(I,J) .LT. 0.0 ) THEN
             Hour(I,J) = Hour(I,J) + 24.0
             Day(I,J) = Day(I,J) - 1
          ELSEIF ( Hour(I,J) .GE. 24.0 ) THEN
             ! print*,'Invalid hour: HOUR(I,J) is ', Hour(I,J)
          ENDIF

          ! Sinbeta now with Beta values from Scadis
          Sinbeta(I,J)   = SIN((90-Beta) / 57.29578)
          IF (Sinbeta(I,J) .EQ. 1) THEN  
             Sinbeta(I,J) = 0.
          ENDIF

          TairK0(I,J)    = Tc((Layers+1))                           !Air temperature above canopy
          Ws0(I,J)       = Wind((Layers+1))                         ! Wind speed above canopy
          PPFD(I,J)      = PPFD0(I,J)
          StomataDI(I,J) = DIstomata(DI(I,J))
          Solar(I,J)     = PPFD(I,J)/ConvertWm2toUmolm2s*2
          Maxsolar(I,J)  = Sinbeta(I,J) * SolarConstant * CalcEccentricity(Day(I,J))

       ENDDO
    ENDDO
    ! Call GaussianIntegration(VPgausWt, VPgausDis, Layers,NCLOS,NROWS)

    !NOTE you might want to see if this( VPgausWt and VPgausDis) is correct. I think it is.

    VPgausWt(1,1, : )  = LADp(1 : Layers)
    VPgausDis(1,1, : ) = midpoint 

    Call SolarFractions(0, Solar, Maxsolar, Transmis, Difffrac, PPFDfrac,NCLOS,NROWS)

    DO I=1, NCLOS
       DO J=1, NROWS

          Qdiffv(I,J) = PPFDfrac(I,J) * Solar(I,J) * Difffrac(I,J)
          Qbeamv(I,J) = PPFDfrac(I,J) * Solar(I,J) * (1 - Difffrac(I,J))
          Qdiffn(I,J) = (1 - PPFDfrac(I,J)) * Solar(I,J) * Difffrac(I,J)
          Qbeamn(I,J) = (1 - PPFDfrac(I,J)) * Solar(I,J) * (1 - Difffrac(I,J))
       ENDDO
    ENDDO

    Call WeightSLW(VPgausDis, VPgausWt, LAI, Layers, VPslwWT,NCLOS,NROWS)


    !Rosa 12.8.2014: this loop is for all different types of vegetation in the canopy. In Hyytiälä we only have trees with needles, so we
    ! don't really need a loop: only K = 2 needs to be considered. Hence it makes the model faster not to include other K. 
    ! HOWEVER, other types of vegetation exist at other stations. So, I added the variables start_ind and end_ind, that define the 
    ! starting and ending indexes of the loop.
    
    IF ( (station == 'HYY') ) THEN
       start_ind = 2
       end_ind   = 2
    ELSE
       start_ind = 1
       end_ind   = NPFT
    END IF    
       
    
    !DO K = 1, NPFT
    !DO K = 2, 2
    
    DO K = start_ind, end_ind

       DO I=1, NCLOS
          DO J=1, NROWS
             PFTP(I,J)=PFTF(K,I,J)
             !Cantype (i,j) = K
             Cantype (i,j) = 1
          ENDDO
       ENDDO

       Call CanopyRad(VPgausDis, VPgausWt, Layers, LAI, Sinbeta, Qbeamv, Qdiffv, Qbeamn, Qdiffn,Cantype ,Canopychar, Sunfrac, &
            QbAbsV, QdAbsV, QsAbsV, QbAbsn, QdAbsn, QsAbsn, SunQv, ShadeQv, SunQn, ShadeQn, SunPPFD,      &
            ShadePPFD,PFTP, NrCha, NrTyp,NCLOS,NROWS)


       DO I=1, NCLOS
          DO J=1, NROWS
             HumidairPa0(I,J) =  WaterVapPres(Humidity((Layers+1)), Pres(I,J), WaterAirRatio)

             Trate(I,J)   =  Stability(Canopychar, Cantype(I,J), Solar(I,J), NrCha, NrTyp)
          ENDDO
       ENDDO

       Call CanopyEB(Trate, Layers, VPgausDis, Canopychar, Cantype, StomataDI, TairK, HumidairPa, Ws, SunPPFD, ShadePPFD, & 
            SunQv, ShadeQv, SunQn, ShadeQn, Sunleaftk, SunleafSH, SunleafLH,SunleafIR, Shadeleaftk, ShadeleafSH, &
            ShadeleafLH, ShadeleafIR, NrCha, NrTyp, Ws0, TairK0, HumidairPa0,NCLOS,NROWS)


       DO I=1, NCLOS
          DO J=1, NROWS
             Ea1tCanopy(I,J)        = 0.
             Ea1pCanopy(I,J)        = 0.
             Ea1Canopy(I,J)         = 0.
             SH(I,J)                = 0.
             !Layers are what we want out. Ea1player = emission activity of light per layer,
             !                             Ea1tlayer = emission activity of temperature per layer,
             !                             Ea1layer  = companied emission activity per layer

             DO II = 1, Layers
             
               ! problems with the leaf temperature from Scadis. Using values for sunlit and shaded leaf temperature from Scadis.
               ! Don't want to risk mixing up the algorithm some how, so commented out these lines where the leaf temperatures are 
               ! used and replacing them with the input values
               
                !Ea1tLayer(I,J,II) = Ea1t99(Sunleaftk(I,J,II), Td)   * Sunfrac(I,J,II) + &
                !     Ea1t99(Shadeleaftk(I,J,II), Td) * (1 - Sunfrac(I,J,II))

                Ea1tLayer(I,J,II) = Ea1t99(Sunleaftk_in(II), Td)   * Sunfrac(I,J,II) + &
                     Ea1t99(Shadeleaftk_in(II), Td) * (1 - Sunfrac(I,J,II))


                Ea1pLayer(I,J,II) = Ea1p99(VPgausDis(I,J,II) * LAI(I,J), SunPPFD(I,J,II))   * Sunfrac(I,J,II) + &
                     Ea1p99(VPgausDis(I,J,II) * LAI(I,J), ShadePPFD(I,J,II)) * (1 - Sunfrac(I,J,II))

                SH(I,J)           = SH(I,J) + (SunleafSH(I,J,II) * Sunfrac(I,J,II) + ShadeleafSH(I,J,II) * &
                     (1 - Sunfrac(I,J,II))) * LAI(I,J) * VPgausWt(I,J,II)

                !Ea1Layer(I,J,II)  = Ea1t99(Sunleaftk(I,J,II), Td)   * Ea1p99(VPgausDis(I,J,II) * &
                !     LAI(I,J), SunPPFD(I,J,II))   * Sunfrac(I,J,II) + Ea1t99(Shadeleaftk(I,J,II), Td) * &
                !     Ea1p99(VPgausDis(I,J,II) * LAI(I,J), ShadePPFD(I,J,II)) * (1 - Sunfrac(I,J,II))

                Ea1Layer(I,J,II)  = Ea1t99(Sunleaftk_in(II), Td)   * Ea1p99(VPgausDis(I,J,II) * &
                     LAI(I,J), SunPPFD(I,J,II))   * Sunfrac(I,J,II) + Ea1t99(Shadeleaftk_in(II), Td) * &
                     Ea1p99(VPgausDis(I,J,II) * LAI(I,J), ShadePPFD(I,J,II)) * (1 - Sunfrac(I,J,II))


                IF (K .eq. 2) then
                   IF (SunPPFD(I,J,1) .GT. 0.) THEN
                      sun_par(II) = SunPPFD(I,J,II) / SunPPFD(I,J,1)
                   ELSE
                      sun_par(II) = 1.
                   ENDIF
                ENDIF

                !What is VPslwWT and should it be used?

                Ea1pCanopy(I,J)   = Ea1pCanopy(I,J) + Ea1pLayer(I,J,II) * VPslwWT(I,J,II) * VPgausWt(I,J,II)
                Ea1tCanopy(I,J)   = Ea1tCanopy(I,J) + Ea1tLayer(I,J,II) * VPslwWT(I,J,II) * VPgausWt(I,J,II)
                Ea1Canopy(I,J)    = Ea1Canopy(I,J)  + Ea1Layer(I,J,II)  * VPslwWT(I,J,II) * VPgausWt(I,J,II)

                gammasout(1, K ,II) = Ea1tLayer(I,J,II)                    ! * VPslwWT(I,J,II)  ! * VPgausWt(I,J,II)
                gammasout(2, K ,II) = Ea1pLayer(I,J,II)                    ! * VPslwWT(I,J,II)  ! * VPgausWt(I,J,II)
                gammasout(3, K ,II) = Ea1Layer(I,J,II) !* Cce * LAI(I,J)    ! * VPgausWt(I,J,II) ! * VPslwWT(I,J,II)

             ENDDO

             Ea1Canopy(I,J) = Ea1Canopy(I,J) !* Cce * LAI(I,J) ! and Cce is a parameter, about what?
             Ealpt(K,I,J)   = Ea1Canopy(I,J)
             Ealt(K,I,J)    = Ea1tCanopy(I,J)
             Ealp(K,I,J)    = Ea1pCanopy(I,J) 
          ENDDO
       ENDDO
       sun_par(1:Layers) = sun_par(Layers : 1 : -1)

    ENDDO
    RETURN
  END SUBROUTINE GAMME_CE

  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   FUNCTION
  !
  !   Calculates the solar zenith angle
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  FUNCTION Calcbeta(Day, Lat, Hour)

    IMPLICIT NONE

    !       INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

    REAL(kind=dp) :: Pi, Rpi, Rpi180, Hour, Lat, SinDelta, &
         CosDelta, A, B, Sinbeta, Calcbeta
    INTEGER :: Day

    Pi       = 3.14159
    Rpi180   = 57.29578

    SinDelta = -SIN(0.40907) * COS(6.28 * (Day + 10) / (365))
    CosDelta = (1 - SinDelta**2.)**0.5

    A = SIN(Lat / Rpi180) * SinDelta
    B = COS(Lat / Rpi180) * Cosdelta
    Sinbeta = A + B * COS(2 * Pi * (Hour - 12) / 24)
    Calcbeta = ASIN(Sinbeta) * 57.29578
  END FUNCTION Calcbeta

  !-------------------------------------------------------------------------
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   FUNCTION DIstomata
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  FUNCTION DIstomata(DI)

    IMPLICIT NONE
    !       INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
    REAL(kind=dp) :: DIhigh, DIlow, DI, DIstomata

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

  Function CalcEccentricity(Day)


    IMPLICIT NONE

    !         INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

    INTEGER :: Day

    REAL(kind=dp) :: CalcEccentricity

    CalcEccentricity = 1 + 0.033 * COS(2 * 3.14 * &
         (Day - 10) / 365)

  END FUNCTION CalcEccentricity

  !--This subroutine (Gaussianintegration) is replaced with values for the layers coming from Sosa--
  !--------------------------------------------------------------------------------
  !        SUBROUTINE GaussianIntegration(Weightgauss, Distgauss, 
  !     &                                   Layers,NCLOS,NROWS)
  !
  !        IMPLICIT NONE
  !        INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
  !        INTEGER :: I,J, II, Layers,NCLOS,NROWS
  !
  !        REAL(kind=dp), DIMENSION(NCLOS,NROWS,Layers) :: 
  !     &                     Weightgauss, Distgauss
  !
  !        DO I=1, NCLOS
  !        DO J=1, NROWS
  !
  !       IF (Layers .EQ. 1) THEN
  !       Weightgauss(I,J,1) = 1
  !       Distgauss(I,J,1)   = 0.5
  !       ELSEIF (Layers .EQ. 3) THEN
  !       Weightgauss(I,J,1) = 0.277778
  !       Weightgauss(I,J,2) = 0.444444
  !       Weightgauss(I,J,3) = 0.277778
  !       Distgauss(I,J,1)   = 0.112702
  !       Distgauss(I,J,2)   = 0.5
  !       Distgauss(I,J,3)   = 0.887298
  !        ELSEIF (Layers .EQ. 5) THEN
  !        Weightgauss(I,J,1) = 0.1184635
  !        Weightgauss(I,J,2) = 0.2393144
  !        Weightgauss(I,J,3) = 0.284444444
  !        Weightgauss(I,J,4) = 0.2393144
  !        Weightgauss(I,J,5) = 0.1184635
  !        Distgauss(I,J,1)   = 0.0469101
  !        Distgauss(I,J,2)   = 0.2307534
  !        Distgauss(I,J,3)   = 0.5
  !        Distgauss(I,J,4)   = 0.7692465
  !        Distgauss(I,J,5)   = 0.9530899
  !        ELSE
  !        DO II = 1, Layers
  !        Weightgauss(I,J,II) = 1. / Layers
  !        Distgauss(I,J,II)   = (II - 0.5) / Layers
  !         ENDDO
  !         ENDIF
  !         ENDDO
  !         ENDDO
  !         RETURN
  !         END SUBROUTINE GaussianIntegration
  
  
  !-------------------------------------------------------------------------------
  !   SUBROUTINE WeightSLW
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  !This should be checked, Anton


  SUBROUTINE WeightSLW(Distgauss, Weightgauss, &
       LAI, Layers, SLW,NCLOS,NROWS)
    IMPLICIT NONE

    !        INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

    INTEGER :: II, Layers,I,J,NCLOS,NROWS

    REAL(kind=dp) ::  SLWsum, RealII

    REAL(kind=dp), DIMENSION(NCLOS,NROWS,Layers) :: &
         Distgauss, Weightgauss, SLW
    REAL(kind=dp), DIMENSION(NCLOS,NROWS) :: LAI

    DO I=1,NCLOS
       DO J=1,NROWS
          SLWsum = 0
          DO II = 1, Layers
             RealII = II
             SLW(I,J,II) = (0.63 + 0.37 * EXP(-(LAI(I,J) * Distgauss(I,J,II)) * ((RealII - 1) / Layers)))
             SLWsum = SLWsum + SLW(I,J,II) * Weightgauss(I,J,II)
          ENDDO
          DO II = 1, Layers
             SLW(I,J,II) = SLW(I,J,II) / SLWsum
          ENDDO
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE WeightSLW

  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   SUBROUTINE SolarFractions
  !
  !   Transmission, fraction of PPFD that is diffuse, fraction of solar rad that is PPFD
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  SUBROUTINE SolarFractions(Timeperiod, Solar, Maxsolar, &
       Transmis, FracDiff, PPFDfrac,NCLOS,NROWS)

    IMPLICIT NONE

    !        INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

    INTEGER :: Timeperiod,NCLOS,NROWS,I,J
    REAL(kind=dp) ,DIMENSION(NCLOS,NROWS) :: Solar, Maxsolar, Transmis, FracDiff, PPFDfrac
    REAL(kind=dp) ::  TransMin, TransSlope
    IF (Timeperiod .EQ. 1) THEN     ! Daily transmission
       TransMin   = 0.26
       TransSlope= 1.655
    ELSE                       ! Hourly transmission
       TransMin   = 0.26
       TransSlope = 1.655
    ENDIF
    DO I=1,NCLOS
       DO J=1,NROWS
          IF (Maxsolar(I,J) <= 0) THEN
             Transmis(I,J) = 0.5
          ELSE
             Transmis(I,J) = Solar(I,J) / Maxsolar(I,J)
          ENDIF

          ! Estimate diffuse fraction based on daily transmission (Roderick 1999, Goudriann and Van Laar 1994- P.33)

          !  Estimate diffuse fraction based on daily transmission (Roderick 1999, Goudriann and Van Laar 1994- P.33)

          IF (Transmis(I,J) > 0.81) THEN
             FracDiff(I,J) = 0.05
          ELSEIF (Transmis(I,J) > TransMin) THEN
             FracDiff(I,J) = 0.96-TransSlope * (Transmis(I,J) - TransMin)
          ELSE
             FracDiff(I,J) = 0.96
          ENDIF

          !The fraction of total solar radiation that is PPFD (43% to 55%) G. and L. 84
          PPFDfrac(I,J) = 0.43 + FracDiff(I,J) * 0.12
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE SolarFractions
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   Subroutine CanopyRad
  !
  !   Canopy light environment model
  !   Code developed by Alex Guenther, based on Spitters et al. (1986), Goudrian and Laar (1994), Leuning (1997)
  !   Initial code 8-99, modified 7-2000 and 12-2001
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  SUBROUTINE CanopyRad(Distgauss, LADp, Layers, LAI, Sinbeta, &
       Qbeamv, Qdiffv, Qbeamn, Qdiffn, Cantype, &
       Canopychar, Sunfrac, QbAbsV, QdAbsV, QsAbsV, & 
       QbAbsn, QdAbsn, QsAbsn, SunQv, &  
       ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD, & 
       PFTP, NrCha, NrTyp,NCLOS,NROWS)



    IMPLICIT NONE

    !        INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

    INTEGER :: II, JJ, Layers, NrCha, NrTyp,NCLOS,NROWS,I,J

    REAL(kind=dp) :: Sb, ScatV, ScatN, RefldV, RefldN, ReflbV, &
         ReflbN,Kb, Kd, KbpV, KbpN, KdpV, KdpN, LAIdepth, &
         Cluster, ConvertPPFD,QdAbsVL, QsAbsVL, QdAbsNL, &
         QsAbsNL, QqV, Qqn

    REAL(kind=dp),DIMENSION(NCLOS,NROWS,Layers) ::  Distgauss, &
         Sunfrac,QdAbsV, QsAbsV, QdAbsn, QsAbsn, SunQv, &
         ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD


    REAL(kind=dp), DIMENSION(NCLOS,NROWS) :: Solar, Maxsolar, PFTP, &
         Difffrac,PPFDfrac, QbAbsn,Qbeamv,Sinbeta,Qdiffv, &
         Qbeamn, Qdiffn, QbAbsV

    INTEGER, DIMENSION(NCLOS,NROWS) :: Cantype

    REAL(kind=dp), DIMENSION(NCLOS,NROWS) :: LAI
    REAL(kind=dp), DIMENSION(NrCha,NrTyp)  :: Canopychar
    
    REAL(kind=dp), DIMENSION(Layers) :: LADp

    DO I=1,NCLOS
       DO J=1,NROWS
          IF ( PFTP(I,J).NE.0.0) THEN
             If (((Qbeamv(I,J) + Qdiffv(I,J)) > 0.001) .AND. (Sinbeta(I,J) > 0.00002) .AND. (LAI(I,J) > 0.001)) Then    ! Daytime
                Sb = 0.0000000567

                ConvertPPFD = 4.766

                ! Scattering coefficients (scatV,scatN), diffuse and beam reflection coefficients (ref..) for visible or near IR
                ScatV   = Canopychar(5,Cantype(I,J))
                ScatN   = Canopychar(6,Cantype(I,J))
                RefldV  = Canopychar(7,Cantype(I,J))
                RefldN  = Canopychar(8,Cantype(I,J))
                Cluster = Canopychar(9,Cantype(I,J))
                !          print*,'cluster',  Cluster
                ! Extinction coefficients for black leaves for beam (kb) or diffuse (kd)
                ! Kb = Cluster * 1 / Sinbeta   
                Kb = Cluster * 0.5 / Sinbeta(I,J)
                ! (0.5 assumes a spherical leaf angle distribution (0.5 = cos (60 deg))
                Kd = 0.8 * Cluster              
                ! (0.8 assumes a spherical leaf angle distribution)

                !       Call CalcExtCoeff(Qbeamv(I,J), ScatV, Kb, Kd, ReflbV,
                !     &              KbpV, KdpV, QbAbsV(I,J))
                Call CalcExtCoeff(Qbeamv(I,J), ScatV, Kb, Kd, ReflbV, &
                     KbpV, KdpV, QqV)
                QbAbsV(I,J) =  QqV
                !       Call CalcExtCoeff(Qbeamn(I,J), ScatN, Kb, Kd, ReflbN, 
                !     &                  KbpN, KdpN, QbAbsn(I,J))
                Call CalcExtCoeff(Qbeamn(I,J), ScatN, Kb, Kd, ReflbN, &
                     KbpN, KdpN, Qqn)
                QbAbsn(I,J) =  Qqn

                DO II = 1, Layers
                   !LAIdepth   = LAI(I,J) * Distgauss(I,J,II)     ! LAI depth at this layer
                   LAIdepth = LAI(I,J) * (SUM(LADp(1:II-1))+LADp(II)*0.5)  ! corrected by Rosa 11.8.2014
                   
                   Sunfrac(I,J,II) = EXP(-Kb * LAIdepth)    !fraction of leaves that are sunlit

                   Call CalcRadComponents(Qdiffv(I,J), Qbeamv(I,J), kdpV, & 
                        kbpV, kb, scatV, refldV, &
                        reflbV, LAIdepth, QbAbsV(I,J), QdAbsVL, QsAbsVL) 

                   Call CalcRadComponents(Qdiffn(I,J), Qbeamn(I,J), kdpN, &
                        kbpN, kb, scatN, refldN, &
                        reflbN, LAIdepth, QbAbsn(I,J), QdAbsNL, QsAbsNL) 

                   ShadePPFD(I,J,II) = (QdAbsVL + QsAbsVL) * &
                        ConvertPPFD / (1 - scatV)
                   SunPPFD(I,J,II) = ShadePPFD(I,J,II) + (QbAbsV(I,J) * & 
                        ConvertPPFD / (1 - scatV))
                   QdAbsV(I,J,II)    = QdAbsVL
                   QsAbsV(I,J,II)    = QsAbsVL
                   QdAbsn(I,J,II)    = QdAbsNL
                   QsAbsn(I,J,II)    = QsAbsNL
                   ShadeQv(I,J,II)   = QdAbsVL + QsAbsVL
                   SunQv(I,J,II)     = ShadeQv(I,J,II) + QbAbsV(I,J)
                   ShadeQn(I,J,II)   = QdAbsNL + QsAbsNL
                   SunQn(I,J,II)     = ShadeQn(I,J,II) + QbAbsn(I,J)
                ENDDO

             Else                                                                            ! Night time
                QbAbsV(I,J) = 0
                QbAbsn (I,J) = 0

                DO II = 1, Layers
                   Sunfrac(I,J,II)   = 0.2
                   SunQn(I,J,II)     = 0
                   ShadeQn(I,J,II)   = 0
                   SunQv(I,J,II)     = 0
                   ShadeQv(I,J,II)   = 0
                   SunPPFD(I,J,II)   = 0
                   ShadePPFD(I,J,II) = 0
                   QdAbsV(I,J,II)    = 0
                   QsAbsV(I,J,II)    = 0
                   QdAbsn(I,J,II)    = 0
                   QsAbsn(I,J,II)    = 0

                ENDDO
             ENDIF

          ELSE
             DO II = 1, Layers
                Sunfrac(I,J,II)   = 0.
                SunQn(I,J,II)     = 0
                ShadeQn(I,J,II)   = 0
                SunQv(I,J,II)     = 0
                ShadeQv(I,J,II)   = 0
                SunPPFD(I,J,II)   = 0
                ShadePPFD(I,J,II) = 0
                QdAbsV(I,J,II)    = 0
                QsAbsV(I,J,II)    = 0
                QdAbsn(I,J,II)    = 0
                QsAbsn(I,J,II)    = 0
             ENDDO
             QbAbsV(I,J) = 0
             QbAbsn (I,J) = 0
          ENDIF

       ENDDO
    ENDDO
    RETURN
  End SUBROUTINE CanopyRad

  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   Subroutine CalcExtCoeff
  !
  ! Edited by Rosa 17.7.2014: intent in/out for variables, commented out the last line (not needed)
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  SUBROUTINE CalcExtCoeff(Qbeam, scat, kb, kd, reflb, &
       kbp, kdp, QbeamAbsorb)

    IMPLICIT NONE

    !        INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

    !REAL(kind=dp) :: Qbeam, scat, Kb, Kd, Reflb, Kbp, Kdp, &
    !     QbeamAbsorb, P
    
    REAL(kind=dp), INTENT(IN) :: Qbeam, scat, Kb, Kd
    REAL(kind=dp), INTENT(OUT) :: Reflb, Kbp, Kdp, QbeamAbsorb
    REAL(kind=dp) :: P

    P     = (1 - scat)**0.5
    Reflb = 1 - Exp((-2 * ((1 - P) / (1 + P)) * kb) / (1 + kb))

    ! Extinction coefficients
    Kbp   = Kb * P
    Kdp   = Kd * P
    QbeamAbsorb = kb * Qbeam * (1 - scat) !calculate Absorbed beam radiation
    RETURN
  END SUBROUTINE CalcExtCoeff
  

  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   Subroutine CalcRadComponents
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  SUBROUTINE CalcRadComponents(Qdiff, Qbeam, kdp, kbp, kb, & 
       scat, refld, reflb, LAIdepth, QbAbs, QdAbs, QsAbs)


    IMPLICIT NONE

    !        INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

    REAL(kind=dp) :: Qdiff, Qbeam, kdp, kbp, kb, scat, refld, reflb, &
         LAIdepth, QbAbs, QdAbs, QsAbs
         
    QdAbs = Qdiff *   Kdp * (1 - Refld) * Exp(-Kdp * LAIdepth)
    QsAbs = Qbeam * ((Kbp * (1 - Reflb) * Exp(-Kbp * LAIdepth)) - & 
         (Kb * (1 - Scat) * Exp(-Kb * LAIdepth)))
               
    RETURN
  END SUBROUTINE CalcRadComponents


  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   FUNCTION
  !
  !   Convert water mixing ratio (kg/kg) to water vapor pressure (Pa or Kpa depending on units of input )
  !   Mixing ratio (kg/kg), temp (C), pressure (KPa)!
  !
  !   Rosa: multiplying Dens with 0.001, so input should be in g/kg!
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  FUNCTION WaterVapPres(Dens, Pres, WaterAirRatio)
    IMPLICIT NONE

    !       INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

    REAL(kind=dp) :: Dens, Dens001, Pres, WaterVapPres, WaterAirRatio
    Dens001 = Dens * 0.001 !uncommented rusan [g kg-1] -> [kg]

    WaterVapPres = (Dens001 / (Dens001 + WaterAirRatio)) * Pres

  END FUNCTION WaterVapPres


  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   FUNCTION
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  FUNCTION Stability(Canopychar, Cantype, Solar, NrCha, NrTyp)
    IMPLICIT NONE

    !       INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

    INTEGER :: Cantype, NrCha, NrTyp

    REAL(kind=dp) :: Solar, Trateboundary, Stability

    REAL(kind=dp), DIMENSION(NrCha, NrTyp)  :: Canopychar


    Trateboundary = 500

    IF (Solar > Trateboundary) THEN
       ! Daytime temperature lapse rate
       Stability = Canopychar(12, Cantype)
    ELSEIF (Solar > 0) THEN
       Stability = Canopychar(12, Cantype) - &
            ((Trateboundary - Solar) / Trateboundary) * &
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
  !
  ! NOTE winnd and temperature are now directly give to this subroutine asa layers, It doesn't calculate them anymore,
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  SUBROUTINE CanopyEB(Trate, Layers, Distgauss, Canopychar, &
       Cantype, StomataDI, TairK, HumidairPa, Ws, &
       SunPPFD, ShadePPFD, SunQv, ShadeQv, SunQn, ShadeQn, &
       Sunleaftk, SunleafSH, SunleafLH, &
       SunleafIR, Shadeleaftk, ShadeleafSH, &
       ShadeleafLH, ShadeleafIR, NrCha, NrTyp, Ws0, &
       TairK0, HumidairPa0,NCLOS,NROWS)


    IMPLICIT NONE

    !       INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

    INTEGER :: II, NrCha, NrTyp, Layers,I,J,NCLOS,NROWS,TAGI1,TAGI2

    REAL(kind=dp) :: Cdepth, Lwidth, Llength, Cheight, Eps, &
         TranspireType,Deltah, Ldepth, Wsh, &
         IRin,IRout

    REAL(kind=dp), DIMENSION(NrCha, NrTyp)  :: Canopychar

    REAL(kind=dp),DIMENSION(NCLOS,NROWS,Layers) ::  Distgauss, &
         SunQv,ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD, &
         TairK, HumidairPa, Ws, Sunleaftk, SunleafSH, SunleafLH, &
         SunleafIR, Shadeleaftk, ShadeleafSH, ShadeleafLH, ShadeleafIR


    REAL(kind=dp), DIMENSION(NCLOS,NROWS) ::  Trate, StomataDI, TairK0, HumidairPa0, Ws0

    INTEGER, DIMENSION(NCLOS,NROWS) :: Cantype


    DO I=1,NCLOS
       DO J=1,NROWS
          Cdepth        = Canopychar(1, Cantype(I,J))
          Lwidth        = Canopychar(2, Cantype(I,J))
          Llength       = Canopychar(3, Cantype(I,J))
          Cheight       = Canopychar(4, Cantype(I,J))
          Eps           = Canopychar(10,Cantype(I,J))
          TranspireType = Canopychar(11,Cantype(I,J))

          DO II = 1, Layers
             IRin                = UnexposedLeafIRin(TairK(I,J,layers+1-II), Eps)
             ShadeleafIR(I,J,II) = 2 * IRin
             SunleafIR(I,J,II)   = 0.5 * ExposedLeafIRin(HumidairPa0(I,J), TairK0(I,J)) + 1.5 * IRin

             ! Sun
             CALL LeafEB(SunPPFD(I,J,II), SunQv(I,J,II) + SunQn(I,J,II), &
                  SunleafIR(I,J,II), Eps, TranspireType, Lwidth, Llength, &
                  TairK(I,J,layers+1-II), HumidairPa(I,J,layers+1-II), Ws(I,J,II), &
                  Sunleaftk(I,J,II), &
                  SunleafSH(I,J,II),SunleafLH(I,J,II), IRout, &
                  StomataDI(I,J))

             SunleafIR(I,J,II) = SunleafIR(I,J,II) - IRout

             ! Shade
             CALL LeafEB(ShadePPFD(I,J,II), ShadeQv(I,J,II) + ShadeQn(I,J,II), &
                  ShadeleafIR(I,J,II), Eps, TranspireType, Lwidth, Llength, &
                  TairK(I,J,layers+1-II), HumidairPa(I,J,layers+1-II), Ws(I,J,II), &
                  Shadeleaftk(I,J,II), &
                  ShadeleafSH(I,J,II),ShadeleafLH(I,J,II), &
                  IRout, StomataDI(I,J))

             ShadeleafIR(I,J,II) = ShadeleafIR(I,J,II) - IRout

          ENDDO
       ENDDO
    ENDDO
    RETURN
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
    !       INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
    REAL(kind=dp) :: Eps, Tk, UnexposedLeafIRin
    REAL(kind=dp), PARAMETER ::  Sb    = 5.67051e-8        ! Boltzman constant (W/m2/K4)

    UnexposedLeafIRin = Eps * Sb * (Tk**4.) !

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

    !       INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
    REAL(kind=dp) :: Sb, Tk, HumidPa, EmissAtm, ExposedLeafIRin
    Sb = 0.0000000567 ! Stefan-boltzman constant  W m-2 K-4

    ! Apparent atmospheric emissivity for clear skies: function of water vapor pressure (Pa) and ambient
    ! Temperature (K) based on Brutsaert(1975) referenced in Leuning 1997

    EmissAtm        = 0.642 * (HumidPa / Tk)**(1./7.)
    ExposedLeafIRin = EmissAtm * Sb * (Tk**4.)

  END FUNCTION ExposedLeafIRin


  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   FUNCTION
  !
  !   Saturation vapor density  (kg/m3)
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  FUNCTION ConvertHumidityPa2kgm3(Pa, Tk)

    IMPLICIT NONE

    !          INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

    REAL(kind=dp) :: ConvertHumidityPa2kgm3
    REAL(kind=dp) :: Pa, Tk

    ConvertHumidityPa2kgm3 = 0.002165 * Pa / Tk


  END FUNCTION ConvertHumidityPa2kgm3




  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   Subroutine LeafEB
  !
  !   Leaf energy balance
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  SUBROUTINE LeafEB(PPFD, Q, IRin, Eps, TranspireType, &
       Lwidth, Llength, TairK, HumidairPa, Ws, Tleaf, &
       SH, LH, IRout, StomataDI)

    IMPLICIT NONE
    !        INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
    INTEGER :: II
    !FIX THIS!!!!
    !REAL(kind=dp) :: PPFD, Q, IRin, Eps, TranspireType, Lwidth, &
    !     Llength, TairK, HumidairPa, &
    !     Ws, Tleaf, SH, LH, IRout, StomataDI, HumidAirKgm3, & 
    !     ConvertHumidityPa2kgm3, &
    !     GHforced, StomRes, ResSC, IRoutairT, LeafIROut, &
    !     LatHv, LHV, LHairT, &  
    !     LeafLE, Tdelt, Balance, E1, GH1, LeafBLC, SH1, &
    !     LeafH, LH1, IRout1, GH

    REAL(kind=dp) :: PPFD, Q, IRin, Eps, TranspireType, Lwidth, &
         Llength, TairK, HumidairPa, &
         Ws, Tleaf, SH, LH, IRout, StomataDI, HumidAirKgm3, & 
         GHforced, StomRes, IRoutairT, &
         LatHv, LHairT, &  
         Tdelt, Balance, E1, GH1,SH1, &
         LH1, IRout1, GH


    IF (Ws <= 0) Then
       Ws = 0.001
    ENDIF


    ! Air vapor density kg m-3
    HumidAirKgm3 = ConvertHumidityPa2kgm3(HumidairPa, TairK)

    ! Heat convection coefficient (W m-2 K-1) for forced convection. Nobel page 366
    GHforced = 0.0259 / (0.004 * ((Llength / Ws)**0.5))

    ! Stomatal resistence s m-1
    StomRes  = ResSC(PPFD, stomataDI)

    IRoutairT = LeafIROut(tairK, eps)

    ! Latent heat of vaporization (J Kg-1)
    LatHv = LHV(TairK)

    ! Latent heat flux
    LHairT = LeafLE(TairK, HumidAirKgm3, LatHv, &
         GHforced, StomRes, TranspireType)

    E1 = (Q + IRin - IRoutairT - LHairT)
    IF (E1 .EQ. 0.) THEN
       E1 = -1.
    ENDIF

    II = 1
    Tdelt = 1
    Balance = 10
    DO II = 1, 10
       IF (ABS(Balance) > 2) THEN
          GH1     = LeafBLC(GHforced, Tdelt, Llength)       ! Boundary layer conductance
          SH1     = LeafH(Tdelt, GH1)                       ! Convective heat flux
          LatHv   = LHV(TairK + Tdelt)                      ! Latent heat of vaporization (J Kg-1)
          LH      = LeafLE(TairK + Tdelt, HumidAirKgm3, &
               LatHv, GH1, StomRes, TranspireType)
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
    LH    = LeafLE(Tleaf, HumidAirKgm3, LatHv, &
         GH, StomRes, TranspireType)
    IRout = LeafIROut(Tleaf, Eps)
    RETURN
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

    !          INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

    REAL(kind=dp) :: Par, StomataDI, SCadj, ResSC


    SCadj = StomataDI * ((0.0027 * 1.066 * Par) / &
         ((1 + 0.0027 * 0.0027 * Par**2.)**0.5))

    IF (SCadj < 0.1) THEN
       ResSC = 2000
    ELSE
       ResSC = 200 / SCadj
    ENDIF

  END FUNCTION ResSC


  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   FUNCTION
  !
  !   IR thermal radiation energy output by leaf
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  FUNCTION LeafIROut(Tleaf, Eps)
    IMPLICIT NONE
    !        INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
    REAL(kind=dp) :: Tleaf, Eps, Sb, LeafIROut

    Sb = 0.0000000567 ! Stefan-boltzman constant  W m-2 K-4
    LeafIROut = Eps * Sb * (2 * (Tleaf**4.))

  END FUNCTION LeafIROut


  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   FUNCTION
  !
  !   Latent Heat of vaporization(J Kg-1)from Stull p641
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  FUNCTION LHV(Tk)
    IMPLICIT NONE

    !        INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
    REAL(kind=dp) :: Tk, LHV
    LHV = 2501000 - (2370 * (Tk - 273))

  END FUNCTION LHV


  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   FUNCTION
  !
  !   Latent energy term in Energy balance
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  FUNCTION LeafLE(Tleaf, Ambvap, LatHv, GH, StomRes, &
       TranspireType)
    IMPLICIT NONE

    !        INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
    REAL(kind=dp) :: Tleaf, Ambvap, LatHv, GH, StomRes, &
         TranspireType, LeafRes, Vapdeficit, LeafLE, LE


    LeafRes    = (1 / (1.075 * (GH / 1231))) + StomRes
    Vapdeficit = (SvdTk(Tleaf) - Ambvap)
    ! Latent heat of vap (J Kg-1) * vap deficit(Kg m-3) / leaf resistence (s m-1)
    LE = TranspireType * (1 / LeafRes) * LatHv * Vapdeficit
    IF (LE < 0) THEN
       LeafLE = 0
    ELSE
       LeafLE = LE
    ENDIF


  END FUNCTION  LeafLE


  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   FUNCTION
  !
  !   Boundary layer conductance
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  FUNCTION LeafBLC(GHforced, Tdelta, Llength)


    IMPLICIT NONE

    !         INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

    REAL(kind=dp) :: GHforced, Tdelta, Llength, Ghfree, LeafBLC


    ! This is based on Leuning 1995 p.1198 except using molecular conductivity (.00253 W m-1 K-1 Stull p 640) instead of molecular
    ! diffusivity so that you end up with a heat convection coefficient (W m-2 K-1) instead of a conductance for free convection

    IF (Tdelta >= 0) THEN
       GhFree = 0.5 * 0.00253 * ((160000000 * Tdelta / & 
            (Llength**3.))**0.25) / Llength
    ELSE
       GhFree = 0
    ENDIF

    LeafBLC = GHforced + GhFree

  END FUNCTION LeafBLC

  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   FUNCTION
  !
  !   Convective energy term in Energy balance (W m-2 heat flux from both sides of leaf)
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  FUNCTION LeafH(Tdelta, GH)

    IMPLICIT NONE

    !        INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

    REAL(kind=dp) :: Tdelta, GH, LeafH
    ! 2 sides X conductance X Temperature gradient
    LeafH = 2 * GH * Tdelta

  END FUNCTION LeafH


  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   FUNCTION
  !
  !   Saturation vapor density  (kg/m3)
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  Function SvdTk(Tk)
    IMPLICIT NONE

    !        INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

    REAL(kind=dp) :: Tk, Svp, SvdTk


    Svp = 10**((-2937.4 / Tk) - (4.9283 * LOG10(Tk)) + 23.5518)  
    ! Saturation vapor pressure (millibars)

    SvdTk = 0.2165 * Svp / Tk

  END FUNCTION  SvdTk


  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   FUNCTION
  !
  !   Temperature dependence activity factor for emission type 1 (e.g. isoprene, MBO)
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  FUNCTION Ea1t99(Tm, Td)


    IMPLICIT NONE

    !         INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)

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
       Ea1t99 = Eopt * Ctm2 * Exp(Ctm1 * X) / &
            (Ctm2 - Ctm1 * (1 - EXP(Ctm2 * X)))
    ENDIF



  END FUNCTION  Ea1t99

  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !
  !   FUNCTION
  !
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  FUNCTION Ea1p99(LAI, PPFD)
    IMPLICIT NONE

    !          INTEGER, PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
    REAL(kind=dp) :: LAI, PPFD, Alpha, C1, Ea1p99


    Alpha  = 0.001 + 0.00085 * LAI
    C1     = 1.42 * EXP(-0.3 * LAI)
    Ea1p99 = (Alpha * C1 * PPFD) / &
         ((1 + Alpha**2. * PPFD**2.)**0.5)

  END FUNCTION  Ea1p99


END MODULE M2_canopy

