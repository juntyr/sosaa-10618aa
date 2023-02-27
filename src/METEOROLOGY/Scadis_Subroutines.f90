! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
!> @defgroup scadis_subroutines 
!> Colection of subroutines and functions used here and there in the meteorology module.   
!
!> For more information see MT_MainMet.f90                                                 
!>                                                                                         
!> Commented and restructured by Rosa Gierens, University of Helsinki                      
!                                                                                               !
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !


MODULE Scadis_Subroutines

  ! variables from other modules
  USE constants_mod, ONLY : dp, T00
  
  ! subroutines and variables part of the meteorology scheme
  ! USE Scadis_Parameters

  IMPLICIT NONE

  PRIVATE

  !Public subroutines:
  PUBLIC ::  New2Old,  New2Old_HalfStep

  !Public variables:
  !PUBLIC ::   
  
  
CONTAINS



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !SUBROUTINE New2Old
   !> copy new 1-values from previous timestep to old 0-values to be used in this time step
   !> @ingroup scadis_subroutines
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE New2Old(n_levels, soil_water, soil_water1, soil_temperature, soil_temperature1, &
                      l, l1,  tke, tke_dissipation, tke1, tke_dissipation1,  &
                      u_wind,u_wind1, v_wind,v_wind1, w_wind,w_wind1, kt, kt1, &
                      temperature, temperature1,alt, alt1, leaf_temp_sunlit, leaf_temp_sunlit1, &
                      leaf_temp_shaded, leaf_temp_shaded1, humidity_abs, humidity_abs1)
      
      ! INPUT
      INTEGER,                            INTENT(IN) :: n_levels                    ! (=kz) number of levels, defined in Sosa_data [1]
      
      
      REAL(kind=dp), DIMENSION(n_levels), INTENT(IN) :: l1                       , & !
                                                        soil_temperature1        , & ! (=tsoil1) soil temperature [K]
                                                        tke1                     , & ! (=bt1) turbulent kinetic energy?
                                                        tke_dissipation1         , & ! (=dbt1) dissipation rate of TKE 
                                                        u_wind1                  , & ! (=u1) u-component of wind [m/s] 
                                                        v_wind1                  , & ! (=v1) v-component of wind [m/s] 
                                                        w_wind1                  , & ! (=w1) w-component of wind [m/s] 
                                                        kt1                      , & 
                                                        temperature1             , & ! (=ta1) air temperature [K]
                                                        leaf_temp_sunlit1        , & ! (=tsn1) sunlit leaf temperature [K]
                                                        leaf_temp_shaded1        , & ! (=tsd1) shaded leaf temperature [K]
                                                        humidity_abs1            , & ! (=qa1) absolute humidity of air [kg/m3]
                                                        alt1                         ! 1/prantdlt number [1]
                                                        
      REAL(kind=dp), DIMENSION(2),        INTENT(IN) :: soil_water1                  ! (=wg1) volumetric water content of the soil [kg/kg]


      ! OUTPUT
      REAL(kind=dp), DIMENSION(n_levels), INTENT(OUT):: l                        , & ! ?
                                                        soil_temperature         , & ! (=tsoil) soil temperature [K]
                                                        tke                      , & ! (=bt) turbulent kinetic energy?
                                                        tke_dissipation          , & ! (=dbt) dissipation rate of TKE 
                                                        u_wind                   , & ! (=u) u-component of wind [m/s] 
                                                        v_wind                   , & ! (=v) v-component of wind [m/s] 
                                                        w_wind                   , & ! (=w) w-component of wind [m/s] 
                                                        kt                       , & 
                                                        temperature              , & ! (=ta) air temperature [K]
                                                        leaf_temp_sunlit         , & ! (=tsn) sunlit leaf temperature [K]
                                                        leaf_temp_shaded         , & ! (=tsd) shaded leaf temperature [K]
                                                        humidity_abs             , & ! (=qa) absolute humidity of air [kg/m3]
                                                        alt                          ! 1/prantdlt number [1]
                                                                     
      REAL(kind=dp), DIMENSION(2),       INTENT(OUT) :: soil_water                   ! (=wg1) volumetric water content of the soil [kg/kg]


      ! LOCAL
      INTEGER :: k                                                            ! foor loop
      
                    
      soil_water = soil_water1
      !154  wg(1)=wg1(1)
      !155  wg(2)=wg1(2)
      
      soil_temperature = soil_temperature1
      !156  tsoil(k)=tsoil1(k)
      
      l = l1  
      !157  l(k)=l1(k)
      
      tke = tke1
      !158  bt(k)=bt1(k)

      tke_dissipation = tke_dissipation1
      !159  dbt(k)=dbt1(k)

      u_wind = u_wind1
      !160  u(k)=u1(k)
      
      v_wind = v_wind1
      !164  v(k)=v1(k)
      
      w_wind = w_wind1
      !163  w(k)=w1(k)
      
      kt = kt1
      !161  kt(k)=kt1(k)
      
      temperature = temperature1
      !162  ta(k)=ta1(k)
      
      alt = alt1
      !165  alt(k)=alt1(k)
      
      leaf_temp_sunlit = leaf_temp_sunlit1
      !166  tsn(k)=tsn1(k)
                              
      leaf_temp_shaded = leaf_temp_shaded1
      !167  tsd(k)=tsd1(k)
      
      humidity_abs = humidity_abs1
      !168  qa(k)=qa1(k)

      do k=1,n_levels
      
        humidity_abs(k) = min(humidity_abs(k), SatAbsoluteHum(temperature1(k)) )
        !169  qa(k)=min(qa(k), (qv0*exp(maga*(ta1(k)-T00)/(magb+ta1(k)-T00)) /ta1(k)) )
        
      enddo  
   

   END SUBROUTINE New2Old

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE New2Old
   !> Get old 0-values for this time step by taking the average of old and new values from previous time step ("half way")
   !> @ingroup scadis_subroutines
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE New2Old_HalfStep(n_levels, soil_water, soil_water1, soil_temperature, soil_temperature1, &
                               tke, tke1,  &
                               u_wind,u_wind1, v_wind,v_wind1, w_wind,w_wind1, kt, kt1, &
                               temperature, temperature1,alt, alt1, leaf_temp_sunlit, leaf_temp_sunlit1, &
                               leaf_temp_shaded, leaf_temp_shaded1, humidity_abs, humidity_abs1)
                               
      
      
      
      ! INPUT
      INTEGER,                            INTENT(IN) :: n_levels                    ! (=kz) number of levels, defined in Sosa_data [1]
      
      
      REAL(kind=dp), DIMENSION(n_levels), INTENT(IN) :: soil_temperature1        , & ! (=tsoil1) soil temperature [K]
                                                        tke1                     , & ! (=bt1) turbulent kinetic energy?
                                                        u_wind1                  , & ! (=u1) u-component of wind [m/s] 
                                                        v_wind1                  , & ! (=v1) v-component of wind [m/s] 
                                                        w_wind1                  , & ! (=w1) w-component of wind [m/s] 
                                                        kt1                      , & 
                                                        temperature1             , & ! (=ta1) air temperature [K]
                                                        leaf_temp_sunlit1        , & ! (=tsn1) sunlit leaf temperature [K]
                                                        leaf_temp_shaded1        , & ! (=tsd1) shaded leaf temperature [K]
                                                        humidity_abs1            , & ! (=qa1) absolute humidity of air [kg/m3]
                                                        alt1                         ! 1/prantdlt number [1]
                                                        
      REAL(kind=dp), DIMENSION(2),        INTENT(IN) :: soil_water1                  ! (=wg1) volumetric water content of the soil [kg/kg]


      ! OUTPUT
      REAL(kind=dp), DIMENSION(n_levels), INTENT(OUT):: soil_temperature         , & ! (=tsoil) soil temperature [K]
                                                        tke                      , & ! (=bt) turbulent kinetic energy?
                                                        u_wind                   , & ! (=u) u-component of wind [m/s] 
                                                        v_wind                   , & ! (=v) v-component of wind [m/s] 
                                                        w_wind                   , & ! (=w) w-component of wind [m/s] 
                                                        kt                       , & 
                                                        temperature              , & ! (=ta) air temperature [K]
                                                        leaf_temp_sunlit         , & ! (=tsn) sunlit leaf temperature [K]
                                                        leaf_temp_shaded         , & ! (=tsd) shaded leaf temperature [K]
                                                        humidity_abs             , & ! (=qa) absolute humidity of air [kg/m3]
                                                        alt                          ! 1/prantdlt number [1]
                                                                     
      REAL(kind=dp), DIMENSION(2),       INTENT(OUT) :: soil_water                   ! (=wg1) volumetric water content of the soil [kg/kg]

   
      tke = (tke+tke1)/2.
      !138  bt(k)=(bt(k)+bt1(k))/2.
      
      soil_water = (soil_water+soil_water1)/2.
      !139  wg(1)=(wg(1)+wg1(1))/2.
      !140  wg(2)=(wg(2)+wg1(2))/2.

      soil_temperature = (soil_temperature+soil_temperature1)/2.
      !141  tsoil(k)=(tsoil(k)+tsoil1(k))/2.          
      
      u_wind = (u_wind+u_wind1)/2.
      !142  u(k)=(u(k)+u1(k))/2.
      
      v_wind = (v_wind+v_wind1)/2.
      !145  v(k)=(v(k)+v1(k))/2.
      
      w_wind = (w_wind+w_wind1)/2.
      !144  w(k)=(w(k)+w1(k))/2.
      
      kt = (kt+kt1)/2.
      !143  kt(k)=(kt(k)+kt1(k))/2.

      temperature = (temperature+temperature1)/2.
      !146  ta(k)=(ta(k)+ta1(k))/2.
      
      humidity_abs = (humidity_abs+humidity_abs1)/2.
      !147  qa(k)=(qa(k)+qa1(k))/2.
      
      leaf_temp_sunlit = (leaf_temp_sunlit+leaf_temp_sunlit1)/2.
      !148  tsn(k)=(tsn(k)+tsn1(k))/2.
      
      leaf_temp_shaded = (leaf_temp_shaded+leaf_temp_shaded1)/2.
      !149  tsd(k)=(tsd(k)+tsd1(k))/2.
      
      alt = (alt+alt1)/2.
      !150  alt(k)=(alt(k)+alt1(k))/2.

   END SUBROUTINE  New2Old_HalfStep


    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION SatAbsoluteHum
   ! calculates saturation absolute huidity (kg/m3) at given temperature
   !
   ! repeats in code
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
   FUNCTION SatAbsoluteHum(temperature)

      ! copied this from Malte (Rosa 5.4.2012)

      ! INPUT
      REAL(kind=dp) :: temperature ! air temperature (K)
      
      ! OUTPUT
      REAL(kind=dp) :: SatAbsoluteHum   !saturated absolute humidity (kg/m3)

      ! LOCAL 
      REAL(kind=dp), PARAMETER ::  qv0 = 1.3318375, &
                                   maga = 17.502, &
                                   magb = 488.97
       ! absolute humidity:
       !   SatAbsoluteHum = rho_w = e_sat/(R_w T)
       !     
       !   saturation vapour pressure (hPa): e_sat = 6.11*exp( (maga*T_celcius)/(magb+T_celcius))
       !   Specific gas constant for water vapor R_w = 461.495 J/(kg K)
       !   qv0 = 6.11 hPa/R_w = 6.11 hPa/461.495 J/(kg K)  = 1.32395801. Close enough.



       SatAbsoluteHum = (qv0*exp(maga*(temperature-T00)/(magb+temperature-T00)) / temperature)

   END FUNCTION SatAbsoluteHum
    
    
END MODULE Scadis_Subroutines
