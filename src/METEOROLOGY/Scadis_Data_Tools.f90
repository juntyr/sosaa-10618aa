! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
! Processing of input data & nudging related parameters                                         !
! For more information see MT_MainMet.f90                                                       !
!                                                                                               !
! Commented and restructured by Rosa Gierens, University of Helsinki                            !
!                                                                                               !
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !


MODULE Scadis_Data_Tools

  ! variables from other modules
  USE constants_mod, ONLY : dp
    
  IMPLICIT NONE

  PRIVATE

  !Public subroutines:
  PUBLIC :: MakeNudgingData, make_nudging_data_new

  !Public variables:
  !PUBLIC ::   


CONTAINS


SUBROUTINE make_nudging_data_new( &
  station, method, &
  mod_kz, mod_z, mod_montime, mod_data, &
  obs_nl, obs_lv, obs_ind, obs_dt, obs_data, &
  nud_data, &
  obs_zbound, obs_databound)
  CHARACTER(LEN=*), INTENT(IN) :: station  ! name of the station, defined in Sosa_data.f90
  CHARACTER(LEN=*), INTENT(IN) :: method   ! interpolation method, 'linear', 'spline'

  INTEGER , INTENT(IN) :: mod_kz  ! (=kz) number of levels, defined in Sosa_data [1]
  REAL(dp), INTENT(IN) :: mod_z(mod_kz)  ! (=z) height of each (vertical) layer [m]
  REAL(dp), INTENT(IN) :: mod_montime  ! (=time_in_month) time [s] from the beginning of the month
  REAL(dp), INTENT(IN) :: mod_data(mod_kz)  ! u, v, ta, qa

  INTEGER , INTENT(IN) :: obs_nl  ! number of levels for observation data
  REAL(dp), INTENT(IN) :: obs_lv(obs_nl)  ! levels for observation data
  INTEGER , INTENT(IN) :: obs_ind  ! (=nxodrad) current index of observation data
  REAL(dp), INTENT(IN) :: obs_dt  ! (=dt_obs) time interval of observation data
  REAL(dp), INTENT(IN) :: obs_data(1489, obs_nl)  ! observation data

  REAL(dp), INTENT(OUT) :: nud_data(mod_kz)

  REAL(dp), OPTIONAL :: obs_zbound(2), obs_databound(2)  ! z and data at 2 boundaries

  REAL(dp) :: obs_data_now(obs_nl)
  REAL(dp) :: obs_data_nonan(obs_nl)
  REAL(dp) :: obs_lv_nonan(obs_nl)
  INTEGER  :: obs_nl_nonan

  INTEGER :: i, j

  ! Interpolate observation data to current time
  obs_nl_nonan = 0
  DO i = 1, obs_nl
    obs_data_now(i) = obs_data(obs_ind, i) + (mod_montime - (obs_ind-1)*obs_dt) * (obs_data(obs_ind+1, i) - obs_data(obs_ind, i)) / obs_dt
    ! If the obs_data_now(i) is not NAN, use it, otherwise, remove it from the data set.
    IF (obs_data_now(i) == obs_data_now(i)) THEN
      obs_nl_nonan = obs_nl_nonan + 1
      obs_lv_nonan(obs_nl_nonan) = obs_lv(i)
      obs_data_nonan(obs_nl_nonan) = obs_data_now(i)
    END IF
  END DO

  IF (obs_nl_nonan == 0) THEN
    nud_data = mod_data
  ELSE
    SELECT CASE ( TRIM(ADJUSTL(method)) )
    CASE ('linear')  ! linear interpolation
      DO i = 1, mod_kz
        ! Interpolate when the height is within the observation range
        IF ( mod_z(i) >= obs_lv_nonan(1) .and. mod_z(i) <= obs_lv_nonan(obs_nl_nonan) ) THEN
          CALL interp_1d(obs_nl_nonan, obs_lv_nonan(1:obs_nl_nonan), obs_data_nonan(1:obs_nl_nonan), 1, mod_z(i), nud_data(i))
          ! If nudging data is NAN, set it to model data
          IF ( ISNAN(nud_data(i)) ) nud_data(i) = mod_data(i)
        ! Set the nudging data to the model data if out of the observation range
        ELSE
          nud_data(i) = mod_data(i)
        END IF
      END DO
    CASE ('spline')  ! cubic spline interpolation
      CALL inter3( (/obs_zbound(1), obs_lv_nonan(1:obs_nl_nonan), obs_zbound(2)/), (/obs_databound(1), obs_data_nonan(1:obs_nl_nonan), obs_databound(2)/), obs_nl_nonan+2, &
        mod_z, mod_kz, nud_data)
      DO i = 1, mod_kz
        ! Only nudge the data within the observation height range
        IF ( mod_z(i) < obs_lv_nonan(1) .or. mod_z(i) > obs_lv_nonan(obs_nl_nonan) ) THEN
          nud_data(i) = mod_data(i)
        END IF
        ! Set nudging data to model data if it is NAN
        IF ( ISNAN(nud_data(i)) ) nud_data(i) = mod_data(i)
      END DO
    END SELECT
  END IF
END SUBROUTINE make_nudging_data_new


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MakeNudgingData
   !! Create variables used for nudging from the input data.
   !! By default assuming data form Hyytiälä (or measurement heights and number of those that correspond to Hyytiälä tower). 
   !! For other stations where numbers and heights of observations is different, edit this subroutine and make new subroutines 
   !! NudgingTemperatureNewStation, NudgingHumidityNewStation, and NudgignWindNewStation that correspond to your observational data.
   ! 
   !! includes rows 96-121 & 181-243 from original Scadis code
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE MakeNudgingData(n_levels, z_height, time_s, ind_now, station, dt_obs, temp_in, humidity_in, uwind_in, vwind_in, &
                              temperature, humidity_abs, u_wind, v_wind, &
                              temp_nudging, hum_nudging, uwind_nudging, vwind_nudging, wind_top)
      
      ! INPUT
      INTEGER,                            INTENT(IN   ) :: n_levels              , & ! (=kz) number of levels, defined in Sosa_data [1]
                                                           ind_now                   ! (=nxodrad) row index for current model time
      CHARACTER(len=3)                                  :: station                   ! name of the station, defined in Sosa_data.f90
      REAL(kind=dp),                      INTENT(IN   ) :: time_s                , & ! (=tmcontr) time [s] from the beginning of the month
                                                           wind_top              , & ! (=windx) top border wind speed [m/s]
                                                           dt_obs                    ! temporal resolution of surface observations (= input data) [s]
      REAL(kind=dp), DIMENSION(6,1489),   INTENT(IN   ) :: temp_in               , & ! (=tem) measured temperature [K] from input file hyy_mix
                                                           humidity_in               ! (=hum) measured humidity [kg/m3] from input file hyy_mix
      REAL(kind=dp), DIMENSION(5,1489),   INTENT(IN   ) :: uwind_in                  ! (=uwind) measured wind speed [m/s] from input file hyy_mix
      REAL(kind=dp), DIMENSION(5,1489),   INTENT(IN   ) :: vwind_in                  ! (=vwind) measured wind speed [m/s] from input file hyy_mix
      REAL(kind=dp), DIMENSION(n_levels), INTENT(IN   ) :: z_height              , & ! (=z) height of each (vertical) layer [m]
                                                           temperature           , & ! (=ta) air temperature [K]
                                                           humidity_abs          , & ! (=qa) absolute humidity of air [kg/m3]
                                                           u_wind                , & ! (=u) u-component of wind [m/s] 
                                                           v_wind                    ! (=v) v-component of wind [m/s] 
      ! OUTPUT
      REAL(kind=dp), DIMENSION(n_levels), INTENT(  OUT) :: temp_nudging         , & ! (=tnud) temperature [K] for nudging towards
                                                           hum_nudging          , & ! (=qnud) humidity [kg/m3] for nudging towards
                                                           uwind_nudging        , & ! (=unud) u-wind speed [m/s] for nudging towards
                                                           vwind_nudging            ! (=vnud) v-wind speed [m/s] for nudging towards

      if (station == 'WLG') then
         CALL NudgingWelgegund(n_levels, z_height, time_s, ind_now, dt_obs, temp_in(1,:), temperature, temp_nudging)
         CALL NudgingWelgegund(n_levels, z_height, time_s, ind_now, dt_obs, humidity_in(1,:), humidity_abs, hum_nudging)
         CALL NudgingWelgegund(n_levels, z_height, time_s, ind_now, dt_obs, uwind_in(1,:), u_wind, uwind_nudging)
         vwind_nudging = 0.
      else ! Hyytiälä
         CALL NudgingTemperature(n_levels, z_height, time_s, ind_now, dt_obs, temp_in, temperature, temp_nudging)
         CALL NudgingHumidity(n_levels, z_height, time_s, ind_now, dt_obs, humidity_in, humidity_abs, hum_nudging)
         CALL NudgignWind(n_levels, z_height, time_s, ind_now, dt_obs, uwind_in, vwind_in, u_wind, v_wind, uwind_nudging, &
                          vwind_nudging, wind_top)
      endif
      

   END SUBROUTINE MakeNudgingData

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE NudgingTemperature
   !! Create variables used for nudging temperature from the input data
   !! Assuming data form Hyytiälä (default station). For other stations where numbers and heights of observations is different, make
   !!  a new subroutine with similar content.
   ! 
   !! includes rows 105-112, 181-183, 190-191, 193-195, 197-199, 201-203, 205-207, 209, 224 from original Scadis code
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE NudgingTemperature(n_levels, z_height, time_s, ind_now, dt_obs, temp_in, temperature, temp_nudging)
      
      ! INPUT
      INTEGER,                            INTENT(IN   ) :: n_levels             , & ! (=kz) number of levels, defined in Sosa_data [1]
                                                           ind_now                  ! (=nxodrad) row index for current model time
      REAL(kind=dp),                      INTENT(IN   ) :: time_s               , & ! (=tmcontr) time [s] from the beginning of the month
                                                           dt_obs                   ! temporal resolution of surface observations (= input data) [s]
      REAL(kind=dp), DIMENSION(6,1488),   INTENT(IN   ) :: temp_in                  ! (=tem) measured temperature [K] from input file hyy_mix
      REAL(kind=dp), DIMENSION(n_levels), INTENT(IN   ) :: temperature          , & ! (=ta) air temperature [K]
                                                           z_height                 ! (=z) height of each (vertical) layer [m]
      ! OUTPUT
      REAL(kind=dp), DIMENSION(n_levels), INTENT(  OUT) :: temp_nudging             ! (=tnud) temperature [K] for nudging towards
      ! LOCAL 
      REAL(kind=dp)                                     :: ta67                  , & ! measured temperature [K] at 67.2 m
                                                           ta50                  , & ! measured temperature [K] at 50.4 m
                                                           ta33                  , & ! measured temperature [K] at 33.6 m
                                                           ta16                  , & ! measured temperature [K] at 16.8 m
                                                           ta8                   , & ! measured temperature [K] at 8.4 m
                                                           ta4                       ! measured temperature [K] at 4.2 m
      INTEGER                                           :: k                         ! foor loop


      ! data read from input file "hyy_mix.txt" columns 6-11
      ! ta67,ta50,ta33,ta16,ta8 and ta4 - temperature measured at 67.2, 50.4, 33.6, 16.8, 8.4 and 4.2 m

      ta67 = temp_in(1,ind_now) + (time_s-(ind_now-1)*dt_obs)*(temp_in(1,ind_now+1)-temp_in(1,ind_now))/dt_obs ! +1.8 !year 2050
      ta50 = temp_in(2,ind_now) + (time_s-(ind_now-1)*dt_obs)*(temp_in(2,ind_now+1)-temp_in(2,ind_now))/dt_obs ! +1.8 !year 2050 
      ta33 = temp_in(3,ind_now) + (time_s-(ind_now-1)*dt_obs)*(temp_in(3,ind_now+1)-temp_in(3,ind_now))/dt_obs ! +1.8 !year 2050
      ta16 = temp_in(4,ind_now) + (time_s-(ind_now-1)*dt_obs)*(temp_in(4,ind_now+1)-temp_in(4,ind_now))/dt_obs ! +1.8 !year 2050
      ta8  = temp_in(5,ind_now) + (time_s-(ind_now-1)*dt_obs)*(temp_in(5,ind_now+1)-temp_in(5,ind_now))/dt_obs  ! +1.8 !year 2050
      ta4  = temp_in(6,ind_now) + (time_s-(ind_now-1)*dt_obs)*(temp_in(6,ind_now+1)-temp_in(6,ind_now))/dt_obs  ! +1.8 !year 2050


      !105         !  INFORMATON FOR NUDGING FROM FILE "TEM" WITH DIMENSION(6,1488)
      !106         !   ta67,ta50,ta33,ta16,ta8 and ta4 - temperature measured at 67.2, 50.4, 33.6, 16.8, 8.4 and 4.2 m
      !
      !107         ta67=tem(1,nxodrad)+(tmcontr-(nxodrad-1)*1800.)*(tem(1,nxodrad+1)-tem(1,nxodrad))/1800.  !+1.8 !year 2050
      !108         ta50=tem(2,nxodrad)+(tmcontr-(nxodrad-1)*1800.)*(tem(2,nxodrad+1)-tem(2,nxodrad))/1800.  !+1.8 !year 2050 
      !109         ta33=tem(3,nxodrad)+(tmcontr-(nxodrad-1)*1800.)*(tem(3,nxodrad+1)-tem(3,nxodrad))/1800.  !+1.8 !year 2050
      !110         ta16=tem(4,nxodrad)+(tmcontr-(nxodrad-1)*1800.)*(tem(4,nxodrad+1)-tem(4,nxodrad))/1800.  !+1.8 !year 2050
      !111         ta8=tem(5,nxodrad)+(tmcontr-(nxodrad-1)*1800.)*(tem(5,nxodrad+1)-tem(5,nxodrad))/1800.   !+1.8 !year 2050
      !112         ta4=tem(6,nxodrad)+(tmcontr-(nxodrad-1)*1800.)*(tem(6,nxodrad+1)-tem(6,nxodrad))/1800.   !+1.8 !year 2050


      DO k=1,n_levels
      !181         do k=1,kz

         temp_nudging(k)=temperature(k)
         !182            !    NO ANY NUDGING YET
         !183            tnud(k)=ta(k)


         !    DEFINE NUDGING VARIABLES BY SIMPLE LINEAR APPROXIMATION

         IF ( z_height(k) > 50.4  .AND.  z_height(k) <= 67.2 ) THEN
            temp_nudging(k) = ta50 + (ta67-ta50) * (z_height(k)-50.4) / (67.2-50.4)
         END IF
         !190            if(z(k).gt.50.4.and.z(k).le.67.2) then
         !191               tnud(k)=ta50+(ta67-ta50)*(z(k)-50.4)/(67.2-50.4)
         !193            endif

         IF ( z_height(k) > 33.6  .AND.  z_height(k) <= 50.4 ) THEN
            temp_nudging(k) = ta33 + (ta50-ta33) * (z_height(k)-33.6) / (50.4-33.6)
         END IF
         !194            if(z(k).gt.33.6.and.z(k).le.50.4) then
         !195               tnud(k)=ta33+(ta50-ta33)*(z(k)-33.6)/(50.4-33.6)
         !197            endif

         IF ( z_height(k) > 16.8  .and.  z_height(k) <= 33.6) THEN
            temp_nudging(k) = ta16 + (ta33-ta16) * (z_height(k)-16.8) / (33.6-16.8)
         END IF
         !198            if(z(k).gt.16.8.and.z(k).le.33.6) then
         !199               tnud(k)=ta16+(ta33-ta16)*(z(k)-16.8)/(33.6-16.8)
         !201            endif

         IF (z_height(k) > 8.4  .AND.  z_height(k) <= 16.8 ) THEN
            temp_nudging(k) = ta8 + (ta16-ta8) * (z_height(k)-8.4) / (16.8-8.4)
         END IF
         !202            if(z(k).gt.8.4.and.z(k).le.16.8) then
         !203               tnud(k)=ta8+(ta16-ta8)*(z(k)-8.4)/(16.8-8.4)
         !205            endif

         IF( z_height(k) > 4.2  .AND.  z_height(k) <= 8.4 ) THEN
            temp_nudging(k) = ta4 + (ta8-ta4) * (z_height(k)-4.2) / (8.4-4.2)
         END IF
         !206            if(z(k).gt.4.2.and.z(k).le.8.4) then
         !207               tnud(k)=ta4+(ta8-ta4)*(z(k)-4.2)/(8.4-4.2)
         !209            endif

      END DO
      !224         enddo



   END SUBROUTINE NudgingTemperature


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE NudgingHumidity
   !! Create variables used for nudging humidity from the input data.
   !! Assuming data form Hyytiälä (default station). For other stations where numbers and heights of observations is different, make
   !!  a new subroutine with similar content.
   ! 
   !! Includes rows 113-120, 181-182, 184, 190, 192-194, 196-198, 200-202, 204-206, 208-209, 224 from original Scadis code
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE NudgingHumidity(n_levels, z_height, time_s, ind_now, dt_obs, humidity_in, humidity_abs, hum_nudging) 

      ! INPUT
      INTEGER,                            INTENT(IN   ) :: n_levels             , & ! (=kz) number of levels, defined in Sosa_data [1]
                                                           ind_now                  ! (=nxodrad) row index for current model time
      REAL(kind=dp),                      INTENT(IN   ) :: time_s               , & ! (=tmcontr) time [s] from the beginning of the month
                                                           dt_obs                   ! temporal resolution of surface observations (= input data) [s]
      REAL(kind=dp), DIMENSION(6,1488),   INTENT(IN   ) :: humidity_in              ! (=hum) measured humidty [kg/m3] from input file hyy_mix
      REAL(kind=dp), DIMENSION(n_levels), INTENT(IN   ) :: humidity_abs         , & ! (=qa) absolute humidity of air [kg/m3]
                                                           z_height                 ! (=z) height of each (vertical) layer [m]
      ! OUTPUT
      REAL(kind=dp), DIMENSION(n_levels), INTENT(  OUT) :: hum_nudging              ! (! (=qnud) humidity [kg/m3] for nudging towards
      ! LOCAL 
      REAL(kind=dp)                                     :: qa67                 , & ! measured humidity [kg/m3] at 67.2 m
                                                           qa50                 , & ! measured humidity [kg/m3] at 50.4 m
                                                           qa33                 , & ! measured humidity [kg/m3] at 33.6 m
                                                           qa16                 , & ! measured humidity [kg/m3] at 16.8 m
                                                           qa8                  , & ! measured humidity [kg/m3] at 8.4 m
                                                           qa4                      ! measured humidity [kg/m3] at 4.2 m
      INTEGER                                           :: k                        ! for loops


      ! data read from input file "hyy_mix.txt" columns 6-11
      ! qa67,qa50,qa33,qa16,qa8 and qa4 - humidity measured at 67.2, 50.4, 33.6, 16.8, 8.4 and 4.2 m

      qa67 = humidity_in(1,ind_now) + (time_s-(ind_now-1)*dt_obs)*(humidity_in(1,ind_now+1)-humidity_in(1,ind_now))/dt_obs
      qa50 = humidity_in(2,ind_now) + (time_s-(ind_now-1)*dt_obs)*(humidity_in(2,ind_now+1)-humidity_in(2,ind_now))/dt_obs
      qa33 = humidity_in(3,ind_now) + (time_s-(ind_now-1)*dt_obs)*(humidity_in(3,ind_now+1)-humidity_in(3,ind_now))/dt_obs
      qa16 = humidity_in(4,ind_now) + (time_s-(ind_now-1)*dt_obs)*(humidity_in(4,ind_now+1)-humidity_in(4,ind_now))/dt_obs
      qa8  = humidity_in(5,ind_now) + (time_s-(ind_now-1)*dt_obs)*(humidity_in(5,ind_now+1)-humidity_in(5,ind_now))/dt_obs
      qa4  = humidity_in(6,ind_now) + (time_s-(ind_now-1)*dt_obs)*(humidity_in(6,ind_now+1)-humidity_in(6,ind_now))/dt_obs

      !113         !  INFORMATON FOR NUDGING FROM FILE "HUM" WITH DIMENSION(6,1488)
      !114         !   qa67,qa50,qa33,qa16,qa8 and qa4 - moisture measured at 67.2, 50.4, 33.6, 16.8, 8.4 and 4.2 m
      !
      !115         qa67=hum(1,nxodrad)+(tmcontr-(nxodrad-1)*1800.)*(hum(1,nxodrad+1)-hum(1,nxodrad))/1800.
      !116         qa50=hum(2,nxodrad)+(tmcontr-(nxodrad-1)*1800.)*(hum(2,nxodrad+1)-hum(2,nxodrad))/1800.
      !117         qa33=hum(3,nxodrad)+(tmcontr-(nxodrad-1)*1800.)*(hum(3,nxodrad+1)-hum(3,nxodrad))/1800.
      !118         qa16=hum(4,nxodrad)+(tmcontr-(nxodrad-1)*1800.)*(hum(4,nxodrad+1)-hum(4,nxodrad))/1800.
      !119         qa8=hum(5,nxodrad)+(tmcontr-(nxodrad-1)*1800.)*(hum(5,nxodrad+1)-hum(5,nxodrad))/1800.
      !120         qa4=hum(6,nxodrad)+(tmcontr-(nxodrad-1)*1800.)*(hum(6,nxodrad+1)-hum(6,nxodrad))/1800.


      DO k=1,n_levels
      !181         do k=1,kz

         hum_nudging(k)=humidity_abs(k)
         !182            !    NO ANY NUDGING YET
         !184            qnud(k)=qa(k)


         !    DEFINE NUDGING VARIABLES BY SIMPLE LINEAR APPROXIMATION

         IF ( z_height(k) > 50.4  .AND.  z_height(k) <= 67.2 ) THEN
            hum_nudging(k) = qa50 + (qa67-qa50) * (z_height(k)-50.4) / (67.2-50.4)
         END IF
         !190            if(z(k).gt.50.4.and.z(k).le.67.2) then
         !192               qnud(k)=qa50+(qa67-qa50)*(z(k)-50.4)/(67.2-50.4)
         !193            endif

         IF ( z_height(k) > 33.6  .AND.  z_height(k) <= 50.4 ) THEN
            hum_nudging(k) = qa33 + (qa50-qa33) * (z_height(k)-33.6) / (50.4-33.6)
         END IF
         !194            if(z(k).gt.33.6.and.z(k).le.50.4) then
         !196               qnud(k)=qa33+(qa50-qa33)*(z(k)-33.6)/(50.4-33.6)
         !197            endif

         IF ( z_height(k) > 16.8  .AND.  z_height(k) <= 33.6 ) THEN
            hum_nudging(k) = qa16 + (qa33-qa16) * (z_height(k)-16.8) / (33.6-16.8)
         END IF
         !198            if(z(k).gt.16.8.and.z(k).le.33.6) then
         !200               qnud(k)=qa16+(qa33-qa16)*(z(k)-16.8)/(33.6-16.8)
         !201            endif

         IF ( z_height(k) > 8.4  .AND.  z_height(k) <= 16.8 ) THEN
            hum_nudging(k) = qa8 + (qa16-qa8) * (z_height(k)-8.4) / (16.8-8.4)
         END IF
         !202            if(z(k).gt.8.4.and.z(k).le.16.8) then
         !204               qnud(k)=qa8+(qa16-qa8)*(z(k)-8.4)/(16.8-8.4)
         !205            endif

         IF ( z_height(k) > 4.2  .AND.  z_height(k) <= 8.4 ) THEN
            hum_nudging(k) = qa4+(qa8-qa4) * (z_height(k)-4.2) / (8.4-4.2)
         END IF
         !206            if(z(k).gt.4.2.and.z(k).le.8.4) then
         !208               qnud(k)=qa4+(qa8-qa4)*(z(k)-4.2)/(8.4-4.2)
         !209            endif

      END DO
      !224         enddo


   END SUBROUTINE NudgingHumidity


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE NudgignWind
   !! Create variables used for nudging wind speed from the input data.
   !! Assuming data form Hyytiälä (default station). For other stations where numbers and heights of observations is different, make
   !!  a new subroutine with similar content.
   ! 
   !! Two methods for interpolating between measurement heights: linear and cubic spline (natural) interpolation. The quick and dirty
   !! way of selecting one of these: comment one of them out, leave the other one.
   ! 
   !! Includes rows 97-104, 181-182, 185-186, 211-243 from original Scadis code.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE NudgignWind(n_levels, z_height, time_s, ind_now, dt_obs, &
                          uwind_in, vwind_in, u_wind, v_wind, uwind_nudging, vwind_nudging, wind_top)

     ! INPUT
      INTEGER,                            INTENT(IN   ) :: n_levels              , & ! (=kz) number of levels, defined in Sosa_data [1]
                                                           ind_now                   ! (=nxodrad) row index for current model time
      REAL(kind=dp),                      INTENT(IN   ) :: time_s                , & ! (=tmcontr) time [s] from the beginning of the month
                                                           wind_top              , & ! (=windx) top border wind speed [m/s]
                                                           dt_obs                    ! temporal resolution of surface observations (= input data) [s]
      REAL(kind=dp), DIMENSION(5,1489),   INTENT(IN   ) :: uwind_in              ! (=speed) measured wind speed [m/s] from input file hyy_mix
      REAL(kind=dp), DIMENSION(5,1489),   INTENT(IN   ) :: vwind_in              ! (=speed) measured wind speed [m/s] from input file hyy_mix
      REAL(kind=dp), DIMENSION(n_levels), INTENT(IN   ) :: z_height              , & ! (=z) height of each (vertical) layer [m]
                                                           u_wind                , & ! (=u) u-component of wind [m/s] 
                                                           v_wind                    ! (=v) v-component of wind [m/s]                            
      ! OUTPUT
      REAL(kind=dp), DIMENSION(n_levels), INTENT(OUT  ) :: uwind_nudging        , & ! (=unud) u-wind speed [m/s] for nudging towards
                                                           vwind_nudging            ! (=vnud) v-wind speed [m/s] for nudging towards
      ! LOCAL 
      REAL(kind=dp)                                     :: u74,v74                  , & ! measured wind speed [m/s] at 74 m
                                                           u33,v33                  , & ! measured wind speed [m/s] at 33.6 m
                                                           u16,v16                  , & ! measured wind speed [m/s] at 16.8 m
                                                           u8,v8                   , & ! measured wind speed [m/s] at 8.4 m
                                            uwind_output(n_levels),vwind_output(n_levels)    
! (=ffn) output from interp3 function, wind speed
                                                                                    ! to nudge toweards for full column
      INTEGER                                           :: k                        ! foor loop


      ! data read from input file "hyy_mix.txt" columns 6-11
      ! u74,u33,u16,u8 - wind speed measured at 74, 33.6, 16.8 and 8.4 m
      u74 = uwind_in(1,ind_now) + (time_s-(ind_now-1)*dt_obs)*(uwind_in(1,ind_now+1)-uwind_in(1,ind_now))/dt_obs
      u33 = uwind_in(2,ind_now) + (time_s-(ind_now-1)*dt_obs)*(uwind_in(2,ind_now+1)-uwind_in(2,ind_now))/dt_obs
      u16 = uwind_in(3,ind_now) + (time_s-(ind_now-1)*dt_obs)*(uwind_in(3,ind_now+1)-uwind_in(3,ind_now))/dt_obs
      u8  = uwind_in(4,ind_now) + (time_s-(ind_now-1)*dt_obs)*(uwind_in(4,ind_now+1)-uwind_in(4,ind_now))/dt_obs

      v74 = vwind_in(1,ind_now) + (time_s-(ind_now-1)*dt_obs)*(vwind_in(1,ind_now+1)-uwind_in(1,ind_now))/dt_obs
      v33 = vwind_in(2,ind_now) + (time_s-(ind_now-1)*dt_obs)*(vwind_in(2,ind_now+1)-uwind_in(2,ind_now))/dt_obs
      v16 = vwind_in(3,ind_now) + (time_s-(ind_now-1)*dt_obs)*(vwind_in(3,ind_now+1)-uwind_in(3,ind_now))/dt_obs
      v8  = vwind_in(4,ind_now) + (time_s-(ind_now-1)*dt_obs)*(vwind_in(4,ind_now+1)-uwind_in(4,ind_now))/dt_obs

      !97          !  INFORMATON FOR NUDGING FROM FILE "SPEED" WITH DIMENSION(5,1488)
      !98          !   u74,u33,u16,u8 - wind speed measured at 74, 33.6, 16.8 and 8.4 m
      !99          !   b23 - TKE recalculated from U* measured at 23.3 m 
      !100         u74=speed(1,nxodrad)+(tmcontr-(nxodrad-1)*1800.)*(speed(1,nxodrad+1)-speed(1,nxodrad))/1800.
      !101         u33=speed(2,nxodrad)+(tmcontr-(nxodrad-1)*1800.)*(speed(2,nxodrad+1)-speed(2,nxodrad))/1800.
      !102         u16=speed(3,nxodrad)+(tmcontr-(nxodrad-1)*1800.)*(speed(3,nxodrad+1)-speed(3,nxodrad))/1800.
      !103         u8=speed(4,nxodrad)+(tmcontr-(nxodrad-1)*1800.)*(speed(4,nxodrad+1)-speed(4,nxodrad))/1800.
      !104         !bt23=( speed(5,nxodrad)+(tmcontr-(nxodrad-1)*1800.)*(speed(5,nxodrad+1)-speed(5,nxodrad))/1800.)**2/(cc2**0.5)   

      DO k=1,n_levels
      !181         do k=1,kz
         uwind_nudging(k) = u_wind(k)
         vwind_nudging(k) = v_wind(k)
         
         !182            !    NO ANY NUDGING YET
         !185            unud(k)=u(k)
         !186            vnud(k)=v(k)


         ! LINEAR INTERPOLATION
         ! comment/uncomment rows between *:s to (not) use linear interpolation
         ! * * * * * * * * * * * * 
 !         
 !        IF ( z_height(k) > 33.6  .AND.  z_height(k) <= 74.0 ) THEN
 !           uwind_nudging(k) = u33 + (u74-u33) * (z_height(k)-33.6) / (74.0-33.6)
 !           vwind_nudging(k) = 0.0
 !        END IF
 !        !211            if(z(k).gt.33.6.and.z(k).le.74.) then
 !        !212               unud(k)=u33+(u74-u33)*(z(k)-33.6)/(74.-33.6)
 !        !213               vnud(k)=0.
 !        !214            endif
 !        
 !        IF ( z_height(k) > 16.8  .AND.  z_height(k) <= 33.6 ) THEN
 !           uwind_nudging(k) = u16 + (u33-u16) * (z_height(k)-16.8) / (33.6-16.8)
 !           vwind_nudging(k) = 0.0
 !        END IF      
 !        !215            if(z(k).gt.16.8.and.z(k).le.33.6) then 
 !        !216               unud(k)=u16+(u33-u16)*(z(k)-16.8)/(33.6-16.8)
 !        !217               vnud(k)=0.
 !        !218            endif
 !        
 !        IF ( z_height(k) > 8.4  .AND.  z_height(k) <= 16.8 ) THEN
 !           uwind_nudging(k) = u8 + (u16-u8) * (z_height(k)-8.4) / (16.8-8.4)
 !           vwind_nudging(k) = 0.0
 !        END IF  
 !        !219            if(z(k).gt.8.4.and.z(k).le.16.8) then
 !        !220               unud(k)=u8+(u16-u8)*(z(k)-8.4)/(16.8-8.4)
 !        !221               vnud(k)=0.
 !        !222            endif
 !        !223 600        continue   
 !               
         ! * * * * * * * * * * * * 
         
      END DO
      !224         enddo  

    
      ! CUBIC SPLINE INTERPOLATION
      ! comment/uncomment rows between *:s to (not) use cubic spline interpolation
      ! * * * * * * * * * * * * 
      
      CALL inter3( dble((/ 0.0, 8.4, 16.8, 36.6, 74.0, 3000.0/)) , &
                   (/0.0d0, u8, u16, u33, u74, wind_top/)      , &
                   6, z_height, n_levels, uwind_output)
      
      CALL inter3( dble((/ 0.0, 8.4, 16.8, 36.6, 74.0, 3000.0/)) , &
                   (/0.0d0, u8, u16, u33, u74, wind_top/)      , &
                   6, z_height, n_levels, vwind_output)
  
      DO k = 1,n_levels
         IF(  z_height(k) > 8.4  .AND.  z_height(k) <= 75.0 ) THEN
            uwind_nudging(k) = uwind_output(k)
            vwind_nudging(k) = vwind_output(k)
         END IF
      END DO
        
      ! * * * * * * * * * * * *  
      
      !225         zini(1)=0.
      !226         zini(2)=8.4
      !227         zini(3)=16.8
      !228         zini(4)=36.6
      !229         zini(5)=74.
      !230         zini(6)=3000.
      !
      !231         fini(1)=0.0
      !232         fini(2)=u8
      !233         fini(3)=u16
      !234         fini(4)=u33
      !235         fini(5)=u74
      !236         fini(6)=windx
      !
      !237         CALL inter3(zini,fini,6,z,kz,ffn)
      !238         do k=1,kz
      !
      !239            if(z(k).gt.8.4.and.z(k).le.75.) then
      !240               unud(k)=ffn(k)
      !241               vnud(k)=0.
      !242            endif
      !
      !243         enddo


   END SUBROUTINE NudgignWind
 
 
   ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   ! SUBROUTINE NudgingWelgegund
   !! Create variables used for nudging temperature/humidity/wind speed from the input data.
   !! This subroutine is for Welgegund (South-Africa) data.
   ! 
   ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   SUBROUTINE NudgingWelgegund(n_levels, z_height, time_s, ind_now, dt_obs, data_in, model_data, data_nudging)

      ! INPUT
      INTEGER,                            INTENT(IN   ) :: n_levels             , & ! (=kz) number of levels, defined in Sosa_data [1]
                                                           ind_now                  ! (=nxodrad) row index for current model time
      REAL(kind=dp),                      INTENT(IN   ) :: time_s               , & ! (=tmcontr) time [s] from the beginning of the month
                                                           dt_obs                   ! temporal resolution of surface observations (= input data) [s]
      REAL(kind=dp), DIMENSION(1488),     INTENT(IN   ) :: data_in                  ! measured values
      REAL(kind=dp), DIMENSION(n_levels), INTENT(IN   ) :: model_data           , & ! model values for parameter corresponding input data
                                                           z_height                 ! (=z) height of each (vertical) layer [m]
      ! OUTPUT
      REAL(kind=dp), DIMENSION(n_levels), INTENT(  OUT) :: data_nudging             ! the real world field: data to nudge towards
      ! LOCAL 
      REAL(kind=dp)                                     :: data_now                 ! observed data interpolated to current model time
      INTEGER                                           :: k                        ! for loops



      ! data read from input file "hyy_mix.txt"
      data_now   = data_in(ind_now) + (time_s-(ind_now-1)*dt_obs)*(data_in(ind_now+1)-data_in(ind_now))/dt_obs

      
      data_nudging = model_data

      DO k=1,n_levels
      
         !! measurement height 5 m -> give the observed value to layers 5 +- 0.5 m
         !! (with 3km domain having 51 layers this includes only one model level)

         IF ( z_height(k) > 4.5  .AND.  z_height(k) <= 5.5 ) THEN
            data_nudging(k) = data_now
         END IF
       
      END DO

      
   END SUBROUTINE NudgingWelgegund
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! These subroutines and functions copied "as is"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    
    
    
 SUBROUTINE inter3(z,f,kz,z1,kz1,f2)
    !c     Subroutine inter3 takes in a 3d-array of data (f), and a 2d-array
    !c     describing the height of terrain (rou). Then it interpolates
    !c     the data in z-dimension and produces a new dataset (f2).
    !c     The new, interpolated, dataset is in absolute coordinate system.
    !c     Cubic spline (natural) interpolation is used.
    !c
    !c     z(1:kz)
    !c       z-coordinates of gridpoints in terrain-following system
    !c     f(1:kx,1:ky,1:kz)
    !c       values of some function in terrain-following system
    !c     rou(1:kx,1:ky)
    !c       ground height in the grid nodes
    !c     kx,ky
    !c       x- and y-dimensions of the 3d-arrays
    !c     kz
    !c       z-dimension of the terrain-following data
    !c     kz2
    !c       z-dimension of the new data that we shall produce by interpolation
    !c     z2(1:kz2)
    !c       z-coordinates of the new gridpoints
    !c     f2(1:kx,1:ky,1:kz2)
    !c       here we will put the new, interpolated, values of the function        
    IMPLICIT NONE
    INTEGER :: kz,kz1
    REAL(kind=dp) :: z(kz),f(kz),z1(kz1),f2(kz1)
    INTEGER :: k
    REAL(kind=dp) :: fvals(kz),fvals2(kz1),h

    DO k=1,kz
       fvals(k) = f(k)
    ENDDO
    h = 0.
    CALL inter1(z,fvals,kz,h,z1,kz1,fvals2)
    DO k=1,kz1
       f2(k) = fvals2(k)
    ENDDO
    RETURN
  END SUBROUTINE inter3
  
    SUBROUTINE inter1(za,fa,n,z0,z2a,n2,f2a)
    IMPLICIT NONE
    INTEGER :: n, n2
    REAL(kind=dp) :: za(n),fa(n),z0,z2a(n2),f2a(n2)
    INTEGER :: i
    REAL(kind=dp) :: derivs(n), hh, resul
    CALL splinh(za,fa,n,derivs)
    DO i=1,n2
       IF (z2a(i) .LT. z0) THEN
          f2a(i) = 0.0
       ELSE
          !            IF (z2a(i) .gt. za(n)) PAUSE 'inter1 does not extrapolate'
          hh = z2a(i)-z0
          CALL splint(za,fa,derivs,n,hh,resul)
          f2a(i) = resul
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE inter1
  
  
  SUBROUTINE splint(xa,ya,y2a,n,x,y)
    IMPLICIT NONE
    INTEGER :: n
    REAL(kind=dp) :: xa(n),ya(n),y2a(n),x,y
    INTEGER :: k, klo, khi
    REAL(kind=dp) :: a,b,hs
    klo = 1
    khi = n
1005 IF (khi-klo .GT. 1) THEN
       k = (khi+klo)/2
       IF (xa(k) .GT. x) THEN
          khi = k
       ELSE
          klo = k
       ENDIF
       GOTO 1005
    ENDIF
    hs = xa(khi)-xa(klo)
    if (hs .eq. 0.0) STOP 'bad xa input in splint'
    a = (xa(khi)-x)/hs
    b = (x-xa(klo))/hs
    y = a*ya(klo) + b*ya(khi) +&
         ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi)) * (hs**2)/6.
    RETURN
  END SUBROUTINE splint

  SUBROUTINE splinh(x,y,n,y2)
    IMPLICIT NONE
    INTEGER :: n
    REAL(kind=dp) :: x(n),y(n),y2(n)
    INTEGER, PARAMETER :: NMAX=1000
    INTEGER :: i,k
    REAL(kind=dp) :: p,sig,uu(NMAX)
    y2(1) = 0.
    uu(1) = 0.
    DO i=2,n-1
       sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
       p = sig*y2(i-1)+2.
       y2(i) = (sig-1.)/p
       uu(i) = (6.*((y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))&
            /(x(i)-x(i-1)))/(x(i+1)-x(i-1)) - sig*uu(i-1))/p
    ENDDO
    y2(n) = 0.
    DO k=n-1,1,-1
       y2(k)=y2(k)*y2(k+1)+uu(k)
    ENDDO
    RETURN
  END SUBROUTINE splinh


!==========================================================================================================!
! Linear interpolation for 1D array
!==========================================================================================================!
SUBROUTINE interp_1d(nd, xd, yd, ni, xi, yi)
  INTEGER, INTENT(IN) :: nd  ! number of data points
  INTEGER, INTENT(IN) :: ni  ! number of interpolation points

  REAL(dp), INTENT(IN) :: xd(nd), yd(nd)  ! input data points and values, should be monotonic
  REAL(dp), INTENT(IN) :: xi(ni)  ! interpolation points
  REAL(dp), INTENT(OUT) :: yi(ni)  ! interpolation values

  INTEGER :: i, k  ! loop integers
  REAL(dp) :: t
  REAL(dp) :: xd_incr(nd)
  REAL(dp) :: yd_incr(nd)

  ! Check number of data points
  IF (nd < 2) THEN
    WRITE(*,*) 'Interpolation can not be done because there is only one available data point.'
    RETURN
  END IF

  ! Make sure x values are monotonically increasing
  IF (xd(2) > xd(1)) THEN
    xd_incr = xd
    yd_incr = yd
  ELSE
    xd_incr = xd(nd:1:-1)
    yd_incr = yd(nd:1:-1)
  END IF

  ! Default values
  yi(1:ni) = 0.0_dp

  ! When there is only one data point, the interpolation values are equal to that only value
  IF (nd == 1) THEN
    yi(1:ni) = yd_incr(1)
    RETURN
  END IF

  DO i = 1, ni
    ! The point is beyond the left boundary
    IF (xi(i) <= xd_incr(1)) THEN
      ! Extrapolation
      t = (xi(i) - xd_incr(1)) / (xd_incr(2) - xd_incr(1))
      yi(i) = (1.0_dp - t)*yd_incr(1) + t*yd_incr(2)
    ! The point is beyond the right boundary
    ELSE IF (xd_incr(nd) <= xi(i)) THEN
      t = (xi(i) - xd_incr(nd-1)) / (xd_incr(nd) - xd(nd-1))
      yi(i) = (1.0_dp - t)*yd_incr(nd-1) + t*yd_incr(nd)
    ! The point is within the domain
    ELSE
      DO k = 2, nd
        IF (xd_incr(k-1) <= xi(i) .and. xi(i) <= xd_incr(k)) THEN
          t = (xi(i) - xd_incr(k-1)) / (xd_incr(k) - xd_incr(k-1))
          yi(i) = (1.0_dp - t)*yd_incr(k-1) + t*yd_incr(k)
          EXIT
        END IF
      END DO
    END IF
  END DO
END SUBROUTINE interp_1d
    
END MODULE Scadis_Data_Tools
