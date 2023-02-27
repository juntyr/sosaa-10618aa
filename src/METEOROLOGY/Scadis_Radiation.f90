!==============================================================================!
!
! Module Scadis_Radiation       
!> @defgroup scadis_radiation                         
!> All parts of SCADIS(1D) related to radiation       
!
!> For more information see MT_MainMet.f90
!>                                                        
!> Commented and restructured by Rosa Gierens, University of Helsinki         
!              
!==============================================================================!

MODULE Scadis_Radiation

  ! variables from other modules

  use constants_mod, only: dp, PI, PIx2, c432, c864, solar_constant, deg2rad, rad2deg
  use Scadis_Parameters

  ! subroutines and variables part of the meteorology scheme
  ! USE Scadis_Tools

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: get_solar_declination, get_cos_zenith, &
    get_sunrise_and_sunset_time_UTC, &
    RadiationIntoATM, TotalSWRAboveCanopy, SWRDirectDiffuse, &
    PenetrationFDirect, PenetrationFDiffuse, IntegralFLeafOrientation


CONTAINS


  !============================================================================!
  !
  ! Solar declination could be positive (spring and summer) or negative (autumn
  ! and winter), and since the declination is the same for the whole Earth,
  ! julian day should use the UTC date.
  !
  ! Method 1
  !
  !   The formula is modified from the equation in Appendix A in Sogachev 2002,
  !   they said it referred to Tjernstrom (o with two dots) 1989, but I did not
  !   find the same formula in Tjernstrom 1989. This formula is also used in
  !   previous SOSAA versions.
  !
  ! Method 2
  !
  !   The formula refers to Eq. 2.5 in Stull (2011), and we assume the summer
  !   solstice on June 21st (172), and do not consider the leap year now. So it
  !   is an approximation but enough for our simulations.
  !
  !   Stull 2011: https://www.eoas.ubc.ca/books/Practical_Meteorology/mse3.html
  !
  ! Method 3
  !
  !   From TM5
  !
  ! Summary
  !
  !   It seems the first method is the most accurate according to the comparison
  !   with the results calculated in:
  !   https://gml.noaa.gov/grad/solcalc/azel.html
  !
  !============================================================================!
  function get_solar_declination(julian, method) result(deltas)
    ! Input
    real(dp) :: julian
    integer, optional :: method

    ! Output
    real(dp) :: deltas  ! [deg]

    ! Local
    real(dp), parameter :: phir = 23.44d0
    real(dp), parameter :: julian_solstice = 172.0d0
    integer :: method_

    ! Set default value of method to 1
    if (present(method)) then
      method_ = method
    else
      method_ = 1
    end if

    select case (method_)
    case (1)
      deltas = 0.409d0 * sin(PIx2 *(julian - 79.8d0)/365.24d0)
    case (2)
      deltas = phir*deg2rad * cos( PIx2 * (julian - julian_solstice)/365.0d0 )
    case (3)
      deltas = asin( sin(phir*deg2rad) * sin(4.88d0 + 2.0d0*PI*real(julian)/365.0d0) )
    case default
      write(*,*) 'Wrong method value.'
    end select
  end function get_solar_declination


  !============================================================================!
  ! Cosine of solar zenith angle is negative at night (after sunset)
  !
  ! The formula refers to Eq. 2.6 in Stull (2011).
  !
  ! Stull 2011: https://www.eoas.ubc.ca/books/Practical_Meteorology/mse3.html
  !============================================================================!
  function get_cos_zenith(time_UTC, julian_UTC, lat, lon) result(cos_zenith)
    ! Input
    real(dp) :: time_UTC  ! [s], UTC time of day, 0...86400
    real(dp) :: julian_UTC    ! [day], day of year, 1...365
    real(dp) :: lat  ! [deg], [-90, 90]
    real(dp) :: lon  ! [deg], [-180, 180]
                                                 
    ! Output
    real(dp) :: cos_zenith  ! [-]

    ! Local
    real(dp) :: declination  ! [rad], solar declination
    real(dp) :: solar_hour   ! [rad], solar hour angle

    declination = get_solar_declination(julian_UTC) * deg2rad
    solar_hour = PIx2 * time_UTC/c864 + lon*deg2rad

    ! Notice here minus sign should be used
    cos_zenith = sin(lat*deg2rad) * sin(declination) &
                 - cos(lat*deg2rad) * cos(declination) * cos(solar_hour)
  end function get_cos_zenith


  !============================================================================!
  ! Calculate sunrise and sunset time in UTC time
  ! All the tables and equations refer to Chapter 2 in Stull (2017, v102b) if
  ! not specified explicitly.
  !
  ! Elevation angles for diurnal events PSI (Table 2-2)
  ! The values are: deg (rad)
  ! Sunrise and sunset
  ! - Geometric:  0     ( 0      )
  ! - Apparent : -0.833 (-0.01454)
  ! Twilight
  ! - Civil       :  -6 (-0.10472)
  ! - Military    : -12 (-0.20944)
  ! - Astronomical: -18 (-0.31416)
  !============================================================================!
  subroutine get_sunrise_and_sunset_time_UTC(julian, lat, lon, &
      sunrise, sunset)
    ! Input
    real(dp), intent(in) :: julian       ! [day], day of year
    real(dp), intent(in) :: lat, lon     ! [deg], latitude and longitude

    ! Output
    real(dp), intent(out) :: sunrise      ! [s], sunrise time in UTC
    real(dp), intent(out) :: sunset       ! [s], sunset time in UTC

    ! Local parameters
    real(dp), parameter :: P     = 365.256363d0  ! [day], orbital period
    real(dp), parameter :: day_p = 4.0d0         ! [day], perihelion day
    real(dp), parameter :: PSI   = -0.833        ! [deg], elevation angle
    real(dp), parameter :: a = 7.659d0*60.0d0    ! [s]
    real(dp), parameter :: b = 9.863d0*60.0d0    ! [s]
    real(dp), parameter :: c = 3.588d0           ! [rad], = 205.58 [deg]

    ! Local variables
    real(dp) :: declination  ! [deg], solar declination
    real(dp) :: M  ! [rad], mean anomaly
    real(dp) :: delta_ta  ! [s], equation of time
    real(dp) :: terma  ! tmp term

    ! Solar declination
    declination = get_solar_declination(julian) * deg2rad

    ! Eq. 2.2
    M = PIx2 * (julian - day_p) / P

    ! Eq. 2.8a
    terma = ( sin(lat*deg2rad)*sin(declination*deg2rad) - sin(PSI*deg2rad) ) &
            / ( cos(lat*deg2rad)*cos(declination*deg2rad) )

    !>> Sunrise
    sunrise = c864 / PIx2 * (-lon*deg2rad + acos(terma))
    if (sunrise < 0) then
      sunrise = sunrise + c864
    end if

    !>> Sunset
    sunset  = c864 / PIx2 * (-lon*deg2rad - acos(terma))
    if (sunset < 0) then
      sunset = sunset + c864
    end if

    ! Correct the time with the tilted, elliptical orbit of the earth (Eq. 2.8b
    ! and Eq. 2.8c)
    delta_ta = - a*sin(M) + b*sin(2.0d0*M + c)

    sunrise = sunrise - delta_ta
    sunset  = sunset  - delta_ta
  end subroutine get_sunrise_and_sunset_time_UTC


  !============================================================================!
  ! Subroutine RadiationIntoATM
  !> Calculating the arrival of total solar radiation at the top of the atmosphere.
  !>
  !> 
  !>    Incoming radiation outside the atmosphere (\f$ Q_{0} \f$), which is also maximum possible radiation inside the model domain,
  !>       \f[  
  !>          Q_{0} = I_{s} \sin{h_{s}}                                                                                \qquad (1)
  !>       \f]
  !>    (equation A1 in Sogachev et al. 2002), where 
  !>       \f[ 
  !>          \sin{h_{s}} = \cos{\theta} = \sin{\phi} \sin{\delta} +\cos{\phi} \cos{\delta} \cos{t_{s}}                \qquad (2)
  !>       \f]
  !>    and \f$ \cos{\theta} \f$ is not allowed to exceed the value 0.01.  \f$ \delta \f$ can be calculated from 
  !>       \f[ 
  !>          \delta = 0.4145 \sin{ \left( \frac{2 \pi}{ 365.24} J  \right)}                                           \qquad (3)
  !>       \f]
  !>    (Tjernstöm, 1989) as it appears in the model description (Sogachev et al.2002), or by
  !>       \f[   
  !>          \delta = 0.409 \sin{ \left(\frac{2 \pi}{ 365.24} (Julian - 79.8) \right)}                                \qquad (4)
  !>       \f]
  !>    as it appears in the model with reference to Henderson-Sellers. The solar hour angle \f$ t_{s} \f$ is calculated by 
  !>       \f[ 
  !>          t_{s} = \pi \left(\frac{t_{in day}}{12 \mathrm{h} \times 3600\mathrm{s/h}} - 1.  \right)   ,             \qquad (5)
  !>       \f]
  !>    don't know why.
  !>
  !> Reference: Appendix A in Sogachev et al.2002
  !>
  !> Includes rows 11-13 and 16-17 from original Scadis code.
  !>
  !>
  !>   variable            | model variable       | description
  !>  -------------------- | -------------------- | -------------------------------------------------
  !> \f$  Q_{0}        \f$ | radiation_into_earth | incoming radiation outside the atmosphere (W/m^2)
  !> \f$  I_{s}        \f$ | SOLAR_CONSTANT       | solar constant (W/m^2)
  !> \f$  h_{s}        \f$ |  -                   | solar height angle
  !> \f$  \theta       \f$ |  -                   | solar zenith angle
  !> \f$  \cos{\theta} \f$ | cos_zenith           | cos(solar zenith angle) = sin(sun height angle)
  !> \f$  \phi         \f$ | latitude             | geographic latitude
  !> \f$  \delta       \f$ | declination          | solar declination
  !> \f$  t_{s}        \f$ | solar_hour           | solar hour angle
  !> \f$  J            \f$ |  -                   | calendar day with J = 1 at March 23 
  !> \f$  Julian       \f$ | julian_day           | Julian day (day of year), 1...365
  !> \f$  t_{in day}   \f$ | time_of_day          | time of the day, 0...86400 (s)
  !> 
  !
  !> @ingroup scadis_radiation
  !
  !============================================================================!
  SUBROUTINE RadiationIntoATM(time_of_day, julian_day, lat_rad, &
      declination, cos_zenith, radiation_into_earth)

      ! INPUT
      REAL(dp), INTENT(IN)  :: time_of_day, &  ! [s], time of day, 0...86400
                               lat_rad         ! [rad], latitude
      INTEGER,  INTENT(IN)  :: julian_day      ! [day], Julian day (day of year), 1...365
                                                   
      ! OUTPUT
      REAL(dp), INTENT(OUT) :: declination, &  ! (=wsun) Solar declination
                               cos_zenith , &  ! (=cos_zenith) cos(zenith angle)
                               radiation_into_earth  ! (=rads) incoming radiation outside the atmosphere [W/m2]
                                                   
      ! LOCAL
      REAL(dp)              :: solar_hour      ! (=zeit) solar hour angle


      ! solar declination, after Tjernström (1989)
      ! r_julian_day - 79.8 should be equal to calendar day (J), with 1 at March 23 (according to Sogachev et al.2002)
      ! (Equation 4)
      declination = 0.409d0  * sin(PIx2 *(julian_day - 79.8d0)/365.24d0) ! Henderson-Sellers
 
      ! solar hour angle, don't know why calculated like this
      ! (Equation 5)
      solar_hour=PI*(time_of_day/c432-1.0d0)

      ! sin(sun height angle) = cos(zenith angle) = ...  (Equation below eq. A1 in Sogachev et al.2002)
      ! (Equation 2)
      cos_zenith = sin(lat_rad) * sin(declination) &
                 + cos(lat_rad) * cos(declination) * cos(solar_hour)

      ! cos_zenith = max(1.0d-2,cos_zenith)

      ! incoming radiation outside the atmosphere [W/m2]
      ! (equation 1)
      radiation_into_earth = SOLAR_CONSTANT * cos_zenith   !max possible radiation
      
   END SUBROUTINE RadiationIntoATM



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine TotalSWRAboveCanopy
   !>    Calculation of arrival of total radiation at the top of the canopy: which part of the radiation at the 
   !>    top of the atosphere penetrates trough the atm and reaches the canopy.
   !
   !>    The amount of solar radiation at the top of the canopy (lowest model level above the canopy, \f$ h \f$) \f$ Q_{h} \f$ is 
   !>    calculated as  
   !>       \f[
   !>          Q_{h} = Q_{0} (0.3 + 0.7 \tau^{ m} )(1 - 0.8n_{Low} - 0.5n_{Mid} - c_{High} n_{High} )                   \qquad (1)   
   !>       \f]
   !>    (equation A2 in Sogachev et al. 2002), where \f$ \tau \f$ is atmospheric transmittance, \f$ m \f$ the optical mass of the
   !>    atmosphere, and \f$ n_{i} \f$ the cloud cover fraction for low, middle and high clouds. The empirical coefficient 
   !>    \f$ c_{High} \f$ depends on zenith angle of the sun:   
   !>       \f{eqnarray*}{
   !>          c_{High} = \left\{ \begin{array}{rclll}
   !>          0.2 & \mbox{for} & h_{s} > 20^{\mathrm{o}}, & \mbox{or} \cos{\theta} < 0.3420                          \qquad (2) \\
   !>          0.4 & \mbox{for} &  h_{s} <= 20^{\mathrm{o}}, & \mbox{or} \cos{\theta} >= 0.3420 \\
   !>          \end{array}\right.
   !>       \f}
   !>    The term \f$ x = (1 - 0.8n_{Low} - 0.5n_{Mid} - c_{High} n_{High} ) \f$ has been set to have a minimum value of 0.2. The  
   !>    optical mass \f$ m \f$ is calculated by 
   !>       \f[ 
   !>          m = 796(\sqrt{\sin^2{h_{s}}  + 0.002514 - \sin h_{s}  })                                                 \qquad (3)
   !>       \f]
   !>    or, since \f$ \sin{h_{s}} = \cos{\theta} \f$, 
   !>       \f[ 
   !>          m = 796(\sqrt{\cos^2{\theta}  + 0.002514 - \cos \theta })	.                                            \qquad (4)
   !>       \f]
   !>    In the model \f$ m \f$ has been set to have a maximum value of 6. \f$\tau \f$, \f$ n_{Low} \f$, \f$ n_{Mid} \f$, 
   !>    \f$ n_{High} \f$ are  constant parameters given in Scadis_Parameters.f90. 
   !> 
   !> Reference: Appendix A in Sogachev et al.2002
   !>
   !> Includes rows 14-15, 22, 26, 28 and 30 from original Scadis code.
   !>
   !>
   !>   variable            | model variable        | description
   !>  -------------------- | --------------------- | -------------------------------------------------
   !> \f$  Q_{h}        \f$ | total_global          | amount of solar radiation above the canopy (W/m^2)
   !> \f$  h            \f$ | canopy_top_level      | height of the canopy: lowest model level above the canopy 
   !> \f$  Q_{0}        \f$ | radiation_into_earth  | incoming radiation outside the atmosphere (W/m^2)
   !> \f$  \tau         \f$ | TRANSMITTANCE         | atmospheric transmittance
   !> \f$  m            \f$ | optical_mass          | optical mass of atmosphere
   !> \f$  n_{Low}      \f$ | CLOUDCOVER_LOW        | cloud cover fraction of low clouds
   !> \f$  n_{Mid}      \f$ | CLOUDCOVER_MIDDLE     | cloud cover fraction of middle clouds
   !> \f$  n_{High}     \f$ | CLOUDCOVER_HIGH       | cloud cover fraction of high clouds   
   !> \f$  c_{High}     \f$ | highcloud_coefficient | empirical constant
   !> \f$  h_{s}        \f$ |  -                    | solar height angle
   !> \f$  x            \f$ | cloud_effect          | help variable
   !> \f$  \theta       \f$ |  -                    | solar zenith angle
   !> \f$  \cos{\theta} \f$ | cos_zenith            | cos(solar zenith angle) = sin(sun height angle)
   !>
   !
   !> @ingroup scadis_radiation
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE TotalSWRAboveCanopy(cos_zenith, radiation_into_earth, total_global)
   
      ! INPUT
      REAL(kind=dp),                INTENT(IN)  :: cos_zenith              , & ! (=cos_zenith) cos(zenith angle) = sin(sun height angle)
                                                   radiation_into_earth       ! (=rads) incoming radiation outside the

      ! OUTPUT
      REAL(kind=dp),                INTENT(OUT) :: total_global               ! (=rads2) total global radiation at the top of the  canopy

      ! LOCAL
      REAL(kind=dp)                             :: optical_mass           , & ! (=optmas) optical mass of the atmosphere
                                                   cloud_effect           , & ! (=trancl) help variable, describes the effect of clouds on radiation penetrating trought the atmosphere
                                                   highcloud_coefficient      ! (=ct) coefficient for high clouds
                                                   
      ! optical mass of the atmosphere
      ! HEY!!! in the paper the last term is also inside the squer root!!
      ! (equation 9) 
      optical_mass=1.0d0*796.0d0*(sqrt(cos_zenith**2+0.0025d0)-cos_zenith)  ! (according to Sogachev et al.2002, the 0.0025 should be 0.002514)
      
      optical_mass=min(6.d0,optical_mass)
      
      highcloud_coefficient=0.4d0

      IF (cos_zenith >= 0.342d0) highcloud_coefficient=0.2d0
      
      ! radiation coming through the atmosphere (Equation A2 in Sogachev et al.2002)
      ! (equation 6) 
      cloud_effect = (1.0d0 - 0.8d0*CLOUDCOVER_LOW - 0.5d0*CLOUDCOVER_MIDDLE - highcloud_coefficient*CLOUDCOVER_HIGH)
      cloud_effect = max(2.d-1,cloud_effect)
      
      total_global = radiation_into_earth * (0.3d0+0.7d0*TRANSMITTANCE**optical_mass) * cloud_effect  !real radiation      
   
   END SUBROUTINE TotalSWRAboveCanopy


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine SWRDirectDiffuse
   !> Separates the total radiation at the top of the canopy into direct and diffuse components.
   !
   !>    The total solar (or short wave) radiation reaching the top of the canopy \f$ Q_{h} \f$, is divided into downward direct  
   !>    solar radiation flux \f$ I_{h} \f$ and downward diffuse sky radiation flux \f$ D_{h} \f$. The ratio \f$ D_{h}/Q_{h} \f$ is 
   !>    calculated depending on the ratio \f$ r = Q_{h}/Q_{0} \f$: 
   !>
   !>      \f{equation}{
   !>                    \frac{D_{h}}{Q_{h}} = f = \left\{ \begin{array}{lll}
   !>                                                       1 & \mbox{for} & r < 0.22 \\
   !>                                                       \alpha_D(1 - 6.4(0.22r)^{2}) & \mbox{for} & 0.22 < r < 0.35 \\
   !>                                                       \alpha_D(1.47 - 1.66r) & \mbox{for} & 0.35 <r < \frac{1.47 - R}{1.66} \\
   !>                                                       \alpha_D R & \mbox{for} &\frac{1.47 - R}{1.66} < r
   !>                                                      \end{array}\right.
   !>       \f}
   !>    (Table 2 in Sogachev et al. 2002) where
   !>       \f[
   !>          R = 0.847 - 1.61\cos{\theta} + 0.9\cos^2{\theta} - 0.05\cos^3{\theta}                                    \qquad (2)  
   !>       \f]
   !>    (don't know why like this, it's different in the paper), and
   !>       \f[
   !>          \alpha_D = \left(  \frac{0.22}{r}  \right)^{0.37}	.                                                     \qquad (3)
   !>       \f]   
   !>    \f$ \alpha_D \f$ is a conversion factor to correct for the overestimation of \f$ D_{h} \f$ under clear-sky conditions.
   !>
   !>    So the calculation of the diffuse solar radiation is simply:
   !>       \f[   
   !>          D_{h} = Q_{h} \times f                                                                                   \qquad (4)
   !>       \f]   
   !>    and direct solar radiation:
   !>       \f[
   !>          I_{h} = Q_{h} - D_{h}                                                                                    \qquad (5)
   !>       \f]      
   !>
   !> Reference: Appendix A in Sogachev et al. 2002
   !>   
   !> Includes rows 70-93 from original Scadis code.
   !>
   !>   variable            | model variable       | description
   !>  -------------------- | -------------------- | ---------------------------------------------------------
   !> \f$  Q_{h}        \f$ | total_global         | amount of solar radiation above the canopy (W/m^2)
   !> \f$  h            \f$ | canopy_top           | height of the canopy: lowest model level above the canopy 
   !> \f$  I_{h}        \f$ | direct_swr_down      | downward direct solar radiation flux at \f$ h \f$ (W/m^2)
   !> \f$  D_{h}        \f$ | diffuse_swr_down     | downward diffuse sky radiation flux at \f$ h \f$ (W/m^2)
   !> \f$  Q_{0}        \f$ | radiation_into_earth | incoming radiation outside the atmosphere (W/m^2)
   !> \f$  r            \f$ | ratio_radiation      | \f$ Q_{h}/Q_{0} \f$
   !> \f$  R            \f$ | R_help               | help variable
   !> \f$ \cos{\theta}  \f$ | cos_zenith           | cos(solar zenith angle)
   !> \f$ \alpha_D      \f$ | conversion_factor    | conversion factor 
   !>
   !> @ingroup scadis_radiation
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE SWRDirectDiffuse(radiation_into_earth, total_global, cos_zenith, &
                                  direct_swr_down, diffuse_swr_down)
   
      ! INPUT
      REAL(kind=dp),                INTENT(IN)  :: radiation_into_earth   , & ! (=rads) incoming radiation outside the atm
                                                   total_global           , & ! (=rads2) total radiation at the top of the  canopy
                                                   cos_zenith                  ! (=cos_zenith) cos(zenith angle) = sin(sun height angle)
      ! OUTPUT
      REAL(kind=dp),                INTENT(OUT) :: direct_swr_down        , & ! (=rsnt) downward direct solar radiation flux at the top of the canopy (W/m2)
                                                   diffuse_swr_down           ! (=rskt) downward diffuse solar radiation flux at the top of the canopy (W/m2)
      
      ! LOCAL
      REAL(kind=dp)                             :: ratio_radiation        , & ! (=relrad) ratio of radiation at the top of the canopy and at the top of the atmosphere
                                                   conversion_factor      , & ! (=apar) conversion factor for ratio_radiation >= 0.22
                                                   R_help                 , & ! (=rhelp) help variable, see Rosa's notes
                                                   boundary_help              ! (=radk) help variable
                                                   
      ! help functions           
      
      ! (equation 2)                               
      R_help = 0.847 - 1.61*cos_zenith + 0.9*cos_zenith**2 - 0.05*cos_zenith**3
      !70 rhelp=0.847-1.61*zenith+0.9*zenith**2-0.05*zenith**3 ! help function
   
      boundary_help = (1.47-R_help) / 1.66  
      !71 radk=(1.47-rhelp)/1.66  

      ! none of these calculations necessary, if there is no incoming solar radiation -> moved this if outside of this subroutine
      !IF (radiation_into_earth > 0.) THEN
      
      ! calculate ratio of radiation at the top of the canopy and at the top of the atmosphere
      ratio_radiation = total_global/radiation_into_earth
   
      ! (equation 3) 
      conversion_factor=(0.22/ratio_radiation)**0.37
         
      !END IF
               
      !72 if(rads.gt.0.) then
      !73    relrad=rads2/rads
      !74    apar=(0.22/relrad)**0.37
      !75 endif

      
      
      ! calculating the direct and diffuse radiation fluxes: the calculation depends on 
      ! the ratio of radiation at the top of the canopy and at the top of the atmosphere
      ! (equation 1 and 4) 

      IF (ratio_radiation < 0.22) THEN
         diffuse_swr_down = total_global
      END IF
      !76 if(relrad.lt.0.22) then
      !77   rskt=rads2
      !78   rsnt=0.
      !79 end if
      
      
      IF (ratio_radiation >= 0.22  .AND.  ratio_radiation < 0.35) THEN
         diffuse_swr_down = conversion_factor * total_global * (1.-6.4*(0.22*ratio_radiation)**2)
      END IF
      !80 if(relrad>=0.22.or.relrad<0.35) then
      !81   rskt=min(rads2,rads2*(1.-6.4*(0.22*relrad)**2)*apar)
      !82   rsnt=rads2-rskt
      !83 end if
      
      
      IF (ratio_radiation >= 0.35  .AND.  ratio_radiation < boundary_help) THEN
         diffuse_swr_down = conversion_factor * total_global*(1.47-1.66*ratio_radiation)
      END IF
      
      !84 if(relrad>=0.35.or.relrad<radk) then
      !85   rskt=min(rads2,rads2*(1.47-1.66*relrad)*apar)
      !86   rsnt=rads2-rskt
      !87 end if
      
         
      IF (ratio_radiation >= boundary_help) THEN
         diffuse_swr_down = conversion_factor * total_global * R_help
      END IF
      
      !88 if(relrad>=radk) then
      !89   rskt=min(rads2,rads2*rhelp*apar)
      !90   rsnt=rads2-rskt
      !91 end if
      
      
      diffuse_swr_down = min(total_global,diffuse_swr_down) ! diffuse radiation cannot exceed the total amount of radiation
      direct_swr_down  = total_global - diffuse_swr_down    !(equation 5)
      
      
      ! these just commented out, because rsktm and rsntm not used anywhere
      !92 rsktm=rskt
      !93 rsntm=rsnt
         
   END SUBROUTINE SWRDirectDiffuse
  

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine PenetrationFDirect
   !> Calculates the penetration function for direct solar radiation.
   !> The variable is used for calculating radiation propagation inside the canopy.
   !
   !>    The penetration function for direct radiation, \f$ \eta \f$, is used to calculate the amount of solar radiation at each  
   !>    model level inside the canopy. It's calculated for each model level inside the canopy:
   !>       \f[
   !>          \eta(z) =\exp  \left( - \frac{G_l LAI(z)}{\cos{\theta}} \right)                                          \qquad (1)
   !>       \f]
   !>    (Equation A3 in Sogachev et al. 2002), where \f$ G_l \f$ is the integral function for leaf orientation, and \f$ LAI(z) \f$ 
   !>    is cumulative leaf area index. When \f$ G_l = \cos{\theta} \f$, which is the case for 
   !>       - horizontally oriented phytoelements (leaf_type = 2) 
   !>       - spherically oriented phytoelements (leaf_type = 1) when \f$ Q_{h} <= 10 \quad W/m^2\f$ ,
   !>
   !>    the equation becomes simply
   !>       \f[
   !>          \eta(z) =\exp  \left( - LAI(z)\right)	.                                                                 \qquad (2)
   !>       \f]
   !>
   !>    Don't yet understand why, but in the model cumulative leaf area index is calculated as 
   !>       \f[
   !>          -LAI(z) = lai(z)-lai0                                                                                    \qquad (3)
   !>       \f]
   !>
   !> Reference: Appendix A in Sogachev et al. 2002
   !>   
   !> Includes rows ??? from original Scadis code.
   !>
   !>   variable            | model variable          | description
   !>  -------------------- | ----------------------- | ---------------------------------------------------------
   !> \f$  \eta         \f$ | penetation_direct       | penetration function for direct radiation
   !> \f$  G_l          \f$ | leaf_orientation        | integral function for leaf orientation
   !> \f$  LAI          \f$ | lai_cumulative          | cumulative leaf area index (m^2 / m^2)
   !> \f$ \cos{\theta}  \f$ | cos_zenith              | cos(solar zenith angle)
   !> \f$  h            \f$ | canopy_top              | height of the canopy: lowest model level above the canopy 
   !> \f$  Q_{h}        \f$ | total_global            | amount of solar radiation above the canopy (W/m^2)
   !> \f$  ??           \f$ | n_levels                | number of levels in model domain [1]
   !>
   !
   !> @todo add reference line numbers from original code
   !> @ingroup scadis_radiation
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE PenetrationFDirect(canopy_top, n_levels, leaf_type, total_global, cos_zenith, leaf_orientation, lai_cumulative, &
                                 penetration_direct)
                                 
      ! INPUT
      REAL(kind=dp),                      INTENT(IN   ) :: total_global          , & ! (=rads2) total radiation at the top of the  canopy
                                                           cos_zenith                ! (=cos_zenith) cos(zenith angle) = sin(sun height angle)
      REAL(kind=dp), DIMENSION(n_levels), INTENT(IN   ) :: leaf_orientation      , & ! (=gl) integral function for leaf orientation
                                                           lai_cumulative            ! (= -(lai(k)-lai0) ), cumulative leaf area index (m^2 / m^2)
      INTEGER,                            INTENT(IN   ) :: canopy_top            , & ! (=nz, znrad) top of the canopy (model level)
                                                           n_levels              , & ! (=kz) number of levels  [1]
                                                           leaf_type                 ! (=su) flag for leaf type. 1 for conifirous trees, 2 for broadleaf trees
      ! OUPUT
      REAL(kind=dp), DIMENSION(n_levels), INTENT(  OUT) :: penetration_direct        ! (=psn) penetration function for diffuse radiation

      ! LOCAL
      INTEGER                                   :: k                         ! for loops for height
      
      
      DO k = canopy_top,1,-1
      
         IF (leaf_type == 2   .OR.  total_global <=  10)   THEN
          
            penetration_direct(k) = exp( -lai_cumulative(k) ) 
         
         ELSE  
         
            penetration_direct(k) = exp( -leaf_orientation(k) * lai_cumulative(k) / cos_zenith) 
            
         END IF
          
      ENDDO 
      
      ! DO k=nz,1,-1
      !   psn(k)=exp(0.5*(lai(k)-lai0)/zenith)  ! penetration function for direct incident radiation
      !   if(rads2.le.10.) psn(k)=exp((lai(k)-lai0))
      !   if(su.eq.2) then  !  horizontally oriented phytoelements
      !      psn(k)=exp((lai(k)-lai0))
      !   endif
      ! ENDDO

   
   END SUBROUTINE PenetrationFDirect


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine PenetrationFDiffuse
   !> Calculates the penetration function for diffuse solar radiation.
   !> The variable is used for calculating radiation propagation inside the canopy.
   !
   !>    The penetration function for diffuse radiation, \f$ \mu \f$, is used to calculate the amount of solar radiation at each 
   !>    model level inside the canopy. It's calculated for each model level inside the canopy:
   !>       \f[
   !>          \mu(z) =2  \int_{0}^{\pi/2}  \exp \left( - \frac{G_l LAI(z)}{\cos{\theta}} \right) \cos{\theta} \sin{\theta}\,d\theta 
   !>                                                                                                                   \qquad (1)
   !>       \f]
   !>    (Equation A4 in Sogachev et al. 2002), where \f$ G_l \f$ is the integral function for leaf orientation, and \f$ LAI(z) \f$ 
   !>    is cumulative leaf area index. When \f$ G_l = \cos{\theta} \f$, which is the case for 
   !>       - horizontally oriented phytoelements (leaf_type = 2) 
   !>       - spherically oriented phytoelements (leaf_type = 1) when \f$ Q_{h} <= 10 \quad W/m^2\f$ ,
   !>
   !>    the solution becomes simple:
   !>       \f[
   !>          \mu(z) =2   \exp \left( - LAI(z) \right) \int_{0}^{\pi/2}  \cos{\theta} \sin{\theta}\,d\theta 
   !>                                                                                                                   \qquad (2)
   !>       \f]
   !>    The solution of the integral is 1/2, so we get for the penetration function
   !>       \f[
   !>          \mu(z) =   \exp \left( - LAI(z) \right)  
   !>                                                                                                                   \qquad (3)
   !>       \f]   
   !>    which is the same result as for the penetration function for direct radiation (Equation 2 in PenetrationFDirect). In this 
   !>    case (\f$ G_l = \cos{\theta} \f$) we can use \f$ \mu  = \eta \f$. In other cases, the integral needs to be solved 
   !>    numerically. For this the trapezoidal method is used.
   !>
   !>    The trapezoidal rule is Newton-Cotes formula to solve integrals of one dimensional functions. The integral of function 
   !>    \f$ f(x) \f$ from \f$ a \f$ to \f$ b \f$ can be approximated as 
   !>       \f[  
   !>          \int_a^b f(x)\,dx \approx (b - a) \frac{f(a) + f(b)}{2}                                                  \qquad (4)
   !>       \f]   
   !>    If the interval \f$ [a,b] \f$ is wide, the method is of course very inaccurate. This can be improved by dividing the 
   !>    interval into \f$ N \f$ subintervals, and using the Equation (2) \f$ N - 1 \f$ times and do the integrating in intervals 
   !>    \f$ [x_1,x_2] \f$, \f$ [x_2,x_3] \f$, ... \f$ [x_{N-1},x_N] \f$. For an uniform spacing, obtained by
   !>       \f[
   !>          x_i = a + \Delta x (i-1)                                                                                 \qquad (5)
   !>       \f]
   !>    where
   !>       \f[
   !>          \Delta x = \frac { (a-b) }{ N }  \quad ,                                                                 \qquad (6)
   !>       \f]
   !>    the integral of the entire interval \f$ [a,b] \f$ can then be calculated by (the composite rule)
   !>       \f[
   !>          \int_a^b f(x)\,dx \approx \frac{(b - a)}{N} \left[  \frac{f(x_1)}{2} + f(x_2) + f(x_3) + f(x_{N - 1}) + ... 
   !>                         \frac{f(x_N)}{2} \right]	                                                               \qquad (7)
   !>       \f]
   !>    or using the summation notation
   !>       \f[
   !>          \int_a^b f(x)\,dx \approx \frac{(b - a)}{N} \left[  \frac{f(x_1)}{2} + \sum_{i=2}^{N-1} f(x_i) + 
   !>                         \frac{f(x_N)}{2} \right]	.                                                              \qquad (8)
   !>       \f]
   !>
   !>    In this case \f$ a = 0 \f$ and \f$ b = \pi/2 \f$ (Equation 1), and \f$ f(0) = f(\pi/2) = 0 \f$, so the equation simplifies 
   !>    to
   !>       \f[
   !>          \int_a^b f(x)\,dx \approx \frac{(b - a)}{N} \left[ \sum_{i=2}^{N-1} f(x_i) \right]	.                    \qquad (9)
   !>       \f]
   !>    Finally, for the solution of Equation 1, we get
   !>       \f[
   !>          \mu(z) \approx 2 \frac{\pi/2}{N} \left[ \sum_{i=2}^{N-1} f(\theta_i) \right]      
   !>                  \approx \frac{\pi}{N} \left[ \sum_{i=2}^{N-1} f(\theta_i) \right]                                \qquad (10)
   !>       \f]   
   !>    where
   !>       \f[   
   !>          f(\theta_i) = \exp \left( - \frac{G_l LAI(z)}{\cos{\theta}} \right) \cos{\theta} \sin{\theta}            \qquad (11)
   !>       \f]
   !>
   !>
   !> References: 
   !> - Appendix A in Sogachev et al. 2002
   !> - Numerical Recipes in Fortran 90, p.123, 127 (http://apps.nrbook.com/fortran/index.html)
   !> - Wikipedia: Trapezoidal rule (http://en.wikipedia.org/wiki/Trapezoidal_rule)
   !>   
   !> Includes rows ??? from original Scadis code.
   !>
   !>   variable             | model variable       | description
   !>  --------------------- | -------------------- | ---------------------------------------------------------
   !> \f$  \mu           \f$ | penetation_diffuse   | penetration function for diffuse radiation
   !> \f$  G_l           \f$ | leaf_orientation     | integral function for leaf orientation
   !> \f$  LAI           \f$ | lai_cumulative       | cumulative leaf area index (m^2 / m^2)
   !> \f$ \cos{\theta}   \f$ | cos_zenith           | cos(solar zenith angle)
   !> \f$  h             \f$ | canopy_top           | height of the canopy: lowest model level above the canopy 
   !> \f$  Q_{h}         \f$ | total_global         | amount of solar radiation above the canopy (W/m^2)
   !> \f$  ??            \f$ | n_levels             | number of levels in model domain [1]
   !> \f$  f             \f$ | function4integration | the function, whose integral needs to be solved
   !> \f$  a             \f$ | -                    | lower limit for integration
   !> \f$  b             \f$ | -                    | upper limit for integration   
   !> \f$  N             \f$ | n_subint             | number of subintervals used for trapezoidal rule
   !> \f$  \theta_i      \f$ | thetas               | array definining the subintervals for the composite trapezoidal rule
   !> \f$  \Delta \theta \f$ | delta_theta          | distance of the points of the subintervals
   !>
   !>
   !
   !> @todo add reference line numbers from original code
   !> @todo read trough description
   !>
   !> @ingroup scadis_radiation
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE PenetrationFDiffuse(n_levels, canopy_top, leaf_orientation, lai_cumulative, penetration_diffuse)


      ! INPUT                                  
      REAL(kind=dp), DIMENSION(n_levels), INTENT(IN   ) :: leaf_orientation      , & ! (=gl) integral function for leaf orientation
                                                           lai_cumulative            ! (= -(lai(k)-lai0) ), cumulative leaf area index (m^2 / m^2)
      INTEGER,                            INTENT(IN   ) :: canopy_top            , & ! (=nz, znrad) top of the canopy (model level)
                                                           n_levels             ! , & ! (=kz) number of levels  [1]
                                                           !leaf_type                 ! (=su) flag for leaf type. 1 for conifirous trees, 2 for broadleaf trees
      ! OUPUT
      REAL(kind=dp), DIMENSION(n_levels), INTENT(  OUT) :: penetration_diffuse        ! (=psk) penetration function for diffuse radiation

      ! LOCAL
      INTEGER                                           :: k, i                  , & ! for loops
                                                           n_subint              , & ! (=ntet) number of subintervals used
                                                           allocation_status         ! used for checking whether variable allocation successfull
      REAL(kind=dp)                                     :: delta_theta               ! (=dtet) width of the subintervals
      REAL(kind=dp), ALLOCATABLE                        :: thetas(:)             , & ! (=tet)
                                                           function4integration(:)   ! (=fl) the function that is integrated
   

      ! The structure and logic of the code is not 100% identical to original scadis, but the principle is the same. The changes 
      ! influenced the result by a difference less than 0.0001%, so knowone will even notice.

      ! defining the number of subintervals used
      n_subint = 30
      
      ! allocating the varianles, which dimensions depend on n_subint
      ALLOCATE(thetas(2:n_subint+1), STAT = allocation_status)
         IF (allocation_status /= 0) THEN ! in case something goes wrong
            WRITE(*,*) 'problems in subroutine PenetrationFDiffuse, Scadis_Radiation.f90' 
            CALL abort()  
         END IF

      ALLOCATE(function4integration(2:n_subint+1), STAT = allocation_status)
         IF (allocation_status /= 0) THEN ! in case something goes wrong
            WRITE(*,*) 'problems in subroutine PenetrationFDiffuse, Scadis_Radiation.f90' 
            CALL abort()  
         END IF
         
         
      ! calculating the spacing of the subintervals (Equation 4)
      delta_theta = ( pi/2. ) /n_subint  ! ( upper limit  - lower limit ) / number of subsections


      ! calcuating thetas, the values of the boundaries of the subintervals (Equation 3)
      ! lower limit = 0,  would have thetas(1) = 0 but ignored since f(0) = 0
      
      ! for some reason the subintervals shifted by -delta_theta/2, so that 
      ! thetas(2) = delta_theta/2
      ! thetas(3) = delta_theta/2 + delta_theta
      ! thetas(4) = delta_theta/2 + delta_theta*2
      ! etc...

      DO i=2,n_subint+1
         thetas(i)=delta_theta/2+delta_theta*(i-2)
      ENDDO

 
      ! penetration function at the top of the canopy 
      penetration_diffuse(canopy_top)=1.         

      ! calculation of the penetration function for each model level inside the canopy
      DO k = canopy_top-1,1,-1
      
         ! Equation 9
         function4integration=exp(leaf_orientation(k)*(- lai_cumulative(k) )/cos(thetas)) *cos(thetas)*sin(thetas)
         
         ! Equation 8
         penetration_diffuse(k) = (pi /(n_subint)) * ( sum (function4integration)) 
 
      END DO
     
     
      DEALLOCATE(thetas)
      DEALLOCATE(function4integration)


     !DO ja=1,ntet
     !   tet(ja)=teto+dtet*(ja-1)
     !ENDDO
     !DO k =nz-1,1,-1
     !   DO ja=1,ntet
     !      cote=cos(tet(ja))
     !      site=sin(tet(ja))
     !      fl(ja)=2.*exp(gl(k)*(lai(k)-lai0)/cote)*cote*site
     !   ENDDO
     !   psk(k)=(fl(1)+fl(ntet))*dtet/2.
     !   DO ja=2,ntet
     !      psk(k)=psk(k )+(fl(ja-1)+fl(ja))*dtet/2.
     !      if(su.eq.2) then
     !         psk(k)=psn(k)
     !      endif
     !   ENDDO
     !sENDDO


  
   END SUBROUTINE PenetrationFDiffuse
      

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine IntegralFLeafOrientation
   !> Calcualtes the integral function for leaf orientation for 1) spherically and 2) horizontally oriented phytoelements.
   !> The variable is used for calculating radiation propagation inside the canopy.
   !
   !>    Calculates the integral function for leaf orientation, \f$ G_l \f$, for 1) spherically 
   !>       \f{equation}{
   !>          G_l  = \left\{ \begin{array}{rcl}
   !>                                           0.5 & \mbox{for} & Q_{h} > 10 \mathrm{W/m}^2 \\
   !>                                           \cos{\theta} & \mbox{for} & Q_{h} <= 10 \mathrm{W/m}^2  \\
   !>                         \end{array}\right.
   !>       \f}
   !>    and 2) horizontally oriented phytoelements
   !>       \f[
   !>          G_l = \cos{\theta}		.                                                                                \qquad (2)
   !>       \f]
   !>    The variable is used for calculating radiation propagation inside the canopy. If the total incident radiation at the 
   !>    top of the atmosphere \f$ Q_{h} <= 5 \mathrm{W/m}^2 \f$, then  \f$ G_l = 0 \f$. Above the canopy \f$ G_l = 1 \f$.
   !>
   !>
   !> Reference: Appendix A in Sogachev et al. 2002
   !>   
   !> Includes rows ??? from original Scadis code.
   !>
   !>   variable             | model variable       | description
   !>  --------------------- | -------------------- | ---------------------------------------------------------
   !> \f$  G_l           \f$ | leaf_orientation     | integral function for leaf orientation
   !> \f$ \cos{\theta}   \f$ | cos_zenith           | cos(solar zenith angle)
   !> \f$  h             \f$ | canopy_top           | height of the canopy: lowest model level above the canopy 
   !> \f$  Q_{h}         \f$ | total_global         | amount of solar radiation above the canopy (W/m^2)
   !> \f$  ??            \f$ | n_levels             | number of levels in model domain [1]
   !>
   !>
   !> @todo check description
   !> @todo add reference line numbers from original code
   !>
   !> @ingroup scadis_radiation
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE IntegralFLeafOrientation(n_levels, canopy_top, radiation_total, cos_zenith, leaf_type, &
                                                 leaf_orientation)
  
      ! INPUT
      REAL(kind=dp),              INTENT(IN)    :: radiation_total        , & ! (=rads2) total radiation at the top of the  canopy
                                                   cos_zenith                  ! (=cos_zenith) cos(zenith angle) = sin(sun height angle)
                                                   
      INTEGER,                    INTENT(IN)    :: n_levels               , & ! (=kz) number of levels used [1]
                                                   canopy_top       , & ! (=nz, znrad) height of the canopy: lowest model level above the canopy
                                                   leaf_type                  ! (=su) flag for leaf type. 1 for conifirous trees, 2 for broadleaf trees
                                                         
      ! OUPUT
      REAL(kind=dp), DIMENSION(n_levels), INTENT(OUT)  :: leaf_orientation  ! (=gl) integral function for leaf orientation


      ! LOCAL
      INTEGER                                   :: k                         ! for loops
            
      
      DO k=2,n_levels
         leaf_orientation(k) = 1.
      ENDDO
      ! DO k=2,kz
      !    gl(k)=1.
      ! ENDDO
     
      IF (radiation_total <= 5) THEN
            leaf_orientation = 0.
      ELSE
   
         DO k = canopy_top,1,-1
         
               IF (leaf_type == 1) THEN ! spherically oriented phytoelements
                  leaf_orientation(k) = 0.5
                  IF(radiation_total <= 10.) leaf_orientation(k) = cos_zenith  
               
               ELSEIF (leaf_type == 2) THEN !  horizontally oriented phytoelements
                  leaf_orientation(k) = cos_zenith
               ENDIF
            
         ENDDO
         
      ENDIF
      
      ! DO k=nz,1,-1
      !   gl(k)=0.5    ! integral function for leaf orientation, spherically oriented phytoelements
      !   if(rads2.le.10.) gl(k)=zenith  !  horizontally oriented phytoelements
      !   if(rads2.le.5.) gl(k)=0.
      !   if(su.eq.2) then
      !      gl(k)=zenith
      !      if(rads2.le.5.) gl(k)=0.
      !   endif
      ! ENDDO

   END SUBROUTINE IntegralFLeafOrientation
   
END MODULE Scadis_Radiation
