MODULE SATCONC_MOD
IMPLICIT NONE

contains


    pure elemental function saturation_conc_m3(A,B, Temperature) result(Vapour_concentration)
      implicit none

      INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(14,300)
      real(dp), PARAMETER   :: kb = 8.3144598/6.022140857d23    ! [J/K] Boltzmann constant
      real(dp), intent(in)  :: A, B, temperature
      real(dp)              :: Vapour_concentration, vapour_pressure

      ! Using antoine equation log_10(p) = A- (B/T)
      vapour_pressure      = 10 ** (A - (B/temperature)) ! in atm
      Vapour_concentration = (vapour_pressure*101325)/(kb * temperature) ! #/m3

    end function saturation_conc_m3


END MODULE SATCONC_MOD
