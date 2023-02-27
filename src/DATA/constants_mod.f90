!-------------------------------------------------------------------------------
!
!  Universal constants and global parameters
!
!-------------------------------------------------------------------------------

module constants_mod

implicit none

public

!-----------------------------------------------------------------------------
! Definition of different levels of accuracy for REAL variables using KIND
! parameterization
!-----------------------------------------------------------------------------

! SP - Single precision kind
integer, parameter :: sp = selected_real_kind(6,30)
! DP - Double precision kind
integer, parameter :: dp = selected_real_kind(14,300)
! QP - Quadruple precision kind
integer, parameter :: qp = selected_real_kind(18,400)

! this data format is used in simbim so we have to keep it
integer, parameter :: dp_sb = selected_real_kind(15,307)


!-----------------------------------------------------------------------------
! Time constants
!-----------------------------------------------------------------------------

integer, parameter :: NMONTH = 12

! Number of days in each month
integer, dimension(NMONTH) :: &
  MONTH_DAYS = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/), &
  MONTH_DAYS_LEAP = (/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)


!-----------------------------------------------------------------------------
! Math and physics constants
!-----------------------------------------------------------------------------

real(dp), parameter :: PI   = 4.0d0*atan(1.0d0)  ! PI
real(dp), parameter :: PIx2 = PI*2               ! PI*2
real(dp), parameter :: deg2rad = PI/180.0d0      ! degree to radian
real(dp), parameter :: rad2deg = 180.0d0/PI      ! degree to radian

REAL(dp), PARAMETER :: T00  = 273.15d0   ! [K], 0 deg to K
REAL(dp), PARAMETER :: Patm = 1.01325d5  ! [Pa], standard air pressure
REAL(dp), PARAMETER :: grav = 9.80665d0  ! [m s-2], standard gravity

! Stefan-Boltzmann constant, [W m-2 K-4]
REAL(dp), PARAMETER :: sigma_sb = 5.670374d-8

! Solar constant [W m-2], varying in fact, but in our simulation cases, it can
! be considered as a constant. Here refer to Stull (2017, version 1.02b)
! Eq. 2.17.
REAL(dp), PARAMETER :: solar_constant = 1361.0d0

! Avogadro's number, [molec mol-1]
REAL(dp), PARAMETER :: Avog    = 6.02214076d23

! Universal gas constant, [J mol-1 K-1]
REAL(dp), PARAMETER :: Rgas    = 8.314462618d0

! Specific gas constant, [J kg-1 K-1]
real(dp), parameter :: xmm_d   = 28.9647d-3   ! [kg mol-1], molar mass of dry air
real(dp), parameter :: xmm_H2O = 18.01528d-3  ! [kg mol-1], molar mass of H2O
REAL(dp), PARAMETER :: Rgas_d    = Rgas/xmm_d    ! dry air
REAL(dp), PARAMETER :: Rgas_H2O  = Rgas/xmm_H2O  ! water vapor


REAL(dp), PARAMETER :: HALF  = 0.5_dp
REAL(dp), PARAMETER :: ONE   = 1.0_dp
REAL(dp), PARAMETER :: THREE = 3.0_dp
REAL(dp), PARAMETER :: FOUR  = 4.0_dp
REAL(dp), PARAMETER :: FIVE  = 5.0_dp
REAL(dp), PARAMETER :: SECONDS_IN_ONE_MINUTE = 60.0d0
REAL(dp), PARAMETER :: SECONDS_IN_HALF_HOUR  = 1800.0d0
REAL(dp), PARAMETER :: SECONDS_IN_ONE_HOUR   = 3600.0d0
REAL(dp), PARAMETER :: SECONDS_IN_HALF_DAY   = 43200.0d0
REAL(dp), PARAMETER :: SECONDS_IN_ONE_DAY    = 86400.0d0


!-----------------------------------------------------------------------------
! Others
!-----------------------------------------------------------------------------

INTEGER, PARAMETER :: LINE_WIDTH   = 300  ! maximum line width for input files
INTEGER, PARAMETER :: SPC_NAME_LEN = 120  ! maximum length of a species name

! Used in the soil emisison calculations, not sure what it means
REAL(dp), PARAMETER :: conversion = 1d-9 * 6.02214129d23 * 1.0d0/60.0d0 * 1d-4

REAL(dp), PARAMETER :: c432=43200.0d0  ! [s], seconds of half day
REAL(dp), PARAMETER :: c864=86400.0d0  ! [s], seconds of one day

REAL(dp):: hh=3000.0d0

!***** Coefficient for calculate of the heat resistence *****!
REAL(dp), PARAMETER :: h(28) = (/ &
  1000., 0.528, 0.5, 200000., 0.267, 0.6, 0.018,            &
  0.84, 0.000143, 0.00143, 714., 1.128, 0.125, 28500000.,     &
  0.433, 0.25, 0.12, 0.33, 0.506, 0.254, 0.017, 0.00016,      &
  0.0016, 800., 1.112, 32000000., 0.48, 1.115 &
  /)


!-----------------------------------------------------------------------------
! Not used or need not sure what used for, need to remove in future
!-----------------------------------------------------------------------------

! Number of gas species
INTEGER, PARAMETER :: NG = 37

!===== Aerosol related variables =====!
! factor controlling reaction products from NO3 oxidation. 2 means magnifying concentration by 2 times.
REAL(dp), PARAMETER :: no3_coe = 1.0d0
! factor controlling reaction products from OH, O3, NO3 oxidation.  
REAL(dp), DIMENSION(3), PARAMETER :: vap_coe = (/1.0d0, 0.2d0, 1.0d0/)
REAL(dp), PARAMETER :: cut_off_diameter=1.*1d-7
INTEGER::cut_off_bin

! Temporary
  real(dp), parameter :: maga=17.572d0, magb=241.9d0
  real(dp), parameter :: qv0=1.3318375d0

  ! Constants for saturation vapour pressure 
  REAL(dp), PARAMETER :: &
    a0 = 6.107799961d0,  &
    a1 = 4.436518524d-1, &
    a2 = 1.428945805d-2, &
    a3 = 2.650648471d-4, &
    a4 = 3.031240396d-6, &
    a5 = 2.034080948d-8, &
    a6 = 6.136820929d-11


CONTAINS


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


! Day of year
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


INTEGER FUNCTION MonthDay(y, m)
  INTEGER, INTENT(IN) :: y, m
  INTEGER, PARAMETER :: md(12) = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
  
  MonthDay = md(m)
  IF (m==2 .AND. IsLeapYear(y)) THEN
    MonthDay = MonthDay + 1
  END IF
END FUNCTION MonthDay

!------------------------------------------------------------------------------!
!
! Saturation vapor pressure
! 
! Different equations can be used to calculate the saturation vapor pressure P
! [Pa] with T [degC]:
!   Magnus equation: P = 610.94 * exp( (17.625*T) / (T + 243.04) )
!   Buck equation: P = 611.21 * exp( (18.678 - T/234.5) * ( T/(257.14+T) ) )
!
! Here Buck equation is more accurate and thus recommended. In this model, the
! absolute humidity [kg m-3] is used, so the SVP should be converted to rhovs
! [kg m-3] by:
!   rhovs = P / (Rgas_H2O*TK)
!
! Reference: https://en.wikipedia.org/wiki/Vapour_pressure_of_water
!
!------------------------------------------------------------------------------!

pure function svp_magnus(T)
  real(dp) :: svp_magnus  ! [Pa], saturation vapor pressure with Magnus method
  real(dp), intent(in) :: T  ! [K], temperature

  ! Coefficients
  real(dp), parameter :: p0 = 610.94d0
  real(dp), parameter :: a=17.625d0, b=243.04d0

  svp_magnus = p0*exp( a*(T-T00) / (b+T-T00) )
end function svp_magnus


pure function rhovs_magnus(T)
  real(dp) :: rhovs_magnus  ! [kg m-3], saturation absolute humidity
  real(dp), intent(in) :: T  ! [K], temperature

  rhovs_magnus = svp_magnus(T) / (Rgas_H2O*T)
end function rhovs_magnus


pure function svp_buck(T)
  real(dp) :: svp_buck  ! [Pa], saturation vapor pressure with Magnus method
  real(dp), intent(in) :: T  ! [K], temperature

  ! Coefficients
  real(dp), parameter :: p0 = 611.21d0
  real(dp), parameter :: a=18.678d0, b=234.5d0, c=257.14d0

  real(dp) :: TC

  TC = T - T00
  svp_buck = p0 * exp( (a - TC/b) * ( TC/(c+TC) ) )
end function svp_buck


pure function rhovs_buck(T)
  real(dp) :: rhovs_buck  ! [Pa], saturation vapor pressure with Magnus method
  real(dp), intent(in) :: T  ! [K], temperature

  rhovs_buck = svp_buck(T) / (Rgas_H2O*T)
end function rhovs_buck


pure function svp(T)
  real(dp) :: svp  ! [Pa], saturation vapor pressure
  real(dp), intent(in) :: T

  real(dp) :: TC

  ! Constants for saturation vapour pressure 
  REAL(dp), PARAMETER :: &
    a0 = 6.107799961d0,  &
    a1 = 4.436518524d-1, &
    a2 = 1.428945805d-2, &
    a3 = 2.650648471d-4, &
    a4 = 3.031240396d-6, &
    a5 = 2.034080948d-8, &
    a6 = 6.136820929d-11

  TC = T - T00
  svp = (a0 + a1*TC + a2*TC**2 + &
    a3*TC**3 + a4*TC**4 + a5*TC**5 + a6*TC**6)*100.0d0
end function svp

end module constants_mod
