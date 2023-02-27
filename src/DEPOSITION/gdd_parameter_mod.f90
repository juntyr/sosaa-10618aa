MODULE gdd_parameter_mod
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)
  REAL(dp), PARAMETER :: Eps = 1.0d-20  ! to avoid divided by 0
  REAL(dp), PARAMETER :: NA = 6.022045d23  ! [# mol-1]
  REAL(dp), PARAMETER :: Mair = 28.970d-3  ! [kg mol-1]


END MODULE gdd_parameter_mod
