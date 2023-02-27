! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
! All parts of SCADIS that deal with the components of the enery balance (except radiation,     !
! those are in the Scadis_Radiation):  turbulent heat fluxes, energy balance closure and the    !
! canopy interactions (leaf temperature etc)                                                    !
!                                                                                               !
! For more information see MT_MainMet.f90   co                                                  !
!                                                                                               !
! Commented and restructured by Rosa Gierens, University of Helsinki                            !
!                                                                                               !
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !


MODULE Scadis_Energy_Canopy

  ! variables from other modules
  USE Sosa_data, ONLY : dp    ! For simplicity, the precision of real numbers is applied to the whole program. L.
  !USE SOSA_DATA !, ONLY : variablename
  
  ! subroutines and variables part of the meteorology scheme
  ! USE Scadis_Parameters
  ! USE Scadis_Subroutines

  IMPLICIT NONE

  ! will make this privata later
  !PRIVATE

  !Public subroutines:
  !PUBLIC :: 

  !Public variables:
  !PUBLIC ::   
  
  
CONTAINS


! executable code here
    
END MODULE Scadis_Energy_Canopy
