! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
! Global parameters for SCADIS (meteorology) in sosa.                                           !
! For more information see MT_MainMet.f90                                                       !
!                                                                                               !
! Commented and restructured by Rosa Gierens, University of Helsinki                            !
!                                                                                               !
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

MODULE Scadis_Parameters

  ! variables from other modules
  USE Sosa_data, ONLY: dp    ! For simplicity, the precision of real numbers is applied to the whole program. L

  IMPLICIT NONE
  
  ! Model parameteres
  REAL(kind=dp), PARAMETER :: TRANSMITTANCE = 0.8                     ! atmospheric transmittance for global radiation  
                              ! MT_MainScadis 23 & 29


  REAL(kind=dp), PARAMETER :: CLOUDCOVER_LOW = 0.                 , & ! Cover fraction of low,       
                              CLOUDCOVER_MIDDLE = 0.              , & ! middle,
                              CLOUDCOVER_HIGH = 0.                    ! and high clouds (0 - 10).    
                              ! from Scadis_Initial.f90         
    
END MODULE Scadis_Parameters
