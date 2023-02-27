! =======================================================================
!     MODULE GAMMA_ETC
! 
!     CONTAINS: 1)GAMMA_LAI
!               2)GAMMA_P
!               3)GAMMA_TISOP
!               4)GAMMA_TNISP
!               5)GAMMA_A -2012
!               6)GAMMA_S -2012
! 
!     History:
!     Made some changes in gamma_TNSIP to make it take temprature in layers
! =======================================================================

      MODULE M2_GAMMA_ETC

      USE constants_mod, only : dp

      IMPLICIT NONE

      CONTAINS

!-----------------------------------------------------------------------
! .....1) Calculate GAM_L (GAMMA_LAI)
!-----------------------------------------------------------------------
!                            0.49[LAI]
!             GAMMA_LAI = ----------------    (non-dimension)
!                         (1+0.2LAI^2)^0.5
! 
!     SUBROUTINE GAMMA_LAI returns the GAMMA_LAI values
!-----------------------------------------------------------------------
      SUBROUTINE GAMMA_LAI( NCOLS, NROWS, &
                           LAI, GAM_L )

      IMPLICIT NONE

      INTEGER NCOLS, NROWS
      REAL(kind=dp) LAI(NCOLS,NROWS), GAM_L(NCOLS,NROWS)

      GAM_L = (0.49*LAI) / ( (1+0.2*(LAI**2))**0.5 )

      RETURN
      END SUBROUTINE GAMMA_LAI
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! .....3) Calculate GAM_T (GAMMA_T) for isoprene
!-----------------------------------------------------------------------
!                          Eopt*CT2*exp(CT1*x)
!             GAMMA_T =  ------------------------
!                        [CT2-CT1*(1-exp(CT2*x))]
!           where x      = [ (1/Topt)-(1/Thr) ] / 0.00831
!                 Eopt   = 1.75*exp(0.08(Tdaily-297)
!                 CT1    = 80
!                 CT2    = 200
!                 Thr    = hourly average air temperature (K)
!                 Tdaily = daily average air temperature (K)
!                 Topt   = 313 + 0.6(Tdaily-297)
! 
!                 Note: AAA = Eopt*CT2*exp(CT1*x)
!                       BBB = [CT2-CT1*(1-exp(CT2*x))]
!                       GAMMA_T = AAA/BBB
! 
!     SUBROUTINE GAMMA_TISOP returns the GAMMA_T value for isoprene
!-----------------------------------------------------------------------
      SUBROUTINE GAMMA_TISOP( NCOLS, NROWS, TEMP, D_TEMP, GAM_T )

      IMPLICIT NONE

      INTEGER NCOLS, NROWS

      REAL(kind=dp)    TEMP(NCOLS,NROWS)              ! hourly surface temperature
      REAL(kind=dp)    D_TEMP(NCOLS,NROWS)            ! daily surface temperature
      REAL(kind=dp)    GAM_T(NCOLS,NROWS)             ! GAMMA_T

! ...  Local parameters
      REAL(kind=dp)    Eopt(NCOLS,NROWS), Topt(NCOLS,NROWS), X(NCOLS,NROWS)
      REAL(kind=dp)    AAA(NCOLS,NROWS), BBB(NCOLS,NROWS)
      REAL(kind=dp), PARAMETER :: CT1 = 80.0
      REAL(kind=dp), PARAMETER :: CT2 = 200.0

      Eopt = 1.75 * exp(0.08*(D_TEMP-297.0))
      Topt = 313.0 + ( 0.6*(D_TEMP-297.0) )
      X = ( (1/Topt)-(1/TEMP) ) / 0.00831

      AAA = Eopt*CT2*exp(CT1*X)
      BBB = (  CT2-CT1*( 1-exp(CT2*X) )  )

      GAM_T = AAA/BBB

      RETURN
      END SUBROUTINE GAMMA_TISOP
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! .....4) Calculate GAM_T (GAMMA_T) for non-isoprene
!-----------------------------------------------------------------------
! 
!             GAMMA_T =  exp[BETA*(T-Ts)]
!           where BETA   = temperature dependent parameter
!                 Ts     = standard temperature (normally 303K, 30C )
! 
!     SUBROUITINE GAMMA_TNISP returns the GAMMA_T value for non-isoprene
!-----------------------------------------------------------------------
      SUBROUTINE GAMMA_TNISP( NCOLS, NROWS, SPCNAM, TEMP, LAYERS, &
                            GAM_T                               )

      USE M2_TEMPD_PRM

      IMPLICIT NONE

      
      !INCLUDE 'TEMPD_PRM.EXT'

      INTEGER       INDEX1
      EXTERNAL      INDEX1
      CHARACTER*16  SPCNAM
      INTEGER       NCOLS, NROWS, VAR, LAYERS
      INTEGER       SPCNUM                             ! Species number
      REAL(kind=dp)           TEMP(NCOLS,NROWS,LAYERS), GAM_T(NCOLS,NROWS,LAYERS)
      REAL(kind=dp), PARAMETER :: Ts = 303.0

      DO VAR = 1 , N_TDF_SPC
	     IF ( (TRIM(SPCNAM)) .EQ. (TRIM(TDF_SPC( VAR )))) THEN
              SPCNUM = VAR
	     end IF
      ENDDO

      !SPCNUM = INDEX1(SPCNAM,N_TDF_SPC,TDF_SPC)
    
      GAM_T = exp( TDF_PRM(SPCNUM)*(TEMP-Ts) )

      RETURN
      END SUBROUTINE GAMMA_TNISP
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! Calculate GAM_A (GAMMA_age) the activity factor for the age of the
! leaves/needles
!-----------------------------------------------------------------------
      SUBROUTINE GAMMA_A(NCOLS, NROWS, SPC_NAME, LAIARp, &
                         LAIARc, TSTLEN, D_TEMP, GAM_A)

      USE M2_REL_EM_ACT !Relative emission activity for leaves of different age
      
      IMPLICIT NONE

      CHARACTER*16 ::  FUNCNAME = 'GAMMA_A'
      INTEGER ::       NCOLS, NROWS
      CHARACTER*8 ::   SPC_NAME
      REAL(kind=dp) :: D_TEMP(NCOLS,NROWS) !this is currently 0K
      REAL(kind=dp) :: LAIARp(NCOLS,NROWS), LAIARc(NCOLS,NROWS)
      INTEGER ::       TSTLEN
      REAL(kind=dp) :: GAM_A(NCOLS,NROWS)
      REAL(kind=dp) :: Fnew, Fgro, Fmat, Fold !Fraction of new, growing, mature, and old foliage
      INTEGER ::       AINDX       ! relative emission acitivity index
      CHARACTER*256 :: MESG        ! message buffer
      REAL(kind=dp) :: LAIp        ! LAI at previous month
      REAL(kind=dp) :: LAIc        ! LAI at current month
      INTEGER ::       t           ! time step [days]
      REAL(kind=dp) :: ti          ! number of days between budbreak and the induction of emission
      REAL(kind=dp) :: tm          ! number of days between budbreak and the initiation of peak emissions rates
      REAL(kind=dp) :: Tt          ! Daily average temperature [K] - this is currently 0K
      INTEGER ::       I, J                

      SELECT CASE ( TRIM(SPC_NAME) )
      CASE ('ACTO','ACTA','FORM','CH4','NO','CO')
         AINDX = 1
      CASE ('MYRC','SABI','LIMO','3CAR','OCIM','BPIN','APIN','OMTP', &
           'CINE','LINA')
         AINDX = 2
      CASE ('FARN','BCAR','OSQT')
         AINDX = 3
      CASE ('MEOH')
         AINDX = 4
      CASE ('ISOP','MBO' )
         AINDX = 5
      CASE DEFAULT
         WRITE(*,*) 'Error: Chemical species, invalid variable: ' ,  TRIM(SPC_NAME)
         STOP
      ENDSELECT

      t = TSTLEN 
      DO I = 1, NCOLS
         DO J = 1, NROWS
            LAIc = LAIARc(I,J)
            LAIp = LAIARp(I,J)
            Tt   = D_TEMP(I,J)
            IF (LAIp .EQ. LAIc) THEN
               Fnew = 0.0
               Fgro = 0.1
               Fmat = 0.8
               Fold = 0.1
            ELSEIF (LAIp .GT. LAIc) THEN
               Fnew = 0.0
               Fgro = 0.0
               Fold = (LAIp-LAIc) / LAIp
               Fmat = 1-Fold
            ELSEIF (LAIp .LT. LAIc) THEN
               IF (Tt .LE. 303.0) THEN
                  ti = 5.0 + (0.7*(300-Tt))
               ELSEIF (Tt .GT. 303.0) THEN
                  ti = 2.9
               ENDIF
               tm = 2.3*ti
               IF (t .LE. ti) THEN
                  Fnew = 1.0 - (LAIp/LAIc)
               ELSEIF (t .GT. ti) THEN
                  Fnew = (ti/t) * ( 1-(LAIp/LAIc) )
               ENDIF
               IF (t .LE. tm) THEN
                  Fmat = LAIp/LAIc
               ELSEIF (t .GT. tm) THEN
                  Fmat = (LAIp/LAIc) + ( (t-tm)/t ) * ( 1-(LAIp/LAIc) )
               ENDIF
               Fgro = 1.0 - Fnew - Fmat
               Fold = 0.0
            ENDIF
         GAM_A(I,J) = Fnew*Anew(AINDX) + Fgro*Agro(AINDX) + &
                      Fmat*Amat(AINDX) + Fold*Aold(AINDX)
         ENDDO    
      ENDDO      

      RETURN
      END SUBROUTINE GAMMA_A

!-----------------------------------------------------------------------
! Calculate GAM_SMT (GAMMA_SM) the activity factor for soil moisture
!-----------------------------------------------------------------------
       SUBROUTINE GAMMA_S(NCOLS, NROWS, SMOIS, ISLTYP, SPC_NAME, GAM_S)

       IMPLICIT NONE

       INTEGER ::       NCOLS, NROWS, I, J
       CHARACTER*8 ::   SPC_NAME
       REAL(kind=dp) :: GAM_S(NCOLS,NROWS), SMOIS(NCOLS,NROWS)
       REAL(kind=dp) :: SOILW(NCOLS,NROWS) !wilting point
       REAL(kind=dp), PARAMETER :: DSOIL = 0.04
       REAL(kind=dp) :: SOIL1
       INTEGER ::       ISLTYP(NCOLS,NROWS)

       SELECT CASE ( TRIM(SPC_NAME) )

       CASE ('MYRC', 'SABI', 'LIMO', '3CAR', 'OCIM', 'BPIN', 'APIN', &
             'OMTP', 'FARN', 'BCAR', 'OSQT', 'MBO', 'MEOH', 'ACTO', &
             'CH4', 'NO', 'ACTA', 'FORM', 'CO', 'CINE', 'LINA')

            DO I = 1,NCOLS
               DO J= 1,NROWS
                  GAM_S(I,J) = 1.0
               ENDDO
            ENDDO

       CASE ('ISOP')

             DO I =1,NCOLS
                DO J= 1,NROWS
                   IF (ISLTYP(I,J) .EQ. 1) THEN !sand
                       SOILW(I,J) = 0.01 
                   ELSEIF (ISLTYP(I,J) .EQ. 14) THEN
                           SOILW(I,J) = 1.0
                   ELSEIF (ISLTYP(I,J) .EQ. 15 .OR. ISLTYP(I,J) .EQ. 16) THEN
                           SOILW(I,J) = 0.05
                   ELSE
                           SOILW(I,J) = 0.138  !clay soils
                   ENDIF
                   SOIL1 =  SOILW(I,J) + DSOIL
                   IF (SMOIS(I,J) > SOIL1) THEN
                       GAM_S(I,J) = 1.0
                   ELSEIF (SMOIS(I,J) < SOIL1 .AND. SMOIS(I,J) > SOILW(I,J)) THEN
                           GAM_S(I,J) = (SMOIS(I,J) - SOILW(I,J)) / DSOIL
                   ELSE
                           GAM_S(I,J) = 0.0
                   ENDIF
                ENDDO
             ENDDO

       CASE DEFAULT
            WRITE(*,*) 'Error: Chemical species, invalid variable: ' , TRIM(SPC_NAME)
            STOP
       ENDSELECT
       
       RETURN
       END SUBROUTINE GAMMA_S

! =======================================================================
      END MODULE M2_GAMMA_ETC
