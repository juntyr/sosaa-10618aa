MODULE simbim
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: emission_sb, neq

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,307)
  ! DLSODA is in f77 DOUBLE PRECISION, so here we try match that with this

  ! for photosynthesis:
  !   constant parameters:
  REAL(dp), PARAMETER :: &
    alpha  = 1.2, &
    beta   = 133.6271, &
    Gmax   = 61.0, &
    kg     = 0.11, &
    phi    = 0.00185, &
    delta  = 3.5, &
    theta  = 0.6, &
    vLight = 6.5, &
    kLight = 463.13, &
    Rd     = 0.2, &
    ca     = 380.0 


  !   variables (that need to be global, so that are visible inside f):
  REAL(dp) :: L   ! Light
  REAL(dp) :: V

  ! for BIM:
  !   constant parameters:
  REAL(dp), PARAMETER :: &
    ! DXP synthesis
    vmDXP30 = 1.9       ,&
    kMDXPGAP = 6.160    ,&
    kMDXPPGA = 6.160    ,&
    ! MEP synthesis
    vmMEP30 = 8.22      ,&
    kMMEPDXP = 3.73333  ,&
    kMMEPNADPH = 9.3333 ,&
    ! IDP synthesis
    vmIDPs30 = 6.34     ,&
    kMIDPs = 7.840      ,&
    ! IDP isomerase
    vmIDPi30 = 1.9467   ,&
    kMIDPiIDP = 1.30667 ,&
    kMIDPiDMDP = 1.30667    ,&
    keqIDPi = 5.78667       ,&
    ! Isoprene synthase
    vmIs30 = 9.893       ,&
    kMIs = 9.33333        ,&
    ! GDP synthase *)
    vmGDPs30 = 81.57  ,&
    kMGDPDMADP = 1.58667    ,&
    kMGDPIDP = 1.04533     ,&
    ! GGDP synthase
    vmGGDPs30 = 8.157 ,&
    kMGGDPs = 1.58667       ,&
    ! Monoterpene synthase
    vmMTs30 = 0.9893      ,&
    kMMTs = 9.3333       ,&
    ! I don't know what the rest are (Sampo)
    fpga = 0.335  ,&
    kMTP = 80.0   ,&
    lR = 0.0035   ,&
    diso = 0.0085 ,&
    dmono = 0.0016

  ! stuff for DLSODA
  INTEGER, PARAMETER ::  & 
    neq   = 14          ,& ! number of differential equations
    itol  = 1           ,& ! so that atol is a scalar, not array
    itask = 1           ,& ! for normal computation from t to tout
    iopt  = 0           ,& ! for no optional inputs
    lrw   = 22 + neq * MAX(16, neq + 9)  ,& ! real workspace size
    liw   = 20 + neq    ,& ! integer workspace size
    jt    = 2              ! for no user-provided jacobian
  REAL(dp), PARAMETER :: &
    rtol = 1e-5         ,& ! relative tolerance
    atol = 1e-6            ! absolute tolerance
  REAL(dp) :: rwork(lrw)   ! real workspace
  INTEGER  :: iwork(liw)   ! integer workspace

  CONTAINS

    SUBROUTINE emission_sb(y,tin,tout, Lin, Vin)
      IMPLICIT NONE
      REAL(dp), INTENT(inout) :: y(neq)
      REAL(dp), INTENT(in) :: tin,tout, Lin, Vin
      
      ! for DLSODA:
      INTEGER  :: istate  ! a flag
      REAL(dp) :: tt
        ! tt cannot be an input for simbim_step, since DLSODA changes tt
        ! call this tt, since T is already used for leaf temperature below

      istate = 1
      tt = tin

      L = Lin 
      V = Vin 


      CALL DLSODA (f, neq, y, tt, tout, itol, rtol, atol, itask, &
                   istate, iopt, rwork, lrw, iwork, liw, dummy, jt)

    END SUBROUTINE emission_sb

    !!! Stomatal target functions

    SUBROUTINE f (neq, t, y, ydot)
      IMPLICIT NONE
      ! this computes the right-hand side of the system y' = f(t,y), where
      ! y and f are vectors of length neq
      INTEGER,  INTENT(in)  :: neq
      REAL(dp), INTENT(in)  :: t
      REAL(dp), INTENT(in)  :: y(neq)
      REAL(dp), INTENT(out) :: ydot(neq)
      REAL(dp) :: vp1,vp2,vp3,vp4,gs,ci,aps
      REAL(dp) :: vb1,vb2
      REAL(dp) :: v1,v2,v3,v4,v5,v6,v7,v8,v9,v10
      REAL(dp) :: GAP,PGA,NADPH,DXP,MEP,IDP,DMADP,Isoprene,GDP,Mono

      ! photosynthesis model
      gs  = y(1)
      ci  = y(2)
      aps = y(3)

      vp1 = kg * ( G(L,V) - gs )
      vp2 = gs * (ca - ci)
      vp3 = r(L) * ci
      vp4 = theta * aps

      ! BIM
      GAP      = y(4)
      PGA      = y(5)
      NADPH    = y(6)
      DXP      = y(7)
      MEP      = y(8)
      IDP      = y(9)
      DMADP    = y(10)
      Isoprene = y(11)
      GDP      = y(12)
      Mono     = y(13)

      vb1 = fpga * (aps/3.0 + Rd)**2 / (kMTP + (aps/3.0 + Rd))
      vb2 = (1.0 - fpga) * (aps/3.0 + Rd)**2 / (kMTP + (aps/3.0 + Rd))

      v1 = vmDXP30 * GAP / (kMDXPGAP + GAP) * PGA / (kMDXPPGA + PGA)
      v2 = vmMEP30 * NADPH * DXP / (kMMEPDXP*DXP + kMMEPNADPH*NADPH + DXP*NADPH)
      v3 = vmIDPs30 * MEP / (kMIDPs + MEP)
      v4 = vmIDPi30 * (IDP - DMADP/keqIDPi) / &
                      (kMIDPiIDP * (1 + DMADP/kMIDPiDMDP) + IDP)

      v5 = vmIs30 * DMADP / (kMIs + DMADP)
      v6 = vmGDPs30 * DMADP / (kMGDPDMADP + DMADP) * IDP / (kMGDPIDP + IDP) 
      v7 = vmMTs30 * GDP / (kMMTs + GDP)
      v8 = vmGGDPs30 * (GDP / (kMGGDPs + GDP))**2

      v9  = diso * Isoprene
      v10 = dmono * Mono

      ydot(1)  = vp1
      ydot(2)  = vp2 - vp3
      ydot(3)  = vp3 - vp4 - Rd
      ydot(4)  = vb2 - v1
      ydot(5)  = vb1 - v1 
      ydot(6)  = 0.0
      ydot(7)  = v1 - v2
      ydot(8)  = v2 - v3
      ydot(9)  = v3 - v4 - v6
      ydot(10) = v4 - v5 -v6
      ydot(11) = v4 - v9
      ydot(12) = v6 - v7 -v8
      ydot(13) = v7 - v10
      ydot(14) = 60*Mono  ! = y(13)
        ! Mono is the production rate of monoterpenes, in units of
        ! (mol or micromol?) per (second or minute?) per m2 of leaf
        ! area We use y(14) to get the accumulated value, time
        ! integral, of monoterpenes produced diring the time step.

    END SUBROUTINE f

    SUBROUTINE dummy
      ! does nothing
    END SUBROUTINE dummy

    FUNCTION GammL(L)
      IMPLICIT NONE
      REAL(dp) :: GammL,L

      GammL = alpha*L / (beta + L)
    END FUNCTION GammL

    FUNCTION GammVPD(V)
      IMPLICIT NONE
      REAL(dp) :: GammVPD,V

      IF (V > 0) THEN
         GammVPD = 1.0 - V/3.0
      ELSE
         GammVPD = 1.0
      END IF
      ! old: GammVPD = delta / V**en
    END FUNCTION GammVPD

    !!! Plant minimizes water loss
    FUNCTION G(L,V)
      IMPLICIT NONE
      REAL(dp) :: G,L,V
     
      G = Gmax * MIN(GammL(L), GammVPD(V))
    END FUNCTION G

    !!! Photosynthetical target function
    FUNCTION r(L) 
      IMPLICIT NONE
      REAL(dp) :: r,L

      r = phi*vLight*L / (kLight + L)
    END FUNCTION r

    !!! Calculate VPD by temperature, but this is not used because
    !!! measurements are read??
!!$    FUNCTION e(T) 
!!$      IMPLICIT NONE
!!$      REAL(dp) :: e,T
!!$
!!$      e = 101.325*EXP(13.3185*(1.0 - tk100/T) - 1.976*(1.0 - tk100/T)**2 &
!!$          - 0.6445*(1.0 - tk100/T)**3 - 0.1229*(1.0 - tk100/T)**4)
!!$    END FUNCTION e

END MODULE simbim
