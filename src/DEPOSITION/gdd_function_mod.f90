MODULE gdd_function_mod
  USE gdd_parameter_mod
  USE gdd_variable_mod

  IMPLICIT NONE

  INTERFACE read_data_from_messy
    MODULE PROCEDURE read_data_from_messy_1
    MODULE PROCEDURE read_data_from_messy_2
  END INTERFACE

  !!!!! Time
  REAL(dp) :: time
  !!!!! Grid
  INTEGER :: nvegl, ntrac  ! canopy level number, trace number
  REAL(dp), ALLOCATABLE :: zz(:)  ! [m], (nvegl)
  !!!!! Meteorology
  REAL(dp), ALLOCATABLE :: PARl(:), PARl_sl(:), PARl_sh(:)  ! [m2 m-2], (nvegl)
  REAL(dp), ALLOCATABLE :: fsl(:)  ! [-], (nvegl)
  REAL(dp), ALLOCATABLE :: u_veg(:)
  REAL(dp), ALLOCATABLE :: ustveg(:)  ! [m s-1], friction velocity
  REAL(dp), ALLOCATABLE :: rho(:), pres(:), Tc(:), RH(:) ! RH ~ [0,1]
  !!!!! Canopy
  REAL(dp), ALLOCATABLE :: LAIl(:), LAIl_sl(:), LAIl_sh(:)  ! [m2 m-2], LAI at each layer (nvegl)
  REAL(dp) :: hc  ! [m], canopy height
  REAL(dp), ALLOCATABLE :: frac_veg(:), frac_ws(:)  ! (nvegl)
  !!!!! Soil parameters
  REAL(dp) :: ws, wsmax, wswilt  ! [m3 m-3], soil moisture, field capacity, wilt point
  !!!!! Chemistry
  CHARACTER(LEN=20), ALLOCATABLE :: trname(:)  ! chemical names, (ntrac)
  REAL(dp), ALLOCATABLE :: molar_mass(:)  ! [g mol-1], (ntrac)
  REAL(dp), ALLOCATABLE :: henry(:)  ! [M atm-1], [M] = [mol -L], effective Henry's law constant
  REAL(dp), ALLOCATABLE :: dryreac(:)  ! [-], reactivity coefficient [0:non-react., 0.1:semi react., 1:react.]  
  INTEGER :: ind_O3, ind_SO2
  !!!!! Other
  REAL(dp) :: stomblock  ! [-], if the stomata block is considered
  LOGICAL :: l_drydep, l_vpd, l_wetskin
  INTEGER :: jk, jt  ! jk for level loop, jt for trace gases loop

  REAL(dp), PARAMETER :: R=8.3144598d0  ! [J mol-1 K-1], gas constant
  REAL(dp), PARAMETER :: R1=0.082057338d0  ! [atm M-1 K-1], gas constant in another unit
  REAL(dp), PARAMETER :: vk=0.41d0  ! [-], von Karman constant
  REAL(dp), PARAMETER :: D_H2O=2.4d-5 ! [m2 s-1], molecular diffusivity of water vapor in air
                                      ! http://webserver.dmt.upm.es/~isidoro/dat1/Mass%20diffusivity%20data.pdf
                                      ! Other D value can be found in Seinfeld and Pandis (P.921, ACP second edition)
  REAL(dp), PARAMETER :: MOLAR_MASS_H2O = 18.02d0  ! [g mol-1], molar mass of water vapor
  REAL(dp), PARAMETER :: MOLAR_MASS_AIR = 28.97d0  ! [g mol-1], molar mass of water vapor
  REAL(dp), PARAMETER :: Rspe_air=R/(MOLAR_MASS_AIR*1.0d-3)  ! [J kg-1 K-1], specific gas constant of air
  REAL(dp), PARAMETER :: r_max = 1.0d10  ! [s m-1], default maximum resistance value

  PRIVATE
  PUBLIC :: get_gas_dry_deposition_velocity
CONTAINS

SUBROUTINE get_gas_dry_deposition_velocity(time_in, &
                                           nvegl_in, ntrac_in, &
                                           l_drydep_in, l_vpd_in, l_wetskin_in, &
                                           frac_veg_in, frac_ws_in, &
                                           zz_in, hc_in, LAIl_in, LAIl_sl_in, LAIl_sh_in, &
                                           PARl_in, PARl_sl_in, PARl_sh_in, fsl_in, &
                                           u_veg_in, ustveg_in, rho_in, pres_in, Tc_in, RH_in, &
                                           ws_in, wsmax_in, wswilt_in, &
                                           stomblock_in, &
                                           gs_h2o_sn_in, gs_h2o_sd_in, gs_h2o_avg_in, &
                                           gstm_h2o_sn_in, gstm_h2o_sd_in, gstm_h2o_avg_in, &
                                           trname_in, molar_mass_in, henry_in, dryreac_in, ind_O3_in, ind_SO2_in, &
                                           vdep, dep_output)
  REAL(dp), INTENT(IN) :: time_in
  INTEGER, INTENT(IN) :: nvegl_in, ntrac_in
  LOGICAL, INTENT(IN) :: l_drydep_in, l_vpd_in, l_wetskin_in
  REAL(dp), INTENT(IN) :: frac_veg_in(nvegl_in), frac_ws_in(nvegl_in)  ! [-], (nvegl)
  REAL(dp), INTENT(IN) :: zz_in(nvegl_in)  ! [m], (nvegl)
  REAL(dp), INTENT(IN) :: hc_in  ! [m]
  REAL(dp), INTENT(IN) :: LAIl_in(nvegl_in), LAIl_sl_in(nvegl_in), LAIl_sh_in(nvegl_in)  ! [m2 m-2], (nvegl)
  REAL(dp), INTENT(IN) :: PARl_in(nvegl_in), PARl_sl_in(nvegl_in), PARl_sh_in(nvegl_in)  ! [W m-2], (nvegl)
  REAL(dp), INTENT(IN) :: fsl_in(nvegl_in)
  REAL(dp), INTENT(IN) :: u_veg_in(nvegl_in), ustveg_in(nvegl_in)  ! [m s-1]
  REAL(dp), DIMENSION(nvegl_in), INTENT(IN) :: rho_in, pres_in, Tc_in, RH_in ! RH ~ [0,1]
  REAL(dp), INTENT(IN) :: ws_in, wsmax_in, wswilt_in  ! [m3 m-3]
  REAL(dp), INTENT(IN) :: stomblock_in
  REAL(dp), INTENT(IN) :: gs_h2o_sn_in(nvegl_in), gs_h2o_sd_in(nvegl_in), gs_h2o_avg_in(nvegl_in)
  REAL(dp), INTENT(IN) :: gstm_h2o_sn_in(nvegl_in), gstm_h2o_sd_in(nvegl_in), gstm_h2o_avg_in(nvegl_in)

  CHARACTER(LEN=20), INTENT(IN) :: trname_in(ntrac_in)
  REAL(dp), INTENT(IN) :: molar_mass_in(ntrac_in)
  REAL(dp), INTENT(IN) :: henry_in(ntrac_in)
  ! REAL(dp), INTENT(IN) :: henry(ntrac, nvegl)
  REAL(dp), INTENT(IN) :: dryreac_in(ntrac_in)
  INTEGER, INTENT(IN) :: ind_O3_in, ind_SO2_in
  REAL(dp) :: dep_output(nvegl_in,ntrac_in,20)

  REAL(dp), INTENT(OUT) :: vdep(nvegl_in,ntrac_in)  ! [m s-1], (nvegl, ntrac)

  LOGICAL, DIMENSION(ntrac_in) :: l_gas, l_atmbios, l_vdbigl

  REAL(dp), DIMENSION(ntrac_in) :: molec_diff  ! molecular diffusivity
  !!!!! Schmidt number, Sc = nu/molec_diff, nu is the kinematic viscosity ([m2 s-1])
  !!!!! Since the viscosity is calculated mostly by interacting with air molecules,
  !!!!! we use nu_air for all the gases in the atmosphere
  REAL(dp), DIMENSION(ntrac_in) :: Sc
  ! REAL(dp), PARAMETER :: 

  REAL(dp), DIMENSION(ntrac_in) :: diff, diffrb
  REAL(dp), DIMENSION(ntrac_in) :: rbs, rsoil, rwater, rws, rsnow, rmes, rcut
  REAL(dp), DIMENSION(nvegl_in, ntrac_in):: rbveg, rleaf, rleafw, rstm, rstm_sl, rstm_sh, ccomp, rs, rsveg, rswet  ! [s m-1], (nvegl, ntrac)
  REAL(dp), DIMENSION(nvegl_in):: rstm_h2o, rstm_h2o_sl, rstm_h2o_sh, fvpd, fvpd_sl, fvpd_sh
  REAL(dp) :: fws, rac

  REAL(dp) :: rtot(nvegl_in,ntrac_in)  ! [s m-1], (nvegl, ntrac)
  REAL(dp) :: rup, rdown

  INTEGER :: method  ! method to calculate rtot
  INTEGER :: m_rstmh2o  ! 1: LAIl=1.0; 2: LAIl=LAIl
  INTEGER :: m_wetskin  ! 1: frac_ws=0; 2: frac_ws=param

  LOGICAL, SAVE :: FIRST_TIME=.TRUE.

  !=============================================================!
  ! Initiation for the gas dry deposition model
  !=============================================================!
  IF (FIRST_TIME) THEN
    !===== Allocate =====!
    nvegl = nvegl_in
    ntrac = ntrac_in
    ALLOCATE(zz(nvegl), LAIl(nvegl), LAIl_sl(nvegl), LAIl_sh(nvegl), PARl(nvegl), PARl_sl(nvegl), PARl_sh(nvegl), fsl(nvegl))
    ALLOCATE(u_veg(nvegl), ustveg(nvegl), rho(nvegl), pres(nvegl), Tc(nvegl), RH(nvegl))
    ALLOCATE(frac_veg(nvegl), frac_ws(nvegl))
    ALLOCATE(trname(ntrac), molar_mass(ntrac), henry(ntrac), dryreac(ntrac))
    !===== Set values only once =====!
    zz = zz_in
    hc = hc_in
    trname = trname_in
    molar_mass = molar_mass_in
    henry = henry_in
    dryreac = dryreac_in
    ind_O3 = ind_O3_in
    ind_SO2 = ind_SO2_in

    FIRST_TIME=.FALSE.
  END IF
  !===== Obtain values every time step =====!
  time = time_in
  DO jk=1,nvegl
    LAIl(jk) = MAX(1.0d-10, LAIl_in(jk))
    LAIl_sl(jk) = MAX(1.0d-10, LAIl_sl_in(jk))
    LAIl_sh(jk) = MAX(1.0d-10, LAIl_sh_in(jk))
  END DO
  PARl = PARl_in
  PARl_sl = PARl_sl_in
  PARl_sh = PARl_sh_in
  u_veg = u_veg_in
  ustveg = ustveg_in
  rho = rho_in
  pres = pres_in
  Tc = Tc_in
  RH = RH_in
  frac_veg = frac_veg_in
  frac_ws = frac_ws_in
  ws = ws_in
  wsmax = wsmax_in
  wswilt = wswilt_in
  stomblock = stomblock_in
  l_drydep = l_drydep_in
  l_vpd = l_vpd_in
  l_wetskin = l_wetskin_in
  !===== Intiate chemical properties =====!
  CALL tracer_init_new(l_gas, l_atmbios, l_vdbigl)

  !=============================================================!
  ! Calculate gas dry deposition
  !=============================================================!
  method = 3
  m_rstmh2o = 1
  m_wetskin = 2
  IF (l_drydep) THEN
    SELECT CASE (method)
    CASE (1)  ! consider total PARl and LAIl
      !===== Get diff, diffrb, rsoil, rwater, rws, rsnow, rmes, rcut =====!
      CALL get_uptake_resistance(diff, diffrb, rsoil, rwater, rws, rsnow, rmes, rcut, Tc)

      !===== Get stomatal resistance rstm =====!
      CALL get_rstm_h2o(LAIl, PARl, gstm_h2o_avg_in, m_rstmh2o, rstm_h2o)

      CALL get_fvpd(ustveg(nvegl), Tc, RH, rstm_h2o, l_vpd, fvpd)
      CALL get_fws(ws, wsmax, wswilt, fws)
      !***** use rstm_h2o from SCADIS
      rstm_h2o = 1.0d0/(2.7d0*gstm_h2o_avg_in)
      fvpd = 1.0d0
      CALL get_rstm(diff, PARl, rstm_h2o, fvpd, fws, rstm)

      !===== Get rac, rbveg =====!
      CALL get_rac(hc, LAIl, ustveg(1), rac)
      CALL get_rbveg(diffrb, u_veg, rbveg)

      !===== Get frac_ws =====!
      CALL get_frac_ws(RH, m_wetskin, frac_ws, frac_veg)
      ! write(*,*) 'RH', RH(1:5)
      ! write(*,*) 'frac_ws', frac_ws(1:5)

      !===== Get rleaf, rleafw =====!
      CALL get_ccomp(frac_veg, rho, Tc, ccomp)
      CALL get_rleaf(rcut, rstm, rmes, ccomp, rleaf)
      CALL get_rleafw(rws, rstm, rmes, stomblock, ccomp, rleafw)

      !===== Get rsveg, rswet =====!
      CALL get_rs(rbveg, rleaf, rleafw, rsoil, rac, LAIl, rsveg, rswet)

      !===== Get rtot =====!
      CALL get_rtot(frac_veg, frac_ws, rsveg, rswet, rtot)
    CASE (2)  ! consider sunlit and shaded part
      !===== Get diff, diffrb, rsoil, rwater, rws, rsnow, rmes, rcut =====!
      CALL get_uptake_resistance(diff, diffrb, rsoil, rwater, rws, rsnow, rmes, rcut, Tc)

      !===== Get stomatal resistance rstm =====!
      CALL get_rstm_h2o(LAIl_sl, PARl_sl, gstm_h2o_avg_in, m_rstmh2o, rstm_h2o_sl)
      CALL get_rstm_h2o(LAIl_sh, PARl_sh, gstm_h2o_avg_in, m_rstmh2o, rstm_h2o_sh)

      CALL get_fvpd(ustveg(nvegl), Tc, RH, rstm_h2o_sl, l_vpd, fvpd_sl)
      CALL get_fvpd(ustveg(nvegl), Tc, RH, rstm_h2o_sh, l_vpd, fvpd_sh)

      CALL get_fws(ws, wsmax, wswilt, fws)

      CALL get_rstm(diff, PARl, rstm_h2o_sl, fvpd_sl, fws, rstm_sl)
      CALL get_rstm(diff, PARl, rstm_h2o_sh, fvpd_sh, fws, rstm_sh)

      !===== Get rac, rbveg =====!
      CALL get_rac(hc, LAIl, ustveg(1), rac)
      CALL get_rbveg(diffrb, u_veg, rbveg)

      !===== Get frac_ws =====!
      CALL get_frac_ws(RH, m_wetskin, frac_ws, frac_veg)

      !===== Get rtot =====!
      DO jt=1,ntrac
        rtot(1,jt) = 1.0d0/( &
                              1.0d0/( rbveg(1,jt)/LAIl(1) + &
                                      1.0d0/( 2.0d0*LAIl(1)*frac_veg(1)/rcut(jt) + 2.0d0*LAIl(1)*frac_ws(1)/rws(jt) + &
                                              LAIl_sl(1)/(rstm_sl(1,jt)+rmes(jt)) + LAIl_sh(1)/(rstm_sh(1,jt)+rmes(jt)) ) &
                                    ) + &
                              1.0d0/( rac + rsoil(jt) ) &
                            )

        DO jk=2,nvegl
          rtot(jk,jt) = rbveg(jk,jt)/LAIl(jk) + &
                        1.0d0/( 2.0d0*LAIl(jk)*frac_veg(jk)/rcut(jt) + 2.0d0*LAIl(jk)*frac_ws(jk)/rws(jt) + &
                                LAIl_sl(jk)/(rstm_sl(jk,jt)+rmes(jt)) + LAIl_sh(jk)/(rstm_sh(jk,jt)+rmes(jt)) &
                              )
        END DO
      END DO
    CASE (3)
      !===== Get diff, diffrb, rsoil, rwater, rws, rsnow, rmes, rcut =====!
      CALL get_uptake_resistance(diff, diffrb, rsoil, rwater, rws, rsnow, rmes, rcut, Tc)

      ! CALL get_fvpd(ustveg(nvegl), Tc, RH, rstm_h2o, l_vpd, fvpd)
      l_vpd = .False.
      CALL get_fvpd(ustveg, Tc, RH, rstm_h2o, l_vpd, fvpd)  ! use ustveg of each layer
      CALL get_fws(ws, wsmax, wswilt, fws)

      !===== Get stomatal resistance rstm =====!
      m_rstmh2o = 3
      CALL get_rstm_h2o(LAIl, PARl, gstm_h2o_avg_in, m_rstmh2o, rstm_h2o)
      fvpd = 1.0d0
      CALL get_rstm(diff, PARl, rstm_h2o, fvpd, fws, rstm)

      !===== Get rac, rbveg =====!
      CALL get_rac(hc, LAIl, ustveg(1), rac)
      CALL get_rbveg(diffrb, u_veg, rbveg)
      CALL get_rbs(ustveg(1), pres(1), Tc(1), molar_mass, rbs)

      !===== Get frac_ws =====!
      CALL get_frac_ws(RH, m_wetskin, frac_ws, frac_veg)

      !===== Get rleaf, rleafw =====!
      CALL get_ccomp(frac_veg, rho, Tc, ccomp)
      CALL get_rleaf(rcut, rstm, rmes, ccomp, rleaf)
      CALL get_rleafw(rws, rstm, rmes, stomblock, ccomp, rleafw)

      DO jt=1,ntrac
        rtot(:, jt) = ( rbveg(:,jt) + 1.0d0/( frac_veg/rleaf(:,jt) + frac_ws/rleafw(:,jt) ) )/LAIl
        rup = ( rbveg(1,jt) + 1.0d0/( frac_veg(1)/rcut(jt) + frac_ws(1)/rws(jt) ) )/(0.5d0*LAIl(1))
        rdown = ( rbveg(1,jt) + &
          1.0d0/( frac_veg(1)/rcut(jt) + frac_ws(1)/rws(jt) + 1.0d0/(rstm(1,jt)+rmes(jt)) ) ) &
          / (0.5d0*LAIl(1))
        rtot(1, jt) = 1.0d0/(1.0d0/rup + 1.0d0/rdown)
      END DO
      rtot(1,:) = 1.0d0/( 1.0d0/(rac+rbs+rsoil) + 1.0d0/rtot(1,:) )
      ! write(*,*) 'rleaf', rleaf(nvegl, ind_O3)
      ! write(*,*) 'rtot', rtot(nvegl, ind_O3)
      ! write(*,*) 'LAIl', LAIl(nvegl)

      !===== Get rsveg, rswet =====!
      ! CALL get_rs(rbveg, rleaf, rleafw, rsoil, rac, LAIl, rsveg, rswet)

      !===== Get rtot =====!
      ! CALL get_rtot(frac_veg, frac_ws, rsveg, rswet, rtot)

    END SELECT

    !===== Get vdep =====!
    vdep = 1.0d0/rtot
  ELSE
    vdep = 0.0d0
  END IF
  dep_output = 0.0d0
  dep_output(:,:,1) = rac
  dep_output(:,:,2) = rstm
  dep_output(:,:,3) = rbveg
  dep_output(:,ind_O3,4) = frac_ws(:)
  dep_output(:,:,5) = rleaf
  dep_output(:,:,6) = rleafw
  dep_output(:,:,7) = rsveg
  dep_output(:,:,8) = rswet
  dep_output(:,:,9) = rtot
  dep_output(:,ind_O3,10) = fvpd
  dep_output(:,ind_O3,11) = rstm_h2o
  dep_output(:,ind_O3,12) = gs_h2o_sn_in
  dep_output(:,ind_O3,13) = gs_h2o_sd_in
  dep_output(:,ind_O3,14) = gs_h2o_avg_in
  ! write(*,*) 'dep_output', RH(1:5)
  dep_output(:,ind_O3,15) = RH
  dep_output(1,:,16) = rbs
  ! write(*,*) 'dep_output 1', dep_output(1:5,ind_O3,15)
  ! write(*,*) 'RH', RH(1:5)
  ! write(*,*) 'frac_ws', frac_ws(1:5)
END SUBROUTINE get_gas_dry_deposition_velocity


SUBROUTINE tracer_init_new(l_gas, l_atmbios, l_vdbigl)
  ! declaration of tracer properties that are used to calculate their surface uptake resistances
  LOGICAL, INTENT(OUT) :: l_gas(ntrac), l_atmbios(ntrac), l_vdbigl(ntrac)

  INTEGER, PARAMETER :: AIR=1, AEROSOL=2
  INTEGER :: jt, jk

  !===========================================!
  ! l_gas
  !===========================================!
  l_gas = .TRUE.

  !===========================================!
  ! l_atmbios
  !===========================================!
  l_atmbios = .FALSE.
  DO jt=1,ntrac
    IF (INDEX(trname(jt), ' ') > 1) l_atmbios(jt)=.TRUE.
  END DO

  !===========================================!
  ! l_vdbigl
  !===========================================!
  l_vdbigl = .FALSE.
  WHERE (l_gas .AND. molar_mass>0.0_dp .AND. henry>0.0_dp)
    l_vdbigl = .TRUE.
  END WHERE
END SUBROUTINE tracer_init_new


SUBROUTINE get_uptake_resistance(diff, diffrb, rsoil, rwater, rws, rsnow, rmes, rcut, Tc)
  !-----------------------------------------------------------------------------
  ! subroutine emdep_xtsurf_calc_rs, to calculate the values of the uptake resistances
  ! required to calculate the trace gas dry deposition velocity. this routine
  ! is based on an approach by wesely, 1989, in which the uptake resistances of
  ! trace gases, for which the dry deposition velocities have not been observed,
  ! are estimated based on the henry coefficient and a reactivity coefficient and
  ! the uptake resistances of so2 and o3, of the "big leaf" dry deposition scheme
  ! by ganzeveld and j. lelieveld j. geophys. res., 100, 20,999-21,012,1995,
  ! ganzeveld et al.,j. geophys. res., 103, 5679-5694, 1998 and ganzeveld et al,
  ! submitted to j. geophys. res., 2001. for more information of the wesely
  ! approach see atmospheric environment vol 23, no 6, 1293-1304.
  !
  ! the program needs as input data the molecular mass of the defined trace
  ! gases, the henry coefficient [m atm-1] and an estimated reactivity
  ! coefficient which has 3 distinct values: 0 for non-reactive species,
  ! (e.g, so2, acetaldehyde), 0.1 for moderately reactive species (e.g., pan),
  ! and 1 for reactive species (e.g., o3, hno3). these values are defined in
  ! the module mo_*_request_tracer
  !-----------------------------------------------------------------------------

  ! mz_lg_20040503+ 
  ! Interface:
  ! ----------
  ! input 
  ! lo_derived : true whenever SO2 and O3 are defined as tracers
  ! lvd_bigl   : true for the tracer when the required parameters are defined
  ! trname     : tracer name
  ! molar_mass : molecular weight 
  ! dryreac    : reactivity coefficient [0:non-react., 0.1:semi react., 1:react.]        
  ! henry      : henry coefficient [M atm-1]
  REAL(dp), INTENT(OUT):: diff(ntrac), diffrb(ntrac), &
                          rsoil(ntrac), rwater(ntrac), rws(ntrac), rsnow(ntrac), rmes(ntrac), rcut(ntrac)

  ! declarations

  INTEGER  :: jt
  REAL(dp) :: diff_so2, diffrb_so2,  &
              rsoil_so2, rwater_so2, rws_so2, rsnow_so2, rmes_so2, rcut_so2,  &
              diff_o3, diffrb_o3,  &
              rsoil_o3, rwater_o3, rws_o3, rsnow_o3, rmes_o3, rcut_o3
  REAL(dp) :: henry_so2
  REAL(dp) :: henry_now  ! H value of current compound

  REAL(dp) :: Tc(nvegl), Tc_avg

  INTEGER :: m_H

  ! mz_lg_20030811+, SO2
  diff_so2=1.9_dp
  diffrb_so2=1.6_dp
  rsoil_so2=250.0_dp
  rwater_so2=1._dp
  rws_so2=100._dp
  rsnow_so2=1._dp
  rmes_so2=1._dp
  rcut_so2=1.e5_dp

  ! mz_lg_20030811+, O3, Ganzeveld 1995 JGR, (Tab. 1)
  diff_o3=1.6_dp
  diffrb_o3=1.2_dp
  rsoil_o3=400._dp
  rwater_o3=2000._dp
  rws_o3=2000.0d0  ! default: 2000._dp, suggested by Laurens later: 750
  ! mz_lg_20030111+ modified for sensitivity analyses
  ! rsnow_o3=750._dp
  rsnow_o3=2000._dp
  rmes_o3=0.0_dp
  rcut_o3=1.0d5  ! 1.0d5

  Tc_avg = SUM(Tc)/nvegl
  ! mz_lg_20020115 the resistances are only calculated for species that
  !     have not been included in the original dry deposition scheme.
  DO jt=1, ntrac
    ! mz_lg_20030811+, modified to deal with the fact that there might
    !     multiple tracers starting with the same basename but with
    !     with different subnames, e.g., SO2 (e.g., SO2, SO2_GM7)
    SELECT CASE(TRIM(trname(jt))) ! mz_pj_20040330

    CASE('SO2')
      diff(jt)=diff_so2
      diffrb(jt)=diffrb_so2
      rsoil(jt)=rsoil_so2
      rwater(jt)=rwater_so2
      rws(jt)=rws_so2
      rsnow(jt)=rsnow_so2
      rmes(jt)=rmes_so2
      rcut(jt)=rcut_so2
    CASE('O3')
      diff(jt)=diff_o3
      diffrb(jt)=diffrb_o3
      rsoil(jt)= rsoil_o3
      rwater(jt)=rwater_o3
      rws(jt)=rws_o3
      rsnow(jt)=rsnow_o3
      rmes(jt)=rmes_o3
      rcut(jt)=rcut_o3
    CASE('SO4')
      diff(jt)=2.7_dp
      diffrb(jt)=1.8_dp
      rsoil(jt)=1.e5_dp
      rwater(jt)=1.e5_dp
      rws(jt)=1.e5_dp
      rsnow(jt)=1.e5_dp
      rmes(jt)=1.e5_dp
      rcut(jt)=1.e5_dp
    ! CASE('HNO3')
    !   diff(jt)=1.9_dp
    !   diffrb(jt)=1.4_dp
    !   rsoil(jt)=1._dp
    !   rwater(jt)=1._dp
    !   rws(jt)=1._dp
    !   rsnow(jt)=1._dp
    !   rmes(jt)=1._dp
    !   rcut(jt)=1._dp
    CASE('NO')
      diff(jt)=1.3_dp
      diffrb(jt)=1.1_dp
      rsoil(jt)=1.e5_dp
      rwater(jt)=1.e5_dp
      rws(jt)=1.e5_dp
      rsnow(jt)=1.e5_dp
      rmes(jt)=500._dp
      rcut(jt)=1.e5_dp
    CASE('NO2')
      diff(jt)=1.6_dp
      diffrb(jt)=1.2_dp
      rsoil(jt)=600._dp
      rwater(jt)=1.e5_dp
      rws(jt)=1.e5_dp
      rsnow(jt)=1.e5_dp
      rmes(jt)=1._dp
      rcut(jt)=1.e5_dp
    ! ESS_lg_20130503+  
    CASE('CO2')
      diff(jt)=1.6_dp
      diffrb(jt)=1.2_dp
      rsoil(jt)=1.e5_dp
      rwater(jt)=1.e5_dp
      rws(jt)=1.e5_dp
      rsnow(jt)=1.e5_dp
      rmes(jt)=1._dp
      rcut(jt)=1.e5_dp
    ! ESS_lg_20130503-  
    CASE DEFAULT
      ! mz_lg_20030411+ modified: only O3 and SO2 need to be defined for estimating
      !     the various resistances, not the other species of the default scheme

      ! mz_lg_20040504+ lo_derived is true whenever O3 is present. SO2 is now
      !     longer required because of the above defined specific resistances for
      !     SO2 and O3

      ! write(*,*) trname(jt)
      ! henry_so2 = henry(ind_SO2)  ! new method
      henry_so2 = 1.0d5  ! old method

      ! IF (TRIM(trname(jt)) == 'HNO3') THEN
      !   henry_now = 1.0d+14  ! Wesely, 1989; Table S4, Nguyen et al., 2015
      ! ELSEIF (TRIM(trname(jt)) == 'H2O2') THEN
      !   henry_now = 5.0d+7  ! Table S4, Nguyen et al., 2015
      ! ELSE
      !   henry_now = henry(jt)
      ! END IF
      henry_now = henry(jt)

      ! calculation of term which is used to correct the stomatal resistance
      ! for differences in the diffusitivy (see also equation 4).
      diff(jt)=SQRT(molar_mass(jt)/18._dp)     ! mz_pj_20040330

      ! calculation of the term to correct for differences in diffusivity
      ! between water vapor and the trace gas. it is calculated from:
      ! diff bl=(v/dx)**2/3, with v=0.189 sm-1 and dx= dh2o/sqrt(mh2o/mx),
      ! with dh2o=0.212
      diffrb(jt)=(0.189_dp/(0.212_dp/(diff(jt)+Eps)))**(2._dp/3._dp)

      m_H = 2
      SELECT CASE(m_H) ! mz_pj_20040330
      CASE(1)
        ! calculation of rmx, the mesophyll resistance
        rmes(jt)=1.0d0/(henry_now/3000._dp+100._dp*dryreac(jt) + Eps)

        ! calculation of rlux, the cuticular resistance, equation 7 of wesely's paper
        rcut(jt)=1.0d0/(henry_now/henry_so2+dryreac(jt)+Eps)*rcut_o3  ! substitute 1.e-5 with 1/henry_SO2

        ! calculation of rgsx, the soil resistance, equation 9 of wesely's paper
        rsoil(jt)=1.0d0 /( henry_now/(henry_so2*rsoil_so2) + dryreac(jt)/rsoil_o3 + Eps )  ! substitute 1.e+5 with henry_SO2

        ! the snow resistance is similar as the soil resistance
        ! mz_lg_20050628+ modified to assure calculations of rsnow for
        ! for species with a Henry and reactivity of that of ozone similar
        ! to the rsnowO3 of 2000 s m-1 (based on feedback from A. Pozzer)
        rsnow(jt)=1.0d0/(henry_now/(henry_so2*rsnow_so2) + dryreac(jt)/rsnow_o3 + Eps)  ! substitute 1.e+5 with henry_SO2

        ! calculation of rlux-wet, the wet skin resistance, equation 14 of wesely's paper
        rws(jt)=1.0d0/(1.0d0/(3.0d0*rws_so2)+henry_now/(100.0d0*henry_so2) + dryreac(jt)/rws_o3 + Eps)  ! substitute 1.e-7 with 1/(100*henry_SO2)

        ! calculation of sea uptake resistance, using equation 9 of wesely's paper
        rwater(jt)=1.0d0/(henry_now/(henry_so2*rwater_so2)+ dryreac(jt)/rwater_o3 + Eps)
      CASE(2)  ! method with intrinsic H values instead of effective H (Nguyen et al., 2015, PNAS)
        ! calculation of rmx, the mesophyll resistance
        rmes(jt)=1.0d0/(henry_now/(50.0d0*R1*Tc_avg)+100.0d0*dryreac(jt) + Eps)

        ! calculation of rlux, the cuticular resistance, equation 7 of wesely's paper
        rcut(jt)=1.0d0/(henry_now*1.0d-4/(R1*Tc_avg)+dryreac(jt)+Eps)*rcut_o3  ! substitute 1.e-5 with 1/henry_SO2

        ! calculation of rgsx, the soil resistance, equation 9 of wesely's paper
        rsoil(jt)=1.0d0 /( henry_now*1.0d-4/(R1*Tc_avg)/rsoil_so2 + dryreac(jt)/rsoil_o3 + Eps )  ! substitute 1.e+5 with henry_SO2

        ! the snow resistance is similar as the soil resistance
        ! mz_lg_20050628+ modified to assure calculations of rsnow for
        ! for species with a Henry and reactivity of that of ozone similar
        ! to the rsnowO3 of 2000 s m-1 (based on feedback from A. Pozzer)
        rsnow(jt)=1.0d0/(henry_now*1.0d-4/(R1*Tc_avg)/rsnow_so2 + dryreac(jt)/rsnow_o3 + Eps)  ! substitute 1.e+5 with henry_SO2

        ! calculation of rlux-wet, the wet skin resistance, equation 14 of wesely's paper
        rws(jt)=1.0d0/(1.0d0/(3.0d0*rws_so2)+henry_now/100.0d0*1.0d-4/(R1*Tc_avg) + dryreac(jt)/rws_o3 + Eps)  ! substitute 1.e-7 with 1/(100*henry_SO2)

        ! calculation of sea uptake resistance, using equation 9 of wesely's paper
        rwater(jt)=1.0d0/(henry_now*1.0d-4/(R1*Tc_avg)/rwater_so2+ dryreac(jt)/rwater_o3 + Eps)
      END SELECT  ! m_H
    END SELECT
  END DO  ! jt=1, ntrac
END SUBROUTINE get_uptake_resistance


SUBROUTINE get_rstm_h2o(LAIl, PARl, gstm_h2o, m_rstmh2o, rstm_h2o)
  REAL(dp), INTENT(IN) :: LAIl(nvegl)
  REAL(dp), INTENT(IN) :: PARl(nvegl)
  REAL(dp), INTENT(IN) :: gstm_h2o(nvegl)
  INTEGER, INTENT(IN) :: m_rstmh2o

  REAL(dp), INTENT(OUT) :: rstm_h2o(nvegl)

  REAL(dp), PARAMETER :: a = 5000.0d0  ! [J m-3], constant to define the stomatal resistance
  REAL(dp), PARAMETER :: b = 10.0d0  ! [W m-2], constant to define the stomatal resistance
  REAL(dp), PARAMETER :: c = 100.0d0  ! [s m-1], minimum stomatal resistance
  REAL(dp), PARAMETER :: k = 0.9d0  ! [-]
  REAL(dp) :: d

  INTEGER :: jk

  !========================================================!
  ! H2O stomatal resistance
  !========================================================!
  IF (m_rstmh2o < 3) THEN
    DO jk=1,nvegl
      ! IF (PARl(jk) > 1.0d-3 .AND. LAIl(jk)>1.0d-2) THEN
      IF (PARl(jk) > 1.0d-3) THEN
        d = (a+b*c)/(c*PARl(jk))

        SELECT CASE (m_rstmh2o)
        CASE (1)  ! LAIl=1.0d0 for each layer
          rstm_h2o(jk) = k*c / &
                         ( b/(d*PARl(jk))*LOG( (d*EXP(k)+1.0_dp)/(d+1.0_dp) ) - LOG( (d+EXP(-k))/(d+1.0_dp) ) )
        CASE (2)  ! multilayer LAIl
          rstm_h2o(jk) = k*c / &
                         ( b/(d*PARl(jk))*LOG( (d*EXP(k*LAIl(jk))+1.0_dp)/(d+1.0_dp) ) - LOG( (d+EXP(-k*LAIl(jk)))/(d+1.0_dp) ) + Eps)
        END SELECT
      ELSE
        rstm_h2o(jk)=1.0d5
      END IF
    END DO
  !== use SCADIS
  ELSEIF (m_rstmh2o == 3) THEN
    rstm_h2o = 1.0d0/gstm_h2o
  END IF
END SUBROUTINE get_rstm_h2o


ELEMENTAL FUNCTION saturation_vapor_pressure(T)
  REAL(dp), INTENT(IN) :: T
  REAL(dp) :: saturation_vapor_pressure

  saturation_vapor_pressure = &
    6.108_dp*10_dp**(7.5_dp*(T-273.15_dp)/(237.3_dp+(T-273.15_dp)))  ! [hPa]
END FUNCTION saturation_vapor_pressure


SUBROUTINE get_fvpd(ustveg, Tc, RH, rstm_h2o, l_vpd, fvpd)
  REAL(dp), INTENT(IN) :: ustveg(nvegl), Tc(nvegl), RH(nvegl)
  REAL(dp), INTENT(IN) :: rstm_h2o(nvegl)
  LOGICAL, INTENT(IN) :: l_vpd

  REAL(dp), INTENT(OUT) :: fvpd(nvegl)

  REAL(dp) :: rveg(nvegl)
  REAL(dp) :: vpair(nvegl), vpleaf(nvegl), vpsfc(nvegl)

  INTEGER :: jk

  IF (l_vpd) THEN
    rveg = 2.0d0/(0.4d0*ustveg)
    vpair = saturation_vapor_pressure(Tc) * RH
    vpleaf = saturation_vapor_pressure(Tc)
    vpsfc = (vpleaf-vpair)*rveg/(rstm_h2o+rveg)+vpair ! Default: ESS_lg_20130516+ using rstm_h2o of top-layer
    ! vpsfc = vpair ! ESS_lg_20130516+ using rstm_h2o of top-layer

    DO jk=1, nvegl
      fvpd(jk) = MIN(MAX(1.0d0 - 0.02d0*(vpleaf(jk)-vpsfc(jk)), 0.5d0), 1.0d0)
    END DO
  ELSE
    fvpd = 1.0d0
  END IF
END SUBROUTINE get_fvpd


!========================================================!
!
!  soil moisture stress attenuation function [0-1]
!
!========================================================!
SUBROUTINE get_fws(ws, wsmax, wswilt, fws)
  REAL(dp), INTENT(IN) :: ws, wsmax, wswilt
  REAL(dp), INTENT(OUT) :: fws
  !===== Original version in Laurens' code: ESS_lg_20130127+ to avoid Fws > 1, soil moisture stress function, from MESSY =====!
  ! fws =MAX(0._dp, MIN(1._dp,  (ws-0.35_dp*wsmax)/(0.75_dp*wsmax-0.35_dp*wsmax)))
  !===== This effect is not considered, the plants have enough water =====!
  fws = 1.0d0
END SUBROUTINE get_fws


!========================================================!
!
!  rstm: stomatal resistance
!
!========================================================!
SUBROUTINE get_rstm(diff, PARl, rstm_h2o, fvpd, fws, rstm)
  REAL(dp), INTENT(IN) :: diff(ntrac), PARl(nvegl), rstm_h2o(nvegl), fvpd(nvegl), fws
  REAL(dp), INTENT(OUT) :: rstm(nvegl, ntrac)

  DO jk=1,nvegl
    IF (PARl(jk) > 1.0d-3) THEN
      rstm(jk,:)=diff*rstm_h2o(jk)/(fvpd(jk)*fws+Eps)
    ELSE
      rstm(jk,:) = 1.0d5
    END IF
  END DO
END SUBROUTINE get_rstm


!========================================================!
!
! rac: In-canopy turbulent resistance needed to consider soil deposition
!
!========================================================!
SUBROUTINE get_rac(hc, LAIl, ustveg, rac)
  REAL(dp), INTENT(IN) :: hc, LAIl(nvegl), ustveg
  REAL(dp), INTENT(OUT) :: rac

  REAL(dp) :: href

  ! href = 0.25d0*hc
  ! rac = MAX(1._dp, 14._dp*SUM(LAIl)*hc/ustveg) * href/hc
  rac = 0.0d0
END SUBROUTINE get_rac


!========================================================!
!
! rbveg: quasi-laminar boundary layer resistance for canopy layers
!
!========================================================!
SUBROUTINE get_rbveg(diffrb, u_veg, rbveg)
  REAL(dp), INTENT(IN) :: diffrb(ntrac), u_veg(nvegl)
  REAL(dp), INTENT(OUT) :: rbveg(nvegl, ntrac)

  DO jk=1,nvegl
    rbveg(jk,:)=diffrb*180._dp * SQRT( 0.07_dp/MAX(1.d-10, u_veg(jk)) )
  END DO
END SUBROUTINE get_rbveg


!========================================================!
!
! rbs: quasi-laminar boundary layer resistance for soil
!
! Refering to Eq. 13 in Launiainen et al. (2013, AFM), here we only consider the O3, others will be added later.
!
!========================================================!
SUBROUTINE get_rbs(ustveg1, pres1, T1, molar_mass, rbs)
  REAL(dp), INTENT(IN) :: ustveg1, pres1, T1
  REAL(dp), INTENT(IN) :: molar_mass(ntrac)
  REAL(dp), INTENT(INOUT) :: rbs(ntrac)

  REAL(dp), PARAMETER :: Sc_O3 = 1.07d0   ! Lamaud et al., 2002, AE
  ! REAL(dp), PARAMETER :: D_O3 = D_H2O / 1.63d0  ! S and P (P.921)
  REAL(dp), PARAMETER :: z1 = 0.1d0  ! [m], below which the wind profile is assumed logarithmic (roughness length???)

  REAL(dp) :: delta0  ! [m], the height at which molecular and turbulent transport equal
  REAL(dp) :: dyn_vis_air  ! [kg m-1 s-1]
  REAL(dp) :: kin_vis_air  ! [m2 s-1]
  REAL(dp) :: rho_air  ! [kg m-3]
  REAL(dp) :: molec_diff(ntrac)  ! [m2 s-1]
  REAL(dp) :: Sc(ntrac)  ! [m2 s-1]

  DO jt=1, ntrac
    IF (molar_mass(jt) > 1.0d-3) THEN  ! when molar mass is larger than 0
      !!!!! molecular diffusivity for gases
      molec_diff(jt) = D_H2O*SQRT(MOLAR_MASS_H2O/molar_mass(jt))

      !!!!! Schmidt number of gases, kinematic viscosity of air is used since gases are diluted in the air
      !! First method
      ! dyn_vis_air = dynamic_viscosity_air(T1)  ! dynamic viscosity of air
      ! rho_air = pres1/(Rspe_air*T1)  ! air density
      ! kin_vis_air = dyn_vis_air / rho_air  ! kinematic viscosity of air
      ! Sc(jt) = kin_vis_air/molec_diff(jt)

      !! Second method
      ! Sc(jt) = Sc_O3*molar_mass(jt)/molar_mass(ind_O3)

      !! Third method
      Sc(jt) = 1.589d-5/molec_diff(jt)  ! niu_air = 1.5571d-5 [m2 s-1] at 25 [degC]

      delta0 = molec_diff(jt) / (vk*ustveg(1))
      rbs(jt) = (Sc(jt) - LOG(delta0/z1))/(vk*ustveg(1))
    ELSE
      rbs(jt) = r_max  ! default maximum resistance values END IF
    END if
  END DO
  ! rbs(ind_O3) = 0.0d0  ! Test RBS0 case
  ! write(*,*) 'O3, CH3COCH3, CH3OH, SO2: ', Sc(ind_O3), Sc(1920), Sc(300), Sc(2025)
  ! write(*,*) 'rbs(O3): ', rbs(ind_O3), rbs(1920), rbs(300), rbs(2025)
END SUBROUTINE get_rbs


!========================================================!
!
! Calculate frac_ws and frac_veg
!
!========================================================!
SUBROUTINE get_frac_ws(RH, m_wetskin, frac_ws, frac_veg)
  REAL(dp), INTENT(IN) :: RH(nvegl)
  INTEGER, INTENT(IN) :: m_wetskin
  REAL(dp), INTENT(OUT) :: frac_ws(nvegl), frac_veg(nvegl)

  INTEGER :: jk
  REAL(dp), PARAMETER :: rhm=0.70d0, rhh=0.90d0  ! default: 0.55, 0.90

  SELECT CASE (m_wetskin)
  CASE (1)  ! frac_ws = 0.0d0
    frac_ws = 0.0d0
  CASE (2)  ! frac_ws is parameterized
    DO jk=1,nvegl
      ! write(*,*) 'RH(jk)', RH(jk)
      IF ( RH(jk)>=rhm .AND. RH(jk)<rhh ) THEN
        ! write(*,*) RH(jk), (RH(jk)-rhm)/(rhh-rhm)
        ! frac_ws(jk)=1.0d0 - (1.0d0-(RH(jk)-rhm)/(rhh-rhm))
        frac_ws(jk)=(RH(jk)-rhm)/(rhh-rhm)
      ELSE IF (RH(jk)>=rhh) THEN
        frac_ws(jk)=1.0d0
      ELSE
        frac_ws(jk)=0.0d0
      END IF
    END DO
  CASE DEFAULT  ! frac_ws = frac_ws_in, use input wetskin fraction
    CONTINUE
  END SELECT

  frac_veg = 1.0d0-frac_ws
END SUBROUTINE get_frac_ws


!========================================================!
!
!  Calculation/definition of the stomatal compenstion
!  point for use in the dry deposition calculations  
!  ----------
!  Interface:
!  ----------
!  input 
!  rho      : air density in kg m-3
!  frac_veg : vegetation fraction
!  Tc       : surface temperature
!  pxtmveg_CO2 : crown layer CO2 concentration
!  idt_CO2  : tracer index CO2
!  idt_NO2  : tracer index NO2
!  idt_NH3  : tracer index NH3
!
!  output
!  ccomp     : tracer compensation points
!
!========================================================!
SUBROUTINE get_ccomp(frac_veg, rho, Tc, ccomp)
  IMPLICIT NONE

  ! I/O
  REAL(dp), INTENT(IN)  :: frac_veg(nvegl)
  REAL(dp), INTENT(IN)  :: rho(nvegl), Tc(nvegl)

  REAL(dp), INTENT(OUT) :: ccomp(nvegl, ntrac)

  ! local parameters
  REAL(dp) :: rho2nc(nvegl)  ! [# m-3], convert rho to number concentration of air
  INTEGER :: jk

  INTEGER :: idt_CO2 = 0
  INTEGER :: idt_NO2 = 0
  INTEGER :: idt_NH3 = 0

  rho2nc=rho*NA/Mair

  ccomp = 0.0_dp
  IF (idt_CO2 > 0) THEN
    ccomp(:, idt_CO2)=365.0_dp*1d-6 * rho2nc * 1d-6  ! air number concentration to CO2 number concentration, 350 ppmv to [# cm-3] ! alternative 0.99*pxtmveg_CO2(jl)   
    ccomp(:, idt_NO2)=0.5_dp*1d-9 * rho2nc * 1d-6  ! 0.5 ppbv to [# cm-3]
  END IF

  ! LG- compensation point of NH3, 10-2001 it is set to a constant value.
  !     However, there is more information in the paper by Asman et al., 1998
  !     (submitted that time, LGP475), showing that the stomatal compensation
  !     point Xs can be calculated as Xs=f(T)[NH4+]/[H+], with f(T) expressing 
  !     the temperature dependence (leaf temperature) and [NH4+] and [H+] being
  !     the apoplastic concentrations, [NH4+] is very sensitive to the leaf   
  !     N status and the external N supply. Typical concentrations, given in
  !     mM range between 0.04 and 2.3 mM, dependent on the external NH4+/NH3 
  !     concentrations. Apoplastic pH (and thus -log[H+]) ranges between 
  !     5.5 and 6.5. For f(T), maybe that papers by Farquhar show some specific
  !     functions, but also the results in Figure 6 of the Asman paper can be
  !     used to define f(T)

  ! ccomp(idt_NH3)=0.5_dp*1d-9 * rho2nc * 1d-6  ! 0.5 ppbv to [# cm-3]
  ! ESS_lg_20070806+ see excercise NH3 compensation model, Cc model approach 
  ! ccomp(:, idt_NH3)=0.0_dp
  IF (idt_NH3 > 0) THEN
    DO jk=1,nvegl
      IF (frac_veg(jk) > 0.0_dp) THEN ! ESS_lg_20070801+ included the vegfrac
         ccomp(jk, idt_NH3)=0.05*                                       & ! scalings factor, typical concentration should be 1 ug m-3
            ((2749644./Tc(jk))*EXP(-10378./Tc(jk))*1e9*2300.) & ! ESS_lg_20070706+ ug m-3
              *(1e-6*1e-6*NA/17.)                                    ! molecules cm-3     
      ENDIF
    END DO
  END IF
END SUBROUTINE get_ccomp


!========================================================!
!
!  rleaf: leaf resistance
!
!========================================================!
SUBROUTINE get_rleaf(rcut, rstm, rmes, ccomp, rleaf)
  REAL(dp), INTENT(IN) :: rcut(ntrac), rstm(nvegl, ntrac), rmes(ntrac), ccomp(nvegl, ntrac)
  REAL(dp), INTENT(OUT) :: rleaf(nvegl, ntrac)

  DO jk=1,nvegl
    rleaf(jk,:)=1._dp/( 1._dp/rcut + 1._dp/(rstm(jk,:)+rmes) )
    WHERE (ccomp(jk,:)>0.0_dp)
      rleaf(jk,:) = -99.99_dp
    END WHERE
  END DO
  !!!!! comment below if do not consider the one-sided stomata
  ! rleaf(1,:)=1._dp/( 1._dp/rcut + 1._dp/(2.0d0*(rstm(1,:)+rmes)) )
END SUBROUTINE get_rleaf


!========================================================!
!
!  rleafw: leaf resistance for wet skin
!
!========================================================!
SUBROUTINE get_rleafw(rws, rstm, rmes, stomblock, ccomp, rleafw)
  REAL(dp), INTENT(IN) :: rws(ntrac), rstm(nvegl, ntrac), rmes(ntrac), stomblock, ccomp(nvegl, ntrac)
  REAL(dp), INTENT(OUT) :: rleafw(nvegl, ntrac)

  DO jk=1,nvegl
    rleafw(jk,:)=1._dp/( 1._dp/rws + &
                 1._dp*stomblock/rws + &
                 1._dp/(rstm(jk,:)+rmes)*(1._dp-stomblock) )
    WHERE (ccomp(jk,:)>0.0_dp)
      rleafw(jk,:) = -99.99_dp
    END WHERE
  END DO
  !!!!! comment below if do not consider the one-sided stomata
  ! rleafw(1,:)=1._dp/( 1._dp/rws + &
  !             1._dp*stomblock/rws + &
  !             1._dp/(2.0d0*(rstm(1,:)+rmes))*(1._dp-stomblock) )
END SUBROUTINE get_rleafw


!========================================================!
!
!  rsveg:
!  rswet:
!
!========================================================!
SUBROUTINE get_rs(rbveg, rleaf, rleafw, rsoil, rac, LAIl, rsveg, rswet)
  REAL(dp), INTENT(IN) :: rbveg(nvegl, ntrac), rleaf(nvegl, ntrac), rleafw(nvegl, ntrac)
  REAL(dp), INTENT(IN) :: rsoil(ntrac), rac
  REAL(dp), INTENT(IN) :: LAIl(nvegl)

  REAL(dp), INTENT(OUT) :: rsveg(nvegl, ntrac), rswet(nvegl, ntrac)

  ! LG-    calculation of the total uptake resistance of the vegetated and
  !        wet skin fraction from the specific land cover type resistances
  !        and the land cover fractions, this gives the total uptake
  !        velocities for these fractions, reflecting solely the within
  !        canopy uptake processes. the turbulent transport is considered
  !        in the calculations of the canopy top fluxes. the bulk dry 
  !        deposition velocities are scaled with the LAD and radiation 
  !        profile for the two layers in order to scale the "big leaf" 
  !        dry deposition velocity accounting for the biomass distribution 
  !        and the larger removal rates for the sunlit leaves

  ! LG-    for the lowest canopy layer, the limiting turbulent transport from
  !        the reference height to the soil surface is also considered.

  !===== Understorey (jk=1) =====!
  jk=1
  DO jt=1,ntrac
    IF (rleaf(jk,jt) == 0.0d0) CYCLE

    !!!!! rsveg
    IF (rleaf(jk,jt) > 0.0d0) THEN  ! dry vegetation
      rsveg(jk,jt) = 1._dp/(                                                            &
                            1._dp/( (rbveg(jk,jt)+rleaf(jk,jt))/MAX(1.d-5,LAIl(jk)) ) + &
                            1._dp/( rsoil(jt)+rac )                                     &
                          )
    ELSE
      rsveg(jk,jt)=1.d10
    END IF

    !!!!! rswet
    IF (rleafw(jk,jt) > 0._dp) THEN   ! ESS_LG_20120721+
      rswet(jk,jt) = 1._dp/(                                                             &
                            1._dp/( (rbveg(jk,jt)+rleafw(jk,jt))/MAX(1.d-5,LAIl(jk)) ) + &
                            1._dp/( rsoil(jt)+rac )                                      &
                           )
    ELSE
      rswet(jk,jt)=1.d10
    ENDIF
  END DO  ! jt=1,ntrac

  !===== Crown layers (jk=2:nvegl) =====!
  DO jt=1,ntrac
    DO jk=2,nvegl
      !!!!! rsveg
      IF (rleaf(jk,jt) > 0._dp) THEN
        rsveg(jk,jt)=(rbveg(jk,jt)+rleaf(jk,jt)) / MAX(1.d-5, LAIl(jk))
      ELSE
        rsveg(jk,jt)=1.d10
      END IF

      !!!!! rswet
      IF (rleafw(jk,jt) > 0._dp) THEN  ! ESS_lg_20120721+
        rswet(jk,jt)=(rbveg(jk,jt)+rleafw(jk,jt)) / MAX(1.d-5, LAIl(jk))
      ELSE
        rswet(jk,jt)=1.d10
      END IF
    END DO
  END DO
END SUBROUTINE get_rs


!========================================================!
!
! rtot: calculation of surface uptake resistance 
!
!========================================================!
SUBROUTINE get_rtot(frac_veg, frac_ws, rsveg, rswet, rtot)
  REAL(dp), INTENT(IN) :: frac_veg(nvegl), frac_ws(nvegl)
  REAL(dp), INTENT(IN) :: rsveg(nvegl, ntrac), rswet(nvegl, ntrac)

  REAL(dp), INTENT(OUT) :: rtot(nvegl, ntrac)

  DO jt=1,ntrac
    rtot(:,jt)=1._dp/(                        &
                      frac_veg/rsveg(:,jt) +  &
                      frac_ws/rswet(:,jt)     &
                     )
  END DO
END SUBROUTINE get_rtot




!**************************************************************************************************
!
!  Procedures used for debugging
!
!**************************************************************************************************
SUBROUTINE read_data_from_messy_1(file_name, skip_row, col, output_data)
  CHARACTER(LEN=*), INTENT(IN) :: file_name
  INTEGER, INTENT(IN) :: skip_row, col
  REAL(dp), INTENT(OUT) :: output_data(:)

  INTEGER :: nd
  INTEGER :: jl
  REAL(dp), ALLOCATABLE :: dummy_data(:,:)  ! save the data left of col
  CHARACTER(LEN=20) :: dummy_string(2)

  nd = SIZE(output_data)
  ALLOCATE(dummy_data(nd, col))

  OPEN(UNIT=101, FILE=TRIM(file_name), STATUS='OLD', ACTION='READ')
  DO jl=1,nd+skip_row
    IF (jl<=skip_row) THEN
      READ(101, *)
    ELSE
      READ(101, *) dummy_data(jl-skip_row,1), dummy_string, dummy_data(jl-skip_row,4:col)
    END IF
  END DO
  CLOSE(UNIT=101)

  output_data = dummy_data(:,col)
END SUBROUTINE read_data_from_messy_1


SUBROUTINE read_data_from_messy_2(file_name, skip_row, cols, output_data)
  CHARACTER(LEN=*), INTENT(IN) :: file_name
  INTEGER, INTENT(IN) :: skip_row, cols(:)
  REAL(dp), INTENT(OUT) :: output_data(:,:)

  INTEGER :: ncols, nd1, nd2, maxcol
  INTEGER :: jl
  REAL(dp), ALLOCATABLE :: dummy_data(:,:)
  CHARACTER(LEN=20) :: dummy_string(2)

  ncols = SIZE(cols)
  nd1 = SIZE(output_data,1)
  nd2 = SIZE(output_data,2)
  maxcol = MAXVAL(cols)
  ALLOCATE(dummy_data(nd1, maxcol))

  OPEN(UNIT=101, FILE=TRIM(file_name), STATUS='OLD', ACTION='READ')
  DO jl=1,nd1+skip_row
    IF (jl<=skip_row) THEN
      READ(101, *)
    ELSE
      READ(101, *) dummy_data(jl-skip_row,1), dummy_string, dummy_data(jl-skip_row,4:maxcol)
    END IF
  END DO
  CLOSE(UNIT=101)

  output_data = dummy_data(:,cols)
END SUBROUTINE read_data_from_messy_2


SUBROUTINE read_data_from_data_file(file_name, skip_row, cols, output_data)
  CHARACTER(LEN=*), INTENT(IN) :: file_name
  INTEGER, INTENT(IN) :: skip_row, cols(:)
  REAL(dp), INTENT(OUT) :: output_data(:,:)

  INTEGER :: ncols, nd1, nd2, maxcol
  INTEGER :: jl
  REAL(dp), ALLOCATABLE :: dummy_data(:,:)
  CHARACTER(LEN=20) :: dummy_string(2)

  ncols = SIZE(cols)
  nd1 = SIZE(output_data,1)
  nd2 = SIZE(output_data,2)
  maxcol = MAXVAL(cols)
  ALLOCATE(dummy_data(nd1, maxcol))

  OPEN(UNIT=101, FILE=TRIM(file_name), STATUS='OLD', ACTION='READ')
  DO jl=1,nd1+skip_row
    IF (jl<=skip_row) THEN
      READ(101, *)
    ELSE
      READ(101, *) dummy_data(jl-skip_row,:)
    END IF
  END DO
  CLOSE(UNIT=101)

  output_data = dummy_data(:,cols)
END SUBROUTINE read_data_from_data_file


SUBROUTINE list_matrix_real_2(mat)
  REAL(dp), INTENT(IN) :: mat(:,:)

  INTEGER :: I, J

  DO I=1, SIZE(mat,1)
    WRITE(*,*) (mat(I,J), J=1,SIZE(mat,2))
  END DO
END SUBROUTINE list_matrix_real_2


ELEMENTAL FUNCTION dynamic_viscosity_air(T)
  REAL(dp), INTENT(IN) :: T
  REAL(dp) :: dynamic_viscosity_air
  REAL(dp), PARAMETER :: b = 1.458d-6  ! [kg m-1 s-1 K-1/2]
  REAL(dp), PARAMETER :: S = 110.4d0  ! [K]

  dynamic_viscosity_air = b*T**1.5/(T+S)  ! [kg m-1 s-1] or [N s m-1]
END FUNCTION dynamic_viscosity_air


! SUBROUTINE write_matrix_real_2(mat)
!   REAL(dp), INTENT(IN) :: mat(:,:)
!   
!   INTEGER :: I, J
! 
!   OPEN(UNIT=101, FILE=TRIM(file_name), STATUS='OLD', ACTION='READ')
! END SUBROUTINE write_matrix_real_2
END MODULE gdd_function_mod
