!==============================================================================!
!
! DATA_FORMAT
!
! Contains global variables, model parameters and constants
!
!==============================================================================!

module Sosa_data
  use mpi
  use constants_mod
  use tool_mod
  use second_Parameters, kpp_dp=>dp  ! NSPEC and chemical indices, ind_xxxx, come from here
  use Chemistry_Mod, ONLY : NKVALUES, NWL
  use second_Global, ONLY : NPHOT
  use psd_constants, a_pi=>pi,a_dp=>dp
  use traj_data_mod

  IMPLICIT NONE


  !============================================================================!
  ! MODEL PARAMETERS
  !============================================================================!

  ! The DIR should be given with absolute paths without the last '/',
  ! or if your program is in default structure you can use relative paths.

  ! Dirs in the name list
  CHARACTER(len=200) :: WORK_DIR   = 'SOSAA'             ! root SOSAA dir
  CHARACTER(len=200) :: CODE_DIR   = 'sosaa'             ! source code dir
  CHARACTER(len=200) :: CASE_DIR   = 'cases'             ! case directory
  CHARACTER(len=200) :: CHEM_DIR   = 'sample'            ! chemistry scheme name
  CHARACTER(len=200) :: INPUT_DIR  = 'sosa_in'           ! input dir
  CHARACTER(len=200) :: OUTPUT_DIR = 'output'            ! output dir
  CHARACTER(len=200) :: STATION    = 'hyytiala'          ! station name

  ! Other dirs
  CHARACTER(LEN=200) :: input_dir_general     , &  ! general input files
                        input_dir_station     , &  ! input data from station
                        input_dir_station_info, &  ! info files from station
                        input_dir_station_data     ! input data in month folders

  NAMELIST /NML_MAIN/ WORK_DIR, CODE_DIR, CASE_DIR, CHEM_DIR, INPUT_DIR, OUTPUT_DIR, STATION

  !===== NML_FLAG: defines different flags =====!

  ! Emission flag: 1 = old megan code, 2 = new megan code, 3 = SIMBIM, 0 = clear cut
  INTEGER :: flag_emis = 2

  ! Gas dry deposition flag: 1 = gas dry deposition included, 0 = not included
  INTEGER :: flag_gasdrydep = 0

  ! 1 = Chemistry included
  INTEGER :: flag_chem = 0

  ! 1 = aerosol dynamics included
  INTEGER :: flag_aero = 0

  ! 1 = measured soil VOC emission included (do not use this at the moment,
  ! cause it needs to be updated).
  ! 2 = include soil VOC parameterisation from Hermanni
  INTEGER :: flag_emis_soil = 0

  ! 1 = pine, 2 = spruce, 3 = birch.
  ! BIOMASS: remember to change this for different tree ages! 4 = clear cut
  INTEGER :: flag_tree = 1

  ! Read data from ECMWF datasets for boundary conditions, need to improve the
  ! flag in future ??
  INTEGER :: flag_ecmwf = 1

  ! 1 = assume monoterpene emission is purely by de novo.
  ! Otherwise it is assumed it is purely by pool
  INTEGER :: flag_syn = 0

  ! 1 = concentration in PPB, 2 = concentration in molec cm-3
  INTEGER :: flag_NH3 = 1

  ! 1: output debug message, 0: no output
  INTEGER :: flag_debug = 1

  INTEGER :: flag_vapor = 0

  INTEGER :: flag_outlist = 0

  INTEGER :: flag_format = 0

  integer, parameter :: STATION_MODE = 1
  integer, parameter :: TRAJECTORY_MODE = 2
  integer :: flag_model_type = STATION_MODE

  integer :: flag_mix_chem = 1
  integer :: flag_mix_aero = 1

  ! if calculating aerosol part in parallel
  logical :: use_parallel_aerosol = .true.
  logical :: FLAG_SKIP_EXISTING = .false.

  NAMELIST /NML_FLAG/ &
    flag_emis, flag_chem, flag_gasdrydep, flag_aero, flag_emis_soil, &
    flag_tree, flag_syn, flag_NH3, flag_ecmwf, flag_debug, &
    flag_vapor, flag_outlist, flag_format, flag_model_type, flag_mix_chem, &
    flag_mix_aero, use_parallel_aerosol, weigh_with_srr,FLAG_SKIP_EXISTING


  !===== AER_FLAG: aerosol options =====!
  NAMELIST /AER_FLAG/ OPTIONS, ACDC_SYSTEMS, ACDC_links, Vap_names


  INTEGER, PARAMETER ::  kz = 50

  !===== NML_GRID: location information =====!

  ! When using Height_flag = 15, kz=51 should be used   !100 layers for output
  ! meters above sea level for the simulation site, Hyytiala: 181.0, Manitou: 2286.0
  REAL(dp) ::  masl = 181.0_dp

  ! latitude and longitude
  ! station: can be set from the namelist, use hyytiala location as default values
  !          (https://www.ingos-infrastructure.eu/smear-ii/), usually constant
  !          in the simulations
  ! trajectory: will be read from the input file, usually changing with time
  real(dp) :: lat_deg = 61.85d0, lon_deg = 24.28d0
  real(dp) :: lat_rad, lon_rad
  NAMELIST /NML_GRID/ masl, lat_deg, lon_deg

  !! Time
  INTEGER, DIMENSION(6) :: start_date,aero_start_date, end_date  ! start date and end date for simulation, in the form of 'yyyy mm dd HH MM SS'
  INTEGER, DIMENSION(6) :: now_date  ! start date and end date for simulation, in the form of 'yyyy mm dd HH MM SS'
  !! Timesteps for different module
  REAL(dp) :: dt_obs = 1800.0d0 ! temporal resolution of surface observations (= input data) [s], added by Rosa
  REAL(dp) :: dt_mete = 10.    ! Time step for Scadis
  REAL(dp) :: dt_emis = 60.    ! Time step for emission
  REAL(dp) :: dt_chem = 60.    ! Time step for Kpp
  REAL(dp) :: dt_depo = 60.    ! Time step for Kpp
  REAL(dp) :: dt_aero = 10.    ! aerosol simulation time
  REAL(dp) :: dt_uhma = 10.    ! Not in use anymore, legacy

  ! UTC time zone, e.g., +2 for Hyytiala in winter. For traj mode, the UTC time
  ! should be used for the simulation, thus the time_zone should be set to 0.
  ! When the local time is needed, the UTC time can be converted according to
  ! the location, now this has not been done since this also relates to the
  ! country regulations. The unit is [hour].
  real(dp) :: time_zone = 2.0d0

  NAMELIST /NML_TIME/ dt_obs, start_date, aero_start_date, end_date, &
    dt_mete, dt_emis, dt_chem, dt_depo, dt_aero, dt_uhma, dt_aero, time_zone


  !============================================================================!
  ! Time variables
  !============================================================================!

  CHARACTER(len=4) :: year_str   ! 1999, 2000, ...
  CHARACTER(len=2) :: year_str2  ! 99, 00, 01, ..., used for dmps data file names
  character(len=2) :: month_str  ! 01, 02, ..., 12
  character(len=2) :: day_str    ! 01, 02, ..., 31
  CHARACTER(len=3) :: month_str_Mxx  ! M01, M02, ..., M12, used for input folder names

  real(dp) :: time           ! time in model, starting from 0
  real(dp) :: time_in_month  ! time in month
  real(dp) :: time_in_day    ! time in day
  INTEGER :: ttime, mon, julian, time_par, &
             date, second, hour, minute,   &  ! for printing time in output files
             nclo_1, nmin_1,               &  ! Starting clock
             nmin_end, nclo_end,           &  ! Ending clock
             nstep                            ! Number of steps
  real(dp) :: sunrise_UTC, sunset_UTC, time_in_day_UTC, julian_UTC

  ! Used for system_clock
  ! call system_clock(count, count_rate, count_maximum)
  ! c: count
  ! cr: count_rate, usually constant
  ! cm: maximum count, usually constant
  ! elapsed time can be calculated as (sc_c1 - sc_c0)/sc_cr
  integer :: sc_c0, sc_c1, sc_cr, sc_cm

  type timer
    integer :: c0, c1, cr, cm
    real(dp) :: used  ! accumulated used time
    real(dp) :: last  ! used time of the last round
  end type timer

  type(timer) :: timer_total_aero
  type(timer) :: root_recv_aero(kz)
  type(timer) :: root_send_aero(kz)


  REAL(dp) :: aero_start
  REAL(dp) :: day, time_chemie
  LOGICAL :: is_newday, is_newmonth
  LOGICAL :: first_time_loop

  !!Grid parameters
  REAL(dp),DIMENSION(kz) :: z, dz

  !!Surface
  REAL(dp) :: z0                     ! Roughness at a the surface z0 (0.0001-1.0 m)
  REAL(dp) :: clob, clom, clot,   &  ! Cloudiness of bottom, middle, and top layer accordingly
              mclob, mclom, mclot    ! Cover fraction of low, middle, and high clouds (0-10)

  REAL(dp) :: rou, & ! Snow cover depth (0. - 2.0 m). 0 means no consideration of snow.
              psnow, abb, sb, szb

  !!Vegetation and canopy
  INTEGER :: su        ! flag for leaf type. 1 for conifirous trees, 2 for broadleaf trees
  REAL(dp) :: hc       ! height of vegetation (0 - 35m)
  REAL(dp) :: pa, pa2  ! Shape of vegetation (1.0-20.0). pa for the meteorology, pa2 only for MEGAN or general emission
  ! activity
  ! pa is alfa-parameter in beta-probability functon, the beta-parameter in beta-probobility function
  ! is fixed as 3.
  !
  !   pa=1 *       pa=3 *      pa=7 *    '
  !        **           **          **** '
  !        ***          ***         ****** '
  !        ****         ****        ***** '
  !        *****        ****        **** '
  !        ******       ***         *** '
  !        ****         **          ** '
  !        *            *           * '
  !
  ! For Hyytiala it is better to take 3.

  REAL(dp)                          :: EM_LAI_year              ! LAI factor for different years
  REAL(dp), DIMENSION(12)            :: EM_LAI_year_Hyy          ! Array of data for different LAI per year (1996-2010)
  REAL(dp) :: LAI_proj, LAI_curv, LAI_tot  ! projected, curved and total LAI

  REAL(dp), DIMENSION(kz) :: s,s1,s2, sl, LAD_P      ! canopy parameters

  !!Met variables
  REAL(dp), DIMENSION(kz) ::               &
       inv_conc_array,                     &  ! temporary array (still deving this- P.C.)
       pres,    &    ! Pressure [Pa]
       air,     &    ! air molecules density ?? [#/cm³] or [#/m³]???
       ta1,tsn1,tsd1,tsoil1, ta1_C,       & ! temperature variables
                                !  EW, ES                            &
       Rh      ! Relative humidity in %

  REAL(dp), DIMENSION(kz, 35) :: Other_out  ! used for output some data, including:
                                                 ! Qabs_sn, Qabs_sd, Labs, Lout_sn, Lout_sd,
                                                 ! SH_sn, SH_sd, LH_sn, LH_sd, Su_out,
                                                 ! Sd_out

  CHARACTER(LEN=*),PARAMETER ::       &
       FORMAT1 = "(6I6, 100E20.8)",    &    ! for aerosol variables (sections), e.g. N_CONC, Rk_Vel
       FORMAT2 = "(6I6, 5E20.8)",     &    ! for aerosol variables (components)  e.g. Sink

       FORMAT3 = "(6I6, 20E20.8)",    &    ! for aerosol variables, one number one layer, e.g. NUC_RATE
       FORMAT4 = "(6I6, 51E20.8)",    &
       FORMAT5 = "(6I6, 110E20.8)",   &

       FORMAT6 = "(6I6, 15E20.8)",    &    ! for gases
       FORMAT7 = "(6I6, 51E20.8)",    &
       FORMAT8 = "(6I6, 110E20.8)"

  !!Input

  !============================================================================!
  ! Input raw data from station input data files
  ! There variables are public ones.
  !============================================================================!

  ! SMEAR II dataset, 1489 = 31*48 + 1, from 00:00 at the first day to 24:00
  ! at the last day every half hour, 31 is the maximum day number
  INTEGER, PARAMETER :: NPMON=1489

  ! ECMWF dataset, 249 = 31*8 + 1, from 00:00 at the first day to 24:00 at the
  ! last day, 31 is the maximum day number
  INTEGER, PARAMETER :: NPECMWF=249

  REAL(dp), DIMENSION(NPMON) :: hyy_ch4
  REAL(dp), DIMENSION(45,90)     :: E_o_E
  REAL(dp), DIMENSION(47,NPMON)     :: swr_hyy
  REAL(dp), DIMENSION(32,NPMON)     :: hyy_mix
  REAL(dp), DIMENSION(10,NPMON)      :: hyy_soil
  real(dp), dimension(2,NPMON)     :: EF_meas
  real(dp), dimension(16,9)     :: canopyfile

  REAL(dp), SAVE                ::  SAMI(366,3)
  REAL(dp), SAVE                ::  LAI_function(366) !change in LAI for birch ([0 1])
  REAL(dp) :: SEP_clear, sep_factor

  ! Time series of input scalar variables at the station
  REAL(dp), DIMENSION(NPMON)        :: stat_cs_H2SO4, stat_cs_HNO3
  real(dp), dimension(NPMON)        :: stat_PAR, stat_glob
  real(dp), dimension(NPMON)        :: stat_sm_1, stat_sm_2, stat_sm_3
  real(dp), dimension(NPMON)        :: stat_albedo
  REAL(dp), DIMENSION(NPMON)        :: stat_gsoil

  REAL(dp), DIMENSION(NPECMWF, 43)     :: input_O3_mixing_ratio, input_temp, input_vwind, &
                                           input_uwind, input_wind, input_RH, input_geop, input_geom, input_rhov, input_qv
  REAL(dp), DIMENSION(NPECMWF, 7)      :: input_strd
  ! REAL(dp), DIMENSION(37)         :: Input_Pressure_levels
  INTEGER :: h1, h2
  CHARACTER(LEN=2) :: next_month

  REAL(dp), DIMENSION(NPECMWF,kz) :: ECMWF_O3_mixing_ratio, &
                                      ECMWF_geop, &
                                      ECMWF_temp, &
                                      ECMWF_rhov, &
                                      ECMWF_qv, &
                                      ECMWF_uwind, ECMWF_vwind, ECMWF_wind, &
                                      ECMWF_RH  !Luxi
  REAL(dp), DIMENSION(NPECMWF) :: tempgrad, ECMWF_dtemp, ECMWF_strd
  REAL(dp), DIMENSION(NPMON,6)  :: Gas_in
  REAL(dp), DIMENSION(NPMON,35) :: Mix_in


  INTEGER, DIMENSION(NPECMWF,6) :: ECMWF_TIME

  ! added by Rosa for using ECMWF data for long wave radation
  REAL(dp), DIMENSION(NPECMWF) :: LWR_in    ! downward long wave radiation at the surface (W/m2)
  REAL(dp)                 :: LWRdown   ! downward long wave radiation at the surface, interpolated to model time (W/m2)
  !!Output

  CHARACTER(LEN=20) :: AR_OUTPUT_FORMAT, CH_OUTPUT_FORMAT


  !!COUNTING VARIABLES
  INTEGER :: I, J, m, n, II, JJ


  !!Miscellaneous
  INTEGER :: layer_canopy = 19     ! for MAN when kz = 51
  REAL(dp), DIMENSION(1809,8)       :: sounding


  !!Du know what they are. L.
  REAL(dp), DIMENSION(kz) ::  &
       da,db,dc,df,                &
       terma,termb,termc,termd,    &   ! subroutine coefficients
       b1,b2,d1,d2,bg,dg,          &
       e, f,ffn,                   &

       alt,                        &
       u,v,w,                      &
       bt,dbt,l,kt,                &
       ta,tsn,tsd,tsoil,           &
       tsn_megan,tsd_megan,        &
       tsn_megan_C,tsd_megan_C,    &
       qa,qsn,qsd,                 &

       unud,vnud,bnud,tnud,qnud,   &

       alt1,                       &
       u1,v1,w1,                   &
       bt1,dbt1,l1,kt1,            &

       qa1,qsn1,qsd1,              &  ! [kg m-3], absolute humidity

       gl,gd,                                        &
       ur,                                           &
       ksnow,                                        &
       psn,psk,                                      &
       fphu,fphd,fniu,fnid,                          &
       iru,ird,firu,fird,                            &
       nusn,nusd,grsn,grsd,re,nure,                  &
       shsd,shsn,shre,                               &
       dhsn,dhsd,dvsd,dvsn,                          &
       rosn,rasn,rasd,rosd,                          &
       fluxle,fluxle3,fluxh,fluxh3,                  &
       asun,asky,airu,aird,aphu,aphd,aniu,anid,      &
       rih,rih1,f2,bt2,det2,prod1,prod2, vtran,      &

       EaP, EaT, EaS, gas_ver, gas_ver_new,   &
       SUN_PAR, EMI_VER, VOC_sink, col_rate,         &
       nuc, nuc_cs, nuc_org

  INTEGER :: var_1, var_2, var_3

  REAL(dp) :: iru0,iruh,nua,kta,kva,lv,rads2,UVA,BETA,cote,dtet,teto,lai0,akt_dvz,akt_duz,pa0,pa1, &
                   btbot, &
                   deep_soil_hum, surf_soil_hum, deep_soil_temp   ! added by rosa: value interpolated from input data to current time step
  REAL(dp), DIMENSION(kz) :: fu0,fu1,fd0,fd1,lai
  REAL(dp), DIMENSION(30) :: check_rh

  ! Horizontal wind velocity
  real(dp) :: hwind(kz)

  !!Added by Sampo, to get rid of implicit variable typing:

  REAL(dp):: abl_max

  INTEGER :: abl
  REAL(dp) :: dirtop1(NPMON),difftop1(NPMON),tem(6,NPMON),hum(6,NPMON), uwind(5, NPMON), vwind(5, NPMON)

  ! Local input data and their corresponding levels
  REAL(dp) :: local_uwind(NPMON, 5)    , loclv_uwind(5)      ! [m s-1], u component
  REAL(dp) :: local_vwind(NPMON, 5)    , loclv_vwind(5)      ! [m s-1], v component
  REAL(dp) :: local_temp(NPMON, 7)     , loclv_temp(7)       ! [degC], air temperature
  REAL(dp) :: local_rhov(NPMON, 7)     , loclv_rhov(7)       ! [kg m-3], absolute humidity
  REAL(dp) :: local_pres(NPMON, 1)     , loclv_pres(1)       ! [Pa], air pressure
  REAL(dp) :: local_soilmoist(NPMON, 5), loclv_soilmoist(5)  ! [m3 m-3], soil volumetric water content
  REAL(dp) :: local_soiltemp(NPMON, 5) , loclv_soiltemp(5)   ! [degC], soil temperature
  REAL(dp) :: local_gsoil(NPMON, 1)    , loclv_gsoil(1)      ! [W m-2], soil heat flux
  REAL(dp) :: local_glob(NPMON, 1)     , loclv_glob(1)       ! [W m-2], global radiation
  REAL(dp) :: local_PAR(NPMON, 1)      , loclv_PAR(1)        ! [umol m-2 s-1], PAR
  REAL(dp) :: local_albedo(NPMON, 1)   , loclv_albedo(1)     ! [-], albedo for global radiation at 125 m
  REAL(dp) :: local_cs_H2SO4(NPMON, 1)   , loclv_cs_H2SO4(1)     ! [s-1],condensation sink without HNO3 correction
  REAL(dp) :: local_cs_HNO3(NPMON, 1)   , loclv_cs_HNO3(1)     ! [s-1],condensation sink with HNO3 correction
  REAL(dp) :: local_coa(NPMON, 1)   , loclv_coa(1)     ! [s-1],condensation sink without HNO3 correction

  REAL(dp) :: local_deep_soiltemp(NPMON)  ! [K], deep soil temperature from measurements

  REAL(dp) :: speed(5,NPMON),border(5,NPMON),border_abl(6,NPECMWF), border_dummy(6,NPECMWF)

  REAL(dp) :: time_end,dis,dtree,dtree2,alai,ztop,zdown,tau3,wind,f1,eps,profile
  REAL(dp) :: a0005,a005, abb1, abb3, ddd1, cc2,c833,c52,alf,face,alph,roa,rh2o,cpa,ch2o,hleaf,alni,dels,delf,asl
  REAL(dp) :: cd,gamma,proatm,wsun,zeit,cos_zenith,optmas,rads,albedo_f

  real(dp) :: zenith_deg
  real(dp) :: declination

  REAL(dp) :: tata,upsn,downsn,ct,teil,rsnt,rskt,tah,qah,fsoil1,pp,temtran,temdif,windx,windy
  REAL(dp) :: udin,sk00,sks00,al1,al100, wg1(2),wg1i,wg2i,wg(2),ts0,smol,fsoil,transp,sumrad
  REAL(dp) :: clo,tran,trancl,emm2,dtmet,shift,tax,temgrad,qax,emm,rhelp,radk,relrad,apar,rsktm
  REAL(dp) :: rsntm,u74,u33,u16,u8,bt23,ta67,ta50,ta33,ta16,ta8,ta4,qa67,qa50,qa33,qa16,qa8,qa4
  REAL(dp) :: cor1,c_nud,c_nud1,zini(6),fini(6),fktt,fktd,dkz1,daz,dkz,dtz1,dqz1,duz1,dvz1,dwz1
  REAL(dp) :: utt,bttop,cc22,c15,c16,cc33,beta0,gf,gf2,sand,snow,dksz,den,site
  REAL(dp) :: tau,wee,tphu,tphd,rphu,rphd,tniu,tnid,rniu,rnid,phsn,phsk,two,tu,td,ru,rd,al,pasn
  REAL(dp) :: pask,trem1,trem2,trem3,trem4,trek1,trek2,trek3,trek4,trom1,trom2,trom3,trom4
  REAL(dp) :: tram1,tram2,tram3,tram4,pgr,pre,rlf,dmax,dmin,topt,tmin,tmax,par50,humg,tc,rhg
  REAL(dp) :: defsoil,wgwilt,wgf,qmg,tmg,gw,g0,tasoil,qasoil,temlim,dk00,dkq,dkt,qmgs,tmgs
  REAL(dp) :: fluxles,fluxhs,haag,haags,ff1,tar,fluxc,tau2,tau4,wgmax,bsoil,cgsat,psoil,asoil
  REAL(dp) :: cw1sat,wps,cw2ref,wgeq,cw1,cw2,cgw,vla,balans,balans1
  REAL(dp) :: Ktt(kz), Ktb(kz)  ! eddy diffusivity for momentum at k+1/2 (Ktt) and k-1/2 (Ktb)
  REAL(dp) :: Kht(kz), Khb(kz)  ! eddy diffusivity for heat and scalar at k+1/2 (Kht) and k-1/2 (Khb)

  ! Planetary boundary layer height and its index
  real(dp) :: pblh
  integer  :: pblh_ind

  REAL(dp), DIMENSION(5) :: wind_mast, q_mast, tem_mast

  INTEGER :: i20,j20,k,kk,nz,kkk,nsmear,kkk5,n233,seq
  INTEGER :: kw,nturb,nxod,nxodrad,nschet,nmetka,krelf,nxodmet,kmix,m2
  INTEGER :: nmix,kf,ks,jk,ja,nnna,k_canopy,nm1, nm

  ! Current albedo, [-]
  real(dp) :: albedo

  ! Incoming photosynthetic active radiation at canopy top, [umol m-2 s-1]
  real(dp) :: PAR

  ! Incoming short wave solar radiation in [W m-2]
  real(dp) :: glob

  ! land sea mask [0-1], 1 means land
  real(dp) :: land_sea_mask

  ! Soil moisture [m3 m-3] in different layers, currently:
  ! 1: organic layer 0 => 5 cm
  ! 2: average of layers -2 => -6 cm and -14 => -25 cm
  ! 3: average of layers -26 => -36 cm and -38 => -61 cm
  ! Direction is upward.
  real(dp) :: soil_moist_1, soil_moist_2, soil_moist_3



  !============================================================================!
  ! Now comes all the parameters that are related to the
  ! chemistry
  !
  ! OBS: local chemistry related MPI parameters and MPI
  ! indexes are NOT denoted with 'CH_'
  !============================================================================!

  ! real(dp) :: Nair(kz)  ! [molec cm-3], air number concentration

  ! Number concentrations of H2O, N2 and O2
  real(dp) :: H2O(kz), N2(kz), O2(kz)

  REAL(dp) :: CH_CONS(NSPEC)                  ! concentration of chemical species in gas phase at one layer
  REAL(dp) :: CH_CONS_all(kz,NSPEC)           ! concentration of chemical species in gas phase at all layers
  REAL(dp) :: CH_CONS_FLUX(kz,NSPEC)          ! flux of chemical species in gas phase
  REAL(dp) :: CH_CONS_VER(kz)                 ! concentration of one chemical species in vertical layers
  REAL(dp) :: CH_CONS_VER_NEW(kz)
  REAL(dp) :: CH_J_values_ALL(kz, NPHOT)       ! J-values   - for MPI
  REAL(dp) :: CH_K_values_ALL(kz, NKVALUES)    ! K-values   - for MPI
  REAL(dp) :: CH_ACF_ALL(kz, NWL)
  REAL(dp) :: CH_RO2(kz)                      ! Peroxyradical concentration
  ! REAL(dp) :: CH_RES1                         ! Condensation sink of H2SO4
  ! REAL(dp) :: CH_RES2                         ! Condensation sink of nitric acid
  REAL(dp) :: CH_TIME_kpp                     ! Time in KPP
  REAL(dp) :: CH_END_kpp                      ! Total time in KPP
  REAL(dp),  DIMENSION(NPHOT)    :: CH_J_values    ! J-value change
  REAL(dp),  DIMENSION(NKVALUES) :: CH_K_values
  REAL(dp),  DIMENSION(NWL)    :: CH_ACF
  REAL(dp),  DIMENSION(12, NPMON)  :: CH_gas_hyy     !Input file with inorganic gas concentrations
  REAL(dp),  DIMENSION(NPMON)  :: NH3_hyy
  REAL(dp),  DIMENSION(kz)     :: NH3_hyy_all
  REAL(dp),  DIMENSION(12, NPMON)  :: CH_VOC_hyy     !Input file with organic gas concentrations
  REAL(dp),  DIMENSION(NPMON)  :: CH_H2SO4_hyy
  REAL(dp),  DIMENSION(kz) :: CH_cs_H2SO4 = 0d0
  REAL(dp),  DIMENSION(kz) :: CH_cs_HNO3 = 0d0

  ! REAL(dp),  DIMENSION(n_cond_tot)::CH_RES_org!,sat_vap,parameter_a,parameter_b
  ! INTEGER,   DIMENSION(n_cond_tot)::index_vapor
  ! character(len=60), dimension(n_cond_tot) ::vapor_names

  ! Reacitivity time step
  INTEGER, PARAMETER :: CH_step2 = 5, CH_step3 = 30
  ! The following is all for the OH-reactivity calculation
  INTEGER :: CH_oh_i,CH_oh_count
  INTEGER, DIMENSION(:), ALLOCATABLE :: CH_oh_indices

  LOGICAL, DIMENSION(NSPEC) :: CH_oh_flag
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: CH_oh_cons3
  ! pseudoconcentrations of OH-reactivities averaged for last 1, 5 (CH_oh_step2), and 30 (CH_oh_step3) minutes, if the
  ! chemistry time step dt_chem is 60 second. If not, the for every 1, 5 (CH_oh_step2) and 30 (CH_oh_step3) dt_chem
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: CH_oh_prev3

  ! The following is all for the O3-reactivity calculation
  INTEGER :: CH_o3_i,CH_o3_count
  INTEGER, DIMENSION(:), ALLOCATABLE :: CH_o3_indices

  LOGICAL, DIMENSION(NSPEC) :: CH_o3_flag
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: CH_o3_cons3
  ! pseudoconcentrations of O3-reactivities averaged for last 1, 5 (CH_o3_step2), and 30 (CH_o3_step3) minutes, if the
  ! chemistry time step dt_chem is 60 second. If not, the for every 1, 5 (CH_o3_step2) and 30 (CH_o3_step3) dt_chem
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: CH_o3_prev3

  ! The following is all for the NO3-reactivity calculation
  INTEGER :: CH_no3_i,CH_no3_count
  INTEGER, DIMENSION(:), ALLOCATABLE :: CH_no3_indices

  LOGICAL, DIMENSION(NSPEC) :: CH_no3_flag
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: CH_no3_cons3
  ! pseudoconcentrations of NO3-reactivities averaged for last 1, 5 (CH_no3_step2), and 30 (CH_no3_step3) minutes, if the
  ! chemistry time step dt_chem is 60 second. If not, the for every 1, 5 (CH_no3_step2) and 30 (CH_no3_step3) dt_chem
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: CH_no3_prev3

  ! The following is all for the reactivity calculation
  INTEGER :: CH_reactivity_i,CH_reactivity_count
  INTEGER, DIMENSION(:), ALLOCATABLE :: CH_reactivity_indices

  LOGICAL, DIMENSION(NSPEC) :: CH_reactivity_flag
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: CH_reactivity_cons3
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: CH_reactivity_prev3


  !============================================================================!
  ! Gas dry deposition variables
  !============================================================================!

  LOGICAL :: l_drydep, l_vpd, l_wetskin
  ! REAL(dp), DIMENSION(:), ALLOCATABLE :: frac_veg, frac_ws
  REAL(dp) :: frac_veg(kz), frac_ws(kz)
  REAL(dp) :: taub(kz), taud(kz)
  REAL(dp) :: LAIl(kz), LAIl_c(kz), LAIl_debug(kz), LAIl_sh(kz), LAIl_sl(kz)
  REAL(dp) :: LAD(kz), LAD_meteo(kz), LAD_megan(kz)
  REAL(dp) :: PARl(kz), PARl_sh(kz), PARl_sl(kz)
  REAL(dp) :: u_veg(kz)
  REAL(dp) :: rho_veg(kz), Tc_veg(kz)
  REAL(dp) :: ws, wsmax
  REAL(dp) :: stomblock
  CHARACTER(LEN=50) :: trname(NSPEC)  ! trace gas names, used for test
  REAL(dp) :: HenrySE(NSPEC), HenryEE(NSPEC), HenryEG(NSPEC), HenryEB(NSPEC), HenryA(NSPEC), HenryB(NSPEC), Henry(NSPEC, kz)
  REAL(dp) :: dryreac(NSPEC)
  REAL(dp) :: molar_mass(NSPEC)
  REAL(dp) :: vdep(kz, NSPEC)  ! [m s-1],

  REAL(dp) :: gasdep_flux(kz, NSPEC)  ! [# m-2 s-1], gas dry deposition flux
  REAL(dp) :: gasdep_flux_acum(kz, NSPEC)  ! accumulated gas dry deposition flux within one output interval

  REAL(dp) :: o3_weight(kz), o3_new
  REAL(dp) :: gstm_h2o_sn(kz)
  REAL(dp) :: gstm_h2o_sd(kz)

  REAL(dp) :: qnud_dep(kz), RH_dep(kz), Tc_dep(kz)

  INTEGER :: dummy_io, dummy_index
  CHARACTER(LEN=100) :: dummy_flag, dummy_MCM, dummy_SMILES, dummy_INCHIkey, dummy_f0Flag, dummy_HFlag, dummy_formula
  REAL(dp) :: dummy_HenrySE, dummy_HenryEE, dummy_HenryEG, dummy_HenryEB, dummy_HenryA, dummy_HenryB
  REAL(dp) :: dummy_heffA, dummy_heffB, dummy_f0, dummy_molar_mass

  REAL(dp) :: dep_output(kz, NSPEC, 20)

  !***** Special output *****!
  REAL(dp) :: Qconc(kz, NSPEC), Qemis(kz, NSPEC), Qchem(kz, NSPEC), Qturb(kz, NSPEC), Qdepo(kz, NSPEC)
  REAL(dp) :: CH_CONS_temp(kz, NSPEC)
  REAL(dp) :: CH_CONS_0(kz, NSPEC), CH_CONS_1(kz, NSPEC), Qfluxdep(kz, NSPEC)
  REAL(dp) :: Qturb_now(kz, NSPEC)


  !============================================================================!
  ! Aerosol variables
  !============================================================================!

  CHARACTER(LEN=6)    :: layerstamp
  CHARACTER(LEN=1500) :: dummy_line

  INTEGER, PARAMETER  :: maxcol=70
  INTEGER             :: cp, csec, time_par_aer, ccount, abl_layer, line

  REAL(dp)       :: Time_aer, End_aer, time_to_write, cs_check

  REAL(dp), DIMENSION(kz)       :: first_run, NUC_RATE, ION_NUC_RATE
  REAL(dp), DIMENSION(maxcol)   :: val  ! used to save values in dmps data files
  character(200) :: dmps_file_name  ! currently used dmps file name

  REAL(dp), DIMENSION(kz, n_bins_par) ::    &
       N_CONC, RADIUS, RDRY,  CORE, MASS, GR, Rk_vel, PAR_FLUX, Brownian, Settling_veloc

REAL(dp), DIMENSION(n_bins_par) :: sea_salt_psd

  !============================================================================!
  ! Emission variables
  !============================================================================!

  ! Emission rate, [molec cm-3 s-1]
  real(dp) :: emis(kz, NSPEC)
  ! Emission rate, [kg m^-1 s^-1]
  real(dp) :: emi_aer(kz, 9)

  INTEGER, PARAMETER :: MEGAN_NSPEC=22
  CHARACTER(LEN=10), PARAMETER :: MEGAN_SPC_NAMES(MEGAN_NSPEC) = (/ &
     'C5H8      ', &  ! 1. C5H8 (Isoprene)
     'Myrcene   ', &  ! 2. Myrcene
     'Sabinene  ', &  ! 3. Sabinene
     'LIMONENE  ', &  ! 4. LIMONENE
     'Carene    ', &  ! 5. Carene
     'Ocimene   ', &  ! 6. Ocimene
     'BPINENE   ', &  ! 7. Bpinene
     'APINENE   ', &  ! 8. Apinene
     'OMT       ', &  ! 9. Other monoterpenes
     'Farnesene ', &  ! 10. Farnesene
     'BCARY     ', &  ! 11. BCARY (Beta-Carophyllene)
     'OSQ       ', &  ! 12. Other sesquiterpenes
     'MBO       ', &  ! 13. MBO (2methyl-3buten-2ol)
     'CH3OH     ', &  ! 14. CH3OH (Methanol)
     'CH3COCH3  ', &  ! 15. CH3COCH3 (Aceton)
     'CH4       ', &  ! 16. CH4 (Methane)
     'NO        ', &  ! 17. NO
     'CH3CHO    ', &  ! 18. CH3CHO (Acetaldehyde)
     'HCHO      ', &  ! 19. HCHO (Formaldehyde)
     'CO        ', &  ! 20. CO
     'Cineole   ', &  ! 21. Cineole (Eucalyptol) (not included in the new megan - used is the same values as Sabinene)
     'Linalool  '  &  ! 22. Linalool
     /)

  ! All variables used in the emissions starting with the three characters 'EM_' to make life easier for implementation
!!$

  REAL                          :: EMI_MEAS_MONO, CHAM_temp

  INTEGER                       :: Can_Lay                           ! Number of vertical layers inside the canopy
  INTEGER                       :: EM_Year                           ! Year in emission module
  INTEGER                       :: EM_run                            ! Run parameter in emission module
  INTEGER                       :: EM_TIME                           ! Time in seconds
  INTEGER                       :: EM_TIME_M2                        ! Time for new Megan in hour minutes second format
  INTEGER                       :: EM_DATE                           ! Date for new Megan in format YYYDDD scalar
  INTEGER                       :: EM_Julian                         ! Julian day of the year
  INTEGER                       :: EM_kz                             ! Amount of layers in the model run (not only canopy)
  INTEGER                       :: EM_Can_Lay                        ! Amount of layers inside the canopy

  REAL(dp), PARAMETER               :: EM_DT = 60               ! Time step for the emission code

  REAL(dp_sb), DIMENSION(kz,14)     :: EM_BIM_S, EM_BIM_C       ! Input (:,1:13) and output (:,14) values for SIMBIM
  REAL(dp_sb)                       :: EM_Time_in, EM_Time_out  ! Start and end time in minutes for SIMBIM
  REAL(dp_sb)                       :: EM_PAR_SB                ! Photosynthetically active radiation for SIMBIM

  REAL(dp)                          :: EM_TIME_M2_R             ! Time in seconds
  REAL(dp)                          :: EM_Lat, EM_Long          ! Latitude and longitude
  REAL(dp)                          :: EM_PAR                   ! Photosynthetically active radiation
  REAL(dp)                          :: EM_PAR_day               ! Daily average photosynthetically active radiation
  REAL(dp)                          :: EM_LAI                   ! Leaf area index
  REAL(dp)                          :: EM_LAI_past              ! LAI from the last month
 ! REAL(dp)                          :: EM_LAI_old
 ! REAL(dp)                          :: EM_LAI_new
 ! REAL(dp)                          :: EM_LAI_past_old
 ! REAL(dp)                          :: EM_LAI_past_new
  REAL(dp)                          :: EM_SRAD                  ! Incoming short wave solar radiation in (W/m²)
  REAL(dp)                          :: EM_SRAD_day              ! Daily average short wave radiation (W/m2)
  REAL(dp)                          :: EM_Beta                  ! Sin of solar angle above horizon
  REAL(dp)                          :: EM_TempK_day             ! Daily average of temperature (K)
  REAL(dp)                          :: EM_SMOIST                ! Soil moisture (%)
 ! REAL(dp)                          :: EM_LAI_year              ! LAI factor for different years

  REAL(dp), DIMENSION(kz)           :: EM_VPD_S, EM_VPD_C       ! Water vapour pressure deficit for sun and shade
  ! leaf temperatures

  REAL(dp), DIMENSION(kz)           :: EM_z                     ! Array of height for each layer (m)
  REAL(dp), DIMENSION(kz)           :: EM_LAD                   ! Array of vertical distribution for the LAI
  REAL(dp), DIMENSION(kz)           :: EM_TempK                 ! Array of temperature (K)
  REAL(dp), DIMENSION(kz)           :: EM_TempC                 ! Array of temperature (C)
  REAL(dp), DIMENSION(kz)           :: EM_WIND                  ! Array of horizontal wind (m/s)
  REAL(dp), DIMENSION(kz)           :: EM_RH                    ! Array of relative humidity (%)
  REAL(dp), DIMENSION(kz)           :: EM_qa1                   ! Array of specific humidity
  REAL(dp), DIMENSION(kz)           :: EM_PRES                  ! Array of Pressure (Pa)
  REAL(dp), DIMENSION(kz)           :: EM_WVM                   ! Array of water vapour mixing ratio (g/g)
  REAL(dp), DIMENSION(kz)           :: EM_EW                    ! Array of water vapour pressure (Pa)
  REAL(dp), DIMENSION(kz)           :: EM_ES, EM_ES_S, EM_ES_C  ! Array of saturation water vapour pressure:
  ! ambient, sun and shade leaf temperatures
  REAL(dp), DIMENSION(8)            :: EM_Mon_Proc              ! Array of data from Jaana chemotypy paper average values for monoterpene distribution
  REAL(dp), DIMENSION(12)           :: EM_LAI_month             ! Array of data for different LAI per month
 ! REAL(dp), DIMENSION(9)            :: EM_LAI_year_Hyy          ! Array of data for different LAI per year (1996-2010)

  ! Output variables:
  REAL(dp), DIMENSION(kz)           :: EM_Sunfrac               ! Array of the fraction of sun leaves. i = 1 for top canopy layer,2 for the next layer, etc.
  REAL(dp), DIMENSION(kz)           :: EM_Sun_Par               ! Array of sun fraction - 1 above the canopy and decreases inside the canopy
  REAL(dp), DIMENSION(kz)           :: EM_SunleafTK, EM_SunleafTC        ! Array of temparture for sun leaf in (K) and (C)
  REAL(dp), DIMENSION(kz)           :: EM_ShadeleafTK, EM_ShadeleafTC    ! Array of temparture for shade leaf in (K) and (C)

  REAL(dp), DIMENSION(kz,22)        :: EM_EMI                   ! Matrix of emission for each layer and each compound
  REAL(dp), DIMENSION(kz)           :: EM_Ea1tL_M1              ! Array of emission activity of temperature per layer
  REAL(dp), DIMENSION(kz)           :: EM_Ea1pL_M1              ! Array of emission activity of light per layer
  REAL(dp), DIMENSION(kz)           :: EM_Ea1NL_M1              ! Array of companied emission activity
  REAL(dp), DIMENSION(4,kz)         :: EM_Ea1tL_M2              ! Array of emission activity of temperature per layer
  REAL(dp), DIMENSION(4,kz)         :: EM_Ea1pL_M2              ! Array of emission activity of light per layer
  REAL(dp), DIMENSION(4,kz)         :: EM_Ea1NL_M2              ! Array of companied emission activity
  REAL(dp), DIMENSION(kz,22)        :: EM_ER_BT                 ! Array of output emission buffer
  REAL(dp), DIMENSION(kz,22)        :: EM_ER_NT                 ! Array of output emission buffer
  REAL(dp), DIMENSION(kz,22)        :: EM_ER_SB                 ! Array of output emission buffer
  REAL(dp), DIMENSION(kz,22)        :: EM_ER_HB                 ! Array of output emission buffer
  REAL(dp), DIMENSION(22, kz)       :: EM_ER                    ! Array of output emission buffer
  REAL(dp), DIMENSION(22,3)         :: EM_GAM_OTHER             ! Array of other gamma factors
  REAL(dp), DIMENSION(22,kz)        :: EM_GAM_TMP               ! Array of temperature response factor
  real(dp), DIMENSION(13,NPMON-1)      :: EM_soil_voc           ! Input file containing soil emission
  real(dp), DIMENSION(12)           :: EM_soil_emi              ! Calculation of soil emission
  real(dp), DIMENSION(kz)       :: M137_soil, M33_soil, M69_soil
  INTEGER :: temp_level

  CHARACTER(16)                           :: EM_VAR(22)                        ! Name of the 22 VOC's

  ! LOGICAL, SAVE                          :: first_call = .TRUE.
  INTEGER, SAVE                          :: first_call = 1  ! 1: true, 0: false

  ! Initial values for simbim for sunny (S) and shadow (C) leaves
  REAL(dp_sb)                       :: BIM_IN_S(kz,14), BIM_IN_C(kz,14)


  !============================================================================!
  ! Variables for the parallel MPI environment
  !============================================================================!

!#ifdef PARALLEL

  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: mpi_status
  INTEGER :: mpi_rc      ! return/error code used by the MPI-subroutine calls
  INTEGER :: mpi_ntasks  ! number of MPI processes
  INTEGER :: mpi_nslaves
  INTEGER :: my_id       ! id number of own process
  INTEGER :: mpi_sendcount
  INTEGER :: mpi_recvcount
  INTEGER :: mpi_dest_id
  INTEGER :: mpi_source_id
  INTEGER :: kz_tag

  REAL(dp) :: kz_tag_dp

  INTEGER :: slave_i

  ! mpi_send_buffer for chemistry
  INTEGER, PARAMETER :: master_id = 0
  INTEGER, PARAMETER :: mpi_task_code_tag = 10
  INTEGER, PARAMETER :: mpi_do_chemistry_code = 3
  INTEGER, PARAMETER :: mpi_do_aerosol_code = 4
  INTEGER, PARAMETER :: mpi_end_task_code = -1
  INTEGER, PARAMETER :: mpi_buffer_tag = 20

  INTEGER, PARAMETER :: n_scalars = 17
  INTEGER, PARAMETER :: index_kz_tag    =  1  ! kz is INTEGER
  INTEGER, PARAMETER :: index_time1     =  2
  INTEGER, PARAMETER :: index_time2     =  3
  INTEGER, PARAMETER :: index_tmcontr   =  4
  INTEGER, PARAMETER :: index_year      =  5
  INTEGER, PARAMETER :: index_month     =  6
  INTEGER, PARAMETER :: index_nxodrad   =  7  ! nxodrad is INTEGER
  INTEGER, PARAMETER :: index_ta1       =  8  ! temperature at atmosphere level k
  INTEGER, PARAMETER :: index_pr1       =  9  ! pressure at atmosphere level k
  INTEGER, PARAMETER :: index_zenith_deg= 10
  INTEGER, PARAMETER :: index_sun_par   = 11
  INTEGER, PARAMETER :: index_res1      = 12
  INTEGER, PARAMETER :: index_res2      = 13
  INTEGER, PARAMETER :: index_qa1       = 14  ! temperature at atmosphere level k
  INTEGER, PARAMETER :: index_air       = 15  ! air number concentratio at level k
  INTEGER, PARAMETER :: index_N2        = 16  ! N2 at level k
  INTEGER, PARAMETER :: index_O2        = 17  ! O2 at level k
  INTEGER, PARAMETER :: index_H2O       = 18  ! water vapor at level k
  INTEGER, PARAMETER :: index_glob      = 19  ! incoming shortwave radiation
  INTEGER, PARAMETER :: index_albedo    = 20
  INTEGER, PARAMETER :: index1_CONS     = index_albedo+1
  INTEGER, PARAMETER :: index2_CONS     = index1_CONS + nspec - 1
  ! INTEGER, PARAMETER :: index1_res_org     = index2_CONS+1
  ! INTEGER, PARAMETER :: index2_res_org     = index1_res_org + n_cond_tot - 1
  ! INTEGER, PARAMETER :: mpi_send_buffer_size = index2_res_org !n_scalars+NSPEC
  INTEGER, PARAMETER :: mpi_send_buffer_size = index2_CONS !n_scalars+NSPEC


  ! mpi receive buffer for chemistry
  ! index_kz_tag         = 1 also here
  INTEGER, PARAMETER :: index_RO2            = 2 ! in the recv buffer
  INTEGER, PARAMETER :: index_H2O_recv       = 3 ! in the recv buffer
  INTEGER, PARAMETER :: index1_CONS_recv     = index_H2O_recv + 1
  INTEGER, PARAMETER :: index2_CONS_recv     = index1_CONS_recv + NSPEC -1
  INTEGER, PARAMETER :: index1_J_values_recv = index2_CONS_recv + 1
  INTEGER, PARAMETER :: index2_J_values_recv = index1_J_values_recv + NPHOT - 1
  INTEGER, PARAMETER :: index1_K_values_recv = index2_J_values_recv + 1
  INTEGER, PARAMETER :: index2_K_values_recv = index1_K_values_recv + NKVALUES - 1
  INTEGER, PARAMETER :: index1_ACF_recv      = index2_K_values_recv + 1
  INTEGER, PARAMETER :: index2_ACF_recv      = index1_ACF_recv + NWL - 1
  INTEGER, PARAMETER :: mpi_recv_buffer_size = index2_ACF_recv

  ! mpi aerosol buffers are in the PSD_INTERFACE

  REAL(dp),DIMENSION(mpi_send_buffer_size)      :: mpi_send_buffer
  REAL(dp),DIMENSION(mpi_recv_buffer_size)      :: mpi_recv_buffer
  ! REAL(dp),DIMENSION(mpi_send_aer_buffer_size)  :: mpi_send_aer_buffer
  ! REAL(dp),DIMENSION(mpi_recv_aer_buffer_size)  :: mpi_recv_aer_buffer

  INTEGER, PARAMETER :: mpi_yearmonth_buffer_size = 7
  CHARACTER(len=mpi_yearmonth_buffer_size) :: mpi_yearmonth_buffer

!#endif
contains


subroutine debug_message(msg, show)
  CHARACTER(*), INTENT(IN) :: msg
  LOGICAL, OPTIONAL, INTENT(IN) :: show
  LOGICAL :: temp_show

  IF (PRESENT(show)) THEN
    temp_show = show
  ELSE
    temp_show = .TRUE.
  END IF

  IF (flag_debug == 1 .AND. temp_show) THEN
    WRITE(*,*) msg
  END IF
end subroutine debug_message


subroutine init_timer(t)
  type(timer) :: t

  t%used = 0.0d0
end subroutine init_timer


subroutine start_timer(t)
  type(timer) :: t

  call system_clock(t%c0, t%cr)
end subroutine start_timer


subroutine end_timer(t)
  type(timer) :: t

  call system_clock(t%c1, t%cr)
  t%used = t%used + real(t%c1 - t%c0)/real(t%cr)
  t%last = real(t%c1 - t%c0)/real(t%cr)
end subroutine end_timer

end module Sosa_data
