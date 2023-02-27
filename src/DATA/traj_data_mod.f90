module traj_data_mod

  use constants_mod, only: dp

  implicit none

  type input_data
    character( 50) :: vname      ! variable name in the input file
    ! character( 50) :: vname_raw  ! variable name in the input file
    character(200) :: fabspath   ! input file absolute path
    real(dp), pointer :: f1d(:)
    real(dp), pointer :: f2d(:, :)
    real(dp), pointer :: f3d(:, :, :)
    ! number of dimensions, now only 1, 2 and 3, meaning f1d, f2d and f3d
    integer :: ndims
  end type input_data

  ! type(input_field), allocatable :: infield(:)

  type input_field
    integer :: ntime
    integer :: nlevel_conc, nlevel_meteo, nlevel_emisant, nlevel_emisbio, nlevel_emisaer
    integer :: nvar
    type(input_data), allocatable :: var(:)

    integer :: &
      ! prescribed concentrations
      id_conc_time       , &
      id_conc_o3         , &
      id_conc_lev        , &
      ! meteorology
      id_t               , &
      id_u               , &
      id_v               , &
      id_q               , &
      id_mla             , &
      id_lp              , &
      id_ssr             , &
      id_lsm             , &
      id_time            , &
      id_lat             , &
      id_lon             , &
      ! anthropogenic emissions
      id_emisant_time    , &
      id_emisant_blh     , &
      id_emisant_mlh     , &
      id_emisant_tlh     , &
      id_emisant_CO      , &
      id_emisant_NOx     , &
      id_emisant_NH3     , &
      id_emisant_SO2     , &
      id_emisant_C5H8    , &
      id_emisant_MT      , &
      ! biogenic emissions
      id_emisbio_time    , &
      id_emisbio_SRRsum  , &
      id_emisbio_CH3CHO  , &
      id_emisbio_HCHO    , &
      id_emisbio_CH3COCH3, &
      id_emisbio_C5H8    , &
      id_emisbio_APINENE , &
      id_emisbio_BPINENE , &
      id_emisbio_CH3OH   , &
      id_emisbio_MBO     , &
      id_emisbio_OMT     , &
      id_emisbio_SQT     , &
      id_emisbio_DMS     , &
      ! aerosol emissions
      id_emisaer_time    , &
      id_emisaer_nconc   , &
      id_emisaer_lev     , &
      id_emisaer_3_10nm  , &
      id_emisaer_10_20nm , &
      id_emisaer_20_30nm , &
      id_emisaer_30_50nm , &
      id_emisaer_50_70nm , &
      id_emisaer_70_100nm   , &
      id_emisaer_100_200nm  , &
      id_emisaer_200_400nm  , &
      id_emisaer_400_1000nm

  end type input_field

  type(input_field) :: infield

  real(dp), allocatable :: input_levels_curr(:)  ! [m], input level heights at current time

  ! Input time pointer, the time index which current time points to, the model time is between
  ! input_time(pointer) and input_time(pointer+1)
  integer :: input_ptr

end module traj_data_mod
