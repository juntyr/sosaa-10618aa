!
! The input and output module.
! Try to avoid using these predefined file units: 0, 5, 6, 100, 101, 102
!   0: stderr, some compilers may use 100
!   5: stdin, some compilers may use 101
!   6: stoutd, some compilers may use 102
!
MODULE SOSA_IO

  USE SOSA_DATA
  USE second_Monitor, ONLY:SPC_NAMES
  USE second_Parameters, kpp_dp=>dp
  USE dios_mod

  IMPLICIT NONE

  PRIVATE

  public :: READ_INPUTS, CONVERT_ECMWF, open_output_files !, vapor_no_uhma
  public :: output_init
  public :: output_step
  public :: output_write
  public :: output_write_old
  public :: output_done
  public :: output_netcdf_done
  public :: linear_interp
  public :: outspc_count, outspc_names, outspc_cheminds
  public :: outemi_count, outemi_names, outemi_cheminds, outemi_meganinds
  public :: outvd_cheminds, outflux_cheminds, outvap_vapinds

  INTEGER, PARAMETER :: OUTPUT_NUMBER_MAX = 400
  REAL(kind=dp), DIMENSION(kz, OUTPUT_NUMBER_MAX) :: output_array  ! save output variable values
  CHARACTER(LEN=15), DIMENSION(OUTPUT_NUMBER_MAX) :: output_name   ! save output variable names
  INTEGER :: output_number = 0
  INTEGER, PARAMETER :: UNIT_START = 100
  LOGICAL :: isopen(OUTPUT_NUMBER_MAX)

  CHARACTER(LEN=LINE_WIDTH) :: rawline
  CHARACTER(LEN=SPC_NAME_LEN) :: gname  ! save the temporary gas or species name

  LOGICAL :: exists
  CHARACTER(LEN=350) :: fname_init

  INTEGER :: ntime_local, ntime_ECMWF  ! number of time points for local data and ECMWF dataset


  !-----------------------------------------------------------------------------
  ! NetCDF related types and variables
  !-----------------------------------------------------------------------------

  type metafields
    character(40) :: vname  ! variable name
    character(80) :: lname  ! variable long name
    character(20) :: unit   ! variable unit
    ! character(150) :: standard_name  ! variable standard name, not used now
  end type metafields

  ! Define field types, here dimids does not include time dimension
  type field0d
    type(metafields)                :: mf
    real(dp)                        :: field
    integer                         :: varid
    integer                         :: dimids  ! not used, set to -1
    integer                         :: ndims   ! not used, set to -1
    logical                         :: change_with_time
  end type field0d

  type field1d
    type(metafields)                    :: mf
    real(dp), dimension(:), allocatable :: field
    integer                             :: varid
    integer, dimension(1)               :: dimids
    integer, dimension(1)               :: ndims
    logical                             :: change_with_time
  end type field1d

  type field2d
    type(metafields)                       :: mf
    real(dp), dimension(:, :), allocatable :: field
    integer                                :: varid
    integer, dimension(2)                  :: dimids
    integer, dimension(2)                  :: ndims
    logical                                :: change_with_time
  end type field2d

  type field3d
    type(metafields)                          :: mf
    real(dp), dimension(:, :, :), allocatable :: field
    integer                                   :: varid
    integer, dimension(3)                     :: dimids
    integer, dimension(3)                     :: ndims
    logical                                   :: change_with_time
  end type field3d

  type field4d
    type(metafields)                             :: mf
    real(dp), dimension(:, :, :, :), allocatable :: field
    integer                                      :: varid
    integer, dimension(4)                        :: dimids
    integer, dimension(4)                        :: ndims
    logical                                      :: change_with_time
  end type field4d

  type mixfile
    ! Here fields do not include time dimension, so only 1d and 2d are included.
    ! Notice: different variables with the same dimensions may have different
    ! shapes, e.g., nconc(outspc_count, kz), nconc_par(n_bins_par, kz)
    type(field0d), dimension(:), pointer :: f0d
    type(field1d), dimension(:), pointer :: f1d
    type(field2d), dimension(:), pointer :: f2d
    ! type(field1d), dimension(:), pointer :: f1d_z   ! (lev)
    ! type(field2d), dimension(:), pointer :: f2d_dz  ! (diameter, lev)
    ! type(field2d), dimension(:), pointer :: f2d_cz  ! (icondv, lev)

    character(200) :: fname
    integer        :: funit

    ! dimenstion id
    integer :: dimid_time, dimid_lev
    integer :: dimid_outspc
    integer :: dimid_dp_dry_fs, dimid_icondv

    ! dimenstion varid
    integer :: varid_time, varid_lev
    integer :: varid_outspc, varid_outflux
    integer :: varid_dp_dry_fs, varid_icondv
  end type mixfile

  type(mixfile) :: mixf

  ! Maximum number of fields
  integer, parameter :: NMAX_0D=100, NMAX_1D=100, NMAX_2D=100

  ! Current number of fields
  integer :: n_0d, n_1d, n_2d

  ! Indices of variables in the mixf fields
  integer :: id0_lat, id0_lon
  integer :: id0_pblh, id0_gsoil, id0_zenith
  integer :: id0_Rdir_down, id0_Rdiff_down, id0_glob
  integer :: id0_NIR_down, id0_NIR_up, id0_PAR_down, id0_PAR_up
  integer :: id0_L_down, id0_L_up

  integer :: id1_lad_meteo, id1_lad_megan
  integer :: id1_ua, id1_va
  integer :: id1_temp, id1_rhov, id1_pres, id1_Nair
  integer :: id1_tke, id1_diffm, id1_diffh, id1_Rich, id1_ustar, id1_mixlen
  integer :: id1_shf, id1_lhf

  integer :: id2_nconc, id2_flux, id2_nconc_par

  INTEGER :: stat  ! record io status

  INTEGER :: itrec  ! current time record
  INTEGER :: varid_lad_meteo, varid_lad_megan
  INTEGER :: varid_temp, varid_rhov, varid_rh, varid_pres, varid_ua, varid_va
  INTEGER :: varid_tke, varid_mixlen, varid_diffm, varid_diffh, varid_ri
  INTEGER :: varid_ustar, varid_shf, varid_lhf
  INTEGER :: varid_nair
  INTEGER :: varid_ceb, varid_gsoil, varid_pblh, varid_albedo, varid_zenith
  INTEGER :: varid_rd_dir, varid_rd_dif, varid_ld, varid_lu, varid_paru, varid_niru
  INTEGER :: varid_nconc_particle, varid_flux_particle, varid_growth_rate, varid_nuc_rate
  INTEGER :: varid_nconc_condv_gas, varid_vconc_condv_par, varid_sink_condv_gas

  type netcdf_general
    CHARACTER(LEN=300) :: fname
    INTEGER :: fid
    INTEGER :: dimid_lev
    INTEGER :: varid_lev
  end type netcdf_general

  type netcdf_meteo
    CHARACTER(LEN=300) :: fname
    INTEGER :: fid
    INTEGER :: dimid_lev, dimid_time
    INTEGER :: dimids_2d(2)  ! (/dimid_lev, dimid_time/)
    INTEGER :: varid_lev, varid_time
  end type netcdf_meteo

  type netcdf_chem
    CHARACTER(LEN=300) :: fname
    INTEGER :: fid
    INTEGER :: dimid_lev, dimid_time
    INTEGER :: dimids_2d(2)
    INTEGER :: varid_lev, varid_time
    ! INTEGER, ALLOCATABLE :: varid_nconc(:)
    ! INTEGER, ALLOCATABLE :: varid_flux(:)
    ! INTEGER, ALLOCATABLE :: varid_vd(:)
    ! INTEGER, ALLOCATABLE :: varid_emi(:)
  end type netcdf_chem

  type netcdf_aerosol
    CHARACTER(LEN=300) :: fname
    INTEGER :: fid
    INTEGER :: dimid_dp_dry_fs, dimid_lev, dimid_time
    INTEGER :: dimids_3d_dp_dry_fs(3)
    INTEGER :: dimids_2d(2)
    INTEGER :: varid_dp_dry_fs, varid_lev, varid_time
  end type netcdf_aerosol

  type netcdf_Tvapor
    CHARACTER(LEN=300) :: fname
    INTEGER :: fid
    INTEGER :: dimid_icondv, dimid_lev, dimid_time
    INTEGER :: dimids_3d_condv(3)
    INTEGER :: dimids_2d_condv(2)
    INTEGER :: varid_icondv, varid_lev, varid_time
  end type netcdf_Tvapor

  type netcdf_SBVapor
    CHARACTER(LEN=300) :: fname
    INTEGER :: fid
    INTEGER :: dimid_icondv, dimid_dp_dry_fs, dimid_lev, dimid_time
    INTEGER :: dimids_4d(4)
    INTEGER :: dimids_3d(3)
    INTEGER :: varid_icondv, varid_dp_dry_fs, varid_lev, varid_time
  end type netcdf_SBVapor

  type(netcdf_general) :: nc_general
  type(netcdf_meteo)   :: nc_meteo
  type(netcdf_chem)    :: nc_chem
  type(netcdf_aerosol) :: nc_aerosol
  type(netcdf_Tvapor)  :: nc_Tvapor
  type(netcdf_SBVapor) :: nc_SBVapor

  ! INTEGER, PARAMETER :: MAX_OUTPUT_SPC_COUNT = 100

  ! Output species list
  INTEGER                                  :: outspc_count       ! number of output species
  CHARACTER(LEN=SPC_NAME_LEN), ALLOCATABLE :: outspc_names(:)     ! name list of output species
  INTEGER                    , ALLOCATABLE :: outspc_cheminds(:)  ! indices of output species in the chemistry scheme
  INTEGER                                  :: outflux_count       ! number of output species
  CHARACTER(LEN=SPC_NAME_LEN), ALLOCATABLE :: outflux_names(:)     ! name list of output species
  INTEGER                    , ALLOCATABLE :: outflux_cheminds(:)  ! indices of output species in the chemistry scheme
  INTEGER                                  :: outvd_count       ! number of output species
  CHARACTER(LEN=SPC_NAME_LEN), ALLOCATABLE :: outvd_names(:)     ! name list of output species
  INTEGER                    , ALLOCATABLE :: outvd_cheminds(:)  ! indices of output species in the chemistry scheme
  INTEGER                                  :: outvap_count       ! number of output species
  CHARACTER(LEN=SPC_NAME_LEN), ALLOCATABLE :: outvap_names(:)     ! name list of output species
  INTEGER                    , ALLOCATABLE :: outvap_cheminds(:)  ! indices of output species in the chemistry scheme
  INTEGER                    , ALLOCATABLE :: outvap_vapinds(:)  ! indices of output species in the vapor list

  ! Emission species list
  INTEGER                                  :: outemi_count        ! number of emitted species
  CHARACTER(LEN=SPC_NAME_LEN), ALLOCATABLE :: outemi_names(:)      ! name list of emitted species
  INTEGER                    , ALLOCATABLE :: outemi_cheminds(:)   ! indices of emitted species in the chemistry scheme
  INTEGER                    , ALLOCATABLE :: outemi_meganinds(:)  ! indices of emitted species in MEGAN

  INTEGER, ALLOCATABLE :: varid_outspc_nconc(:), varid_outspc_flux(:), varid_outspc_vd(:),varid_nconc_condv_gas_one(:), varid_vconc_condv_par_one(:),varid_sink_condv_gas_one(:)
  INTEGER, ALLOCATABLE :: varid_outemi(:)

  ! Reactivity list
  INTEGER, ALLOCATABLE :: varid_outspc_rOH(:)
  INTEGER, ALLOCATABLE :: varid_outspc_rO3(:)
  INTEGER, ALLOCATABLE :: varid_outspc_rNO3(:)
  INTEGER, ALLOCATABLE :: varid_outspc_reactivity(:)

  ! Photolysis rate list
  INTEGER, ALLOCATABLE :: varid_J(:)

  ! HOM10 species list
  INTEGER                                  :: outspc_HOM10_count    ! number of output species
  CHARACTER(LEN=SPC_NAME_LEN), ALLOCATABLE :: outspc_HOM10_name(:)  ! name list of output species
  INTEGER                    , ALLOCATABLE :: outspc_HOM10_ind(:)   ! indices of output species in the chemistry scheme

  ! HOM20 species list
  INTEGER                                  :: outspc_HOM20_count    ! number of output species
  CHARACTER(LEN=SPC_NAME_LEN), ALLOCATABLE :: outspc_HOM20_name(:)  ! name list of output species
  INTEGER                    , ALLOCATABLE :: outspc_HOM20_ind(:)   ! indices of output species in the chemistry scheme

  ! Temporary array to save the position of each species in the name list
  INTEGER, ALLOCATABLE :: tmp_pos(:, :)

  ! A string of a name list separated by comma (space is ignored)
  ! Example:
  !   output_list_spc = "OH, O3, SO2"
  !   output_list_emi = "APINENE, BPINENE"
  CHARACTER(LEN=1000) :: output_layer_id
  CHARACTER(LEN=1000) :: output_list_spc  ! number concentration, flux, dry deposition velocity
  CHARACTER(LEN=1000) :: output_list_emi  ! emissions
  CHARACTER(LEN=1000) :: output_list_vd
  CHARACTER(LEN=1000) :: output_list_vap

  NAMELIST /NML_OUTPUT/ output_layer_id, &
    output_list_spc, output_list_emi, output_list_vd, output_list_vap

  CHARACTER(50) :: start_date_string, first_day_of_month_string


CONTAINS

!
!===============================================================================================================================
!

!> Read namelists from initiation file fname_init. Default values will be used if not specified in the namelist.
!! The namelists are defined in DATA/Sosa_data.f90.
subroutine read_input_init()
  ! Get the init file name from the command line
  CALL GETARG(1, fname_init)  ! Get the first argument of the executable command

  ! Print debug information
  WRITE(*,*) 'Process ', my_id, 'is reading namelists from ', &
    TRIM(ADJUSTL(fname_init))

  ! Open the file read and then close it
  OPEN(UNIT=99, FILE=TRIM(ADJUSTL(fname_init)), STATUS='OLD')
  READ(UNIT=99, NML=NML_MAIN  , IOSTAT=stat)  ! directories and station
  READ(UNIT=99, NML=NML_FLAG  , IOSTAT=stat)  ! scheme flags
  READ(UNIT=99, NML=NML_MISC  , IOSTAT=stat)  ! misc options
  READ(UNIT=99, NML=AER_FLAG  , IOSTAT=stat)  ! aerosol flags
  READ(UNIT=99, NML=NML_GRID  , IOSTAT=stat)  ! e.g., masl
  READ(UNIT=99, NML=NML_TIME  , IOSTAT=stat)  ! start and end time, time steps
  READ(UNIT=99, NML=NML_OUTPUT, IOSTAT=stat)  ! start and end time, time steps
  CLOSE(99)

end subroutine read_input_init


!==============================================================================!
! Read the input from different files
!==============================================================================!
subroutine READ_INPUTS(stage)
  integer :: stage
  ! Read init file
  if (stage == 1) THEN

  call read_input_init()

  !--------------------------------------------------------!
  ! Get date infomation
  !--------------------------------------------------------!

  ! Set current date
  now_date = start_date

  ! Set current month number
  mon = now_date(2)

  ! Display debug information
  if (my_id == master_id) then
    WRITE(*,*) 'Simulation starts on ', start_date
    WRITE(*,*) 'Simulation ends on ', end_date
  end if

  ! Prepare the strings of time variables
  WRITE(year_str     , '(I4)'      ) now_date(1)            ! '2000', '2001', ...
  WRITE(year_str2    , '(I2.2)'    ) MOD(now_date(1), 100)  ! '00', '01', '18', ...
  WRITE(month_str_Mxx, '(A1, I2.2)') 'M', now_date(2)       ! 'M01', 'M02', ..., 'M12'
  WRITE(month_str    , '(I2.2)'    ) now_date(2)            ! '01', '02', ..., '12'
  WRITE(day_str      , '(I2.2)'    ) now_date(3)            ! '01', '02', ..., '31'

  !--------------------------------------------------------!
  ! Set output and input dirs.
  ! The previous filename1 is substituted by the INPUT_DIR.
  !--------------------------------------------------------!

  input_dir_general = TRIM(ADJUSTL(INPUT_DIR)) // '/general'
  input_dir_station = TRIM(ADJUSTL(INPUT_DIR)) // '/station' // '/' // STATION
  input_dir_station_info = TRIM(ADJUSTL(input_dir_station)) // '/info'
  input_dir_station_data = TRIM(ADJUSTL(input_dir_station)) //  '/' // year_str // '/' // month_str_Mxx ! data input

  ! Check if the output directory path exists, if not, create folders needed.
  INQUIRE(FILE=TRIM(ADJUSTL(output_dir)) // '/.', EXIST=exists)
  IF (.NOT. exists) THEN
    CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(output_dir)))
  ELSE
      IF (FLAG_SKIP_EXISTING) then
          inquire( file=TRIM(ADJUSTL(output_dir))//'/output.nc', exist=exists)
          if (exists) &
            STOP 'Skipping this case, existing already, FLAG_SKIP_EXISTING = .true.'
      end if
  END IF

  ! Display debug information
  if (my_id == master_id) then
    write(*,*) 'WORK_DIR = ', TRIM(WORK_DIR)  ! SOSAA project folder
    write(*,*) 'CODE_DIR = ', TRIM(CODE_DIR)  ! code folder
    write(*,*) 'CASE_DIR = ', TRIM(CASE_DIR)  ! case folder
    write(*,*) 'CHEM_DIR = ', TRIM(CHEM_DIR)  ! chemistry scheme folder
    write(*,*) 'INPUT_DIR = ', TRIM(INPUT_DIR)
    write(*,*) 'input_dir_general = ', TRIM(input_dir_general)
    write(*,*) 'input_dir_station = ', TRIM(input_dir_station)
    write(*,*) 'input_dir_station_info = ', TRIM(input_dir_station_info)
    write(*,*) 'input_dir_station_data = ', TRIM(input_dir_station_data)
    write(*,*) 'output_dir = ', TRIM(output_dir)
  end if

  ! Read names and properties of condensable vapors
  ! IF (flag_aero==0) &
    ! NOTE XXX replace this BS with VAPOUR INIT
    ! CALL vapor_no_uhma(trim(adjustl(chem_dir)) // '/INPUT/condensable_vapors', VAPOUR_PROP%vapour_names,VAPOUR_PROP%psat_a,VAPOUR_PROP%psat_b)

elseif (stage==2) THEN
  ! Get outspc_count, outspc_names, outspc_cheminds and outemi_count, outemi_names, outemi_cheminds, outemi_meganinds
  IF (flag_outlist==0) THEN
    CALL get_output_spc_list_from_string()
  ELSE
    CALL get_output_spc_list_from_file()
  END IF

  ! Print debug information
  if (my_id == master_id) then
    write(*,*) 'Output concentration list with chemistry index and name:'
    do i = 1, outspc_count
      write(*,*) outspc_cheminds(i), trim(adjustl(outspc_names(i)))
    end do

    write(*,*) 'Output emission list with chemistry index, megan index and name:'
    do i = 1, outemi_count
      write(*,*) outemi_cheminds(i), outemi_meganinds(i), trim(adjustl(outemi_names(i)))
    end do

    write(*,*) 'Output vd list with chemistry index and name:'
    do i = 1, outvd_count
      write(*,*) outvd_cheminds(i), trim(adjustl(outvd_names(i)))
    end do

    write(*,*) 'Output vapor list with chemistry index, vap index and name (this has not been really added to the output file):'
    do i = 1, outvap_count
      write(*,*) outvap_cheminds(i), outvap_vapinds(i), trim(adjustl(outvap_names(i)))
    end do
  end if

  !--------------------------------------------------------!
  ! Read input of measurement data
  !--------------------------------------------------------!
  call read_input_data()
END IF

CONTAINS


  subroutine read_input_data()
    select case (flag_model_type)
      case (STATION_MODE)
        ! For station hyytiala with text file input
        call read_input_data_hyytiala_ascii()
        ! For station hyytiala with netcdf file input
        ! call read_input_data_hyytiala_netcdf()
        ! For very old version hyytiala data files, not updated and not recommended
        ! call READ_DATA_OLD()
      case (TRAJECTORY_MODE)
        call read_input_data_traj_netcdf()
      case default
        if (my_id == master_id) then
          write(*,*) 'Station ', trim(adjustl(STATION)), ' not found.'
        end if
        call MPI_Finalize(mpi_rc)
        stop
    end select
  end subroutine read_input_data


  subroutine get_output_spc_list_from_string()
    !----------------------------------------------------------------------------!
    ! Get the output species list
    !----------------------------------------------------------------------------!
    ! Allocate temporary position array

    ALLOCATE(tmp_pos(len_trim(adjustl(output_list_spc)), 2))
    ! Get outspc_count
    CALL string_split_trim2(trim(adjustl(output_list_spc)), ',', outspc_count, tmp_pos)


    ! Get outspc_names
    ALLOCATE(outspc_names(outspc_count))

    CALL get_substring_list(trim(adjustl(output_list_spc)), outspc_count, tmp_pos, outspc_names)

    ! Get outspc_cheminds
    ALLOCATE(outspc_cheminds(outspc_count))
    CALL get_inds_outnames_in_rawnames(outspc_names, SPC_NAMES, outspc_cheminds)

    DEALLOCATE(tmp_pos)

    !----------------------------------------------------------------------------!
    ! Get the output Vd species list
    !----------------------------------------------------------------------------!
    ! Allocate temporary position array
    ALLOCATE(tmp_pos(len_trim(adjustl(output_list_Vd)), 2))

    ! Get outvd_count
    CALL string_split_trim2(trim(adjustl(output_list_Vd)), ',', outvd_count, tmp_pos)

    ! Get outvd_names
    ALLOCATE(outvd_names(outvd_count))
    CALL get_substring_list(trim(adjustl(output_list_Vd)), outvd_count, tmp_pos, outvd_names)

    ! Get outvd_cheminds
    ALLOCATE(outvd_cheminds(outvd_count))
    CALL get_inds_outnames_in_rawnames(outvd_names, SPC_NAMES, outvd_cheminds)

    DEALLOCATE(tmp_pos)


    !----------------------------------------------------------------------------!
    ! Get the output emission species list
    !----------------------------------------------------------------------------!
    ! Allocate temporary position array
    ALLOCATE(tmp_pos(len_trim(adjustl(output_list_emi)), 2))

    ! Get outemi_count
    CALL string_split_trim2(trim(adjustl(output_list_emi)), ',', outemi_count, tmp_pos)

    ! Get outemi_names
    ALLOCATE(outemi_names(outemi_count))
    CALL get_substring_list(trim(adjustl(output_list_emi)), outemi_count, tmp_pos, outemi_names)

    ! Get outemi_cheminds, outemi_meganinds
    ALLOCATE(outemi_cheminds(outemi_count))
    ALLOCATE(outemi_meganinds(outemi_count))
    CALL get_inds_outnames_in_rawnames(outemi_names, SPC_NAMES, outemi_cheminds)
    CALL get_inds_outnames_in_rawnames(outemi_names, MEGAN_SPC_NAMES, outemi_meganinds)

    DEALLOCATE(tmp_pos)


    !----------------------------------------------------------------------------!
    ! Get the output Vapor species list
    !----------------------------------------------------------------------------!
    ! Allocate temporary position array
    ALLOCATE(tmp_pos(len_trim(adjustl(output_list_Vap)), 2))

    ! Get outvap_count
    CALL string_split_trim2(trim(adjustl(output_list_Vap)), ',', outvap_count, tmp_pos)

    ! Get outvap_names
    ALLOCATE(outvap_names(outvap_count))
    CALL get_substring_list(trim(adjustl(output_list_Vap)), outvap_count, tmp_pos, outvap_names)

    ! Get outvap_cheminds
    ALLOCATE(outvap_cheminds(outvap_count))
    CALL get_inds_outnames_in_rawnames(outvap_names, SPC_NAMES, outvap_cheminds)

    ! Get outspc_vapinds
    ALLOCATE(outvap_vapinds(outvap_count))
    CALL get_inds_outnames_in_rawnames(outvap_names, VAPOUR_PROP%vapour_names, outvap_vapinds)

    DEALLOCATE(tmp_pos)

  end subroutine get_output_spc_list_from_string


  subroutine get_output_spc_list_from_file()
    !--------------------------------------------------------------------------!
    ! Get the information of species for output from outspc.txt
    !--------------------------------------------------------------------------!
    ! Read the line count of outspc.txt
    CALL GET_FILE_LINE_COUNT(TRIM(ADJUSTL(CHEM_DIR)) // '/INPUT/outspc.txt', outspc_count)

    ! Read the file again to get outspc_name and outspc_cheminds.
    ! Here we assume that the names are unique.
    ALLOCATE(outspc_names(outspc_count), outspc_cheminds(outspc_count))
    CALL GET_SPC_LIST_FROM_FILE(TRIM(ADJUSTL(CHEM_DIR)) // '/INPUT/outspc.txt', &
      outspc_names, outspc_cheminds, spc_names)

    !--------------------------------------------------------------------------!
    ! Get the information of HOM species for output from HOM10.txt and HOM20.txt
    !--------------------------------------------------------------------------!
    IF (flag_aero == 1) THEN
      ! Read the line count of HOM10.txt
!      CALL GET_FILE_LINE_COUNT(TRIM(ADJUSTL(CHEM_DIR)) // '/INPUT/HOM10.txt', outspc_HOM10_count)

      ! Read the file again to get outspc_name and outspc_ind.
      ! Here we assume that the names are unique.
!      ALLOCATE(outspc_HOM10_name(outspc_HOM10_count), outspc_HOM10_ind(outspc_HOM10_count))
!      CALL GET_SPC_LIST_FROM_FILE(TRIM(ADJUSTL(CHEM_DIR)) // '/INPUT/HOM10.txt', outspc_HOM10_name, outspc_HOM10_ind)

      ! Read the line count of HOM20.txt
!      CALL GET_FILE_LINE_COUNT(TRIM(ADJUSTL(CHEM_DIR)) // '/INPUT/HOM20.txt', outspc_HOM20_count)

      ! Read the file again to get outspc_name and outspc_ind.
      ! Here we assume that the names are unique.
!      ALLOCATE(outspc_HOM20_name(outspc_HOM20_count), outspc_HOM20_ind(outspc_HOM20_count))
!      CALL GET_SPC_LIST_FROM_FILE(TRIM(ADJUSTL(CHEM_DIR)) // '/INPUT/HOM20.txt', outspc_HOM20_name, outspc_HOM20_ind)
      ! write(*,*) 'outspc_count ', outspc_count
      ! write(*,*) 'outspc_name ', outspc_name
      ! write(*,*) 'outspc_ind ', outspc_ind
    END IF

    !--------------------------------------------------------------------------!
    ! Get the information of species for emission output from outemi.txt
    !--------------------------------------------------------------------------!
    ! Read the line count of outemi.txt
    CALL GET_FILE_LINE_COUNT(TRIM(ADJUSTL(CHEM_DIR)) // '/INPUT/outemi.txt', outemi_count)

    ! Read the file again to get emispc_name and emispc_ind.
    ! Here we assume that the names are unique.
    ALLOCATE(outemi_names(outemi_count), outemi_cheminds(outemi_count), outemi_meganinds(outemi_count))
    CALL GET_SPC_LIST_FROM_FILE(TRIM(ADJUSTL(CHEM_DIR)) // '/INPUT/outemi.txt', &
      outemi_names, outemi_cheminds, spc_names)
    CALL GET_SPC_MEGAN_IND(outemi_names, outemi_meganinds)
    !--------------------------------------------------------------------------!
    ! Get the information of species for vapors from outflux.txt
    !--------------------------------------------------------------------------!
    ! Read the line count of outflux.txt
    CALL GET_FILE_LINE_COUNT(TRIM(ADJUSTL(CHEM_DIR)) // '/INPUT/outflux.txt', outflux_count)

    ! Read the file again to get outflux_name and outflux_cheminds.
    ! Here we assume that the names are unique.
    ALLOCATE(outflux_names(outflux_count), outflux_cheminds(outflux_count))
    CALL GET_SPC_LIST_FROM_FILE(TRIM(ADJUSTL(CHEM_DIR)) // '/INPUT/outflux.txt', &
      outflux_names, outflux_cheminds, VAPOUR_PROP%vapour_names)
    !--------------------------------------------------------------------------!
    ! Get the information of species for vapors from outvd.txt
    !--------------------------------------------------------------------------!
    ! Read the line count of outvd.txt
    CALL GET_FILE_LINE_COUNT(TRIM(ADJUSTL(CHEM_DIR)) // '/INPUT/outvd.txt', outvd_count)

    ! Read the file again to get outflux_name and outflux_cheminds.
    ! Here we assume that the names are unique.
    ALLOCATE(outvd_names(outvd_count), outvd_cheminds(outvd_count))
    CALL GET_SPC_LIST_FROM_FILE(TRIM(ADJUSTL(CHEM_DIR)) // '/INPUT/outvd.txt', &
      outvd_names, outvd_cheminds, VAPOUR_PROP%vapour_names)
    !--------------------------------------------------------------------------!
    ! Get the information of species for vapors from outvap.txt
    !--------------------------------------------------------------------------!
    ! Read the line count of outvap.txt
    CALL GET_FILE_LINE_COUNT(TRIM(ADJUSTL(CHEM_DIR)) // '/INPUT/outvap.txt', outvap_count)

    ! Read the file again to get outvap_name and outvap_vapinds.
    ! Here we assume that the names are unique.
    ALLOCATE(outvap_names(outvap_count), outvap_vapinds(outvap_count))
    CALL GET_SPC_LIST_FROM_FILE(TRIM(ADJUSTL(CHEM_DIR)) // '/INPUT/outvap.txt', &
      outvap_names, outvap_vapinds, VAPOUR_PROP%vapour_names)
  end subroutine get_output_spc_list_from_file

end subroutine READ_INPUTS


subroutine GET_FILE_LINE_COUNT(fname, n)
  ! Count the non-empty lines

  CHARACTER(LEN=*), INTENT(IN) :: fname
  INTEGER, INTENT(OUT) :: n

  ! Initial value
  n = 0

  ! Open the file
  OPEN(UNIT=99, FILE=TRIM(ADJUSTL(fname)))
  REWIND(99)

  ! Count the line number
  DO
    READ(99, '(a)', IOSTAT=stat) rawline  ! read every line including empty lines
    IF (stat /= 0) EXIT                   ! exit if eof
    IF (LEN(TRIM(rawline)) > 0) n=n+1          ! add one line count if the line is not empty
  END DO

  ! Close the file
  CLOSE(99)
end subroutine GET_FILE_LINE_COUNT


subroutine GET_SPC_LIST_FROM_FILE(fname, spc_name, spc_ind, full_NAMES)
  CHARACTER(LEN=*           ), INTENT(IN   ) :: fname
  CHARACTER(LEN=*           ), INTENT(IN   ) :: full_names(:)
  CHARACTER(LEN=SPC_NAME_LEN), INTENT(  OUT) :: spc_name(:)
  INTEGER                    , INTENT(  OUT) :: spc_ind(:)

  INTEGER :: spc_count, ispc

  ! Number of output species
  spc_count = SIZE(spc_name)

  ! Read spc names one by one
  ispc = 0

  ! Open the file and go back the beginning then read the spc names line by line
  OPEN(UNIT=99, FILE=TRIM(ADJUSTL(fname)))
  REWIND(99)
  DO
    READ(99, '(a)', IOSTAT=stat) rawline  ! read every line including empty lines
    IF (stat /= 0) EXIT                   ! exit if eof
    ! Add one species if the line is not empty
    IF (LEN(TRIM(rawline)) > 0) THEN
      ispc=ispc+1
      spc_name(ispc) = TRIM(rawline)
    END IF
  END DO
  CLOSE(99)

  ! Obtain indices of output species (case insensative)
  ! full_NAMES is from chemistry module or vapor list
  spc_ind(:) = -1  ! default value
  DO I=1, spc_count
    DO J=1, size(full_names)
      IF ( upper( TRIM(ADJUSTL( full_NAMES(J) )) ) == upper( TRIM(ADJUSTL(spc_name(I))) ) ) THEN
        spc_ind(I) = J
        EXIT
      END IF
    END DO

    ! Show error message if spc not found
    if (spc_ind(I) == -1) then
      write(*,*) 'Error: The species ', trim(adjustl(spc_name(I))), ' is not found.'
      stop
    end if
  END DO
end subroutine geT_SPC_LIST_FROM_FILE


subroutine get_spC_MEGAN_IND(spc_name, spc_mind)
  CHARACTER(LEN=SPC_NAME_LEN), INTENT(IN   ) :: spc_name(:)
  INTEGER                    , INTENT(  OUT) :: spc_mind(:)

  INTEGER :: spc_count

  spc_count = SIZE(spc_name)

  spc_mind(:) = -1  ! default value
  DO I=1, spc_count
    DO J=1, MEGAN_NSPEC
      IF ( upper( TRIM(ADJUSTL( MEGAN_SPC_NAMES(J) )) ) == upper( TRIM(ADJUSTL(spc_name(I))) ) ) THEN
        spc_mind(I) = J
        EXIT
      END IF
    END DO
  END DO
end subroutine geT_SPC_MEGAN_IND


subroutine read_input_data_hyytiala_ascii()
  ! Save original data from input files, first 6 columns are date, then for different height levels.
  ! First line is header, first 6 are 0, then height levels.
  ! The maximum time points are 48*31 + 1, then add another 1 for header. So totally 1490.
  ! Now 15 is enough for at most 7 height levels.
  ! For ECMWF data, the maximum time points are 8*31+1 (249).
  REAL(dp) :: raw_data(1490, 15)

  INTEGER :: i, j

  ! Number of time points for local dataset.
  ! yyyymm01-00:00:00 to yyyy(mm+1)01-00:00:00, so 1 is added.
  ntime_local = MonthDay(now_date(1), now_date(2))*48 + 1

  !
  ! Input of measured gases O3, SO2, NO, NO2, CO and air (every half-hour data)
  !

  ! O3
  OPEN(UNIT=101, FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_O3avg.txt', STATUS='old')
  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
  CH_gas_hyy(2, 1:ntime_local) = raw_data(2:ntime_local+1, 7)
  CLOSE(UNIT=101)

  ! SO2
  OPEN(UNIT=101, FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_SO2avg.txt', STATUS='old')
  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
  CH_gas_hyy(3, 1:ntime_local) = raw_data(2:ntime_local+1, 7)
  CLOSE(UNIT=101)

  ! NO
  OPEN(UNIT=101, FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_NOavg.txt', STATUS='old')
  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
  CH_gas_hyy(4, 1:ntime_local) = raw_data(2:ntime_local+1, 7)
  CLOSE(UNIT=101)

  ! NO2
  OPEN(UNIT=101, FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_NO2avg.txt', STATUS='old')
  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
  CH_gas_hyy(5, 1:ntime_local) = raw_data(2:ntime_local+1, 7)
  CLOSE(UNIT=101)

  ! CO
  OPEN(UNIT=101, FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_COavg.txt', STATUS='old')
  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
  CH_gas_hyy(6, 1:ntime_local) = raw_data(2:ntime_local+1, 7)
  CLOSE(UNIT=101)

  ! NH3, if ACDC is used
  IF (flag_aero==1) THEN
    IF (flag_NH3==1) THEN
      OPEN(UNIT=101, FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_NH3avg.txt', STATUS='old')
    ELSE IF (flag_NH3==2) THEN  !if number concentration is used
      OPEN(UNIT=101, FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_NH3avg_NC.txt', STATUS='old')
    END IF
    READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
    NH3_hyy(1:ntime_local) = raw_data(2:ntime_local+1, 7)
    CLOSE(UNIT=101)
  END IF

  if(flag_emis==0)then

  OPEN(UNIT=101, FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_mono_4_2.txt', STATUS='old')
  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
  CH_VOC_hyy(1,1:ntime_local) = raw_data(2:ntime_local+1, 7)
  CLOSE(UNIT=101)

  OPEN(UNIT=101, FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_C5H8_MBO_4_2.txt', STATUS='old')
  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
  CH_VOC_hyy(2,1:ntime_local) = raw_data(2:ntime_local+1, 7)
  CLOSE(UNIT=101)

  OPEN(UNIT=101, FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_MVK_4_2.txt', STATUS='old')
  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
  CH_VOC_hyy(3,1:ntime_local) = raw_data(2:ntime_local+1, 7)
  CLOSE(UNIT=101)

!  OPEN(UNIT=101, FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_MEK_4_2.txt', STATUS='old')
!  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
!  CH_VOC_hyy(4,1:ntime_local) = raw_data(2:ntime_local+1, 7)
!  CLOSE(UNIT=101)

  OPEN(UNIT=101, FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_CH3OH_4_2.txt', STATUS='old')
  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
  CH_VOC_hyy(5,1:ntime_local) = raw_data(2:ntime_local+1, 7)
  CLOSE(UNIT=101)

  OPEN(UNIT=101, FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_CH3CHO_4_2.txt', STATUS='old')
  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
  CH_VOC_hyy(6,1:ntime_local) = raw_data(2:ntime_local+1, 7)
  CLOSE(UNIT=101)

!  OPEN(UNIT=101, FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_C2H5OH_4_2.txt', STATUS='old')
!  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
!  CH_VOC_hyy(7,1:ntime_local) = raw_data(2:ntime_local+1, 7)
!  CLOSE(UNIT=101)

  OPEN(UNIT=101, FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_CH3COCH3_4_2.txt', STATUS='old')
  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
  CH_VOC_hyy(8,1:ntime_local) = raw_data(2:ntime_local+1, 7)
  CLOSE(UNIT=101)

  OPEN(UNIT=101, FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_CH3CHO2H_4_2.txt', STATUS='old')
  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
  CH_VOC_hyy(9,1:ntime_local) = raw_data(2:ntime_local+1, 7)
  CLOSE(UNIT=101)

  OPEN(UNIT=101, FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_H2SO4.txt', STATUS='old')
  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
  CH_H2SO4_hyy(1:ntime_local) = raw_data(2:ntime_local+1, 7)
  CLOSE(UNIT=101)
  endif
  !
  ! Input of measured uwind, vwind, temperature, absolute humidity, pressure
  !

  ! uwind
  ! 8.4, 16.8, 33.6, 67.2, 74.0
  OPEN(UNIT=101,FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_uwind.txt',status='old')
  READ(101, *) ( (raw_data(i, j), j=1, 11), i=1, ntime_local+1 )
  local_uwind(1:ntime_local, 1:5) = raw_data(2:ntime_local+1, 7:11)  ! data
  loclv_uwind(1:5) = raw_data(1, 7:11)  ! height levels
  CLOSE(UNIT=101)

  ! vwind
  ! 8.4, 16.8, 33.6, 67.2, 74.0
  OPEN(UNIT=101,FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_vwind.txt',status='old')
  READ(101, *) ( (raw_data(i, j), j=1, 11), i=1, ntime_local+1 )
  local_vwind(1:ntime_local, 1:5) = raw_data(2:ntime_local+1, 7:11)  ! data
  loclv_vwind(1:5) = raw_data(1, 7:11)  ! height levels
  CLOSE(UNIT=101)

  ! air temperature, [degC]
  ! 4.2, 8.4, 16.8, 33.6, 50.4, 67.2, 125.0
  OPEN(UNIT=101,FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_temp.txt',status='old')
  READ(101, *) ( (raw_data(i, j), j=1, 13), i=1, ntime_local+1 )
  local_temp(1:ntime_local, 1:7) = raw_data(2:ntime_local+1, 7:13) + 273.15_dp  ! data, [degC] --> [K]
  loclv_temp(1:7) = raw_data(1, 7:13)  ! height levels
  CLOSE(UNIT=101)

  ! absolute humidity, [kg m-3]
  ! 4.2, 8.4, 16.8, 33.6, 50.4, 67.2, 125.0
  OPEN(UNIT=101,FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_rhov.txt',status='old')
  READ(101, *) ( (raw_data(i, j), j=1, 13), i=1, ntime_local+1 )
  local_rhov(1:ntime_local, 1:7) = raw_data(2:ntime_local+1, 7:13)  ! data
  loclv_rhov(1:7) = raw_data(1, 7:13)  ! height levels
  CLOSE(UNIT=101)

  ! air pressure, [hPa]
  ! 0.0
  OPEN(UNIT=101,FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_pres.txt',status='old')
  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
  local_pres(1:ntime_local, 1) = raw_data(2:ntime_local+1, 7) * 1e2_dp  ! [hPa] --> [Pa]
  loclv_pres(1) = raw_data(1, 7)  ! height levels
  CLOSE(UNIT=101)

  !
  ! Input of measured soil moist, soil temperature, soil heat flux
  !

  ! soil moist, [m3 m-3]
  ! -5 - 0 cm for organic layer
  ! 2-6 cm in mineral soil
  ! 14-25 cm in mineral soil
  ! 26-36 cm in mineral soil
  ! 38-61 cm in mineral soil
  ! into the soil: -2.5, 4.0, 19.5, 31.0, 49.5
  OPEN(UNIT=101,FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_soilmoist.txt',status='old')
  READ(101, *) ( (raw_data(i, j), j=1, 11), i=1, ntime_local+1 )
  local_soilmoist(1:ntime_local, 1:5) = raw_data(2:ntime_local+1, 7:11)
  loclv_soilmoist(1:5) = raw_data(1, 7:11)
  CLOSE(UNIT=101)

  ! soil temperature, [degC]
  ! -5 - 0 cm for organic layer
  ! 2-5 cm in mineral soil
  ! 9-14 cm in mineral soil
  ! 22-29 cm in mineral soil
  ! 42-58 cm in mineral soil
  ! into the soil: -2.5, 3.5, 11.5, 25.5, 50.0
  OPEN(UNIT=101,FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_soiltemp.txt',status='old')
  READ(101, *) ( (raw_data(i, j), j=1, 11), i=1, ntime_local+1 )
  local_soiltemp(1:ntime_local, 1:5) = raw_data(2:ntime_local+1, 7:11) + 273.15_dp  ! [degC] --> [K]
  loclv_soiltemp(1:5) = raw_data(1, 7:11)
  CLOSE(UNIT=101)
  ! Calculate deep soil temperature from the average of soiltemp at 22-29 cm and 42-58 cm
  local_deep_soiltemp(1:ntime_local) = 0.5_dp * (local_soiltemp(1:ntime_local, 4) + local_soiltemp(1:ntime_local, 5))

  ! soil heat flux, [W m-2]
  ! 0.0
  OPEN(UNIT=101,FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_gsoil.txt',status='old')
  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
  local_gsoil(1:ntime_local, 1) = raw_data(2:ntime_local+1, 7)
  loclv_gsoil(1) = raw_data(1, 7)
  CLOSE(UNIT=101)

  !
  ! Input of measured global radiation, PAR, albedo
  !

  ! downward global radiation
  ! 18.0
  OPEN(UNIT=101,FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_glob.txt',status='old')
  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
  local_glob(1:ntime_local, 1) = raw_data(2:ntime_local+1, 7)
  loclv_glob(1) = raw_data(1, 7)
  CLOSE(UNIT=101)

  ! downward PAR
  ! 18.0
  OPEN(UNIT=101,FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_PAR.txt',status='old')
  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
  local_PAR(1:ntime_local, 1) = raw_data(2:ntime_local+1, 7)
  loclv_PAR(1) = raw_data(1, 7)
  CLOSE(UNIT=101)

  ! albedo
  ! 125.0
  OPEN(UNIT=101,FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_albedo.txt',status='old')
  READ(101, *) ( (raw_data(i, j), j=1, 7), i=1, ntime_local+1 )
  local_albedo(1:ntime_local, 1) = raw_data(2:ntime_local+1, 7)
  loclv_albedo(1) = raw_data(1, 7)
  CLOSE(UNIT=101)
  !Input of CS and CoagS
  !
  !Condensation Sink: sum of dmps and aps with rh-correction
  OPEN(UNIT=101,FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_cs.txt',status='old')
  READ(101,*)((raw_data(i,j),j=1,7),i=1,ntime_local+1)
  local_cs_H2SO4(1:ntime_local,1)=raw_data(2:ntime_local+1,7)
  loclv_cs_H2SO4(1)=raw_data(1,7)
  CLOSE(UNIT=101)
  !
  !Condensation Sink: sum of dmps and aps with HNO3-rh-correction
  OPEN(UNIT=101,FILE=TRIM(ADJUSTL(input_dir_station_data))//'/SMEAR_ch.txt',status='old')
  READ(101,*)((raw_data(i,j),j=1,7),i=1,ntime_local+1)
  local_cs_HNO3(1:ntime_local,1)=raw_data(2:ntime_local+1,7)
  loclv_cs_HNO3(1)=raw_data(1,7)
  CLOSE(UNIT=101)

  ! Check if there are any NANs in the local datasets
  IF ( ANY( ISNAN(CH_gas_hyy(2:6, 1:ntime_local     )) ) ) WRITE(*,*) 'NANs in CH_gas_hyy.'
  IF ( ANY( ISNAN(local_uwind(1:ntime_local, 1:5    )) ) ) WRITE(*,*) 'NANs in local_uwind.'
  IF ( ANY( ISNAN(local_vwind(1:ntime_local, 1:5    )) ) ) WRITE(*,*) 'NANs in local_vwind.'
  IF ( ANY( ISNAN(local_temp(1:ntime_local, 1:7     )) ) ) WRITE(*,*) 'NANs in local_temp.'
  IF ( ANY( ISNAN(local_rhov(1:ntime_local, 1:7     )) ) ) WRITE(*,*) 'NANs in local_rhov.'
  IF ( ANY( ISNAN(local_pres(1:ntime_local, 1       )) ) ) WRITE(*,*) 'NANs in local_pres.'
  IF ( ANY( ISNAN(local_soilmoist(1:ntime_local, 1:5)) ) ) WRITE(*,*) 'NANs in local_soilmoist.'
  IF ( ANY( ISNAN(local_soiltemp(1:ntime_local, 1:5 )) ) ) WRITE(*,*) 'NANs in local_soiltemp.'
  IF ( ANY( ISNAN(local_gsoil(1:ntime_local, 1      )) ) ) WRITE(*,*) 'NANs in local_gsoil.'
  IF ( ANY( ISNAN(local_glob(1:ntime_local, 1       )) ) ) WRITE(*,*) 'NANs in local_glob.'
  IF ( ANY( ISNAN(local_PAR(1:ntime_local, 1        )) ) ) WRITE(*,*) 'NANs in local_PAR.'
  IF ( ANY( ISNAN(local_albedo(1:ntime_local, 1     )) ) ) WRITE(*,*) 'NANs in local_albedo.'
  IF ( ANY( ISNAN(local_cs_H2SO4(1:ntime_local, 1     )) ) ) WRITE(*,*) 'NANs in local_hyy.'
  IF ( ANY( ISNAN(local_cs_HNO3(1:ntime_local, 1     )) ) ) WRITE(*,*) 'NANs in local_hyy2.'
  !
  ! Canopy information
  !
  OPEN(UNIT=101,FILE=TRIM(ADJUSTL(input_dir_station_info)) // '/canopy.txt',status='old')
  canopyfile = 0.0
  ! There are some texts after 9 columns of number. So a do-loop is needed.
  DO i = 1, 16
    READ(101, *) (canopyfile(i, j), j=1, 9)
  END DO
  CLOSE(UNIT=101)
  ! hc = canopyfile(4,1) !This is the 2012 value

  ! Canopy height is growing ~18 cm/year in Hyytiala, but I now used the
  ! measured values - check Megan_version2.f90 for more info
  IF (now_date(1) == 2018) hc = 20.04_dp
  IF (now_date(1) == 2017) hc = 19.86_dp
  IF (now_date(1) == 2016) hc = 19.68_dp  ! hc(2014) + 0.18*2
  IF (now_date(1) == 2015) hc = 19.50_dp  ! hc(2014) + 0.18*1
  IF (now_date(1) == 2014) hc = 19.32_dp
  IF (now_date(1) == 2013) hc = 18.97_dp
  IF (now_date(1) == 2012) hc = 18.63_dp
  IF (now_date(1) == 2011) hc = 18.29_dp
  IF (now_date(1) == 2010) hc = 17.98_dp
  IF (now_date(1) == 2009) hc = 17.58_dp
  IF (now_date(1) == 2008) hc = 17.16_dp
  IF (now_date(1) == 2007) hc = 16.86_dp
  IF (now_date(1) == 2006) hc = 16.57_dp
  IF (now_date(1) == 2005) hc = 16.27_dp
  IF (now_date(1) == 2004) hc = 15.97_dp
  IF (now_date(1) == 2003) hc = 15.68_dp
  IF (flag_tree == 4) THEN
    hc = 4.82
  END IF
  WRITE(*,*) 'The updated canopy height (hc) is ', hc

  !The SEP for all years are based on year 2010 measurements
  !C1: measured chamber temperature (celcius): C2: measured monoterpene emission rate (ng g-1 s-1)
  !'J' refers to the fact that this data comes from Juho
  OPEN(UNIT=101,FILE=TRIM(ADJUSTL(input_dir_station_data))//'/EFJ_mono.txt',status='old')
  READ(101, *) ( (EF_meas(i, j),i=1, 2),j=1, 1488 )
  CLOSE(UNIT=101)

  ! The concentration of methane for all years are based on year 2014 meaurements
  OPEN(UNIT=101,FILE=TRIM(ADJUSTL(input_dir_station)) // '/2014/' // month_str_Mxx // '/hyy_conc_CH4.txt', status='old')
  READ(101, *) ( hyy_ch4(j), j=1, ntime_local-1 )  ! last time point is not included
  CLOSE(UNIT=101)

  DO j = 1, ntime_local
    ! Wind: the order is from high level to low level in nudging subroutines
    ! Used height levles: 8.4, 16.8, 33.6, 74.0
    uwind(1:4, j) = local_uwind(j, (/5, 3, 2, 1/))
    vwind(1:4, j) = local_vwind(j, (/5, 3, 2, 1/))

    ! Temperature: the order is from high level to low level in nudging subroutines
    ! Used height levles: 4.2, 8.4, 16.8, 33.6, 50.4, 67.2
    tem(1:6, j) = local_temp(j, 6:1:-1)

    ! Absolute humidity: the order is from high level to low level in nudging subroutines
    ! Used height levles: 4.2, 8.4, 16.8, 33.6, 50.4, 67.2
    hum(1:6, j) = local_rhov(j, 6:1:-1)

    ! Incoming diffuse and direct radiation
    ! Current measurement is not reliable, so they are obtained from global radiation
    difftop1(j) = MAX(0.5_dp*local_glob(j, 1), 0.0_dp)
    dirtop1(j)  = MAX(0.5_dp*local_glob(j, 1), 0.0_dp)

    stat_glob(j) = MAX(local_glob(j, 1),0.0_dp)
    stat_PAR(j) = MAX(local_PAR(j, 1),0.0_dp)
    stat_albedo(j)  = MIN(MAX(local_albedo(j, 1), 0.1_dp), 0.5_dp)

    stat_sm_3(j) = 0.5_dp * ( local_soilmoist(j, 4) + local_soilmoist(j, 5) )  ! average of 26-36 cm and 38-61 cm
    stat_sm_2(j)   = 0.5_dp * ( local_soilmoist(j, 2) + local_soilmoist(j, 3) )  ! average of 2-6 cm and 14-25 cm
    stat_sm_1(j)   = local_soilmoist(j, 1)  ! Volumetric water content (m3/m3) in organic layer in -5-0 cm
    stat_gsoil(j) = local_gsoil(j, 1)

    ! Condensation sinks, commented out because it's being updated below
    ! stat_cs_H2SO4(j)     = 0.001_dp  ! hyy_mix(4,j20)
    ! stat_cs_HNO3(j)    = 0.001_dp  ! hyy_mix(24,j20)

    !new Condensation sinks
    stat_cs_H2SO4(j) = local_cs_H2SO4(j,1)
    stat_cs_HNO3(j)= local_cs_HNO3(j,1)
  END DO

  IF (flag_tree == 4) THEN !CLEARCUT
    !This is used for the clear cut and includes a modification of the
    !measured data from Haapanala et al., BGS, 2012.
    !C1: daily monoterpene flux assuming constant flux from day
    !116-177. Unit: ug m^-2 h^-1
    !C2: daily monoterpene flux assuming decreasing flux from day
    !116-177 with a start value of around 7400. Unit: ug m^-2 h^-1
    !C3:  daily monoterpene flux assuming decreasing flux from day
    !116-177 with a start value of around 10000. Unit: ug m^-2 h^-1
    OPEN(2715,FILE=TRIM(ADJUSTL(input_dir_station_info)) // '/SAMI_IN.txt')
    DO J = 1,366
       READ(2715,*) (SAMI(J,I), I=1,3)
    ENDDO
    CLOSE(2715)
  ENDIF

  IF (flag_tree == 3) THEN !BIRCH
      open(2716,file=TRIM(ADJUSTL(input_dir_station_info)) // '/lai_function.txt')
      do J = 1,366
         read(2716,*) LAI_function(J)
      enddo
      close(2716)
  ENDIF

  abl = 2

  IF (flag_ecmwf == 1) THEN
    ! time: 0, 3, 6, 9, 12, 15, 18, 21 every day, with 24:00 in last day
    ! space: 6 date numbers + 37 pressure layers, or 6 date numbers + 1 surface value
    ntime_ECMWF = MonthDay(now_date(1), now_date(2))*8 + 1

    OPEN(unit=102,file = TRIM(ADJUSTL(input_dir_station_data))//'/ECMWF_uwind.txt',status='old')
    READ(102,*) ( (input_uwind(i, j), j=1, 43), i=1, ntime_ECMWF )
    CLOSE(102)

    OPEN(unit=102,file = TRIM(ADJUSTL(input_dir_station_data))//'/ECMWF_vwind.txt',status='old')
    READ(102,*) ( (input_vwind(i, j), j=1, 43),i=1, ntime_ECMWF )
    CLOSE(102)

    OPEN(unit=102,file = TRIM(ADJUSTL(input_dir_station_data))//'/ECMWF_temp.txt',status='old')
    READ(102,*) ( (input_temp(i, j), j=1, 43), i=1, ntime_ECMWF )
    CLOSE(102)

    OPEN(unit=102,file = TRIM(ADJUSTL(input_dir_station_data))//'/ECMWF_rhov.txt',status='old')
    READ(102,*) ( (input_rhov(i, j), j=1, 43), i=1, ntime_ECMWF )
    CLOSE(102)

    OPEN(UNIT=102,FILE = TRIM(ADJUSTL(input_dir_station_data))//'/ECMWF_geop.txt',status='old')
    READ(102,*) ( (input_geop(i, j), j=1, 43), i=1, ntime_ECMWF)
    CLOSE(102)
    input_geom = input_geop / grav  ! geopotential height to geometric height

    OPEN(UNIT=102,FILE = TRIM(ADJUSTL(input_dir_station_data))//'/ECMWF_strd.txt',status='old')
    READ(102,*) ( (input_strd(i, j), j=1, 7), i=1, ntime_ECMWF )
    CLOSE(UNIT=102)

    ! Check if there are any NANs in the input ECMWF datasets.
    IF ( ANY(ISNAN(input_geom )) ) WRITE(*,*) 'NANs in the ECMWF geometric heights.'
    IF ( ANY(ISNAN(input_uwind)) ) WRITE(*,*) 'NANs in the ECMWF u wind.'
    IF ( ANY(ISNAN(input_vwind)) ) WRITE(*,*) 'NANs in the ECMWF v wind.'
    IF ( ANY(ISNAN(input_temp )) ) WRITE(*,*) 'NANs in the ECMWF temperature.'
    IF ( ANY(ISNAN(input_qv   )) ) WRITE(*,*) 'NANs in the ECMWF specific humidity.'
    IF ( ANY(ISNAN(input_rhov )) ) WRITE(*,*) 'NANs in the ECMWF absolute humidity.'

  END IF  ! flag_ecmwf == 1

end subroutine read_input_data_hyytiala_ascii


subroutine read_input_data_hyytiala_netcdf()
  integer :: fid, stat
  character(len=300) :: fname_input
  integer :: dimid_lev, dimid_time, dimid_timeval

  ! Read input netcdf file
  fname_input = trim(adjustl(input_dir_station))//'/input.nc'

  call dios_open( fname_input, DIOS_NETCDF, DIOS_READ, fid, stat )

  call dios_inq_dimid(fid, 'lev'    , dimid_lev    , stat)
  call dios_inq_dimid(fid, 'time'   , dimid_time   , stat)
  call dios_inq_dimid(fid, 'timeval', dimid_timeval, stat)

  ! call dios_inq_varid(fid, 'SO2',

  ! call dios_show( fname_input, stat )

  call MPI_Finalize(mpi_rc)
  stop
end subroutine read_input_data_hyytiala_netcdf


subroutine read_input_data_traj_netcdf()
  integer :: stat

  character(len=300) :: fabspath_meteo
  integer :: fid_meteo
  integer :: dimid_meteo_time, dimid_meteo_lev

  character(len=300) :: fabspath_conc
  integer :: fid_conc
  integer :: dimid_conc_time, dimid_conc_lev

  character(len=300) :: fabspath_emisant
  integer :: fid_emisant
  integer :: dimid_emisant_time, dimid_emisant_layer

  character(len=300) :: fabspath_emisbio
  integer :: fid_emisbio
  integer :: dimid_emisbio_time

  character(len=300) :: fabspath_emisaer
  integer :: fid_emisaer
  integer :: dimid_emisaer_time

  character(len=8)   :: traj_datestring
  character(len=2)   :: release

  ! ?? Molar masses of species [kg mol-1], used for converting emission units,
  ! this will be set in a general way for all the chemical species
  real(dp), parameter :: xmm_CO       = 28.01d-3
  real(dp), parameter :: xmm_NOx      = 46.0055d-3  ! ?? assume NO2 now
  real(dp), parameter :: xmm_NH3      = 17.031d-3
  real(dp), parameter :: xmm_SO2      = 64.066d-3
  real(dp), parameter :: xmm_C5H8     = 68.12d-3
  real(dp), parameter :: xmm_MT       = 136.23d-3
  real(dp), parameter :: xmm_CH3CHO   = 44.05d-3
  real(dp), parameter :: xmm_HCHO     = 30.031d-3
  real(dp), parameter :: xmm_CH3COCH3 = 58.08d-3
  real(dp), parameter :: xmm_APINENE  = 136.23d-3
  real(dp), parameter :: xmm_BPINENE  = 136.23d-3
  real(dp), parameter :: xmm_CH3OH    = 32.04d-3
  real(dp), parameter :: xmm_MBO      = 86.1323d-3
  real(dp), parameter :: xmm_OMT      = 136.23d-3
  real(dp), parameter :: xmm_SQT      = 246.3016d-3
  real(dp), parameter :: xmm_DMS      = 62.13d-3

  ! real(dp), allocatable :: emisant_layer_thickness(:)
  integer, allocatable :: tmp_int_1d(:)
  real(dp), allocatable :: tmp_real_1d(:)

  write(traj_datestring,'(i0.4,i0.2,i0.2)') end_date(1),end_date(2),end_date(3)
  write(release,'(i0.2)') 24-end_date(4)

  !--------------------------------------------------------!
  ! Set input netcdf file paths
  !--------------------------------------------------------!
  ! fabspath_conc   = trim(adjustl(INPUT_DIR)) // &
  !   '/OUTPUT_bwd_'//traj_datestring//'/CONC/CONC_'//traj_datestring//'_R'//release//'.nc'
  fabspath_conc   = '/scratch/project_2004394/PETRI/FP_OUT/HYDE_BASE_Y2018' // &
    '/OUTPUT_bwd_'//traj_datestring//'/CONC/CONC_'//traj_datestring//'_R'//release//'.nc'
  fabspath_meteo   = trim(adjustl(INPUT_DIR)) // &
    '/OUTPUT_bwd_'//traj_datestring//'/METEO/METEO_'//traj_datestring//'_R'//release//'.nc'
  fabspath_emisant = trim(adjustl(INPUT_DIR)) // &
    '/OUTPUT_bwd_'//traj_datestring//'/EMISSIONS_0422/'//traj_datestring//'_7daybwd_Hyde_traj_ANT_'//release//'_L3.nc'
  fabspath_emisbio = trim(adjustl(INPUT_DIR)) // &
    '/OUTPUT_bwd_'//traj_datestring//'/EMISSIONS_0422/'//traj_datestring//'_7daybwd_Hyde_traj_BIO_'//release//'_L3.nc'
  fabspath_emisaer = trim(adjustl(INPUT_DIR)) // &
    '/OUTPUT_bwd_'//traj_datestring//'/EMISSIONS_0422/'//traj_datestring//'_7daybwd_Hyde_traj_AER_'//release//'_L3.nc'

  if (my_id == master_id) print*, 'Path to concentr : ', fabspath_conc
  if (my_id == master_id) print*, 'Path to meteo    : ', fabspath_meteo
  if (my_id == master_id) print*, 'Path to emi ant  : ', fabspath_emisant
  if (my_id == master_id) print*, 'Path to emi bio  : ', fabspath_emisbio
  if (my_id == master_id) print*, 'Path to emi aer  : ', fabspath_emisaer

  !--------------------------------------------------------!
  ! Read concentration file
  !--------------------------------------------------------!
  call dios_open(fabspath_conc, DIOS_NETCDF, DIOS_READ, fid_conc, stat)

  ! Get dimension ids and lengths
  call dios_inq_dimid(fid_conc, 'time', dimid_conc_time, stat)
  call dios_inq_dimid(fid_conc, 'lev' , dimid_conc_lev , stat)

  call dios_inquire_dimension(fid_conc, dimid_conc_time, stat, &
    length=infield%ntime)
  call dios_inquire_dimension(fid_conc, dimid_conc_lev, stat, &
    length=infield%nlevel_conc)


  !--------------------------------------------------------!
  ! Read meteorology file
  !--------------------------------------------------------!
  call dios_open(fabspath_meteo, DIOS_NETCDF, DIOS_READ, fid_meteo, stat)

  ! Get dimension ids and lengths
  call dios_inq_dimid(fid_meteo, 'time', dimid_meteo_time, stat)
  call dios_inq_dimid(fid_meteo, 'lev' , dimid_meteo_lev , stat)

  call dios_inquire_dimension(fid_meteo, dimid_meteo_time, stat, &
    length=infield%ntime)
  call dios_inquire_dimension(fid_meteo, dimid_meteo_lev, stat, &
    length=infield%nlevel_meteo)


  !--------------------------------------------------------!
  ! Read anthropogenic emission file
  !--------------------------------------------------------!
  call dios_open(fabspath_emisant, DIOS_NETCDF, DIOS_READ, fid_emisant, stat)

  ! Get dimension ids and lengths
  call dios_inq_dimid(fid_emisant, 'time' , dimid_emisant_time , stat)
  call dios_inquire_dimension(fid_emisant, dimid_emisant_time , stat, &
    length=infield%ntime)

  call dios_inq_dimid(fid_emisant, 'layer', dimid_emisant_layer, stat)
  call dios_inquire_dimension(fid_emisant, dimid_emisant_layer, stat, &
    length=infield%nlevel_emisant)

  write(*,*) 'nlevel, ntime: ', infield%nlevel_emisant, infield%ntime


  !--------------------------------------------------------!
  ! Read aerosol emission file
  !--------------------------------------------------------!
  call dios_open(fabspath_emisaer, DIOS_NETCDF, DIOS_READ, fid_emisaer, stat)

  ! Get dimension ids and lengths
  call dios_inq_dimid(fid_emisaer, 'time' , dimid_emisaer_time , stat)
  call dios_inquire_dimension(fid_emisaer, dimid_emisaer_time , stat, &
    length=infield%ntime)

  call dios_inq_dimid(fid_emisaer, 'layer', dimid_emisant_layer, stat)
  call dios_inquire_dimension(fid_emisaer, dimid_emisant_layer, stat, &
    length=infield%nlevel_emisaer)

  write(*,*) 'nlevel, ntime: ', infield%nlevel_emisaer, infield%ntime


  !--------------------------------------------------------!
  ! Read biogenic emission file
  !--------------------------------------------------------!

  call dios_open(fabspath_emisbio, DIOS_NETCDF, DIOS_READ, fid_emisbio, stat)

  ! Get dimension ids and lengths
  ! call dios_inq_dimid(fid_emisbio, 'time' , dimid_emisbio_time , stat)
  ! call dios_inquire_dimension(fid_emisbio, dimid_emisbio_time , stat, &
  !   length=infield%ntime)


  !------------------------------------------------------!
  ! Put input vairables to infield
  !------------------------------------------------------!

  !===== Set unique ids for variables and count the total number of them =====!
  call add_input_field_id(infield%id_conc_time       , infield%nvar)
  call add_input_field_id(infield%id_conc_o3         , infield%nvar)
  call add_input_field_id(infield%id_conc_lev        , infield%nvar)

  call add_input_field_id(infield%id_t               , infield%nvar)
  call add_input_field_id(infield%id_u               , infield%nvar)
  call add_input_field_id(infield%id_v               , infield%nvar)
  call add_input_field_id(infield%id_q               , infield%nvar)
  call add_input_field_id(infield%id_mla             , infield%nvar)
  call add_input_field_id(infield%id_lp              , infield%nvar)
  call add_input_field_id(infield%id_ssr             , infield%nvar)
  call add_input_field_id(infield%id_lsm             , infield%nvar)
  call add_input_field_id(infield%id_time            , infield%nvar)
  call add_input_field_id(infield%id_lat             , infield%nvar)
  call add_input_field_id(infield%id_lon             , infield%nvar)

  call add_input_field_id(infield%id_emisant_time    , infield%nvar)
  call add_input_field_id(infield%id_emisant_blh     , infield%nvar)
  call add_input_field_id(infield%id_emisant_mlh     , infield%nvar)
  call add_input_field_id(infield%id_emisant_tlh     , infield%nvar)
  call add_input_field_id(infield%id_emisant_CO      , infield%nvar)
  call add_input_field_id(infield%id_emisant_NOx     , infield%nvar)
  call add_input_field_id(infield%id_emisant_NH3     , infield%nvar)
  call add_input_field_id(infield%id_emisant_SO2     , infield%nvar)
  call add_input_field_id(infield%id_emisant_C5H8    , infield%nvar)
  call add_input_field_id(infield%id_emisant_MT      , infield%nvar)

  call add_input_field_id(infield%id_emisbio_time    , infield%nvar)
  call add_input_field_id(infield%id_emisbio_SRRsum  , infield%nvar)
  call add_input_field_id(infield%id_emisbio_CH3CHO  , infield%nvar)
  call add_input_field_id(infield%id_emisbio_HCHO    , infield%nvar)
  call add_input_field_id(infield%id_emisbio_CH3COCH3, infield%nvar)
  call add_input_field_id(infield%id_emisbio_C5H8    , infield%nvar)
  call add_input_field_id(infield%id_emisbio_APINENE , infield%nvar)
  call add_input_field_id(infield%id_emisbio_BPINENE , infield%nvar)
  call add_input_field_id(infield%id_emisbio_CH3OH   , infield%nvar)
  call add_input_field_id(infield%id_emisbio_MBO     , infield%nvar)
  call add_input_field_id(infield%id_emisbio_OMT     , infield%nvar)
  call add_input_field_id(infield%id_emisbio_SQT     , infield%nvar)
  call add_input_field_id(infield%id_emisbio_DMS     , infield%nvar)

  call add_input_field_id(infield%id_emisaer_time       , infield%nvar)
  call add_input_field_id(infield%id_emisaer_3_10nm     , infield%nvar)
  call add_input_field_id(infield%id_emisaer_10_20nm    , infield%nvar)
  call add_input_field_id(infield%id_emisaer_20_30nm    , infield%nvar)
  call add_input_field_id(infield%id_emisaer_30_50nm    , infield%nvar)
  call add_input_field_id(infield%id_emisaer_50_70nm    , infield%nvar)
  call add_input_field_id(infield%id_emisaer_70_100nm   , infield%nvar)
  call add_input_field_id(infield%id_emisaer_100_200nm  , infield%nvar)
  call add_input_field_id(infield%id_emisaer_200_400nm  , infield%nvar)
  call add_input_field_id(infield%id_emisaer_400_1000nm , infield%nvar)

  ! Allocate input field variables according to number of variables
  allocate(infield%var(infield%nvar))

  !===== Read the data and save them to input field =====!

  !----- concentrations

  call add_input_field( infield, infield%id_conc_time, 'time', &
    fabspath_conc, fid_conc, &
    (/infield%ntime/) )

  call add_input_field( infield, infield%id_conc_o3, 'go3', &
    fabspath_conc, fid_conc, &
    (/infield%nlevel_conc, infield%ntime/) )

  call add_input_field( infield, infield%id_conc_lev, 'lev', &
    fabspath_conc, fid_conc, &
    (/infield%nlevel_conc/) )


  !----- Meteorological variables

  call add_input_field( infield, infield%id_time, 'time', &
    fabspath_meteo, fid_meteo, &
    (/infield%ntime/) )

  call add_input_field( infield, infield%id_t, 't', &
    fabspath_meteo, fid_meteo, &
    (/infield%nlevel_meteo, infield%ntime/) )

  call add_input_field( infield, infield%id_u, 'u', &
    fabspath_meteo, fid_meteo, &
    (/infield%nlevel_meteo, infield%ntime/) )

  call add_input_field( infield, infield%id_v, 'v', &
    fabspath_meteo, fid_meteo, &
    (/infield%nlevel_meteo, infield%ntime/) )

  call add_input_field( infield, infield%id_q, 'q', &
    fabspath_meteo, fid_meteo, &
    (/infield%nlevel_meteo, infield%ntime/) )

  call add_input_field( infield, infield%id_mla, 'mla', &
    fabspath_meteo, fid_meteo, &
    (/infield%nlevel_meteo, infield%ntime/) )

  call add_input_field( infield, infield%id_lp, 'lp', &
    fabspath_meteo, fid_meteo, &
    (/infield%nlevel_meteo, infield%ntime/) )

  call add_input_field( infield, infield%id_ssr, 'ssr', &
    fabspath_meteo, fid_meteo, &
    (/infield%ntime/) )

  call add_input_field( infield, infield%id_lsm, 'lsm', &
    fabspath_meteo, fid_meteo, &
    (/infield%ntime/) )

  call add_input_field( infield, infield%id_lat, 'lat', &
    fabspath_meteo, fid_meteo, &
    (/infield%ntime/) )

  call add_input_field( infield, infield%id_lon, 'lon', &
    fabspath_meteo, fid_meteo, &
    (/infield%ntime/) )

  !----- Anthropogenic emission variables
  ! 2d variables are saved as (layer, time) in the netcdf file, so they should
  ! be read as (time, layer), need to modify the input files later.

  call add_input_field( infield, infield%id_emisant_time, 'time', &
    fabspath_emisant, fid_emisant, &
    (/infield%ntime/) )

  call add_input_field( infield, infield%id_emisant_blh, 'bottom_layer_height', &
    fabspath_emisant, fid_emisant, &
    (/infield%nlevel_emisant/) )

  call add_input_field( infield, infield%id_emisant_mlh, 'mid_layer_height', &
    fabspath_emisant, fid_emisant, &
    (/infield%nlevel_emisant/) )

  call add_input_field( infield, infield%id_emisant_tlh, 'top_layer_height', &
    fabspath_emisant, fid_emisant, &
    (/infield%nlevel_emisant/) )

  call add_input_field( infield, infield%id_emisant_CO, 'co', &
    fabspath_emisant, fid_emisant, &
    (/infield%ntime, infield%nlevel_emisant/) , infield%id_emisant_time)
  ! write(*,*) 'co', infield%var(infield%id_emisant_co)%f2d

  call add_input_field( infield, infield%id_emisant_NOx, 'nox', &
    fabspath_emisant, fid_emisant, &
    (/infield%ntime, infield%nlevel_emisant/) , infield%id_emisant_time)

  call add_input_field( infield, infield%id_emisant_NH3, 'nh3', &
    fabspath_emisant, fid_emisant, &
    (/infield%ntime, infield%nlevel_emisant/) , infield%id_emisant_time)

  call add_input_field( infield, infield%id_emisant_SO2, 'so2', &
    fabspath_emisant, fid_emisant, &
    (/infield%ntime, infield%nlevel_emisant/) , infield%id_emisant_time)

  call add_input_field( infield, infield%id_emisant_C5H8, 'isoprene', &
    fabspath_emisant, fid_emisant, &
    (/infield%ntime, infield%nlevel_emisant/) , infield%id_emisant_time)

  call add_input_field( infield, infield%id_emisant_MT, 'monoterpenes', &
    fabspath_emisant, fid_emisant, &
    (/infield%ntime, infield%nlevel_emisant/) , infield%id_emisant_time)

!----- Aerosol emission variables

  ! call add_input_field( infield, infield%id_emisaer_lev, 'lev', &
  !   fabspath_emisaer, fid_emisaer, &
  !   (/infield%nlevel_emisant/) )

  call add_input_field( infield, infield%id_emisaer_time, 'time', &
    fabspath_emisaer, fid_emisaer, &
    (/infield%ntime/) )

  call add_input_field( infield, infield%id_emisaer_3_10nm, '3-10nm', &
    fabspath_emisaer, fid_emisaer, &
    (/infield%ntime, infield%nlevel_emisant/) , infield%id_emisaer_time)

  call add_input_field( infield, infield%id_emisaer_10_20nm, '10-20nm', &
    fabspath_emisaer, fid_emisaer, &
    (/infield%ntime, infield%nlevel_emisant/) , infield%id_emisaer_time)

  call add_input_field( infield, infield%id_emisaer_20_30nm, '20-30nm', &
    fabspath_emisaer, fid_emisaer, &
    (/infield%ntime, infield%nlevel_emisant/) , infield%id_emisaer_time)

  call add_input_field( infield, infield%id_emisaer_30_50nm, '30-50nm', &
    fabspath_emisaer, fid_emisaer, &
    (/infield%ntime, infield%nlevel_emisant/) , infield%id_emisaer_time)

  call add_input_field( infield, infield%id_emisaer_50_70nm, '50-70nm', &
    fabspath_emisaer, fid_emisaer, &
    (/infield%ntime, infield%nlevel_emisant/) , infield%id_emisaer_time)

  call add_input_field( infield, infield%id_emisaer_70_100nm, '70-100nm', &
    fabspath_emisaer, fid_emisaer, &
    (/infield%ntime, infield%nlevel_emisant/) , infield%id_emisaer_time)

  call add_input_field( infield, infield%id_emisaer_100_200nm, '100-200nm', &
    fabspath_emisaer, fid_emisaer, &
    (/infield%ntime, infield%nlevel_emisant/) , infield%id_emisaer_time)

  call add_input_field( infield, infield%id_emisaer_200_400nm, '200-400nm', &
    fabspath_emisaer, fid_emisaer, &
    (/infield%ntime, infield%nlevel_emisant/) , infield%id_emisaer_time)

  call add_input_field( infield, infield%id_emisaer_400_1000nm, '400-1000nm', &
    fabspath_emisaer, fid_emisaer, &
    (/infield%ntime, infield%nlevel_emisant/) , infield%id_emisaer_time)



  !----- Biogenic emission variables

  call add_input_field( infield, infield%id_emisbio_time, 'time', &
    fabspath_emisbio, fid_emisbio, &
    (/infield%ntime/) )

  call add_input_field( infield, infield%id_emisbio_SRRsum, 'SRRsum', &
    fabspath_emisbio, fid_emisbio, &
    (/infield%ntime/) , infield%id_emisbio_time)

  call add_input_field( infield, infield%id_emisbio_CH3CHO, 'acetaldehyde', &
    fabspath_emisbio, fid_emisbio, &
    (/infield%ntime/) , infield%id_emisbio_time)

  call add_input_field( infield, infield%id_emisbio_HCHO, 'formaldehyde', &
    fabspath_emisbio, fid_emisbio, &
    (/infield%ntime/) , infield%id_emisbio_time)

  call add_input_field( infield, infield%id_emisbio_CH3COCH3, 'acetone', &
    fabspath_emisbio, fid_emisbio, &
    (/infield%ntime/) , infield%id_emisbio_time)

  call add_input_field( infield, infield%id_emisbio_C5H8, 'isoprene', &
    fabspath_emisbio, fid_emisbio, &
    (/infield%ntime/) , infield%id_emisbio_time)

  call add_input_field( infield, infield%id_emisbio_APINENE, 'pinene-a', &
    fabspath_emisbio, fid_emisbio, &
    (/infield%ntime/) , infield%id_emisbio_time)

  call add_input_field( infield, infield%id_emisbio_BPINENE, 'pinene-b', &
    fabspath_emisbio, fid_emisbio, &
    (/infield%ntime/) , infield%id_emisbio_time)

  call add_input_field( infield, infield%id_emisbio_CH3OH, 'methanol', &
    fabspath_emisbio, fid_emisbio, &
    (/infield%ntime/) , infield%id_emisbio_time)

  call add_input_field( infield, infield%id_emisbio_MBO, 'MBO', &
    fabspath_emisbio, fid_emisbio, &
    (/infield%ntime/) , infield%id_emisbio_time)

  call add_input_field( infield, infield%id_emisbio_OMT, 'other-monoterpenes', &
    fabspath_emisbio, fid_emisbio, &
    (/infield%ntime/) , infield%id_emisbio_time)

  call add_input_field( infield, infield%id_emisbio_SQT, 'sesquiterpenes', &
    fabspath_emisbio, fid_emisbio, &
    (/infield%ntime/) , infield%id_emisbio_time)

  call add_input_field( infield, infield%id_emisbio_DMS, 'DMS', &
    fabspath_emisbio, fid_emisbio, &
    (/infield%ntime/) , infield%id_emisbio_time)


  !------------------------------------------------------!
  ! Set other variables
  !------------------------------------------------------!

  ! ?? temporary solution, need to be improved later
  hc = 10.0d0  ! [m], canopy height, this should be read from input file
  dirtop1 = 0d0!200.0d0
  difftop1 = 0d0! 200.0d0
  glob = 10.0d0
  ! PAR = 500.0d0
  albedo = 0.2d0

  ! Current input levels for meteorology data should be obtained from interpolation
  ! every time step. But for anthropogenic emissions the level heights do not
  ! change with time, so there is no need to save it to an additional array.
  allocate(input_levels_curr(infield%nlevel_meteo))

  ! The traj time starts from negative values, so convert them to time elapsed
  ! instead of time backwards.
  if (infield%var(infield%id_time)%f1d(1)>infield%var(infield%id_time)%f1d(infield%ntime)) &
    STOP 'TIME VECTOR IS REVERSED'

  IF (weigh_with_srr) print*, 'NOTE! ----------- Weighing emissions with SRR'

  infield%var(infield%id_time)%f1d = infield%var(infield%id_time)%f1d &
                                     - infield%var(infield%id_time)%f1d(1)

  ! Convert the unit of anthropogenic emission rates from [kg m-2 s-1] to
  ! [molec cm-3 s-1]:
  ! nconc = (mconc / xmm) * NA / layer_thickness
  ! The biogenic emission rates also have the unit of [kg m-2 s-1], but there is
  ! no information about layer thickness in the file, so I assume it is the same
  ! as the lowest layer in anthropogenic emissions. NOTE waat? P.C.
  allocate(tmp_int_1d(6), tmp_real_1d(6))

  tmp_int_1d = (/ &
    infield%id_emisant_CO      , &
    infield%id_emisant_NOx     , &
    infield%id_emisant_NH3     , &
    infield%id_emisant_SO2     , &
    infield%id_emisant_C5H8    , &
    infield%id_emisant_MT        &
    /)
  tmp_real_1d = (/ xmm_CO, xmm_NOx, xmm_NH3, xmm_SO2, xmm_C5H8, xmm_MT /)

  do i = 1, infield%ntime
    do j = 1, 6
      infield%var(tmp_int_1d(j))%f2d(i, :) = &
        infield%var(tmp_int_1d(j))%f2d(i, :) &
        / tmp_real_1d(j) * Avog * 1d-6 / & ! 1d-6 is unit conversion from m3 to cm3
        ( infield%var(infield%id_emisant_tlh)%f1d -  &
        infield%var(infield%id_emisant_blh)%f1d )

        IF (weigh_with_srr) infield%var(tmp_int_1d(j))%f2d(i, :) = &
            infield%var(tmp_int_1d(j))%f2d(i, :) * ( infield%var(infield%id_emisbio_SRRsum)%f1d(i) / 3600d0 ) ! normalisation factor in SRR
    end do
  end do

  deallocate(tmp_int_1d, tmp_real_1d)

  allocate(tmp_int_1d(11), tmp_real_1d(11))

  tmp_int_1d = (/ &
    infield%id_emisbio_CH3CHO  , &
    infield%id_emisbio_HCHO    , &
    infield%id_emisbio_CH3COCH3, &
    infield%id_emisbio_C5H8    , &
    infield%id_emisbio_APINENE , &
    infield%id_emisbio_BPINENE , &
    infield%id_emisbio_CH3OH   , &
    infield%id_emisbio_MBO     , &
    infield%id_emisbio_OMT     , &
    infield%id_emisbio_SQT     , &
    infield%id_emisbio_DMS       &
    /)
  tmp_real_1d = (/ &
    xmm_CH3CHO, xmm_HCHO, xmm_CH3COCH3, xmm_C5H8, xmm_APINENE, &
    xmm_BPINENE, xmm_CH3OH, xmm_MBO, xmm_OMT, xmm_SQT, &
    xmm_DMS &
    /)

  ! NOTE .
  do i = 1, infield%ntime
    do j = 1, size(tmp_int_1d,1)
    infield%var(tmp_int_1d(j))%f1d(i) = &
      ! Line below coverts from kg/m2/s -> molecules/m2/s -> molecules/cm2/s
      infield%var(tmp_int_1d(j))%f1d(i) / tmp_real_1d(j) * Avog * 1d-6
      ! NOTE biogenic emissions are later divided by bio_dz

      IF (weigh_with_srr) infield%var(tmp_int_1d(j))%f1d(i) = &
        infield%var(tmp_int_1d(j))%f1d(i) * ( infield%var(infield%id_emisbio_SRRsum)%f1d(i) / 3600d0 ) ! normalisation factor in SRR
      ! / ( infield%var(infield%id_emisant_tlh)%f1d(1) - infield%var(infield%id_emisant_blh)%f1d(1) )
    end do
  end do

  deallocate(tmp_int_1d, tmp_real_1d)


! XXX NOTE normalize particle emissions here
  allocate(tmp_int_1d(9), tmp_real_1d(9))

  tmp_int_1d = (/ &
    infield%id_emisaer_3_10nm    , &
    infield%id_emisaer_10_20nm   , &
    infield%id_emisaer_20_30nm   , &
    infield%id_emisaer_30_50nm   , &
    infield%id_emisaer_50_70nm   , &
    infield%id_emisaer_70_100nm  , &
    infield%id_emisaer_100_200nm , &
    infield%id_emisaer_200_400nm , &
    infield%id_emisaer_400_1000nm  &
    /)
  tmp_real_1d = (/ &
    log10(10d0/3d0), &
    log10(20d0/10d0), &
    log10(30d0/20d0), &
    log10(50d0/30d0), &
    log10(70d0/50d0), &
    log10(100d0/70d0), &
    log10(200d0/100d0), &
    log10(400d0/200d0), &
    log10(1000d0/400d0)  &
    /)

  do i = 1, infield%ntime
    do j = 1, 9
      infield%var(tmp_int_1d(j))%f2d(i, :) = &
        infield%var(tmp_int_1d(j))%f2d(i, :) &

        ! normalize with Delta(log(dp)) and convert to particle numbers / m² (emissions are divided by 10^21)
        / tmp_real_1d(j) * 1D21 &

        ! change from m2 to m3
        / ( infield%var(infield%id_emisant_tlh)%f1d -  &
        infield%var(infield%id_emisant_blh)%f1d )

        IF (weigh_with_srr) infield%var(tmp_int_1d(j))%f2d(i, :) = &
            infield%var(tmp_int_1d(j))%f2d(i, :) * ( infield%var(infield%id_emisbio_SRRsum)%f1d(i) / 3600d0 ) ! normalisation factor in SRR
    end do
  end do

  deallocate(tmp_int_1d, tmp_real_1d)

  ! write(*,*) 'ta: ', infield%var(infield%id_t)%f2d(:, 1)
  ! write(*,*) 'ta: ', infield%var(infield%id_t)%f2d(1, :)
  ! stop

end subroutine read_input_data_traj_netcdf


! ?? Put the emission indices to several arrays in infield variable
! subroutine update_nconc_from_emission_traj(nconc, )
! end subroutine update_nconc_from_emission_traj

!
!===============================================================================================================================
!

SUBROUTINE READ_DATA_OLD()
  IF (TRIM(ADJUSTL(STATION)) == 'hyytiala') THEN

    ! Input of mixed data from Hyytiala (half-hour data):
    ! 1  = date and time
    ! 2  = diffuse global radiation [W/m2]
    ! 3  = direct global radiation [W/m2]       !r for Welgegund this the TOTAL global radiation
    ! 4  = cs - sum of dmps and aps with rh-correction
    ! 5  = u* in m/s
    ! 6  = temperature in K at 67.2 m
    ! 7  = temperature in K  at 50.4 m
    ! 8  = temperature in K  at 33.6 m
    ! 9  = temperature in K  at 16.8 m
    ! 10 = temperature in K  at 8.4 m
    ! 11 = temperature in K  at 4.2 m
    ! 12 = absolute humidity in kg/m3 at 67.2 m
    ! 13 = absolute humidity in kg/m3 at 50.4 m
    ! 14 = absolute humidity in kg/m3 at 33.6 m
    ! 15 = absolute humidity in kg/m3 at 16.8 m
    ! 16 = absolute humidity in kg/m3 at 8.4 m
    ! 17 = absolute humidity in kg/m3 at 4.2 m
    ! 18 = ws in m/s at 74 m
    ! 19 = ws in m/s  at 33.6 m
    ! 20 = ws in m/s  at 16.8 m
    ! 21 = ws in m/s  at 8.4 m
    ! 22 = wd at 33.6 m
    ! 23 = wd at 16.8 m
    ! 24 = cs - sum of dmps and aps with rh-correction for HNO3
    ! 25 = global radiation
    ! 26 = solar par radiation
    ! 27 = coa for 1 nm - sum of dmps and aps with rh-correction for HNO3
    ! 28 = albedo
    ! 29-32 wind data from sonic anemometers
    ntime_local = MonthDay(now_date(1), now_date(2))*48 + 1  ! yyyymm01-00:00:00 to yyyy(mm+1)01-00:00:00, so 1 is added

    ! Input of measured spectral radiation data (every half-hour data), not used now due to measurement problems
    ! OPEN(unit=4,file = TRIM(ADJUSTL(input_dir_station_data))//'/hyy_swr.txt',status='old')
    ! READ(4,*)((swr_hyy(i20,j20),i20=1,47),j20=1,ntime_local)
    ! CLOSE(unit=4)

    ! Input of measured gases O3, SO2, NO, NO2, CO and air (every half-hour data)
    OPEN(unit=4,file = TRIM(ADJUSTL(input_dir_station_data))//'/hyy_gas_new.txt',status='old')
    READ(4,*)((CH_gas_hyy(i20,j20),i20=1,12),j20=1,ntime_local)
    CLOSE(unit=4)

    if (now_date(1)>=2013) then
       OPEN(unit=2,file = TRIM(ADJUSTL(input_dir_station_data))//'/hyy_mix.txt',status='old')
       READ(2,*)((hyy_mix(i20,j20),i20=1,32),j20=1,1488)
       CLOSE(unit=2)
    else

       OPEN(unit=2,file = TRIM(ADJUSTL(input_dir_station_data))//'/hyy_mix.txt',status='old')
       READ(2,*)((hyy_mix(i20,j20),i20=1,28),j20=1,ntime_local)
       CLOSE(unit=2)

    endif

    !r content of the hyy_soil file:
    !r 1  date and time
    !r 2  Volumetric water content (m3/m3) in organic layer in -5-0 cm
    !r 3  Volumetric water content (m3/m3) in 0-5 cm
    !r 4  Volumetric water content (m3/m3) in 5-23 cm
    !r 5  Volumetric water content (m3/m3) in 23-60 cm
    !r 6  Temperature (C) in organic layer in -5-0 cm
    !r 7  Temperature (C) in 0-5 cm
    !r 8  Temperature (C) in 5-23 cm
    !r 9  Temperature (C) in 23-60 cm

    OPEN(unit=5,file = TRIM(ADJUSTL(input_dir_station_data))//'/hyy_soil.txt',status='old')
    READ(5,*)((hyy_soil(i20,j20),i20=1,9),j20=1,ntime_local)
    CLOSE(unit=5)

    open(unit=2037,file=TRIM(ADJUSTL(input_dir_station_info)) // '/canopy.txt', status='old')
    canopyfile = 0.0
    DO j20 = 1,16
       READ(2037,*)(canopyfile(j20,i20),i20=1,9)
    ENDDO
    CLOSE(unit=2037)
    !hc = canopyfile(4,1) !This is the 2012 value

    !Canopy height is growing ~18 cm/year in HyytiÃlÃ, but I now used the
    !measured values - check Megan_version2.f90 for more info
    IF (year_str == '2018') hc = 20.61
    IF (year_str == '2017') hc = 20.36
    IF (year_str == '2016') hc = 20.01
    IF (year_str == '2015') hc = 19.66
    IF (year_str == '2014') hc = 19.32
    IF (year_str == '2013') hc = 18.97
    IF (year_str == '2012') hc = 18.63
    IF (year_str == '2011') hc = 18.29
    IF (year_str == '2010') hc = 17.98
    IF (year_str == '2009') hc = 17.58
    IF (year_str == '2008') hc = 17.16
    IF (year_str == '2007') hc = 16.86
    IF (year_str == '2006') hc = 16.57
    IF (year_str == '2005') hc = 16.27
    IF (year_str == '2004') hc = 15.97
    IF (year_str == '2003') hc = 15.68
    IF (flag_tree == 4) THEN
       hc = 4.82
    ENDIF
    write(*,*) 'this is the updated hc: ', hc

    !The SEP for all years are based on year 2010 measurements
    !C1: measured chamber temperature (celcius): C2: measured monoterpene emission rate (ng g-1 s-1)
    !'J' refers to the fact that this data comes from Juho
    OPEN(unit=8890,file = TRIM(ADJUSTL(input_dir_station_data))//'/EFJ_mono.txt',status='old')
    READ(8890,*)((EF_meas(i20,j20),i20=1,2),j20=1,1488)
    CLOSE(unit=8890)

    ! The concentration of methane for all years are based on year 2014 meaurements
    OPEN(unit=8891,file = TRIM(ADJUSTL(input_dir_station)) // '/2014/' // month_str_Mxx // '/hyy_conc_CH4.txt', status='old')
    READ(8891,*) (hyy_ch4(j20),j20=1,ntime_local)
    CLOSE(unit=8891)

    do j20 = 1,ntime_local
      difftop1(j20) = hyy_mix(2,j20)
      dirtop1(j20)  = hyy_mix(3,j20)           !r for Welgegund this the total global radiation

      do i20 = 1,6
         tem(i20,j20) = hyy_mix(i20+5,j20)
      enddo

      do i20 = 1,6
         hum(i20,j20) = hyy_mix(i20+11,j20)
      enddo
      do i20 = 1,4
         speed(i20,j20) = hyy_mix(i20+17,j20)
      enddo
      stat_cs_H2SO4(j20)     = hyy_mix(4,j20)
      stat_cs_HNO3(j20)    = hyy_mix(24,j20)
      stat_glob(j20)    = hyy_mix(25,j20)
      stat_PAR(j20)    = hyy_mix(26,j20)
      !speed(5,j20)    = hyy_mix(5,j20)       !r this data not used
      stat_albedo(j20)     = hyy_mix(28,j20)
      stat_sm_3(j20)     = hyy_soil(5,j20)       !r 5  Volumetric water content (m3/m3) in 23-60 cm
      stat_sm_2(j20)       = hyy_soil(4,j20)       !r 4  Volumetric water content (m3/m3) in 5-23 cm
      stat_sm_1(j20)       = hyy_soil(2,j20)       !r 2  Volumetric water content (m3/m3) in organic layer in -5-0 cm

    enddo

    do i=1,ntime_local
      dirtop1(i)  = max(0.,dirtop1(i))
      difftop1(i) = max(0.,difftop1(i))
    enddo

    IF (flag_tree == 4) THEN !CLEARCUT
        !This is used for the clear cut and includes a modification of the
        !measured data from Haapanala et al., BGS, 2012.
        !C1: daily monoterpene flux assuming constant flux from day
        !116-177. Unit: ug m^-2 h^-1
        !C2: daily monoterpene flux assuming decreasing flux from day
        !116-177 with a start value of around 7400. Unit: ug m^-2 h^-1
        !C3:  daily monoterpene flux assuming decreasing flux from day
        !116-177 with a start value of around 10000. Unit: ug m^-2 h^-1
        open(2715,file=TRIM(ADJUSTL(input_dir_station_info)) // '/SAMI_IN.txt')
        do J = 1,366
           read(2715,*) (SAMI(J,I), I=1,3)
        enddo
        close(2715)
    ENDIF

      IF (flag_tree == 3) THEN !BIRCH
          open(2716,file=TRIM(ADJUSTL(input_dir_station_info)) // '/lai_function.txt')
          do J = 1,366
             read(2716,*) LAI_function(J)
          enddo
          close(2716)
      ENDIF

  ELSEIF (STATION .EQ. 'MAN') THEN

     OPEN(unit=4,file = TRIM(ADJUSTL(input_dir_station_data))//'/Gas_filled.txt',status='old')
     READ(4,*) ((Gas_in(i20,j20),j20=1,6),i20=1,1488)

     CLOSE(unit=4)

     OPEN(unit=2,file = TRIM(ADJUSTL(input_dir_station_data))//'/Mix_filled.txt',status='old')

     READ(2,*) ((Mix_in(i20,j20),j20=1,35),i20=1,1488)
     CLOSE(unit=2)


     OPEN(unit=5,file = TRIM(ADJUSTL(input_dir_station_data))//'/Soil.txt',status='old')
     READ(5,*)((hyy_soil(i20,j20),i20=1,10),j20=1,1488)
     CLOSE(unit=5)


     OPEN(unit=5,file = TRIM(ADJUSTL(input_dir_station_data))//'/sounding.txt',status='old')
     READ(5,*)((sounding(i20,j20),j20=1,8),i20=1,1809)
     CLOSE(unit=5)

! Columns in Mix_filled.txt
!! 1          % day of year (DOY)
!! 2          % Diffuse solar radiation (W/m2)
!! 3          % Direct solar radiation (W/m2)
!! 4          % Global radiation/ incoming shortwave radiation (W/m2)
!! 5          % incoming PAR (micromols/m2/s) at 3.7m
!! 6:10       % Temperature (K) at 2 7 16 30 43m
!! 11:15      % Humidity (kg/m3) at 2 7 16 30 43m
!! 16:19      % Pressure (Pa) at 2 7 16 43m
!! 20:24      % Wind speed (m/s) at 2 7 16 30 43m
!! 25         % condensation sink for h2so4
!! 26         % condensation sink for hno3
!! 27:31      % in %.
!! 32:35      % Rlw.in, Rlw.out, Rsw.in, Rsw.out at 22m


       do j20 = 1,1488

          difftop1(j20) = Mix_in(j20,2)
          dirtop1(j20)  = Mix_in(j20,3)

          do i20 = 1,5
             tem(i20,j20) = Mix_in(j20,i20+5)
          enddo

          do i20 = 1,5
             hum(i20,j20) = Mix_in(j20,i20+10)
          enddo

          do i20 = 1,5
             speed(i20,j20) = Mix_in(j20, i20+19)
          enddo

          stat_cs_H2SO4(j20)     = Mix_in(j20,25)
          stat_cs_HNO3(j20)    = Mix_in(j20,26)
          stat_glob(j20)    = Mix_in(j20,4)
          stat_PAR(j20)    = Mix_in(j20,5)
          !speed(5,j20)    = hyy_mix(5,j20)
          stat_sm_3(j20)     = hyy_soil(2,j20)

       enddo

       do i=1,1488
          dirtop1(i)  = max(0.,dirtop1(i))
          difftop1(i) = max(0.,difftop1(i))
       enddo

    ENDIF

!!$
!!$    open(unit=12222,file = TRIM(ADJUSTL(output_dir))//'/QQQ1.dat')        ! air moisture, kg/m3
!!$    do j20 = 1, 1488
!!$       write(12222,'(5E20.8)') (hum(i20,j20),i20=1,5)
!!$    enddo
!!$    close(unit=12222)

    abl = 2
    if (TRIM(ADJUSTL(STATION)) .eq. 'hyytiala') then
       if (abl .eq. 1) then

          open(unit=3,file = TRIM(ADJUSTL(input_dir_station_data))//'/border.txt',status='old')
          read(3,*)((border(i20,j20),i20=1,5),j20=1,1488)
          close(unit=3)

       else

          !open(unit=3,file = TRIM(ADJUSTL(input_dir_station_data))//'/border_new.txt',status='old')
          !read(3,*)((border_abl(i20,j20),i20=1,6),j20=1,125)
          !close(unit=3)


       endif
    endif

    IF (flag_ecmwf == 1) THEN
      ntime_ECMWF = 4 * MONTH_DAYS(mon)  ! 0, 6, 12, 18 every day

      ! OPEN(unit=11,file = TRIM(ADJUSTL(input_dir_station_data))//'/O3',status='old')
      ! READ(11,*)((INPUT_O3_mixing_ratio(i20,j20),j20=1,43),i20=1,ntime_ECMWF)

      OPEN(unit=12,file = TRIM(ADJUSTL(input_dir_station_data))//'/ECMWF_qv.txt',status='old')
      READ(12,*)((input_qv(i20,j20),j20=1,43),i20=1,ntime_ECMWF)

      OPEN(unit=13,file = TRIM(ADJUSTL(input_dir_station_data))//'/ECMWF_uwind.txt',status='old')
      READ(13,*)((input_uwind(i20,j20),j20=1,43),i20=1,ntime_ECMWF)

      OPEN(unit=14,file = TRIM(ADJUSTL(input_dir_station_data))//'/ECMWF_vwind.txt',status='old')
      READ(14,*)((input_vwind(i20,j20),j20=1,43),i20=1,ntime_ECMWF)

      OPEN(unit=15,file = TRIM(ADJUSTL(input_dir_station_data))//'/ECMWF_temp.txt',status='old')
      READ(15,*)((input_temp(i20,j20),j20=1,43),i20=1,ntime_ECMWF)

      OPEN(unit=16,file = TRIM(ADJUSTL(input_dir_station_data))//'/ECMWF_geop.txt',status='old')
      READ(16,*)((input_geop(i20,j20),j20=1,43),i20=1,ntime_ECMWF)

!      OPEN(unit=17,file = TRIM(ADJUSTL(input_dir_station_data))//'/RH',status='old')
!      READ(17,*)((INPUT_RH(i20,j20),j20=1,43),i20=1,ntime_ECMWF)

      ! READ ONE MORE DATAPOINT IN THE FOWLLING DAY
      ntime_ECMWF = ntime_ECMWF+1
      !READ(11,*,END=8888)(INPUT_O3_mixing_ratio(ntime_ECMWF,j20),j20=1,43)
      !CLOSE(unit=11)

      READ(12,*,END=8888)(input_qv(ntime_ECMWF,j20),j20=1,43)
      CLOSE(unit=12)

      READ(13,*)(input_uwind(ntime_ECMWF,j20),j20=1,43)
      CLOSE(unit=13)

      READ(14,*)(input_vwind(ntime_ECMWF,j20),j20=1,43)
      CLOSE(unit=14)

      READ(15,*)(input_temp(ntime_ECMWF,j20),j20=1,43)
      CLOSE(unit=15)

      READ(16,*)(input_geop(ntime_ECMWF,j20),j20=1,43)
      CLOSE(unit=16)

      !READ(17,*)(INPUT_RH(ntime_ECMWF,j20),j20=1,43)
      !CLOSE(unit=17)

      GOTO 9999

8888   CONTINUE


       ! if (mon .eq. 12) then
       if (.TRUE.) then
          ! CLOSE (unit=11)
          CLOSE (unit=12)
          CLOSE (unit=13)
          CLOSE (unit=14)
          CLOSE (unit=15)
          CLOSE (unit=16)
          ! INPUT_O3_mixing_ratio(ntime_ECMWF,:) = INPUT_O3_mixing_ratio(ntime_ECMWF-1,:)
          input_qv(ntime_ECMWF,:) = input_qv(ntime_ECMWF-1,:)
          input_uwind(ntime_ECMWF,:) = input_uwind(ntime_ECMWF-1,:)
          input_vwind(ntime_ECMWF,:) = input_vwind(ntime_ECMWF-1,:)
          input_temp(ntime_ECMWF,:) = input_temp(ntime_ECMWF-1,:)
          input_geop(ntime_ECMWF,:) = input_geop(ntime_ECMWF-1,:)

       else

          ! CLOSE (unit=11)
          CLOSE (unit=12)
          CLOSE (unit=13)
          CLOSE (unit=14)
          CLOSE (unit=15)
          CLOSE (unit=16)
          !CLOSE (unit=17)
          WRITE(next_month, '(I2.2)') mon+1
          month_str_Mxx = 'M'//next_month

          !! read in the first data point in the following month
          ! OPEN(unit=11,file = TRIM(ADJUSTL(input_dir_station_data))//'/O3',status='old')
          ! READ(11,*)(INPUT_O3_mixing_ratio(ntime_ECMWF,j20),j20=1,43)
          ! CLOSE(unit=11)

          OPEN(unit=12,file = TRIM(ADJUSTL(input_dir_station_data))//'/ECMWF_qv.txt',status='old')
          READ(12,*)(input_qv(ntime_ECMWF,j20),j20=1,43)
          CLOSE(unit=12)

          OPEN(unit=13,file = TRIM(ADJUSTL(input_dir_station_data))//'/ECMWF_uwind.txt',status='old')
          READ(13,*)(input_uwind(ntime_ECMWF,j20),j20=1,43)
          CLOSE(unit=13)

          OPEN(unit=14,file = TRIM(ADJUSTL(input_dir_station_data))//'/ECMWF_vwind.txt',status='old')
          READ(14,*)(input_vwind(ntime_ECMWF,j20),j20=1,43)
          CLOSE(unit=14)

          OPEN(unit=15,file = TRIM(ADJUSTL(input_dir_station_data))//'/ECMWF_temp.txt',status='old')
          READ(15,*)(input_temp(ntime_ECMWF,j20),j20=1,43)
          CLOSE(unit=15)

          OPEN(unit=16,file = TRIM(ADJUSTL(input_dir_station_data))//'/ECMWF_geop.txt',status='old')
          READ(16,*)(input_geop(ntime_ECMWF,j20),j20=1,43)
          CLOSE(unit=16)

          !OPEN(unit=17,file = TRIM(ADJUSTL(input_dir_station_data))//'/RH',status='old')
          !READ(17,*)(INPUT_RH(ntime_ECMWF,j20),j20=1,43)
          !CLOSE(unit=17)

          WRITE(next_month, '(I2.2)') mon
          month_str_Mxx = 'M'//next_month

       endif

9999   CONTINUE

       input_wind = input_uwind   ! this done for getting the time stamp (first 6 columns)
       DO i = 7, 43
          input_wind(:,i) = sqrt(input_uwind(:,i)**2+input_vwind(:,i)**2)
       ENDDO

       !!       ntime_ECMWF = ntime_ECMWF - 1

       input_geom = input_geop / grav  ! geopotential height to geometric height

    ENDIF

   !r added by rosa
   ! long wave radiation (downwelling) from reanalysis
   OPEN(unit=11,file = TRIM(ADJUSTL(input_dir_station_data))//'/LWR.txt',status='old')
   READ(11,*) (LWR_in(j20),j20=1, (ntime_local/6+2))
   CLOSE(unit=11)
   !r
   stat_gsoil= -999
   OPEN(unit=11,file = TRIM(ADJUSTL(input_dir_station_data))//'/soil_heat_flux.txt',status='old')
   READ(11,*) (stat_gsoil(j20),j20=1,ntime_local)  ! heat flux to the soil based on measurements
   CLOSE(unit=11)
END SUBROUTINE READ_DATA_OLD


SUBROUTINE CONVERT_ECMWF

  ECMWF_uwind = 0.0_dp
  ECMWF_vwind = 0.0_dp
  ECMWF_temp  = 0.0_dp
  ECMWF_qv    = 0.0_dp
  ECMWF_rhov  = 0.0_dp
  ECMWF_strd  = input_strd(:, 7)

  ECMWF_TIME = input_uwind(:, 1:6)

  ! Interpolate input ECMWF data to model levels from ECMWF pressure levels (geometric heights)
  DO j = 1, ntime_ECMWF
    IF (TRIM(ADJUSTL(STATION)) == 'hyytiala') THEN
      CALL interp_1d( 37, input_geom(j, 7:), input_uwind(j, 7:), kz, z+masl, ECMWF_uwind(j, :) )
      CALL interp_1d( 37, input_geom(j, 7:), input_vwind(j, 7:), kz, z+masl, ECMWF_vwind(j, :) )
      CALL interp_1d( 37, input_geom(j, 7:), input_temp (j, 7:), kz, z+masl, ECMWF_temp (j, :) )
      CALL interp_1d( 37, input_geom(j, 7:), input_qv   (j, 7:), kz, z+masl, ECMWF_qv   (j, :) )
      CALL interp_1d( 37, input_geom(j, 7:), input_rhov (j, 7:), kz, z+masl, ECMWF_rhov (j, :) )

      ECMWF_dtemp(j) = (input_temp(j,31)-input_temp(j,33))/(input_geom(j,31)-input_geom(j,33))  ! 31: 850 hPa, 33: 900 hPa

    ELSEIF (TRIM(ADJUSTL(STATION)) == 'Manitou') THEN
      CALL interp_1d( 37, input_geom(j, 7:), input_temp(j, 7:), kz, z+masl, ECMWF_temp(j, :) )
      CALL interp_1d( 37, input_geom(j, 7:), input_uwind(j, 7:), kz, z+masl, ECMWF_uwind(j, :) )
      CALL interp_1d( 37, input_geom(j, 7:), input_vwind(j, 7:), kz, z+masl, ECMWF_vwind(j, :) )
      CALL interp_1d( 37, input_geom(j, 7:), input_qv(j, 7:), kz, z+masl, ECMWF_qv(j, :) )

      ECMWF_dtemp(j) = (input_temp(j,28)-input_temp(j,30))/(input_geom(j,28)-input_geom(j,30)) ! temp gradient 3600 and 2200m above see level
    ENDIF
  ENDDO

  ! Set top boundary layer values to ECMWF data
  border_abl(1, 1:249) = ECMWF_temp (1:249, kz)
  border_abl(2, 1:249) = ECMWF_uwind(1:249, kz)
  border_abl(3, 1:249) = ECMWF_vwind(1:249, kz)
  border_abl(4, 1:249) = ECMWF_rhov (1:249, kz)
  border_abl(5, 1:249) = ECMWF_dtemp(1:249)

END SUBROUTINE CONVERT_ECMWF


SUBROUTINE output_ascii_write()
  CHARACTER(LEN=10) :: str_prefix
  CHARACTER(LEN=2)    :: Jstamp
  output_number = 0

  !================================================================================================!
  ! Meteorology
  !================================================================================================!
  CALL output_ascii_write_one('TEMP', ta1)
  CALL output_ascii_write_one('QQQ', qa1)        ! [kg m-3], absolute humidity
  CALL output_ascii_write_one('UUU', u1)        ! streamwise velocity, m/s
  CALL output_ascii_write_one('VVV', v1)        ! lateral velocity, m/s
  CALL output_ascii_write_one('WWW', w1)        ! vertical velocity, m/s
  CALL output_ascii_write_one('TKE', bt1)        ! turbul kinetic energy, m2/s2
  CALL output_ascii_write_one('DIS', bt1*dbt1)        ! dissipation rate TKE, m2/s3
  CALL output_ascii_write_one('EDD', kt1)        ! eddy diff momentum, m2/s
  CALL output_ascii_write_one('ALTxEDD', alt1*kt1)    ! eddy diff scalar, m2/s
  CALL output_ascii_write_one('LLL', l1)        ! mixing length, m
  CALL output_ascii_write_one('ALT', alt1)        ! converse Schmidt number
  CALL output_ascii_write_one('RIH', rih)        ! Richardson number
  CALL output_ascii_write_one('TSN', tsn)        ! sunlit leaf temperature, K
  CALL output_ascii_write_one('TSD', tsd)        ! shaded leaf temperature, K
  CALL output_ascii_write_one('QSN', qsn)        ! sunlit leaf moisture, kg/m3
  CALL output_ascii_write_one('QSD', qsd)        ! shaded leaf moisture, kg/m3
  CALL output_ascii_write_one('Ustar', ur)      ! local friction velocity, m/s
  CALL output_ascii_write_one('TSOIL', tsoil1)      ! soil temperature, K
  CALL output_ascii_write_one('WIND', sqrt(u1*u1+v1*v1))       ! mean wind, m/s
  CALL output_ascii_write_one('H_TURB', fluxh3)     ! turbulent heat flux Wt/m2
  CALL output_ascii_write_one('H_INT', fluxh)      ! canopy heat flux, Wt/m2
  CALL output_ascii_write_one('LE_TURB', fluxle3)    ! turb latent flux, Wt/m2
  CALL output_ascii_write_one('LE_INT', fluxle)     ! canopy latent flux, Wt/m2
  CALL output_ascii_write_one('Pres', pres)       ! pressure in Pa
  CALL output_ascii_write_one('Rh', RH)         ! realtive humidity in %
  CALL output_ascii_write_one('H2O', H2O)        ! H2O in molecules / cm3 (text file)
  CALL output_ascii_write_one('TSN_M', tsn_megan)      ! sunlit leaf temperature, K from Megan
  CALL output_ascii_write_one('TSD_M', tsd_megan)      ! shaded leaf temperature, K from Megan
  CALL output_ascii_write_one('VPD_S', EM_VPD_S)      ! VPD for sunny leafs in kPa
  CALL output_ascii_write_one('VPD_C', EM_VPD_C)      ! VPD for shaded leafs in kPa
  CALL output_ascii_write_one('t_nud', tnud)      !
  CALL output_ascii_write_one('q_nud', qnud)      !
  CALL output_ascii_write_one('u_nud', unud)      !
  CALL output_ascii_write_one('BALANCE', (/balans/), '(6I6, E20.8)')    ! radiation closure, Wt/m2
  CALL output_ascii_write_one('SKY EMIS', (/emm/), '(6I6, E20.8)')   ! sky emissivity, approx
  CALL output_ascii_write_one('ALBEDO', (/albedo_f/), '(6I6, E20.8)')     ! cover albedo
  CALL output_ascii_write_one('SOIL_FLUX', (/pp/), '(6I6, E20.8)')  ! flux into soil, Wt/m2
  CALL output_ascii_write_one('ABL_HEIGHT', (/pblh/), '(6I6, E20.8)') ! ABL height, m
  CALL output_ascii_write_one('Zenith', (/zenith_deg/), '(6I6, E20.8)')     ! Solar zenith angle
  CALL output_ascii_write_one('fktt_fktd', (/fktt, fktd/), '(6I6, 2E20.8)')
  CALL output_ascii_write_one('DIR_RAD', (/rsnt/), '(6I6, E20.8)')    ! downward direct solar radiation (W/m2)
  CALL output_ascii_write_one('DIFF_RAD', (/rskt/), '(6I6, E20.8)')   ! diffuse direct solar radiation (W/m2)
  CALL output_ascii_write_one('LWR_DOWN', (/fird(k_canopy)/), '(6I6, E20.8)')   ! LWR down at the top of the canopy (W/m2)
  CALL output_ascii_write_one('LWR_UP', (/firu(k_canopy)/), '(6I6, E20.8)')     ! LWR up at the top of the canopy (W/m2)
  CALL output_ascii_write_one('PAR_UP', (/fphu(k_canopy)/), '(6I6, E20.8)')     ! PAR up at the top of the canopy (W/m2)
  CALL output_ascii_write_one('NIR_UP', (/fniu(k_canopy)/), '(6I6, E20.8)')     ! NIR up at the top of the canopy (W/m2)
  CALL output_ascii_write_one('LH_FLUX', (/fluxle3(k_canopy)/), '(6I6, E20.8)')    ! latent heat flux at the top of the canopy (W/m2)
  CALL output_ascii_write_one('SH_FLUX', (/fluxh3(k_canopy)/), '(6I6, E20.8)')     ! sensible heat flux at the top of the canopy (W/m2)
  CALL output_ascii_write_one('SOIL_Q', (/wg1/), '(6I6, 2E20.8)')     ! soil water content
  CALL output_ascii_write_one('Rd_dir', (/rsnt/), '(6I6, E20.8)')     ! downward direction radiation
  CALL output_ascii_write_one('Rd_dif', (/rskt/), '(6I6, E20.8)')     ! downward diffuse radiation
  CALL output_ascii_write_one('air', air)     ! air number concentration

  !================================================================================================!
  ! Gas and gas flux
  !================================================================================================!
  ! str_prefix = 'Gas_'
  DO I=1,outspc_count
    IF (outspc_cheminds(I) > 0) THEN
      gname = SPC_NAMES(outspc_cheminds(I))
      CALL output_ascii_write_one('Gas_'//TRIM(ADJUSTL(gname)), CH_CONS_ALL(:, outspc_cheminds(I)))
      CALL output_ascii_write_one('Flux_'//TRIM(ADJUSTL(gname)), CH_CONS_FLUX(:, outspc_cheminds(I)))
      ! CALL output_ascii_write_one('Qconc_'//TRIM(ADJUSTL(gname)), Qconc(:, outspc_cheminds(I)))
      ! CALL output_ascii_write_one('Qemis_'//TRIM(ADJUSTL(gname)), Qemis(:, outspc_cheminds(I)))
      ! CALL output_ascii_write_one('Qchem_'//TRIM(ADJUSTL(gname)), Qchem(:, outspc_cheminds(I)) - Qemis(:, outspc_cheminds(I)))
      ! CALL output_ascii_write_one('Qdepo_'//TRIM(ADJUSTL(gname)), Qdepo(:, outspc_cheminds(I)))
      ! CALL output_ascii_write_one('Qturb_'//TRIM(ADJUSTL(gname)), Qturb(:, outspc_cheminds(I)))
      ! CALL output_ascii_write_one('Qturbnow_'//TRIM(ADJUSTL(gname)), Qturb_now(:, outspc_cheminds(I)))
      ! CALL output_ascii_write_one('Qfluxdep_'//TRIM(ADJUSTL(gname)), Qfluxdep(:, outspc_cheminds(I)))
    END IF
  END DO
  !================================================================================================!
  ! Aerosol
  !================================================================================================!
  IF (flag_aero == 1) THEN
  !================================================================================================!
  ! HOMs
!    if ( any(outspc_HOM10_ind > 0) ) then
!      CALL output_ascii_write_one('HOM10', sum(CH_CONS_ALL(:, pack(outspc_HOM10_ind, outspc_HOM10_ind>0)),2))
!      CALL output_ascii_write_one('HOM20', sum(CH_CONS_ALL(:, pack(outspc_HOM20_ind, outspc_HOM20_ind>0)),2))
!    end if
  !================================================================================================!
!    DO I=1, n_bins_par
!      CALL output_ascii_write_one('nconc_aerosol', n_conc(:, I), isappend_in=.TRUE.)

!      ! Do not advance the file id
!      output_number = output_number - 1
!    END DO
     CALL output_ascii_write_one('CH_cs_H2SO4', (/CH_cs_H2SO4/))
!     CALL output_ascii_write_one('RADIUS', (/RADIUS(13,:)/), AR_OUTPUT_FORMAT)

   CALL output_ascii_write_one('NUC_RATE', (/NUC_RATE(1:kz)/))
   CALL output_ascii_write_one('ION_NUC_RATE', (/ION_NUC_RATE(1:kz)/))
   DO k = 1, kz
     write(layerstamp,'(F6.1)') z(k)
     CALL output_ascii_write_one('N_CONC'//layerstamp, (/N_CONC(k,:)/), AR_OUTPUT_FORMAT)
     CALL output_ascii_write_one('PAR_FLUX'//layerstamp, (/PAR_FLUX(k,:)/), AR_OUTPUT_FORMAT)
     CALL output_ascii_write_one('GR'//layerstamp, (/GR(k,:)/), AR_OUTPUT_FORMAT)
   ENDDO
  END IF

  !================================================================================================!
  ! OH, O3 and NO3 reactivity
  !================================================================================================!
    IF (CH_oh_count > 0) THEN
       DO CH_oh_i = 1, CH_oh_count
       CALL output_ascii_write_one('GAS_'//TRIM(SPC_NAMES(CH_oh_indices(CH_oh_i))), CH_oh_cons3(:,CH_oh_i,3))
       ENDDO
    ENDIF
    IF (CH_o3_count > 0) THEN
       DO CH_o3_i = 1, CH_o3_count
       CALL output_ascii_write_one('GAS_'//TRIM(SPC_NAMES(CH_o3_indices(CH_o3_i))), CH_o3_cons3(:,CH_o3_i,3))
       ENDDO
    ENDIF
    IF (CH_no3_count > 0) THEN
       DO CH_no3_i = 1, CH_no3_count
       CALL output_ascii_write_one('GAS_'//TRIM(SPC_NAMES(CH_no3_indices(CH_no3_i))), CH_no3_cons3(:,CH_no3_i,3))
       ENDDO
    ENDIF
  !================================================================================================!
  ! Photolysis rate
  !================================================================================================!
  DO I=1, NPHOT
     write(Jstamp,'(I2)') I
  CALL output_ascii_write_one('J_'//Jstamp,(/CH_J_values_ALL(:,I)/))
  END DO
  !================================================================================================!
  ! Emission fluxes
  !================================================================================================!
  ! DO I=1,SIZE(emi_ind)
  !   CALL output_ascii_write_one('Emi_'//SPC_NAMES(emi_ind(I)), EM_EMI(:, I))
  ! END DO

  DO I=1, outemi_count
    CALL output_ascii_write_one('Emi_'//SPC_NAMES(outemi_cheminds(I)), EM_EMI(:, outemi_meganinds(I)))
  END DO

  !================================================================================================!
  ! Output by Putian
  !================================================================================================!
  !***** Vd *****!
  DO I=1,outspc_count
    IF (outspc_cheminds(I) > 0) THEN
      CALL output_ascii_write_one('Vd_'//outspc_names(I), vdep(:, outspc_cheminds(I)))
    END IF
  END DO

  CALL output_ascii_write_one('eta', psn)
  CALL output_ascii_write_one('mu', psk)
  CALL output_ascii_write_one('gl', gl)
  CALL output_ascii_write_one('gd', gd)
  CALL output_ascii_write_one('RSu_PAR', fphu)
  CALL output_ascii_write_one('RSd_PAR', fphd)
  CALL output_ascii_write_one('rac', dep_output(:, ind_O3, 1))
  CALL output_ascii_write_one('rstm', dep_output(:, ind_O3, 2))
  CALL output_ascii_write_one('rbveg', dep_output(:, ind_O3, 3))
  CALL output_ascii_write_one('frac_ws', dep_output(:, ind_O3, 4))
  CALL output_ascii_write_one('rleaf', dep_output(:, ind_O3, 5))
  CALL output_ascii_write_one('rleafw', dep_output(:, ind_O3, 6))
  CALL output_ascii_write_one('rsveg', dep_output(:, ind_O3, 7))
  CALL output_ascii_write_one('rswet', dep_output(:, ind_O3, 8))
  CALL output_ascii_write_one('rtot', dep_output(:, ind_O3, 9))
  CALL output_ascii_write_one('fvpd', dep_output(:, ind_O3, 10))
  CALL output_ascii_write_one('rstm_h2o', dep_output(:, ind_O3, 11))
  CALL output_ascii_write_one('gs_h2o_sn', dep_output(:, ind_O3, 12))
  CALL output_ascii_write_one('gs_h2o_sd', dep_output(:, ind_O3, 13))
  CALL output_ascii_write_one('gs_h2o_avg', dep_output(:, ind_O3, 14))
  CALL output_ascii_write_one('gstm_h2o_sn', gstm_h2o_sn)
  CALL output_ascii_write_one('gstm_h2o_sd', gstm_h2o_sd)

END SUBROUTINE output_ascii_write


SUBROUTINE output_ascii_write_one(var_name, var, output_format, isappend_in)
  CHARACTER(LEN=*) :: var_name
  REAL(dp) :: var(:)
  CHARACTER(LEN=*), OPTIONAL :: output_format
  LOGICAL, OPTIONAL :: isappend_in

  CHARACTER(LEN=300) :: current_format
  INTEGER :: file_id
  LOGICAL :: isappend

  IF (PRESENT(output_format)) THEN
    current_format = output_format
  ELSE
    current_format = CH_OUTPUT_FORMAT
  END IF

  IF (PRESENT(isappend_in)) THEN
    isappend = isappend_in
  ELSE
    isappend = .FALSE.
  END IF

  output_number = output_number + 1
  file_id = UNIT_START+output_number
  IF (.not. isopen(output_number)) THEN
    OPEN( UNIT=file_id, FILE=TRIM(ADJUSTL(output_dir))//'/'//TRIM(var_name)//'.dat' )
    IF (isappend) THEN
      WRITE(file_id, *) ''
      CLOSE(file_id)
      OPEN(UNIT=file_id, FILE=TRIM(ADJUSTL(output_dir))//'/'//TRIM(var_name)//'.dat', POSITION='APPEND', ACTION='WRITE', STATUS='OLD')
    END IF
    isopen(output_number) = .TRUE.
  END IF
  WRITE(file_id, TRIM(current_format)) now_date(1:6), var
END SUBROUTINE output_ascii_write_one


SUBROUTINE open_output_files
  IF (kz .EQ. 51) THEN
    AR_OUTPUT_FORMAT = FORMAT1
    CH_OUTPUT_FORMAT = FORMAT7
  ELSE
    AR_OUTPUT_FORMAT = FORMAT1
    CH_OUTPUT_FORMAT = FORMAT8
  ENDIF

  ! LUXI BE AWARE THAT YOU ARE ALREADY USING UNIT=101 FOR THE GAS OUTPUT!
  open(unit=101,file = TRIM(ADJUSTL(output_dir))//'/border1.dat')
  write(101,'(249E12.3)') border_abl(1,1:249)
  write(101,'(249E12.3)') border_abl(2,1:249)
  write(101,'(249E12.3)') border_abl(3,1:249)
  write(101,'(249E12.3)') border_abl(4,1:249)
  write(101,'(249E12.3)') border_abl(5,1:249)
  close(unit=101)

  open(unit=102,file = TRIM(ADJUSTL(output_dir))//'/border2.dat')
  write(102,'(249E12.3)') border_abl(1,1:249)
  write(102,'(249E12.3)') border_abl(2,1:249)
  write(102,'(249E12.3)') border_abl(3,1:249)
  write(102,'(249E12.3)') border_abl(4,1:249)
  write(102,'(249E12.3)') border_abl(5,1:249)
  close(unit=102)

  !----------------- open file for aerosols----------------------------------------------------------
  IF (flag_aero .EQ. 2) THEN
    open(unit=1001, file = TRIM(ADJUSTL(output_dir))//'/NUC_RATE.dat' )
    write(1001, CH_OUTPUT_FORMAT) 0, 0, 0, 0, 0, 0, z

    open(unit=1002, file = TRIM(ADJUSTL(output_dir))//'/ION_NUC_RATE.dat' )
    write(1002, CH_OUTPUT_FORMAT) 0, 0, 0, 0, 0, 0, z

    open(unit=1003, file = TRIM(ADJUSTL(output_dir))//'/CS_Measured.dat' )
    !write(1003, CH_OUTPUT_FORMAT) 0, 0, 0, 0, 0, 0, z

    open(unit=1011, file = TRIM(ADJUSTL(output_dir))//'/RADIUS.dat' )

    DO k = 1, kz
      write(layerstamp,'(F6.1)') z(k)

      open(unit=1500+k, file=TRIM(ADJUSTL(output_dir))//'/GR'//layerstamp//'.dat')
      open(unit=1600+k, file=TRIM(ADJUSTL(output_dir))//'/PAR_num'//layerstamp//'.dat')
      open(unit=1700+k, file=TRIM(ADJUSTL(output_dir))//'/Sink'//layerstamp//'.dat')
      open(unit=1900+k, file=TRIM(ADJUSTL(output_dir))//'/Vap'//layerstamp//'.dat')

      !          open(unit=2000+k, file=TRIM(ADJUSTL(output_dir))//'/Par_flux'//layerstamp//'.dat')
    ENDDO
  ENDIF

  IF (flag_ecmwf == 1) THEN
    ! write converted ecmwf data to files
    OPEN(unit=700,file = TRIM(ADJUSTL(output_dir))//'/ECMWF_temp.dat')
    write(700,CH_OUTPUT_FORMAT) 0, 0, 0, 0, 0, 0, z

    OPEN(unit=701,file = TRIM(ADJUSTL(output_dir))//'/ECMWF_Specific_Humidity.dat')
    write(701,CH_OUTPUT_FORMAT) 0, 0, 0, 0, 0, 0, z

    OPEN(unit=702,file = TRIM(ADJUSTL(output_dir))//'/ECMWF_wind.dat')
    write(702,CH_OUTPUT_FORMAT) 0, 0, 0, 0, 0, 0, z

    OPEN(unit=703,file = TRIM(ADJUSTL(output_dir))//'/ECMWF_RH.dat')
    write(703,CH_OUTPUT_FORMAT) 0, 0, 0, 0, 0, 0, z

    DO k = 1, ntime_ECMWF
      WRITE (700,FORMAT7) (ECMWF_TIME(k,1:6)), (ECMWF_temp(k,i),i=1,kz)
      WRITE (701,FORMAT7) (ECMWF_TIME(k,1:6)), (ECMWF_qv(k,i),i=1,kz)
      WRITE (702,FORMAT7) (ECMWF_TIME(k,1:6)), (ECMWF_wind(k,i),i=1,kz)
      WRITE (703,FORMAT7) (ECMWF_TIME(k,1:6)), (ECMWF_RH(k,i),i=1,kz)
    ENDDO

    CLOSE(UNIT = 700)
    CLOSE(UNIT = 701)
    CLOSE(UNIT = 702)
    CLOSE(UNIT = 703)
  ENDIF

  !==============================================================================================================!
  ! Initialize stuff for output
  !==============================================================================================================!
  isopen = .FALSE.

END SUBROUTINE open_output_files


subroutine output_init()

  !-----------------------------------------------------------!
  ! Initialize parameters
  !-----------------------------------------------------------!


  !-----------------------------------------------------------!
  ! Initiate DIOS
  !-----------------------------------------------------------!

  call dios_init(stat)


  !-----------------------------------------------------------!
  ! Create nc file
  !-----------------------------------------------------------!

  ! Obtain file name
  mixf%fname = trim(adjustl(output_dir))//'/output.nc'

  ! Create file
  call dios_create(trim(mixf%fname), DIOS_NETCDF, DIOS_REPLACE, mixf%funit, stat)


  !-----------------------------------------------------------!
  ! Define global attributes
  !-----------------------------------------------------------!

  call dios_put_att(mixf%funit, DIOS_GLOBAL, 'flag_emis'      , flag_emis      , stat)
  call dios_put_att(mixf%funit, DIOS_GLOBAL, 'flag_chem'     , flag_chem     , stat)
  call dios_put_att(mixf%funit, DIOS_GLOBAL, 'flag_gasdrydep', flag_gasdrydep, stat)
  call dios_put_att(mixf%funit, DIOS_GLOBAL, 'flag_aero'     , flag_aero     , stat)
  call dios_put_att(mixf%funit, DIOS_GLOBAL, 'flag_emis_soil'     , flag_emis_soil     , stat)
  call dios_put_att(mixf%funit, DIOS_GLOBAL, 'flag_ecmwf'   , flag_ecmwf, stat)
  call dios_put_att(mixf%funit, DIOS_GLOBAL, 'nclo_end'     , nclo_end     , stat)
  call dios_put_att(mixf%funit, DIOS_GLOBAL, 'dt_mete'      , dt_mete      , stat)
  call dios_put_att(mixf%funit, DIOS_GLOBAL, 'dt_chem'      , dt_chem      , stat)
  call dios_put_att(mixf%funit, DIOS_GLOBAL, 'dt_aero'      , dt_aero      , stat)
  call dios_put_att(mixf%funit, DIOS_GLOBAL, 'canopy_height', hc           , stat)
  call dios_put_att(mixf%funit, DIOS_GLOBAL, 'abl'          , abl          , stat)
  call dios_put_att(mixf%funit, DIOS_GLOBAL, 'LAI_tot'      , LAI_tot      , stat)
  call dios_put_att(mixf%funit, DIOS_GLOBAL, 'LAI_curv'     , LAI_curv     , stat)
  call dios_put_att(mixf%funit, DIOS_GLOBAL, 'LAI_proj'     , LAI_proj     , stat)


  !-----------------------------------------------------------!
  ! Define dimensions
  !-----------------------------------------------------------!

  !----- Define dimensions
  call dios_def_dim(mixf%funit, 'time'  , DIOS_UNLIMITED, mixf%dimid_time  , stat)
  call dios_def_dim(mixf%funit, 'lev'   , kz            , mixf%dimid_lev   , stat)
  call dios_def_dim(mixf%funit, 'dp_dry_fs'  , n_bins_par , mixf%dimid_dp_dry_fs  , stat)
  ! call dios_def_dim(mixf%funit, 'icondv', kz            , mixf%dimid_icondv, stat)
  call dios_def_dim(mixf%funit, 'outspc', outspc_count  , mixf%dimid_outspc, stat)

  !----- Define dimension variables

  ! time

  !>> Get strings for start date and first day of month
  write(start_date_string, &
    '(I0.4,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2)') &
    start_date(1), '-', start_date(2), '-', start_date(3), ' ', &
    start_date(4), ':', start_date(5), ':', start_date(6)
  write(first_day_of_month_string, &
    '(I0.4,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2)') &
    start_date(1), '-', start_date(2), '-', 1, ' ', 0, ':', 0, ':', 0

  call dios_def_var(mixf%funit, 'time', DIOS_DOUBLE, &
                    (/mixf%dimid_time/), mixf%varid_time, stat)
  call dios_put_att(mixf%funit, mixf%varid_time, &
                    'unit', 's', stat)
  call dios_put_att(mixf%funit, mixf%varid_time, &
                    'long_name', 'time since the beginning of month', stat)
  call dios_put_att(mixf%funit, mixf%varid_time, &
                    'first_day_of_month', trim(first_day_of_month_string), stat)
  call dios_put_att(mixf%funit, mixf%varid_time, &
                    'start_date', trim(start_date_string), stat)

  ! lev
  call dios_def_var(mixf%funit, 'lev', DIOS_DOUBLE, &
                    (/mixf%dimid_lev/), mixf%varid_lev, stat)
  call dios_put_att(mixf%funit, mixf%varid_lev, &
                    'unit', 'm', stat)
  call dios_put_att(mixf%funit, mixf%varid_lev, &
                    'long_name', 'height above the ground', stat)

  ! dp_dry_fs
  call dios_def_var(mixf%funit, 'dp_dry_fs', DIOS_DOUBLE, &
                    (/mixf%dimid_dp_dry_fs/), mixf%varid_dp_dry_fs, stat)
  call dios_put_att(mixf%funit, mixf%varid_dp_dry_fs, &
                    'unit', 'm', stat)
  call dios_put_att(mixf%funit, mixf%varid_dp_dry_fs, &
                    'long_name', 'dry radius of aerosol particles in each size bin', stat)

  ! outspc
  call dios_def_var(mixf%funit, 'outspc', DIOS_INT, &
                    (/mixf%dimid_outspc/), mixf%varid_outspc, stat)
  call dios_put_att(mixf%funit, mixf%varid_outspc, &
                    'unit', '-', stat)
  call dios_put_att(mixf%funit, mixf%varid_outspc, &
                    'long_name', 'indices of output species in chemistry', stat)

  !-----------------------------------------------------------!
  ! Add variables
  !-----------------------------------------------------------!

  ! Allocate fields pointers
  allocate(mixf%f0d(NMAX_0D))
  allocate(mixf%f1d(NMAX_1D))
  allocate(mixf%f2d(NMAX_2D))

  ! Set counters for 1d and 2d variables
  n_0d = 0
  n_1d = 0
  n_2d = 0

  !----- Add variables to fields

  ! Invariant variables
  call add_variable(mixf, 'lad_meteo', 'LAD used in meteorology module', 'm2 m-3', &
                    (/mixf%dimid_lev/), .false., id1_lad_meteo)
  call add_variable(mixf, 'lad_megan', 'LAD used in megan module', 'm2 m-3', &
                    (/mixf%dimid_lev/), .false., id1_lad_megan)

  ! Time variant scalar variables
  call add_variable(mixf, 'lat', 'latitude', 'deg', &
                    (/-1/), .true., id0_lat)
  call add_variable(mixf, 'lon', 'longitude', 'deg', &
                    (/-1/), .true., id0_lon)
  call add_variable(mixf, 'pblh', 'PBL height', 'm', &
                    (/-1/), .true., id0_pblh)
  call add_variable(mixf, 'gsoil', 'soil heat flux', 'W m-2', &
                    (/-1/), .true., id0_gsoil)
  call add_variable(mixf, 'zenith', 'solar zenith angle', 'deg', &
                    (/-1/), .true., id0_zenith)
  call add_variable(mixf, 'Rdir_down', &
                    'downward direct solar radiation at canopy top', &
                    'W m-2', (/-1/), .true., id0_Rdir_down)
  call add_variable(mixf, 'Rdiff_down', &
                    'downward diffuse solar radiation at canopy top', &
                    'W m-2', (/-1/), .true., id0_Rdiff_down)
  call add_variable(mixf, 'glob', &
                    'global solar radiation', &
                    'W m-2', (/-1/), .true., id0_glob)
  call add_variable(mixf, 'NIR_down', &
                    'downward near-infrared radiation at canopy top', &
                    'W m-2', (/-1/), .true., id0_NIR_down)
  call add_variable(mixf, 'NIR_up', &
                    'upward near-infrared radiation at canopy top', &
                    'W m-2', (/-1/), .true., id0_NIR_up)
  call add_variable(mixf, 'PAR_down', &
                    'downward PAR at canopy top', &
                    'W m-2', (/-1/), .true., id0_PAR_down)
  call add_variable(mixf, 'PAR_up', &
                    'upward PAR at canopy top', &
                    'W m-2', (/-1/), .true., id0_PAR_up)
  call add_variable(mixf, 'L_down', &
                    'downward long wave radiation at canopy top', &
                    'W m-2', (/-1/), .true., id0_L_down)
  call add_variable(mixf, 'L_up', &
                    'upward long wave radiation at canopy top', &
                    'W m-2', (/-1/), .true., id0_L_up)

  ! Time variant 1d variables
  call add_variable(mixf, 'ua', 'u wind', 'm s-1', &
                    (/mixf%dimid_lev/), .true., id1_ua)
  call add_variable(mixf, 'va', 'v wind', 'm s-1', &
                    (/mixf%dimid_lev/), .true., id1_va)
  call add_variable(mixf, 'temp', 'air temperature', 'K', &
                    (/mixf%dimid_lev/), .true., id1_temp)
  call add_variable(mixf, 'rhov', 'absolute humidity', 'kg m-3', &
                    (/mixf%dimid_lev/), .true., id1_rhov)
  call add_variable(mixf, 'pres', 'air pressure', 'Pa', &
                    (/mixf%dimid_lev/), .true., id1_pres)
  call add_variable(mixf, 'tke', 'turbulent kinetic energy', 'm2 s-2', &
                    (/mixf%dimid_lev/), .true., id1_tke)
  call add_variable(mixf, 'diffm', 'eddy diffusivity of momentum', 'm2 s-1', &
                    (/mixf%dimid_lev/), .true., id1_diffm)
  call add_variable(mixf, 'diffh', 'eddy diffusivity of scalar', 'm2 s-1', &
                    (/mixf%dimid_lev/), .true., id1_diffh)
  call add_variable(mixf, 'Rich', 'Richardson number', '-', &
                    (/mixf%dimid_lev/), .true., id1_Rich)
  call add_variable(mixf, 'ustar', 'friction velocity', 'm s-1', &
                    (/mixf%dimid_lev/), .true., id1_ustar)
  call add_variable(mixf, 'mixlen', 'mixing length', 'm', &
                    (/mixf%dimid_lev/), .true., id1_mixlen)
  call add_variable(mixf, 'shf', 'sensible heat flux', 'W m-2', &
                    (/mixf%dimid_lev/), .true., id1_shf)
  call add_variable(mixf, 'lhf', 'latent heat flux', 'W m-2', &
                    (/mixf%dimid_lev/), .true., id1_lhf)
  call add_variable(mixf, 'Nair', 'number concentration of air', 'molec cm-3', &
                    (/mixf%dimid_lev/), .true., id1_Nair)

  ! Time variant 2d variables
  call add_variable(mixf, 'nconc', 'number concentration', 'molec cm-3', &
                    (/mixf%dimid_lev, mixf%dimid_outspc/), .true., id2_nconc)
  call add_variable(mixf, 'flux', 'vertical flux of species', 'molec cm-2 s-1', &
                    (/mixf%dimid_lev, mixf%dimid_outspc/), .true., id2_flux)
  call add_variable(mixf, 'nconc_par', 'particle number concentration', '# m-3', &
                    (/mixf%dimid_lev, mixf%dimid_dp_dry_fs/), .true., id2_nconc_par)


  !-----------------------------------------------------------!
  ! Define variables in the file
  !-----------------------------------------------------------!
  do i = 1, n_0d
    ! The scalar variables must be time variant, otherwise they need to be set
    ! as constant global attributes.
    call dios_def_var(mixf%funit, trim(mixf%f0d(i)%mf%vname), DIOS_DOUBLE, &
      (/mixf%dimid_time/), mixf%f0d(i)%varid, stat)

    ! Put variable attributes
    call dios_put_att(mixf%funit, mixf%f0d(i)%varid, 'unit', &
                      trim(mixf%f0d(i)%mf%unit), stat)
    call dios_put_att(mixf%funit, mixf%f0d(i)%varid, 'long_name', &
                      trim(mixf%f0d(i)%mf%lname), stat)
  end do


  do i = 1, n_1d
    ! Define the variable depending on if it changes with time
    if (mixf%f1d(i)%change_with_time) then
      call dios_def_var(mixf%funit, trim(mixf%f1d(i)%mf%vname), DIOS_FLOAT, &
        (/mixf%f1d(i)%dimids, mixf%dimid_time/), mixf%f1d(i)%varid, stat)
    else
      call dios_def_var(mixf%funit, trim(mixf%f1d(i)%mf%vname), DIOS_FLOAT, &
        mixf%f1d(i)%dimids, mixf%f1d(i)%varid, stat)
    end if

    ! Put variable attributes
    call dios_put_att(mixf%funit, mixf%f1d(i)%varid, 'unit', &
                      trim(mixf%f1d(i)%mf%unit), stat)
    call dios_put_att(mixf%funit, mixf%f1d(i)%varid, 'long_name', &
                      trim(mixf%f1d(i)%mf%lname), stat)
  end do


  do i = 1, n_2d
    write(*,*) mixf%f2d(i)%mf%vname
    ! Define the variable depending on if it changes with time
    if (mixf%f2d(i)%change_with_time) then
      call dios_def_var(mixf%funit, trim(mixf%f2d(i)%mf%vname), DIOS_FLOAT, &
        (/mixf%f2d(i)%dimids, mixf%dimid_time/), mixf%f2d(i)%varid, stat)
    else
      call dios_def_var(mixf%funit, trim(mixf%f2d(i)%mf%vname), DIOS_FLOAT, &
        mixf%f2d(i)%dimids, mixf%f2d(i)%varid, stat)
    end if

    ! Put variable attributes
    call dios_put_att(mixf%funit, mixf%f2d(i)%varid, 'unit', &
                      trim(mixf%f2d(i)%mf%unit), stat)
    call dios_put_att(mixf%funit, mixf%f2d(i)%varid, 'long_name', &
                      trim(mixf%f2d(i)%mf%lname), stat)
  end do


  ! End of definition
  call dios_enddef(mixf%funit, stat)

  ! Put the values of the dimensions
  call dios_put_var(mixf%funit, mixf%varid_lev, z, stat, start=(/1/), count=(/kz/))
  call dios_put_var(mixf%funit, mixf%varid_dp_dry_fs, CURRENT_PSD%dp_dry_fs, stat, &
                    start=(/1/), count=(/n_bins_par/))
  call dios_put_var(mixf%funit, mixf%varid_outspc, outspc_cheminds, stat, &
                    start=(/1/), count=(/outspc_count/))


  ! Set the values for invariant variables
  mixf%f1d(id1_lad_meteo)%field = s1
  mixf%f1d(id1_lad_megan)%field = s2

  ! Write the invariant variables to the output file
  call dios_put_var( mixf%funit, mixf%f1d(id1_lad_meteo)%varid, &
    mixf%f1d(id1_lad_meteo)%field, stat, &
    start=(/1/), count=mixf%f1d(id1_lad_meteo)%ndims )
  call dios_put_var( mixf%funit, mixf%f1d(id1_lad_megan)%varid, &
    mixf%f1d(id1_lad_megan)%field, stat, &
    start=(/1/), count=mixf%f1d(id1_lad_megan)%ndims )

  ! Initilize time record counter
  itrec = 0
end subroutine output_init


subroutine output_step()
  mixf%f0d(id0_lat   )%field = lat_deg
  mixf%f0d(id0_lon   )%field = lon_deg
  mixf%f0d(id0_pblh  )%field = pblh
  mixf%f0d(id0_gsoil )%field = pp
  mixf%f0d(id0_zenith)%field = zenith_deg

  mixf%f0d(id0_Rdir_down )%field = rsnt
  mixf%f0d(id0_Rdiff_down)%field = rskt
  mixf%f0d(id0_glob      )%field = glob
  mixf%f0d(id0_NIR_down  )%field = fnid(k_canopy)
  mixf%f0d(id0_NIR_up    )%field = fniu(k_canopy)
  mixf%f0d(id0_PAR_down  )%field = fphd(k_canopy)
  mixf%f0d(id0_PAR_up    )%field = fphu(k_canopy)
  mixf%f0d(id0_L_down    )%field = fird(k_canopy)
  mixf%f0d(id0_L_up      )%field = firu(k_canopy)

  mixf%f1d(id1_ua    )%field = u1
  mixf%f1d(id1_va    )%field = v1
  mixf%f1d(id1_temp  )%field = ta1
  mixf%f1d(id1_rhov  )%field = qa1
  mixf%f1d(id1_pres  )%field = pres
  mixf%f1d(id1_tke   )%field = bt1
  mixf%f1d(id1_diffm )%field = kt1
  mixf%f1d(id1_diffh )%field = alt1*kt1
  mixf%f1d(id1_Rich  )%field = rih
  mixf%f1d(id1_ustar )%field = ur
  mixf%f1d(id1_mixlen)%field = l1
  mixf%f1d(id1_shf   )%field = fluxh3
  mixf%f1d(id1_lhf   )%field = fluxle3
  mixf%f1d(id1_Nair  )%field = air

  ! write(*,*) 'outspc_cheminds', outspc_cheminds
  ! write(*,*) shape(CH_CONS_ALL), shape(CH_CONS_FLUX)
  ! write(*,*) shape(mixf%f2d(id2_nconc)%field)
  ! write(*,*) 'id2: ', id2_nconc, id2_flux, id2_nconc_par
  ! write(*,*) shape(mixf%f2d)
  ! write(*,*) shape(mixf%f2d(id2_nconc_par)%field)
  ! write(*,*) shape(N_CONC)
  mixf%f2d(id2_nconc)%field     = CH_CONS_ALL(:, outspc_cheminds)
  mixf%f2d(id2_flux)%field      = CH_CONS_FLUX(:, outspc_cheminds)
  mixf%f2d(id2_nconc_par)%field = N_CONC(:, :)
end subroutine output_step


subroutine output_write()
  ! Update time record counter
  itrec = itrec + 1

  ! Write time
  call dios_put_var( mixf%funit, mixf%varid_time, (/time_in_month/), &
    stat, start=(/itrec/), count=(/1/) )

  ! Write 0d variables
  do i = 1, n_0d
    if (mixf%f0d(i)%change_with_time) then
      call dios_put_var( mixf%funit, mixf%f0d(i)%varid, (/mixf%f0d(i)%field/), &
        stat, start=(/itrec/), count=(/1/) )
    end if
  end do

  ! Write 1d variables
  do i = 1, n_1d
    if (mixf%f1d(i)%change_with_time) then
      call dios_put_var(mixf%funit, mixf%f1d(i)%varid, mixf%f1d(i)%field, &
        stat, start=(/1, itrec/), count=(/mixf%f1d(i)%ndims(1), 1/))
    end if
  end do

  ! Write 2d variables
  do i = 1, n_2d
    if (mixf%f2d(i)%change_with_time) then
      call dios_put_var(mixf%funit, mixf%f2d(i)%varid, mixf%f2d(i)%field, &
        stat, &
        start=(/1, 1, itrec/), &
        count=(/mixf%f2d(i)%ndims(1), mixf%f2d(i)%ndims(2), 1/))
    end if
  end do
end subroutine output_write


!==============================================================================!
! New output write subroutines using netcdf format
!==============================================================================!
SUBROUTINE output_write_old()
  INTEGER :: I
  INTEGER, PARAMETER :: n_other=35
  LOGICAL, SAVE :: FIRST_TIME=.TRUE.

  hour = INT(time_in_day/3600.)
  minute = INT((time_in_day-hour*3600.)/60.)
  second = INT(time_in_day-hour*3600.-minute*60.)

  ! DATA AT ABOUT SMEAR LEVEL (time, variable)

  albedo_f=max(0.,(fniu(nz)+fphu(nz))/(rads2+1.))

  ! UVA and swr_hyy are used for old version input and output
  ! UVA = 0.
  ! DO J = 9,24
  !    UVA = UVA + swr_hyy(j,nxodrad)
  ! ENDDO

  ! IF (FIRST_TIME) THEN
  !   FIRST_TIME = .FALSE.
  !   !
  !   ! Write some critical information in the model.
  !   !
  !   OPEN(UNIT=10, FILE = TRIM(ADJUSTL(OUTPUT_DIR))//'/Info.dat')

  !   WRITE(10,'(A20, I10)') 'kz', kz
  !   WRITE(10,'(A20, I10)') 'flag_emis', flag_emis
  !   WRITE(10,'(A20, I10)') 'flag_chem', flag_chem
  !   WRITE(10,'(A20, I10)') 'flag_gasdrydep', flag_gasdrydep
  !   WRITE(10,'(A20, I10)') 'flag_aero', flag_aero
  !   WRITE(10,'(A20, I10)') 'flag_emis_soil', flag_emis_soil
  !   WRITE(10,'(A20, I10)') 'flag_ecmwf', flag_ecmwf
  !   WRITE(10,'(A20, I10)') 'nclo_end', nclo_end
  !   WRITE(10,'(A20, ES25.14E3)') 'dt_mete', dt_mete
  !   WRITE(10,'(A20, ES25.14E3)') 'dt_chem', dt_chem
  !   WRITE(10,'(A20, ES25.14E3)') 'dt_aero', dt_aero
  !   WRITE(10,'(A20, ES25.14E3)') 'hc', hc
  !   WRITE(10,'(A20, I10)') 'abl', abl
  !   WRITE(10,'(A20, ES25.14E3)') 'LAI_tot', LAI_tot
  !   WRITE(10,'(A20, ES25.14E3)') 'LAI_curv', LAI_curv
  !   WRITE(10,'(A20, ES25.14E3)') 'LAI_proj', LAI_proj
  !   WRITE(10,'(A20, 51ES25.14E3)') 'z', z    ! height levels
  !   WRITE(10,'(A20, 51ES25.14E3)') 's1', s1    ! LAD for meteorology
  !   WRITE(10,'(A20, 51ES25.14E3)') 's2', s2    ! LAD for MEGAN

  !   IF (flag_aero == 1) THEN
  !     WRITE(10,*) 'Nucleation scheme:', options%nuc_number
  !     WRITE(10,*) 'Nucleation rate:', ambient%nuc_coeff
  !     WRITE(10,*) 'Concentration of reaction products from NO3 oxidation is first increased by ', &
  !                no3_coe, ' times.'
  !     WRITE(10,*) 'Concentration of reaction products from OH, O3 and NO3 oxidation is then all increased by ', &
  !                vap_coe, ' times.'
  !   ENDIF

  !   CLOSE(10)
  ! END IF

  if (flag_format==0) then
   CALL output_ascii_write()

  ELSE

  IF (FIRST_TIME) THEN
    ! write(start_date_string, '(I0.4,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2)') &
    !   start_date(1), '-', start_date(2), '-', start_date(3), ' ', start_date(4), ':', start_date(5), ':', start_date(6)
    ! write(first_day_of_month_string, '(I0.4,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2)') &
    !   start_date(1), '-', start_date(2), '-', 1, ' ', 0, ':', 0, ':', 0
    !
    ! ! Initiate netcdf output
    call output_netcdf_init()
    !
    ! ! Initiate netcdf output and close since it is only written once
    call output_netcdf_general_init()
    call output_netcdf_general_done()
    !
    ! ! Initiate netcdf for meteorology output
    ! call output_netcdf_meteo_init()
    !
    ! ! Initiate netcdf for chemistry output
    ! call output_netcdf_chem_init()

    ! Initiate netcdf for aerosol output

    IF (flag_aero==1)then
      call output_netcdf_Tvapor_init()
      call output_netcdf_aerosol_init()
      call output_netcdf_SB_vapor_init()
    ! ELSEIF (flag_aero==0)then
    !   call output_netcdf_vapor_no_uhma_init()
    ENDIF
    ! Set current time record to 0 (not written yet)
    ! itrec = 0

    FIRST_TIME = .FALSE.
  END IF

  ! Update time record counter
  ! itrec = itrec + 1

  ! ! Write meteo data to the netcdf file
  ! write(*,*) 'Output netcdf meteo write'
  ! call output_netcdf_meteo_write()
  !
  ! ! Write chemistry data to the netcdf file
  ! write(*,*) 'Output netcdf chem write'
  ! call output_netcdf_chem_write()

  ! Write aerosol data to the netcdf file
  IF (flag_aero==1)then
    ! write(*,*) 'Output aerosol data with UHMA'
    call output_netcdf_aerosol_write()
  ! ELSEIF (flag_aero==0)then
  !   write(*,*) 'Output aerosol data without UHMA'
  !   call output_netcdf_no_uhma_write()
  ENDIF

  ENDIF
END SUBROUTINE output_write_old


!==============================================================================!
! output netcdf init subroutines
!==============================================================================!

subroutine output_netcdf_init()
  ! Initiate DIOS
  call dios_init(stat)
end subroutine output_netcdf_init


subroutine output_netcdf_general_init()
  !-----------------------------------------------------------!
  ! Init general nc file
  ! General information is saved
  !-----------------------------------------------------------!

  ! General file name
  nc_general%fname = trim(adjustl(output_dir))//'/general.nc'

  ! Create file
  call dios_create(trim(nc_general%fname), DIOS_NETCDF, DIOS_REPLACE, nc_general%fid, stat)

  ! Global attributes
  call dios_put_att(nc_general%fid, DIOS_GLOBAL, 'title', 'SOSAA general file', stat)
  call dios_put_att(nc_general%fid, DIOS_GLOBAL, 'number_of_height_levels', kz, stat)
  call dios_put_att(nc_general%fid, DIOS_GLOBAL, 'flag_emis', flag_emis , stat)
  call dios_put_att(nc_general%fid, DIOS_GLOBAL, 'flag_chem', flag_chem, stat)
  call dios_put_att(nc_general%fid, DIOS_GLOBAL, 'flag_gasdrydep', flag_gasdrydep, stat)
  call dios_put_att(nc_general%fid, DIOS_GLOBAL, 'flag_aero', flag_aero, stat)
  call dios_put_att(nc_general%fid, DIOS_GLOBAL, 'flag_emis_soil', flag_emis_soil, stat)
  call dios_put_att(nc_general%fid, DIOS_GLOBAL, 'nclo_end', nclo_end, stat)
  call dios_put_att(nc_general%fid, DIOS_GLOBAL, 'dt_mete', dt_mete, stat)
  call dios_put_att(nc_general%fid, DIOS_GLOBAL, 'dt_chem', dt_chem, stat)
  call dios_put_att(nc_general%fid, DIOS_GLOBAL, 'dt_aero', dt_aero, stat)
  call dios_put_att(nc_general%fid, DIOS_GLOBAL, 'dt_uhma', dt_uhma, stat)
  call dios_put_att(nc_general%fid, DIOS_GLOBAL, 'canopy_height', hc, stat)
  call dios_put_att(nc_general%fid, DIOS_GLOBAL, 'abl', abl, stat)
  call dios_put_att(nc_general%fid, DIOS_GLOBAL, 'LAI_tot', LAI_tot, stat)
  call dios_put_att(nc_general%fid, DIOS_GLOBAL, 'LAI_curv', LAI_curv, stat)
  call dios_put_att(nc_general%fid, DIOS_GLOBAL, 'LAI_proj', LAI_proj, stat)

  ! Define dimensions
  call dios_def_dim(nc_general%fid, 'lev', kz, nc_general%dimid_lev, stat)

  ! Define variables
  call dios_def_var(nc_general%fid, 'lev', DIOS_DOUBLE, (/nc_general%dimid_lev/), &
                    nc_general%varid_lev, stat)
  call dios_put_att(nc_general%fid, nc_general%varid_lev, 'unit', 'm', stat)
  call dios_put_att(nc_general%fid, nc_general%varid_lev, 'long_name', 'height above the ground', stat)

  call dios_def_var(nc_general%fid, 'lad_meteo', DIOS_DOUBLE, (/nc_general%dimid_lev/), &
                    varid_lad_meteo, stat)
  call dios_put_att(nc_general%fid, varid_lad_meteo, 'unit', 'm2 m-3', stat)
  call dios_put_att(nc_general%fid, varid_lad_meteo, &
    'long_name', 'leaf area density for meteorology module', stat)

  call dios_def_var(nc_general%fid, 'lad_megan', DIOS_DOUBLE, (/nc_general%dimid_lev/), &
                    varid_lad_megan, stat)
  call dios_put_att(nc_general%fid, varid_lad_megan, 'unit', 'm2 m-3', stat)
  call dios_put_att(nc_general%fid, varid_lad_megan, &
    'long_name', 'leaf area density for megan module', stat)

  ! End of definition
  call dios_enddef(nc_general%fid, stat)

  ! Put variables
  call dios_put_var(nc_general%fid, nc_general%varid_lev, z, stat, start=(/1/), count=(/kz/))
  call dios_put_var(nc_general%fid, varid_lad_meteo, s1, stat, start=(/1/), count=(/kz/))
  call dios_put_var(nc_general%fid, varid_lad_megan, s2, stat, start=(/1/), count=(/kz/))

end subroutine output_netcdf_general_init


subroutine output_netcdf_meteo_init()

  !----- Init meteo nc file ------------------------------------------------------!

  ! Meteo file name
  nc_meteo%fname = trim(adjustl(output_dir))//'/meteo.nc'

  ! Create file
  call dios_create(trim(nc_meteo%fname), DIOS_NETCDF, DIOS_REPLACE, nc_meteo%fid, stat)

  ! Global attributes
  call dios_put_att(nc_meteo%fid, DIOS_GLOBAL, 'title', 'SOSAA meteorology file', stat)

  ! Definition of dimensions
  call dios_def_dim(nc_meteo%fid, 'time', DIOS_UNLIMITED, nc_meteo%dimid_time, stat)
  call dios_def_dim(nc_meteo%fid, 'lev', kz, nc_meteo%dimid_lev, stat)
  nc_meteo%dimids_2d = (/nc_meteo%dimid_lev, nc_meteo%dimid_time/)

  !----- Definition of variables -----!
  ! Unlimited dimension should be defined in the last dimid (for classic NetCDF),
  ! notice that the dims are inversed when using 'ncdump' to check the data.

  ! Dimension variables
  call dios_def_var(nc_meteo%fid, 'lev', DIOS_DOUBLE, (/nc_meteo%dimid_lev/), nc_meteo%varid_lev, stat)
  call dios_put_att(nc_meteo%fid, nc_meteo%varid_lev, 'unit', 'm', stat)
  call dios_put_att(nc_meteo%fid, nc_meteo%varid_lev, 'long_name', 'height above the ground', stat)

  call dios_def_var(nc_meteo%fid, 'time', DIOS_DOUBLE, (/nc_meteo%dimid_time/), nc_meteo%varid_time, stat)
  call dios_put_att(nc_meteo%fid, nc_meteo%varid_time, 'unit', 's', stat)
  call dios_put_att(nc_meteo%fid, nc_meteo%varid_time, 'long_name', 'time since the beginning of month', stat)
  call dios_put_att(nc_meteo%fid, nc_meteo%varid_time, 'first_day_of_month', trim(first_day_of_month_string), stat)
  call dios_put_att(nc_meteo%fid, nc_meteo%varid_time, 'start_date', trim(start_date_string), stat)

  ! 2D variables
  call dios_def_var(nc_meteo%fid, 'temp', DIOS_DOUBLE, nc_meteo%dimids_2d, varid_temp, stat)
  call dios_put_att(nc_meteo%fid, varid_temp, 'unit', 'K', stat)
  call dios_put_att(nc_meteo%fid, varid_temp, 'long_name', 'air temperature', stat)

  call dios_def_var(nc_meteo%fid, 'rhov', DIOS_DOUBLE, nc_meteo%dimids_2d, varid_rhov, stat)
  call dios_put_att(nc_meteo%fid, varid_rhov, 'unit', 'kg m-3', stat)
  call dios_put_att(nc_meteo%fid, varid_rhov, 'long_name', 'absolute humidity', stat)

  call dios_def_var(nc_meteo%fid, 'rh', DIOS_DOUBLE, nc_meteo%dimids_2d, varid_rh, stat)
  call dios_put_att(nc_meteo%fid, varid_rh, 'unit', '-', stat)
  call dios_put_att(nc_meteo%fid, varid_rh, 'long_name', 'relative humidity', stat)

  call dios_def_var(nc_meteo%fid, 'pres', DIOS_DOUBLE, nc_meteo%dimids_2d, varid_pres, stat)
  call dios_put_att(nc_meteo%fid, varid_pres, 'unit', 'Pa', stat)
  call dios_put_att(nc_meteo%fid, varid_pres, 'long_name', 'air pressure', stat)

  call dios_def_var(nc_meteo%fid, 'ua', DIOS_DOUBLE, nc_meteo%dimids_2d, varid_ua, stat)
  call dios_put_att(nc_meteo%fid, varid_ua, 'unit', 'm s-1', stat)
  call dios_put_att(nc_meteo%fid, varid_ua, 'long_name', 'zonal wind', stat)

  call dios_def_var(nc_meteo%fid, 'va', DIOS_DOUBLE, nc_meteo%dimids_2d, varid_va, stat)
  call dios_put_att(nc_meteo%fid, varid_va, 'unit', 'm s-1', stat)
  call dios_put_att(nc_meteo%fid, varid_va, 'long_name', 'meridional wind', stat)

  call dios_def_var(nc_meteo%fid, 'tke', DIOS_DOUBLE, nc_meteo%dimids_2d, varid_tke, stat)
  call dios_put_att(nc_meteo%fid, varid_tke, 'unit', 'm2 s-2', stat)
  call dios_put_att(nc_meteo%fid, varid_tke, 'long_name', 'turbulent kinetic energy', stat)

  call dios_def_var(nc_meteo%fid, 'diffm', DIOS_DOUBLE, nc_meteo%dimids_2d, varid_diffm, stat)
  call dios_put_att(nc_meteo%fid, varid_diffm, 'unit', 'm2 s-1', stat)
  call dios_put_att(nc_meteo%fid, varid_diffm, 'long_name', 'eddy diffusivity of momentum', stat)

  call dios_def_var(nc_meteo%fid, 'diffh', DIOS_DOUBLE, nc_meteo%dimids_2d, varid_diffh, stat)
  call dios_put_att(nc_meteo%fid, varid_diffh, 'unit', 'm2 s-1', stat)
  call dios_put_att(nc_meteo%fid, varid_diffh, 'long_name', 'eddy diffusivity of scalar', stat)

  call dios_def_var(nc_meteo%fid, 'mixlen', DIOS_DOUBLE, nc_meteo%dimids_2d, varid_mixlen, stat)
  call dios_put_att(nc_meteo%fid, varid_mixlen, 'unit', 'm', stat)
  call dios_put_att(nc_meteo%fid, varid_mixlen, 'long_name', 'mixing length', stat)

  call dios_def_var(nc_meteo%fid, 'ri', DIOS_DOUBLE, nc_meteo%dimids_2d, varid_ri, stat)
  call dios_put_att(nc_meteo%fid, varid_ri, 'unit', '-', stat)
  call dios_put_att(nc_meteo%fid, varid_ri, 'long_name', 'Richardson number', stat)

  call dios_def_var(nc_meteo%fid, 'ustar', DIOS_DOUBLE, nc_meteo%dimids_2d, varid_ustar, stat)
  call dios_put_att(nc_meteo%fid, varid_ustar, 'unit', 'm s-1', stat)
  call dios_put_att(nc_meteo%fid, varid_ustar, 'long_name', 'friction velocity', stat)

  call dios_def_var(nc_meteo%fid, 'shf', DIOS_DOUBLE, nc_meteo%dimids_2d, varid_shf, stat)
  call dios_put_att(nc_meteo%fid, varid_shf, 'unit', 'W m-2', stat)
  call dios_put_att(nc_meteo%fid, varid_shf, 'long_name', 'sensible heat flux', stat)

  call dios_def_var(nc_meteo%fid, 'lhf', DIOS_DOUBLE, nc_meteo%dimids_2d, varid_lhf, stat)
  call dios_put_att(nc_meteo%fid, varid_lhf, 'unit', 'W m-2', stat)
  call dios_put_att(nc_meteo%fid, varid_lhf, 'long_name', 'latent heat flux', stat)

  call dios_def_var(nc_meteo%fid, 'nair', DIOS_DOUBLE, nc_meteo%dimids_2d, varid_nair, stat)
  call dios_put_att(nc_meteo%fid, varid_nair, 'unit', 'molec cm-3', stat)
  call dios_put_att(nc_meteo%fid, varid_nair, 'long_name', 'number concentration of air', stat)

  ! 1D variables
  call dios_def_var(nc_meteo%fid, 'ceb', DIOS_DOUBLE, (/nc_meteo%dimid_time/), varid_ceb, stat)
  call dios_put_att(nc_meteo%fid, varid_ceb, 'unit', 'W m-2', stat)
  call dios_put_att(nc_meteo%fid, varid_ceb, 'long_name', 'energy balance at canopy top', stat)

  call dios_def_var(nc_meteo%fid, 'gsoil', DIOS_DOUBLE, (/nc_meteo%dimid_time/), varid_gsoil, stat)
  call dios_put_att(nc_meteo%fid, varid_gsoil, 'unit', 'W m-2', stat)
  call dios_put_att(nc_meteo%fid, varid_gsoil, 'long_name', 'soil heat flux plus heat storage change in the top soil', stat)

  call dios_def_var(nc_meteo%fid, 'pblh', DIOS_DOUBLE, (/nc_meteo%dimid_time/), varid_pblh, stat)
  call dios_put_att(nc_meteo%fid, varid_pblh, 'unit', 'm', stat)
  call dios_put_att(nc_meteo%fid, varid_pblh, 'long_name', 'PBL height', stat)

  call dios_def_var(nc_meteo%fid, 'albedo', DIOS_DOUBLE, (/nc_meteo%dimid_time/), varid_albedo, stat)
  call dios_put_att(nc_meteo%fid, varid_albedo, 'unit', '-', stat)
  call dios_put_att(nc_meteo%fid, varid_albedo, 'long_name', 'cover albedo', stat)

  call dios_def_var(nc_meteo%fid, 'zenith', DIOS_DOUBLE, (/nc_meteo%dimid_time/), varid_zenith, stat)
  call dios_put_att(nc_meteo%fid, varid_zenith, 'unit', 'deg', stat)
  call dios_put_att(nc_meteo%fid, varid_zenith, 'long_name', 'solar zenith angle', stat)

  call dios_def_var(nc_meteo%fid, 'rd_dir', DIOS_DOUBLE, (/nc_meteo%dimid_time/), varid_rd_dir, stat)
  call dios_put_att(nc_meteo%fid, varid_rd_dir, 'unit', 'W m-2', stat)
  call dios_put_att(nc_meteo%fid, varid_rd_dir, 'long_name', 'downward direct solar radiation', stat)

  call dios_def_var(nc_meteo%fid, 'rd_dif', DIOS_DOUBLE, (/nc_meteo%dimid_time/), varid_rd_dif, stat)
  call dios_put_att(nc_meteo%fid, varid_rd_dif, 'unit', 'W m-2', stat)
  call dios_put_att(nc_meteo%fid, varid_rd_dif, 'long_name', 'downward diffuse solar radiation', stat)

  call dios_def_var(nc_meteo%fid, 'ld', DIOS_DOUBLE, (/nc_meteo%dimid_time/), varid_ld, stat)
  call dios_put_att(nc_meteo%fid, varid_ld, 'unit', 'W m-2', stat)
  call dios_put_att(nc_meteo%fid, varid_ld, 'long_name', 'downward long wave radiation', stat)

  call dios_def_var(nc_meteo%fid, 'lu', DIOS_DOUBLE, (/nc_meteo%dimid_time/), varid_lu, stat)
  call dios_put_att(nc_meteo%fid, varid_lu, 'unit', 'W m-2', stat)
  call dios_put_att(nc_meteo%fid, varid_lu, 'long_name', 'upward long wave radiation', stat)

  call dios_def_var(nc_meteo%fid, 'paru', DIOS_DOUBLE, (/nc_meteo%dimid_time/), varid_paru, stat)
  call dios_put_att(nc_meteo%fid, varid_paru, 'unit', 'W m-2', stat)
  call dios_put_att(nc_meteo%fid, varid_paru, 'long_name', 'upward PAR radiation', stat)

  call dios_def_var(nc_meteo%fid, 'niru', DIOS_DOUBLE, (/nc_meteo%dimid_time/), varid_niru, stat)
  call dios_put_att(nc_meteo%fid, varid_niru, 'unit', 'W m-2', stat)
  call dios_put_att(nc_meteo%fid, varid_niru, 'long_name', 'upward near-infrared radiation', stat)

  ! End of definition
  call dios_enddef(nc_meteo%fid, stat)

  ! Write limited dimension variables
  call dios_put_var(nc_meteo%fid, nc_meteo%varid_lev, z, stat, start=(/1/), count=(/kz/))
end subroutine output_netcdf_meteo_init


subroutine output_netcdf_chem_init
  CHARACTER(LEN=2)    :: Jstamp
  ! Chemistry file name
  nc_chem%fname = trim(adjustl(output_dir))//'/chem.nc'

  ! Create file
  call dios_create(trim(nc_chem%fname), DIOS_NETCDF, DIOS_REPLACE, nc_chem%fid, stat)

  ! Global attributes
  call dios_put_att(nc_chem%fid, DIOS_GLOBAL, 'title', 'SOSAA chemisty file', stat)

  ! Definition of dimensions
  call dios_def_dim(nc_chem%fid, 'time', DIOS_UNLIMITED, nc_chem%dimid_time, stat)
  call dios_def_dim(nc_chem%fid, 'lev', kz, nc_chem%dimid_lev, stat)
  ! call dios_def_dim(nc_chem%fid, 'ispc', outspc_count, nc_chem%dimid_ispc, stat)
  ! call dios_def_dim(nc_chem%fid, 'iemi', outemi_count, nc_chem%dimid_iemi, stat)
  ! nc_chem%dimids_emi2d = (/nc_chem%dimid_iemi, nc_chem%dimid_time/)
  nc_chem%dimids_2d = (/nc_chem%dimid_lev, nc_chem%dimid_time/)

  !----- Definition of variables -----!
  ! Unlimited dimension should be defined in the last dimid (for classic NetCDF),
  ! notice that the dims are inversed when using 'ncdump' to check the data.

  ! Dimension variables
  call dios_def_var(nc_chem%fid, 'lev', DIOS_DOUBLE, (/nc_chem%dimid_lev/), nc_chem%varid_lev, stat)
  call dios_put_att(nc_chem%fid, nc_chem%varid_lev, 'unit', 'm', stat)
  call dios_put_att(nc_chem%fid, nc_chem%varid_lev, 'long_name', 'height above the ground', stat)

  call dios_def_var(nc_chem%fid, 'time', DIOS_DOUBLE, (/nc_chem%dimid_time/), nc_chem%varid_time, stat)
  call dios_put_att(nc_chem%fid, nc_chem%varid_time, 'unit', 's', stat)
  call dios_put_att(nc_chem%fid, nc_chem%varid_time, 'long_name', 'time since the beginning of month', stat)
  call dios_put_att(nc_chem%fid, nc_chem%varid_time, 'first_day_of_month', trim(first_day_of_month_string), stat)
  call dios_put_att(nc_chem%fid, nc_chem%varid_time, 'start_date', trim(start_date_string), stat)

  ! 2D variables

  ! Number concentration
  ALLOCATE(varid_outspc_nconc(outspc_count))
  DO i=1, outspc_count
    if (outspc_cheminds(i) > 0) then  ! only for existed species
      call dios_def_var(nc_chem%fid, 'nconc_'//trim(adjustl(outspc_names(i))), &
        DIOS_DOUBLE, nc_chem%dimids_2d, varid_outspc_nconc(i), stat)
      call dios_put_att(nc_chem%fid, varid_outspc_nconc(i), 'unit', 'molec cm-3', stat)
      call dios_put_att(nc_chem%fid, varid_outspc_nconc(i), 'long_name', 'number concentration', stat)
    end if
  END DO

  ! Vertical flux
  ALLOCATE(varid_outspc_flux(outflux_count))
  DO i=1, outflux_count
    if (outflux_cheminds(i) > 0) then  ! only for existed species
      call dios_def_var(nc_chem%fid, 'flux_'//trim(adjustl(outflux_names(i))), &
        DIOS_DOUBLE, nc_chem%dimids_2d, varid_outspc_flux(i), stat)
      call dios_put_att(nc_chem%fid, varid_outspc_flux(i), 'unit', 'molec cm-2 s-1', stat)
      call dios_put_att(nc_chem%fid, varid_outspc_flux(i), 'long_name', 'vertical flux of molecules', stat)
    end if
  END DO

  ! Gas dry deposition velocity
  ALLOCATE(varid_outspc_vd(outvd_count))
  DO i=1, outvd_count
    if (outvd_cheminds(i) > 0) then  ! only for existed species
      call dios_def_var(nc_chem%fid, 'vd_'//trim(adjustl(outvd_names(i))), &
        DIOS_DOUBLE, nc_chem%dimids_2d, varid_outspc_vd(i), stat)
      call dios_put_att(nc_chem%fid, varid_outspc_vd(i), 'unit', 'm s-1', stat)
      call dios_put_att(nc_chem%fid, varid_outspc_vd(i), 'long_name', 'dry deposition velocity', stat)
    end if
  END DO

  ! Emission rates
  ALLOCATE(varid_outemi(outemi_count))
  DO i=1, outemi_count
    if (outemi_meganinds(i) > 0) then  ! only for existed emission species
      call dios_def_var(nc_chem%fid, 'emirate_'//trim(adjustl(outemi_names(i))), &
        DIOS_DOUBLE, nc_chem%dimids_2d, varid_outemi(i), stat)
      call dios_put_att(nc_chem%fid, varid_outemi(i), 'unit', 'molec cm-3 s-1', stat)
      call dios_put_att(nc_chem%fid, varid_outemi(i), 'long_name', 'emission rate', stat)
    end if
  END DO
  ! OH reactivity
  ALLOCATE(varid_outspc_rOH(CH_oh_count))
  DO i=1,CH_oh_count
    if (CH_oh_indices(i) > 0) then  ! only for existed species
      call dios_def_var(nc_chem%fid, trim(adjustl(SPC_names(CH_oh_indices(i)))), &
        DIOS_DOUBLE, nc_chem%dimids_2d, varid_outspc_rOH(i), stat)
      call dios_put_att(nc_chem%fid, varid_outspc_rOH(i), 'unit', 's-1', stat)
      call dios_put_att(nc_chem%fid, varid_outspc_rOH(i), 'long_name', 'NO3 reactivity', stat)
    end if
  END DO
  ! O3 reactivity
  ALLOCATE(varid_outspc_rO3(CH_o3_count))
  DO i=1,CH_o3_count
    if (CH_o3_indices(i) > 0) then  ! only for existed species
      call dios_def_var(nc_chem%fid, trim(adjustl(SPC_names(CH_o3_indices(i)))), &
        DIOS_DOUBLE, nc_chem%dimids_2d, varid_outspc_rO3(i), stat)
      call dios_put_att(nc_chem%fid, varid_outspc_rO3(i), 'unit', 's-1', stat)
      call dios_put_att(nc_chem%fid, varid_outspc_rO3(i), 'long_name', 'O3 reactivity', stat)
    end if
  END DO
  !NO3 reactivity
  ALLOCATE(varid_outspc_rNO3(CH_no3_count))
  DO i=1,CH_no3_count
    if (CH_no3_indices(i) > 0) then  ! only for existed species
      call dios_def_var(nc_chem%fid, trim(adjustl(SPC_names(CH_no3_indices(i)))), &
        DIOS_DOUBLE, nc_chem%dimids_2d, varid_outspc_rNO3(i), stat)
      call dios_put_att(nc_chem%fid, varid_outspc_rNO3(i), 'unit', 's-1', stat)
      call dios_put_att(nc_chem%fid, varid_outspc_rNO3(i), 'long_name', 'NO3 reactivity', stat)
    end if
  END DO
  !Reactivity of gases with significant atmospheric implication
  ALLOCATE(varid_outspc_reactivity(CH_reactivity_count))
  DO i=1,CH_reactivity_count
    if (CH_reactivity_indices(i) > 0) then  ! only for existed species
      call dios_def_var(nc_chem%fid, trim(adjustl(SPC_names(CH_reactivity_indices(i)))), &
        DIOS_DOUBLE, nc_chem%dimids_2d, varid_outspc_reactivity(i), stat)
      call dios_put_att(nc_chem%fid, varid_outspc_reactivity(i), 'unit', 's-1', stat)
      call dios_put_att(nc_chem%fid, varid_outspc_reactivity(i), 'long_name', trim(adjustl(SPC_names(CH_reactivity_indices(i))))//' reactivity', stat)
    end if
  END DO
  ! Photolysis rate of all photochemical reactions
  ALLOCATE(varid_J(NPHOT))
  DO i=1,NPHOT
     write(Jstamp,'(I2)') I
      call dios_def_var(nc_chem%fid, 'J'//Jstamp, DIOS_DOUBLE, nc_chem%dimids_2d, varid_J(i), stat)
      call dios_put_att(nc_chem%fid, varid_J(i), 'unit', 's-1', stat)
      call dios_put_att(nc_chem%fid, varid_J(i), 'long_name', 'Photolysis rate of reaction NO.'//Jstamp, stat)
  END DO
  ! End of definition
  call dios_enddef(nc_chem%fid, stat)

  ! Write limited dimension variables
  call dios_put_var(nc_chem%fid, nc_chem%varid_lev, z, stat, start=(/1/), count=(/kz/))
end subroutine output_netcdf_chem_init


subroutine output_netcdf_aerosol_init
  ! Aerosol by size bin file name
  nc_aerosol%fname = trim(adjustl(output_dir))//'/aerosol.nc'

  ! Create file
  call dios_create(trim(nc_aerosol%fname), DIOS_NETCDF, DIOS_REPLACE, nc_aerosol%fid, stat)

  ! Global attributes
  call dios_put_att(nc_aerosol%fid, DIOS_GLOBAL, 'title', 'SOSAA aerosol file', stat)

  ! Definition of dimensions
  call dios_def_dim(nc_aerosol%fid, 'dp_dry_fs'  , n_bins_par , nc_aerosol%dimid_dp_dry_fs  , stat)
  call dios_def_dim(nc_aerosol%fid, 'lev'   , kz            , nc_aerosol%dimid_lev   , stat)
  call dios_def_dim(nc_aerosol%fid, 'time'  , DIOS_UNLIMITED, nc_aerosol%dimid_time  , stat)
  nc_aerosol%dimids_3d_dp_dry_fs  = (/nc_aerosol%dimid_lev ,nc_aerosol%dimid_dp_dry_fs  ,  nc_aerosol%dimid_time/)
  nc_aerosol%dimids_2d       = (/nc_aerosol%dimid_lev ,nc_aerosol%dimid_time/)
  !----- Definition of variables -----!
  ! Unlimited dimension should be defined in the last dimid (for classic NetCDF),
  ! notice that the dims are inversed when using 'ncdump' to check the data.

  ! Dimension variables

  call dios_def_var(nc_aerosol%fid, 'dp_dry_fs', DIOS_DOUBLE, (/nc_aerosol%dimid_dp_dry_fs/), nc_aerosol%varid_dp_dry_fs, stat)
  call dios_put_att(nc_aerosol%fid, nc_aerosol%varid_dp_dry_fs, 'unit', 'm', stat)
  call dios_put_att(nc_aerosol%fid, nc_aerosol%varid_dp_dry_fs, 'long_name', 'dry radius of aerosol particles in each size bin', stat)

  call dios_def_var(nc_aerosol%fid, 'lev', DIOS_DOUBLE, (/nc_aerosol%dimid_lev/), nc_aerosol%varid_lev, stat)
  call dios_put_att(nc_aerosol%fid, nc_aerosol%varid_lev, 'unit', 'm', stat)
  call dios_put_att(nc_aerosol%fid, nc_aerosol%varid_lev, 'long_name', 'height above the ground', stat)

  call dios_def_var(nc_aerosol%fid, 'time', DIOS_DOUBLE, (/nc_aerosol%dimid_time/), nc_aerosol%varid_time, stat)
  call dios_put_att(nc_aerosol%fid, nc_aerosol%varid_time, 'unit', 's', stat)
  call dios_put_att(nc_aerosol%fid, nc_aerosol%varid_time, 'long_name', 'time since the beginning of month', stat)
  call dios_put_att(nc_aerosol%fid, nc_aerosol%varid_time, 'first_day_of_month', trim(first_day_of_month_string), stat)
  call dios_put_att(nc_aerosol%fid, nc_aerosol%varid_time, 'start_date', trim(start_date_string), stat)

  ! 3D variables

  ! Number concentration of particles
  call dios_def_var(nc_aerosol%fid, 'nconc_particle', DIOS_DOUBLE, nc_aerosol%dimids_3d_dp_dry_fs, varid_nconc_particle, stat)
  call dios_put_att(nc_aerosol%fid, varid_nconc_particle, 'unit', 'particles m-3', stat)
  call dios_put_att(nc_aerosol%fid, varid_nconc_particle, 'long_name', 'number concentration of particles', stat)

  ! Vertical fluxes of particles
  call dios_def_var(nc_aerosol%fid, 'flux_particle', DIOS_DOUBLE, nc_aerosol%dimids_3d_dp_dry_fs, varid_flux_particle, stat)
  call dios_put_att(nc_aerosol%fid, varid_flux_particle, 'unit', 'particles m-2 s-1', stat)
  call dios_put_att(nc_aerosol%fid, varid_flux_particle, 'long_name', 'vertical flux of particles', stat)

  ! Growth rate of particles
  call dios_def_var(nc_aerosol%fid, 'growth_rate', DIOS_DOUBLE, nc_aerosol%dimids_3d_dp_dry_fs, varid_growth_rate, stat)
  call dios_put_att(nc_aerosol%fid, varid_growth_rate, 'unit', 'm s-1', stat)  ! ?? check the unit later
  call dios_put_att(nc_aerosol%fid, varid_growth_rate, 'long_name', 'growth rate of particles', stat)

  ! Nucleation rate of particles
  call dios_def_var(nc_aerosol%fid, 'nuc_rate', DIOS_DOUBLE, nc_aerosol%dimids_2d, varid_nuc_rate, stat)
  call dios_put_att(nc_aerosol%fid, varid_nuc_rate, 'unit', 'm s-1', stat)  ! ?? check the unit later
  call dios_put_att(nc_aerosol%fid, varid_nuc_rate, 'long_name', 'nucleation rate of particles', stat)

  ! End of definition
  call dios_enddef(nc_aerosol%fid, stat)

  ! Write limited dimension variables
  call dios_put_var(nc_aerosol%fid, nc_aerosol%varid_dp_dry_fs, CURRENT_PSD%dp_dry_fs, stat, start=(/1/), count=(/n_bins_par/))
  call dios_put_var(nc_aerosol%fid, nc_aerosol%varid_lev, z, stat, start=(/1/), count=(/kz/))
end subroutine output_netcdf_aerosol_init

subroutine output_netcdf_TVapor_init
  ! Total vapor file name
  nc_Tvapor%fname = trim(adjustl(output_dir))//'/Vapor_total.nc'

  ! Create file
  call dios_create(trim(nc_Tvapor%fname), DIOS_NETCDF, DIOS_REPLACE, nc_Tvapor%fid, stat)

  ! Global attributes
  call dios_put_att(nc_Tvapor%fid, DIOS_GLOBAL, 'title', 'SOSAA Vapor total concentration and sink file', stat)

  ! Definition of dimensions

  call dios_def_dim(nc_Tvapor%fid, 'lev'   , kz            , nc_Tvapor%dimid_lev   , stat)

  call dios_def_dim(nc_Tvapor%fid, 'time'  , DIOS_UNLIMITED, nc_Tvapor%dimid_time  , stat)

  call dios_def_var(nc_Tvapor%fid, 'lev', DIOS_DOUBLE, (/nc_Tvapor%dimid_lev/), nc_Tvapor%varid_lev, stat)
  call dios_put_att(nc_Tvapor%fid, nc_Tvapor%varid_lev, 'unit', 'm', stat)
  call dios_put_att(nc_Tvapor%fid, nc_Tvapor%varid_lev, 'long_name', 'height above the ground', stat)

  call dios_def_var(nc_Tvapor%fid, 'time', DIOS_DOUBLE, (/nc_Tvapor%dimid_time/), nc_Tvapor%varid_time, stat)
  call dios_put_att(nc_Tvapor%fid, nc_Tvapor%varid_time, 'unit', 's', stat)
  call dios_put_att(nc_Tvapor%fid, nc_Tvapor%varid_time, 'long_name', 'time since the beginning of month', stat)
  call dios_put_att(nc_Tvapor%fid, nc_Tvapor%varid_time, 'first_day_of_month', trim(first_day_of_month_string), stat)
  call dios_put_att(nc_Tvapor%fid, nc_Tvapor%varid_time, 'start_date', trim(start_date_string), stat)


  ! 2D variables
  nc_Tvapor%dimids_2d_condv = (/nc_Tvapor%dimid_lev , nc_Tvapor%dimid_time/)

IF (flag_vapor==0)THEN

  call dios_def_dim(nc_Tvapor%fid, 'icondv', outvap_count     , nc_Tvapor%dimid_icondv, stat)
  ! 3D variables
  nc_Tvapor%dimids_3d_condv = (/nc_Tvapor%dimid_lev ,nc_Tvapor%dimid_icondv,  nc_Tvapor%dimid_time/)
  call dios_def_var(nc_Tvapor%fid, 'icondv', DIOS_DOUBLE, (/nc_Tvapor%dimid_icondv/), nc_Tvapor%varid_icondv, stat)
  call dios_put_att(nc_Tvapor%fid, nc_Tvapor%varid_icondv, 'unit', '-', stat)
  call dios_put_att(nc_Tvapor%fid, nc_Tvapor%varid_icondv, 'long_name', 'indices of condensable vapors', stat)

  call dios_def_var(nc_Tvapor%fid, 'nconc_condv_gas', DIOS_DOUBLE, nc_Tvapor%dimids_3d_condv, varid_nconc_condv_gas, stat)

  call dios_put_att(nc_Tvapor%fid, varid_nconc_condv_gas, 'unit', 'molec cm-3', stat)  ! ?? check the unit later
  call dios_put_att(nc_Tvapor%fid, varid_nconc_condv_gas, 'long_name', 'number concentrations of condensable vapors in gas phase', stat)

  call dios_def_var(nc_Tvapor%fid, 'sink_condv_gas', DIOS_DOUBLE, nc_Tvapor%dimids_3d_condv, varid_sink_condv_gas, stat)

  call dios_put_att(nc_Tvapor%fid, varid_sink_condv_gas, 'unit', 'molec cm-3', stat)  ! ?? check the unit later
  call dios_put_att(nc_Tvapor%fid, varid_sink_condv_gas, 'long_name', 'sink of condensable vapors in gas phase', stat)
  ! End of definition
  call dios_enddef(nc_Tvapor%fid, stat)

  ! Write limited dimension variables
  call dios_put_var(nc_Tvapor%fid, nc_Tvapor%varid_icondv, (/ (I, I = 1, outvap_count) /), stat, start=(/1/), count=(/outvap_count/))
ELSEIF(flag_vapor==1)THEN

  ALLOCATE(varid_nconc_condv_gas_one(outvap_count))
  ALLOCATE(varid_sink_condv_gas_one(outvap_count))
  DO i=1, outvap_count
    if (outvap_vapinds(i) > 0) then  ! only for existed species
  ! 2D Variables
  ! Number concentrations of condensable vapors in particle phase in all layers
  call dios_def_var(nc_Tvapor%fid, 'nconc_condv_gas'//trim(adjustl(outvap_names(i))), DIOS_DOUBLE, nc_Tvapor%dimids_2d_condv, varid_nconc_condv_gas_one(i), stat)

  call dios_put_att(nc_Tvapor%fid, varid_nconc_condv_gas_one(i), 'unit', 'molec cm-3', stat)  ! ?? check the unit later
  call dios_put_att(nc_Tvapor%fid, varid_nconc_condv_gas_one(i), 'long_name', 'number concentrations of condensable vapors in gas phase', stat)
  !
  call dios_def_var(nc_Tvapor%fid, 'sink_condv_gas'//trim(adjustl(outvap_names(i))), DIOS_DOUBLE, nc_Tvapor%dimids_2d_condv, varid_sink_condv_gas_one(i), stat)

  call dios_put_att(nc_Tvapor%fid, varid_sink_condv_gas_one(i), 'unit', 'molec cm-3', stat)  ! ?? check the unit later
  call dios_put_att(nc_Tvapor%fid, varid_sink_condv_gas_one(i), 'long_name', 'sink of condensable vapors in gas phase', stat)
    end if
  END DO
  ! End of definition
  call dios_enddef(nc_Tvapor%fid, stat)

END IF
  ! Write limited dimension variables
  call dios_put_var(nc_Tvapor%fid, nc_Tvapor%varid_lev, z, stat, start=(/1/), count=(/kz/))

end subroutine output_netcdf_TVapor_init

subroutine output_netcdf_SB_Vapor_init
  ! Vapor by size bin file name
  nc_SBvapor%fname = trim(adjustl(output_dir))//'/Vapor_SB.nc'

  ! Create file
  call dios_create(trim(nc_SBvapor%fname), DIOS_NETCDF, DIOS_REPLACE, nc_SBvapor%fid, stat)

  ! Global attributes
  call dios_put_att(nc_SBvapor%fid, DIOS_GLOBAL, 'title', 'SOSAA Vapor concentration by size distribution file', stat)

  ! Definition of dimensions
  call dios_def_dim(nc_SBvapor%fid, 'dp_dry_fs'  , n_bins_par , nc_SBvapor%dimid_dp_dry_fs  , stat)
  call dios_def_dim(nc_SBvapor%fid, 'lev'   , kz            , nc_SBvapor%dimid_lev   , stat)
  call dios_def_dim(nc_SBvapor%fid, 'time'  , DIOS_UNLIMITED, nc_SBvapor%dimid_time  , stat)



  !----- Definition of variables -----!
  ! Unlimited dimension should be defined in the last dimid (for classic NetCDF),
  ! notice that the dims are inversed when using 'ncdump' to check the data.

  ! Dimension variables

  call dios_def_var(nc_SBvapor%fid, 'dp_dry_fs', DIOS_DOUBLE, (/nc_SBvapor%dimid_dp_dry_fs/), nc_SBvapor%varid_dp_dry_fs, stat)
  call dios_put_att(nc_SBvapor%fid, nc_SBvapor%varid_dp_dry_fs, 'unit', 'm', stat)
  call dios_put_att(nc_SBvapor%fid, nc_SBvapor%varid_dp_dry_fs, 'long_name', 'dry radius of aerosol particles in each size bin', stat)

  call dios_def_var(nc_SBvapor%fid, 'lev', DIOS_DOUBLE, (/nc_SBvapor%dimid_lev/), nc_SBvapor%varid_lev, stat)
  call dios_put_att(nc_SBvapor%fid, nc_SBvapor%varid_lev, 'unit', 'm', stat)
  call dios_put_att(nc_SBvapor%fid, nc_SBvapor%varid_lev, 'long_name', 'height above the ground', stat)

  call dios_def_var(nc_SBvapor%fid, 'time', DIOS_DOUBLE, (/nc_SBvapor%dimid_time/), nc_SBvapor%varid_time, stat)
  call dios_put_att(nc_SBvapor%fid, nc_SBvapor%varid_time, 'unit', 's', stat)
  call dios_put_att(nc_SBvapor%fid, nc_SBvapor%varid_time, 'long_name', 'time since the beginning of month', stat)
  call dios_put_att(nc_SBvapor%fid, nc_SBvapor%varid_time, 'first_day_of_month', trim(first_day_of_month_string), stat)
  call dios_put_att(nc_SBvapor%fid, nc_SBvapor%varid_time, 'start_date', trim(start_date_string), stat)


  call dios_put_var(nc_SBvapor%fid, nc_SBvapor%varid_dp_dry_fs, CURRENT_PSD%dp_dry_fs, stat, start=(/1/), count=(/n_bins_par/))
  call dios_put_var(nc_SBvapor%fid, nc_SBvapor%varid_lev, z, stat, start=(/1/), count=(/kz/))
IF (flag_vapor==0)THEN
  call dios_def_dim(nc_SBvapor%fid, 'icondv', outvap_count     , nc_SBvapor%dimid_icondv, stat)
  ! 4D variables
  ! Dimension variables
  call dios_def_var(nc_SBvapor%fid, 'icondv', DIOS_DOUBLE, (/nc_SBvapor%dimid_icondv/), nc_SBvapor%varid_icondv, stat)
  call dios_put_att(nc_SBvapor%fid, nc_SBvapor%varid_icondv, 'unit', '-', stat)
  call dios_put_att(nc_SBvapor%fid, nc_SBvapor%varid_icondv, 'long_name', 'indices of condensable vapors', stat)

  nc_SBvapor%dimids_4d       = (/nc_SBvapor%dimid_lev, nc_SBvapor%dimid_dp_dry_fs, nc_SBvapor%dimid_icondv, nc_SBvapor%dimid_time/)

  ! Number concentrations of condensable vapors in particle phase in all layers
  call dios_def_var(nc_SBvapor%fid, 'vconc_condv_gas', DIOS_DOUBLE, nc_SBvapor%dimids_4d, varid_vconc_condv_par, stat)
  call dios_put_att(nc_SBvapor%fid, varid_vconc_condv_par, 'unit', 'um3 m-3', stat)  ! ?? check the unit later
  call dios_put_att(nc_SBvapor%fid, varid_vconc_condv_par, 'long_name', 'volume concentrations of condensable vapors in particle phase', stat)
  ! End of definition
  call dios_enddef(nc_SBvapor%fid, stat)
  ! Write limited dimension variables
  call dios_put_var(nc_SBvapor%fid, nc_SBvapor%varid_icondv, (/ (I, I = 1, outvap_count) /), stat, start=(/1/), count=(/outvap_count/))

ELSEIF(flag_vapor==1)THEN
  ! 3D variables
  ! Definition of dimensions
  nc_SBvapor%dimids_3d       = (/nc_SBvapor%dimid_lev, nc_SBvapor%dimid_dp_dry_fs, nc_SBvapor%dimid_time/)

  ALLOCATE(varid_vconc_condv_par_one(outvap_count))

  DO i=1, outvap_count
    if (outvap_vapinds(i) > 0) then  ! only for existed species
  ! Number concentrations of condensable vapors in particle phase in all layers
  call dios_def_var(nc_SBvapor%fid, 'vconc_condv_par'//trim(adjustl(outvap_names(i))), DIOS_DOUBLE, nc_SBvapor%dimids_3d, varid_vconc_condv_par_one(i), stat)
  call dios_put_att(nc_SBvapor%fid, varid_vconc_condv_par_one(i), 'unit', 'molec cm-3', stat)  ! ?? check the unit later
  call dios_put_att(nc_SBvapor%fid, varid_vconc_condv_par_one(i), 'long_name', 'number concentrations of condensable vapors in gas phase', stat)
    end if
  END DO
  ! End of definition
  call dios_enddef(nc_SBvapor%fid, stat)

END IF
  ! Write limited dimension variables
  call dios_put_var(nc_SBvapor%fid, nc_SBvapor%varid_lev, z, stat, start=(/1/), count=(/kz/))

end subroutine output_netcdf_SB_Vapor_init

subroutine output_netcdf_vapor_no_uhma_init
  ! Total vapor file name
  nc_Tvapor%fname = trim(adjustl(output_dir))//'/Vapor_total.nc'

  ! Create file
  call dios_create(trim(nc_Tvapor%fname), DIOS_NETCDF, DIOS_REPLACE, nc_Tvapor%fid, stat)

  ! Global attributes
  call dios_put_att(nc_Tvapor%fid, DIOS_GLOBAL, 'title', 'SOSAA Vapor total concentration and sink file with no UHMA', stat)

  ! Definition of dimensions

  call dios_def_dim(nc_Tvapor%fid, 'lev'   , kz            , nc_Tvapor%dimid_lev   , stat)

  call dios_def_dim(nc_Tvapor%fid, 'time'  , DIOS_UNLIMITED, nc_Tvapor%dimid_time  , stat)


  call dios_def_var(nc_Tvapor%fid, 'lev', DIOS_DOUBLE, (/nc_Tvapor%dimid_lev/), nc_Tvapor%varid_lev, stat)
  call dios_put_att(nc_Tvapor%fid, nc_Tvapor%varid_lev, 'unit', 'm', stat)
  call dios_put_att(nc_Tvapor%fid, nc_Tvapor%varid_lev, 'long_name', 'height above the ground', stat)

  call dios_def_var(nc_Tvapor%fid, 'time', DIOS_DOUBLE, (/nc_Tvapor%dimid_time/), nc_Tvapor%varid_time, stat)
  call dios_put_att(nc_Tvapor%fid, nc_Tvapor%varid_time, 'unit', 's', stat)
  call dios_put_att(nc_Tvapor%fid, nc_Tvapor%varid_time, 'long_name', 'time since the beginning of month', stat)
  call dios_put_att(nc_Tvapor%fid, nc_Tvapor%varid_time, 'first_day_of_month', trim(first_day_of_month_string), stat)
  call dios_put_att(nc_Tvapor%fid, nc_Tvapor%varid_time, 'start_date', trim(start_date_string), stat)


  ! 2D variables
  nc_Tvapor%dimids_2d_condv = (/nc_Tvapor%dimid_lev , nc_Tvapor%dimid_time/)

IF (flag_vapor==0)THEN

  call dios_def_dim(nc_Tvapor%fid, 'icondv', outvap_count     , nc_Tvapor%dimid_icondv, stat)
  ! 3D variables
  nc_Tvapor%dimids_3d_condv = (/nc_Tvapor%dimid_lev ,nc_Tvapor%dimid_icondv,  nc_Tvapor%dimid_time/)
  call dios_def_var(nc_Tvapor%fid, 'icondv', DIOS_DOUBLE, (/nc_Tvapor%dimid_icondv/), nc_Tvapor%varid_icondv, stat)
  call dios_put_att(nc_Tvapor%fid, nc_Tvapor%varid_icondv, 'unit', '-', stat)
  call dios_put_att(nc_Tvapor%fid, nc_Tvapor%varid_icondv, 'long_name', 'indices of condensable vapors', stat)

  call dios_def_var(nc_Tvapor%fid, 'nconc_condv_gas', DIOS_DOUBLE, nc_Tvapor%dimids_3d_condv, varid_nconc_condv_gas, stat)

  call dios_put_att(nc_Tvapor%fid, varid_nconc_condv_gas, 'unit', 'molec cm-3', stat)  ! ?? check the unit later
  call dios_put_att(nc_Tvapor%fid, varid_nconc_condv_gas, 'long_name', 'number concentrations of condensables in gas phase', stat)
  ! Write limited dimension variables
  call dios_put_var(nc_Tvapor%fid, nc_Tvapor%varid_icondv, (/ (I, I = 1, outvap_count) /), stat, start=(/1/), count=(/outvap_count/))
ELSEIF(flag_vapor==1)THEN

  ! 2D variables
  ALLOCATE(varid_nconc_condv_gas_one(outvap_count))
  DO i=1, outvap_count
    if (outvap_vapinds(i) > 0) then  ! only for existed species
  ! Number concentrations of condensable vapors in particle phase in all layers
  call dios_def_var(nc_Tvapor%fid, 'nconc_condv_gas'//trim(adjustl(outvap_names(i))), DIOS_DOUBLE, nc_Tvapor%dimids_2d_condv, varid_nconc_condv_gas_one(i), stat)

  call dios_put_att(nc_Tvapor%fid, varid_nconc_condv_gas_one(i), 'unit', 'molec cm-3', stat)  ! ?? check the unit later
  call dios_put_att(nc_Tvapor%fid, varid_nconc_condv_gas_one(i), 'long_name', 'number concentrations of '//trim(adjustl(outvap_names(i)))//' in gas phase', stat)
    end if
  ENDDO

ENDIF
  ALLOCATE(varid_sink_condv_gas_one(1))

  call dios_def_var(nc_Tvapor%fid, 'CS_H2SO4', DIOS_DOUBLE, nc_Tvapor%dimids_2d_condv, varid_sink_condv_gas_one(1), stat)

  call dios_put_att(nc_Tvapor%fid, varid_sink_condv_gas_one(1), 'unit', 'molec cm-3', stat)  ! ?? check the unit later
  call dios_put_att(nc_Tvapor%fid, varid_sink_condv_gas_one(1), 'long_name', 'sink of H2SO4', stat)
  ! End of definition
  call dios_enddef(nc_Tvapor%fid, stat)
  ! Write limited dimension variables
  call dios_put_var(nc_Tvapor%fid, nc_Tvapor%varid_lev, z, stat, start=(/1/), count=(/kz/))
end subroutine output_netcdf_vapor_no_uhma_init
!==============================================================================!
! output netcdf write subroutines
!==============================================================================!

subroutine output_netcdf_meteo_write()
  !----- Write meteo nc file ------------------------------------------------------!

  call dios_put_var(nc_meteo%fid, nc_meteo%varid_time, (/time_in_month/), stat, start=(/itrec/), count=(/1/))

  ! 2D variables
  call dios_put_var(nc_meteo%fid, varid_temp, ta1, stat, start=(/1, itrec/), count=(/kz, 1/))
  call dios_put_var(nc_meteo%fid, varid_rhov, qa1, stat, start=(/1, itrec/), count=(/kz, 1/))
  call dios_put_var(nc_meteo%fid, varid_rh, RH, stat, start=(/1, itrec/), count=(/kz, 1/))
  call dios_put_var(nc_meteo%fid, varid_pres, pres, stat, start=(/1, itrec/), count=(/kz, 1/))
  call dios_put_var(nc_meteo%fid, varid_ua, u1, stat, start=(/1, itrec/), count=(/kz, 1/))
  call dios_put_var(nc_meteo%fid, varid_va, v1, stat, start=(/1, itrec/), count=(/kz, 1/))
  call dios_put_var(nc_meteo%fid, varid_tke, bt1, stat, start=(/1, itrec/), count=(/kz, 1/))
  call dios_put_var(nc_meteo%fid, varid_mixlen, l1, stat, start=(/1, itrec/), count=(/kz, 1/))
  call dios_put_var(nc_meteo%fid, varid_diffm, kt1, stat, start=(/1, itrec/), count=(/kz, 1/))
  call dios_put_var(nc_meteo%fid, varid_diffh, alt1*kt1, stat, start=(/1, itrec/), count=(/kz, 1/))
  call dios_put_var(nc_meteo%fid, varid_ri, rih, stat, start=(/1, itrec/), count=(/kz, 1/))
  call dios_put_var(nc_meteo%fid, varid_ustar, ur, stat, start=(/1, itrec/), count=(/kz, 1/))
  call dios_put_var(nc_meteo%fid, varid_shf, fluxh3, stat, start=(/1, itrec/), count=(/kz, 1/))
  call dios_put_var(nc_meteo%fid, varid_lhf, fluxle3, stat, start=(/1, itrec/), count=(/kz, 1/))
  call dios_put_var(nc_meteo%fid, varid_nair, air, stat, start=(/1, itrec/), count=(/kz, 1/))

  ! 1D variables
  call dios_put_var(nc_meteo%fid, varid_ceb, (/balans/), stat, start=(/itrec/), count=(/1/))
  call dios_put_var(nc_meteo%fid, varid_gsoil, (/pp/), stat, start=(/itrec/), count=(/1/))
  call dios_put_var(nc_meteo%fid, varid_pblh, (/pblh/), stat, start=(/itrec/), count=(/1/))
  call dios_put_var(nc_meteo%fid, varid_albedo, (/albedo_f/), stat, start=(/itrec/), count=(/1/))
  call dios_put_var(nc_meteo%fid, varid_zenith, (/zenith_deg/), stat, start=(/itrec/), count=(/1/))
  call dios_put_var(nc_meteo%fid, varid_rd_dir, (/rsnt/), stat, start=(/itrec/), count=(/1/))
  call dios_put_var(nc_meteo%fid, varid_rd_dif, (/rskt/), stat, start=(/itrec/), count=(/1/))
  call dios_put_var(nc_meteo%fid, varid_ld, (/fird(k_canopy)/), stat, start=(/itrec/), count=(/1/))
  call dios_put_var(nc_meteo%fid, varid_lu, (/firu(k_canopy)/), stat, start=(/itrec/), count=(/1/))
  call dios_put_var(nc_meteo%fid, varid_paru, (/fphu(k_canopy)/), stat, start=(/itrec/), count=(/1/))
  call dios_put_var(nc_meteo%fid, varid_niru, (/fniu(k_canopy)/), stat, start=(/itrec/), count=(/1/))
end subroutine output_netcdf_meteo_write


subroutine output_netcdf_chem_write()
  ! Time
  write(*,*) 'time'
  call dios_put_var(nc_chem%fid, nc_chem%varid_time, (/time_in_month/), stat, start=(/itrec/), count=(/1/))

  write(*,*) 'outspc', outspc_count, nc_chem%fid, varid_outspc_nconc
  DO i=1, outspc_count
      ! Number concentration
      call dios_put_var(nc_chem%fid, varid_outspc_nconc(i), CH_CONS_ALL(:, outspc_cheminds(i)), stat, &
        start=(/1, itrec/), count=(/kz, 1/))
  END DO
  write(*,*) 'outflux'
  DO i=1, outflux_count
      ! Vertical flux
      call dios_put_var(nc_chem%fid, varid_outspc_flux(i), CH_CONS_FLUX(:, outflux_cheminds(i)), stat, &
        start=(/1, itrec/), count=(/kz, 1/))
  END DO
  write(*,*) 'outvd'
  DO i=1, outvd_count
      ! Gas dry deposition
      call dios_put_var(nc_chem%fid, varid_outspc_vd(i), vdep(:, outvd_cheminds(i)), stat, &
        start=(/1, itrec/), count=(/kz, 1/))
  END DO

  ! Emission rate
  write(*,*) 'outemi'
  DO i=1, outemi_count
    if (outemi_meganinds(i) > 0) then  ! only for existed emission species
      call dios_put_var(nc_chem%fid, varid_outemi(i), EM_EMI(:, outemi_meganinds(i)), stat, &
        start=(/1, itrec/), count=(/kz, 1/))
    end if
  END DO

  ! Reactivity
  write(*,*) 'reactivity'
  DO i=1, CH_oh_count
    if (CH_oh_indices(i) > 0) then  ! only for existed emission species
      call dios_put_var(nc_chem%fid, varid_outspc_rOH(i), CH_oh_cons3(:,i,3), stat, &
        start=(/1, itrec/), count=(/kz, 1/))
    end if
  END DO

  DO i=1, CH_o3_count
    if (CH_o3_indices(i) > 0) then  ! only for existed emission species
      call dios_put_var(nc_chem%fid, varid_outspc_rO3(i), CH_o3_cons3(:,i,3), stat, &
        start=(/1, itrec/), count=(/kz, 1/))
    end if
  END DO

  DO i=1, CH_no3_count
    if (CH_no3_indices(i) > 0) then  ! only for existed emission species
      call dios_put_var(nc_chem%fid, varid_outspc_rNO3(i), CH_no3_cons3(:,i,3), stat, &
        start=(/1, itrec/), count=(/kz, 1/))
    end if
  END DO

  DO i=1, CH_reactivity_count
    if (CH_reactivity_indices(i) > 0) then  ! only for existed emission species
      call dios_put_var(nc_chem%fid, varid_outspc_reactivity(i), CH_reactivity_cons3(:,i,3), stat, &
        start=(/1, itrec/), count=(/kz, 1/))
    end if
  END DO

  DO i=1, NPHOT
       call dios_put_var(nc_chem%fid, varid_J(i), CH_J_values_ALL(:,i), stat, &
        start=(/1, itrec/), count=(/kz, 1/))
  END DO
end subroutine output_netcdf_chem_write

! subroutine output_netcdf_no_uhma_write()
!
!   ! Time
!   call dios_put_var(nc_Tvapor%fid, nc_Tvapor%varid_time, (/time_in_month/), stat, start=(/itrec/), count=(/1/))
!
! IF (flag_vapor==0)THEN
!   ! Number concentration of condensable vapors in gas phase
!   call dios_put_var(nc_Tvapor%fid, varid_nconc_condv_gas, VAPOR(:, outvap_vapinds), stat, &
!     start=(/1, 1, itrec/), count=(/ kz, outvap_count, 1/))
! ELSEIF (flag_vapor==1)THEN
!   DO i=1, outvap_count
!   call dios_put_var(nc_Tvapor%fid, varid_nconc_condv_gas_one(i), VAPOR(:, outvap_vapinds(i)), stat, &
!     start=(/1, itrec/), count=(/ kz, 1/))
!   END DO
! ENDIF
!   call dios_put_var(nc_Tvapor%fid, varid_sink_condv_gas_one(1), CH_cs_H2SO4, stat, &
!     start=(/1, itrec/), count=(/ kz, 1/))
! end subroutine output_netcdf_no_uhma_write

subroutine output_netcdf_aerosol_write()

  ! Time
  call dios_put_var(nc_aerosol%fid, nc_aerosol%varid_time, (/time_in_month/), stat, start=(/itrec/), count=(/1/))
  call dios_put_var(nc_SBvapor%fid, nc_SBvapor%varid_time, (/time_in_month/), stat, start=(/itrec/), count=(/1/))
  call dios_put_var(nc_Tvapor%fid, nc_Tvapor%varid_time, (/time_in_month/), stat, start=(/itrec/), count=(/1/))
  ! Number concentration
  call dios_put_var(nc_aerosol%fid, varid_nconc_particle, N_CONC(:, :), stat, &
    start=(/1, 1, itrec/), count=(/ kz, n_bins_par, 1/))

  ! Vertical flux
  call dios_put_var(nc_aerosol%fid, varid_flux_particle, PAR_FLUX(:, :), stat, &
    start=(/1, 1, itrec/), count=(/ kz, n_bins_par, 1/))

  ! Growth rate
  call dios_put_var(nc_aerosol%fid, varid_growth_rate, GR(:, :), stat, &
    start=(/1, 1, itrec/), count=(/ kz,n_bins_par, 1/))

  ! Nucleation rate
  call dios_put_var(nc_aerosol%fid, varid_nuc_rate, NUC_RATE( :), stat, &
    start=(/1, itrec/), count=(/ kz, 1/))

IF (flag_vapor==0)THEN
  ! mass concentration of condensable vapors in particle phase [kg/particle]
  call dios_put_var(nc_SBvapor%fid, varid_vconc_condv_par, MASS_COMPO(:, :, outvap_vapinds), stat, &
    start=(/1, 1, 1, itrec/), count=(/kz, n_bins_par, outvap_count, 1/))
ELSEIF (flag_vapor==1)THEN
  DO i=1, outvap_count
  call dios_put_var(nc_SBvapor%fid, varid_vconc_condv_par_one(i), MASS_COMPO(:, :, outvap_vapinds(i)), stat, &
    start=(/1, 1, itrec/), count=(/kz, n_bins_par, 1/))
  END DO
END IF

IF (flag_vapor==0)THEN
  ! Number concentration of condensable vapors in gas phase
  call dios_put_var(nc_Tvapor%fid, varid_nconc_condv_gas, VAPOR(:, outvap_vapinds), stat, &
    start=(/1, 1, itrec/), count=(/ kz, outvap_count, 1/))
  call dios_put_var(nc_Tvapor%fid, varid_sink_condv_gas, SINK(:, outvap_vapinds), stat, &
    start=(/1, 1, itrec/), count=(/ kz, outvap_count, 1/))
ELSEIF (flag_vapor==1)THEN
  DO i=1, outvap_count
  call dios_put_var(nc_Tvapor%fid, varid_nconc_condv_gas_one(i), VAPOR(:, outvap_vapinds(i)), stat, &
    start=(/1, itrec/), count=(/ kz, 1/))
  call dios_put_var(nc_Tvapor%fid, varid_sink_condv_gas_one(i), SINK(:, outvap_vapinds(i)), stat, &
    start=(/1, itrec/), count=(/ kz, 1/))
  END DO
ENDIF
end subroutine output_netcdf_aerosol_write


!==============================================================================!
! output netcdf done subroutines
!==============================================================================!

subroutine output_netcdf_general_done()
  call dios_close(nc_general%fid, stat)
end subroutine output_netcdf_general_done


subroutine output_netcdf_meteo_done()
  call dios_close(nc_meteo%fid, stat)
end subroutine output_netcdf_meteo_done


subroutine output_netcdf_chem_done()
  call dios_close(nc_chem%fid, stat)
end subroutine output_netcdf_chem_done


subroutine output_netcdf_aerosol_done()
  call dios_close(nc_aerosol%fid, stat)
end subroutine output_netcdf_aerosol_done

subroutine output_netcdf_TVapor_done()
  call dios_close(nc_TVapor%fid, stat)
end subroutine output_netcdf_TVapor_done

subroutine output_netcdf_SBVapor_done()
  call dios_close(nc_SBVapor%fid, stat)
end subroutine output_netcdf_SBVapor_done

subroutine output_netcdf_done()
  ! call output_netcdf_meteo_done()
  ! call output_netcdf_chem_done()
IF (flag_aero==1)then
  call output_netcdf_aerosol_done()
  call output_netcdf_SBVapor_done()
ENDIF
  call output_netcdf_TVapor_done()
  ! call dios_done(stat)
end subroutine output_netcdf_done


subroutine output_done()
  ! call output_netcdf_done()
  call dios_close(mixf%funit, stat)
  call dios_done(stat)
end subroutine output_done


!------------------------------------------------------------------------------!
!
! String subprograms
!
!------------------------------------------------------------------------------!


!==============================================================================!
! Convert all the letters to upper case
!==============================================================================!
FUNCTION upper(s1)  RESULT (s2)
  CHARACTER(*)       :: s1
  CHARACTER(LEN(s1)) :: s2
  CHARACTER          :: ch
  INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
  INTEGER            :: i

  DO i = 1,LEN(s1)
    ch = s1(i:i)
    IF (ch >= 'a'.AND.ch <= 'z') ch = CHAR(ICHAR(ch)+DUC)
    s2(i:i) = ch
  END DO
END FUNCTION upper


!==============================================================================!
! Convert all the letters to lower case
!==============================================================================!
FUNCTION lower(s1)  RESULT (s2)
  CHARACTER(*)       :: s1
  CHARACTER(LEN(s1)) :: s2
  CHARACTER          :: ch
  INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
  INTEGER            :: i

  DO i = 1,LEN(s1)
    ch = s1(i:i)
    IF (ch >= 'A'.AND.ch <= 'Z') ch = CHAR(ICHAR(ch)-DUC)
    s2(i:i) = ch
  END DO
END FUNCTION Lower


!==============================================================================!
! Split the string to substrings separated by wall.
! It does not include the rooms with only spaces.
!==============================================================================!
subroutine string_split_trim2(str, wall, nrooms, irooms)
  ! The rooms are not trimmed

  character(*), intent(in) :: str
  character   , intent(in) :: wall

  integer, intent(out) :: nrooms
  integer, intent(out) :: irooms(:, :)

  integer :: nwall
  integer :: str_len
  integer, allocatable :: wall_pos(:)

  integer :: i, ir

  str_len = len(str)

  ! maximum rooms separated by walls
  allocate(wall_pos(str_len))
  call count_substring(str, wall, .false., nwall, wall_pos)

  ir = 0

  ! the room before first wall
  if (len_trim( adjustl( str(1:wall_pos(1)-1) ) ) /= 0) then
    ir = ir + 1
    irooms(ir, 1) = 1
    irooms(ir, 2) = wall_pos(1)-1
  end if

  ! rooms between wall 1 and nwall-1
  do i=1, nwall-1
    if (len_trim( adjustl( str(wall_pos(i)+1:wall_pos(i+1)-1) ) ) /= 0) then
      ir = ir + 1
      irooms(ir, 1) = wall_pos(i)+1
      irooms(ir, 2) = wall_pos(i+1)-1
    end if
  end do

  ! the room after nwall
  if (len_trim( adjustl( str(wall_pos(nwall)+1:) ) ) /= 0) then
    ir = ir + 1
    irooms(ir, 1) = wall_pos(nwall)+1
    irooms(ir, 2) = str_len
  end if

  nrooms = ir
end subroutine string_split_trim2


!==============================================================================!
! Split the string to substrings separated by wall.
! It includes the rooms with only spaces.
!==============================================================================!
subroutine string_split(str, wall, nrooms, irooms)
  character(*), intent(in) :: str
  character   , intent(in) :: wall

  integer, intent(out) :: nrooms
  integer, intent(out) :: irooms(:, :)

  integer :: nwall
  integer :: str_len
  integer, allocatable :: wall_pos(:)

  integer :: i, ir

  str_len = len(str)

  ! maximum rooms separated by walls
  allocate(wall_pos(str_len))
  call count_substring(str, wall, .false., nwall, wall_pos)

  ir = 0
  ! the room before first wall
  if (len(str(1:wall_pos(1)-1)) /= 0) then
    ir = ir + 1
    irooms(ir, 1) = 1
    irooms(ir, 2) = wall_pos(1)-1
  end if

  ! rooms between wall 1 and nwall-1
  do i=1, nwall-1
    if (len(str(wall_pos(i)+1:wall_pos(i+1)-1)) /= 0) then
      ir = ir + 1
      irooms(ir, 1) = wall_pos(i)+1
      irooms(ir, 2) = wall_pos(i+1)-1
    end if
  end do

  ! the room after nwall
  if (len(str(wall_pos(nwall)+1:)) /= 0) then
    ir = ir + 1
    irooms(ir, 1) = wall_pos(nwall)+1
    irooms(ir, 2) = str_len
  end if

  nrooms = ir
end subroutine string_split


!==============================================================================!
! Count how many substrings are in the string.
!==============================================================================!
subroutine count_substring(s0, ss, greedy, c, pa)
  ! input raw string s0 and substring ss
  character(*), intent(in) :: s0, ss

  ! greedy, e.g. 'oo' in 'oooo'
  ! true : c = 3, pa = [1, 2, 3]
  ! false: c = 2, pa = [1, 3]
  logical, intent(in) :: greedy

  ! count of ss in s0
  integer, intent(out) :: c

  ! position array: starting indices of every ss in s0
  ! dimension size should be the same as the length of s0, but can be reduced if necessary.
  ! it depends on the input array
  ! set to zero when the position is not used
  integer, intent(out) :: pa(:)

  ! current position of s0
  integer :: p

  ! the first occurence of substring
  integer :: posn

  ! initiate
  c = 0   ! count is 0
  pa = 0  ! indices of ss are 0 (does not exist)

  ! raw string is empty, return 0
  if (len(s0) == 0) return

  ! current position
  p = 1

  ! search until not found ss in s0's substrings
  do
    ! search the first occurence of ss in the rest of s0
    posn = index(s0(p:), ss)

    ! return current condition if not found any more
    if (posn == 0) return

    ! add one to counter if found
    c = c + 1

    ! save the position
    pa(c) = p + posn - 1

    ! advance 1 if greedy
    if (greedy) then
      p = p + posn - 1 + 1
    ! skip all ss if not greedy
    else
      p = p + posn - 1 + len(ss)
    end if
  end do
end subroutine count_substring


!==============================================================================!
! Get substring from the position array
!==============================================================================!
SUBROUTINE get_substring_list(s0, nrooms, irooms, rooms)
  character(*), intent(in) :: s0
  integer, intent(in) :: nrooms
  integer, intent(in) :: irooms(:, :)

  character(*), intent(out) :: rooms(:)

  integer :: i

  do i = 1, nrooms
    rooms(i) = trim(adjustl( s0(irooms(i, 1):irooms(i, 2)) ))
  end do
END SUBROUTINE get_substring_list


!==============================================================================!
! Get the indices of names in the original name list
!==============================================================================!
subroutine get_inds_outnames_in_rawnames(outnames, rawnames, outinds)
  character(*), intent(in) :: outnames(:)
  character(*), intent(in) :: rawnames(:)

  integer, intent(out)     :: outinds(:)

  integer :: i

  ! Obtain indices of output species (case insensative) from raw name list
  outinds = -1  ! default value
  do i = 1, size(outnames)
    do j=1, size(rawnames)
      if ( upper( trim(adjustl( rawnames(j) )) ) == upper( trim(adjustl( outnames(i) )) ) ) then
        outinds(i) = j
        exit
      end if
    end do
  end do
end subroutine get_inds_outnames_in_rawnames

    ! subroutine vapor_no_uhma(filename, vapor_names,parameter_A, parameter_B)
    !   character(len=*), intent(in) :: filename  ! for input_flag == 0, it is the folder of condensable_vapors
    !   character(len=60) :: chem_name
    !   character(len=60), dimension(VAPOUR_PROP%n_cond_tot)  :: vapor_names
    !   real(dp):: molecular_mass
    !   real(dp),  DIMENSION(VAPOUR_PROP%n_cond_tot)::parameter_A, parameter_B
    !   !List of all chemical properties is available at start, thus we can
    !   !allocate everything based on the file (and we don't need to change the
    !   ! array size afterwards).
    !   integer :: i, io, io2,  j, n
    !   character(len=80) :: line
    !   if(.true.) then
    !     ! Qi:input_flag = 0: use Pontus Model's method and file
    !     ! 1: use old version UHMA's file. change it in uhma_datatypes.f90
    !     open(100, file=trim(adjustl(filename)) // '/Vapour_properties.dat',action='read')
    !     open(101, file=trim(adjustl(filename)) // '/Vapour_names.dat',action='read')
    !     i=0
    !     do
    !       read(100,*, iostat=io) line
    !       if(io < 0) then
    !         exit
    !       end if
    !       i=i+1
    !     end do
    !
    !     if(io < 0) then
    !       rewind(100)
    !     else
    !       write (*,*) io
    !       stop 'File read error'
    !     end if
    !
    !     do j=1,i-1
    !       read(100,*,iostat=io)  molecular_mass, parameter_A(j), parameter_B(j)
    !       read(101,*,iostat=io2) chem_name
    !       vapor_names(j)=TRIM(chem_name)
    !     enddo
    !
    !     close(100)
    !     close(100)
    !     vapor_names(i) = 'dummy'
    !     vapor_names(i+1) = 'H2SO4'
    !     VAPOUR_PROP%ind_H2SO4 = i+1  !index for sulfuric acid
    !   endif
    ! end subroutine vapor_no_uhma


  ! This will modify the public variable infield
  subroutine add_input_field(infield, id, vname, fabspath, fid, dimlens, time_id)
    type(input_field), intent(inout) :: infield

    integer     , intent(in) :: id
    character(*), intent(in) :: vname
    character(*), intent(in) :: fabspath

    integer, intent(in) :: fid
    integer, intent(in) :: dimlens(:)
    integer, intent(in), optional :: time_id
    logical :: FLIP
    integer :: varid, k

    ! Nullify all the data fields
    nullify(infield%var(id)%f1d, infield%var(id)%f2d, infield%var(id)%f3d)

    ! Set meta data
    infield%var(id)%vname     = vname
    infield%var(id)%fabspath  = fabspath

    ! Get the varid from the iput file
    call dios_inq_varid(fid, trim(adjustl(vname)), varid, stat)

    infield%var(id)%ndims = size(dimlens)

    FLIP = .false.
    if (PRESENT(time_id)) Then
        ! check that array is not reversed, in SOSAA we excpect time to move forward
        k = size(infield%var(time_id)%f1d,1)
        if (infield%var(time_id)%f1d(1)>infield%var(time_id)%f1d(k)) THEN
            do j=1, infield%var(id)%ndims
                if (k==dimlens(j)) THEN
                    FLIP = .true.
                    exit
                end if
            end do
            if (.not.FLIP) CALL debug_message(TRIM(vname)//' TERRRRIBLE MISTAKE IS ABOUT TO HAPPEN!')
        END IF
    end if

    ! Allocate the needed field
    select case (infield%var(id)%ndims)
    case (1)
      allocate( infield%var(id)%f1d(dimlens(1)) )
      call dios_get_var(fid, varid, infield%var(id)%f1d(:), stat)
      if (FLIP) THEN
          CALL debug_message("I'm flipping 1d variable "//TRIM(vname))
          infield%var(id)%f1d(:) = infield%var(id)%f1d(dimlens(1):1:-1)
      END IF
    case (2)
      allocate( infield%var(id)%f2d(dimlens(1), dimlens(2)) )
      call dios_get_var(fid, varid, infield%var(id)%f2d(:, :), stat)
      if (FLIP) THEN
          CALL debug_message("I'm flipping 2d variable "//TRIM(vname))
          IF (dimlens(1) == k) infield%var(id)%f2d(:, :) = infield%var(id)%f2d(dimlens(1):1:-1, :)
          IF (dimlens(2) == k) infield%var(id)%f2d(:, :) = infield%var(id)%f2d(:, dimlens(2):1:-1)
      END IF
    case (3)
      allocate( infield%var(id)%f3d(dimlens(1), dimlens(2), dimlens(3)) )
      call dios_get_var(fid, varid, infield%var(id)%f3d(:, :, :), stat)
      if (FLIP) THEN
          CALL debug_message("I'm flipping 3d variable "//TRIM(vname))
          IF (dimlens(1) == k) infield%var(id)%f3d(:, :, :) = infield%var(id)%f3d(dimlens(1):1:-1, :,:)
          IF (dimlens(2) == k) infield%var(id)%f3d(:, :, :) = infield%var(id)%f3d(:, dimlens(2):1:-1,:)
          IF (dimlens(3) == k) infield%var(id)%f3d(:, :, :) = infield%var(id)%f3d(:, :, dimlens(3):1:-1)
      END IF
    end select


  end subroutine add_input_field


  subroutine add_input_field_id(id, nid)
    integer, intent(out) :: id
    integer, intent(out) :: nid  ! number of ids

    integer, save :: current_id
    logical, save :: first_call = .true.

    if (first_call) then
      current_id = 0
      first_call = .false.
    end if

    current_id = current_id + 1
    id  = current_id
    nid = current_id
  end subroutine add_input_field_id


  subroutine add_variable(mixf, var_name, long_name, var_unit, dimids, cwt, n_out)
    ! Which file
    type(mixfile), intent(in) :: mixf

    ! Input metafields
    character(*), intent(in) :: var_name, long_name, var_unit

    ! The dimension shape, (ndim1, -1) for 1d variables and (ndim1, ndim2) for
    ! 2d variables
    ! integer, intent(in) :: dim_shape(2)

    ! Dimension ids, including time dimension
    integer, dimension(:), intent(in) :: dimids

    ! If the variable changes with time
    logical, intent(in) :: cwt  ! change with time

    ! Current counter of 1d or 2d variables, the return value will be used as
    ! the variable index
    integer, intent(out) :: n_out

    integer :: ndim1, ndim2
    integer :: stat

    ! write(*,*) 'Shape of dimids is: ', size(dimids)

    if (size(dimids) == 1) then
      if (dimids(1) < 0) then  ! scalar variable
        ! Variable counter plus 1
        n_0d = n_0d + 1
        if (n_0d > NMAX_0D) then
          write(*,*) 'Out of maximum number of output variables.'
          stat = 1
        end if

        ! Set metafields
        mixf%f0d(n_0d)%mf = metafields(var_name, long_name, var_unit)

        ! Set dimension ids
        mixf%f0d(n_0d)%dimids = -1
        mixf%f0d(n_0d)%ndims = -1

        ! Initate the field value
        mixf%f0d(n_0d)%field = 0.0d0  ! initiate values to zero

        ! Set the change with time parameter
        mixf%f0d(n_0d)%change_with_time = cwt

        ! Return the counter as the index
        n_out = n_0d
      else  ! 1d variable
        ! Variable counter plus 1
        n_1d = n_1d + 1
        if (n_1d > NMAX_1D) then
          write(*,*) 'Out of maximum number of output variables.'
          stat = 1
        end if

        ! Set metafields
        mixf%f1d(n_1d)%mf = metafields(var_name, long_name, var_unit)

        ! Set dimension ids
        mixf%f1d(n_1d)%dimids = dimids

        ! Get the length of dimension
        call dios_inquire_dimension(mixf%funit, dimids(1), stat, length=ndim1)
        mixf%f1d(n_1d)%ndims = (/ndim1/)

        ! Allocate the field
        allocate( mixf%f1d(n_1d)%field(1:ndim1) )

        ! Initate the field value
        mixf%f1d(n_1d)%field(:) = 0.0d0  ! initiate values to zero

        ! Set the change with time parameter
        mixf%f1d(n_1d)%change_with_time = cwt

        ! Return the counter as the index
        n_out = n_1d
      end if
    else
      ! Variable counter plus 1
      n_2d = n_2d + 1

      ! Set metafields
      mixf%f2d(n_2d)%mf = metafields(var_name, long_name, var_unit)

      ! Set dimension ids
      mixf%f2d(n_2d)%dimids = dimids

      ! Get the length of dimension
      call dios_inquire_dimension(mixf%funit, dimids(1), stat, length=ndim1)
      call dios_inquire_dimension(mixf%funit, dimids(2), stat, length=ndim2)
      mixf%f2d(n_2d)%ndims = (/ndim1, ndim2/)

      ! Allocate the field
      allocate( mixf%f2d(n_2d)%field(1:ndim1, 1:ndim2) )

      ! Initate the field value
      mixf%f2d(n_2d)%field(:, :) = 0.0d0  ! initiate values to zero

      ! Set the change with time parameter
      mixf%f2d(n_2d)%change_with_time = cwt

      ! Return the counter as the index
      n_out = n_2d
    end if
  end subroutine add_variable

END MODULE SOSA_IO
