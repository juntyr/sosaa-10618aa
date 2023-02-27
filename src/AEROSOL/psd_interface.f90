module psd_interface

  use psd_constants
  use psd_scheme
  use psd_aerosol_dynamics
  use indexfromname_mod
  use filecount_mod
  use formatting
  use satconc_mod
  use multimodal_mod



  implicit none

  REAL(dp), ALLOCATABLE :: nominal_dp(:)        ! [m] array with nominal diameters. Stays constant independent of PSD_style
  REAL(dp), ALLOCATABLE :: J_distr(:)           ! Transient variable for distributing ACDC particles
  REAL(dp), ALLOCATABLE :: d_dpar(:)            ! array reporting relative changes to the diameter array within a single timestep
  REAL(dp), ALLOCATABLE :: d_npar(:)            ! array reporting relative changes to the particle number array within a single timestep
  REAL(dp), ALLOCATABLE :: d_npdep(:)           ! array reporting relative changes to the particle number array within a single timestep
  REAL(dp), ALLOCATABLE :: d_vap(:)             ! array reporting relative changes to the vapour concentration array within a single timestep
  REAL(dp)              :: refdvap  = 0d0       ! holds the reference maximum of d_vap, which is the average of three largest or three
  integer               :: irefdvap = 1         !   :... smallest (if negative) in d_vap, whichever is larger in magnitude. irefvap is the index of
  logical               :: logdvap  = .false.   !   :... the chosen maximum or minimum value. logdvap is true if abs(refdvap) is over the precision limits
  logical               :: update_K = .true.    ! is coagulation coefficient to be updated? In the main loop TRUE only if Temp or Press differs enough from the prev time step.
  REAL(dp)              :: A_vert                              ! chamber vertical wall area [m²]
  REAL(dp)              :: Vol_chamber                         ! chamber volume [m³]
  REAL(dp)              :: E_field = 0d0                       ! chamber Electric field (not implemented yet) [V/m]
  REAL(dp), ALLOCATABLE :: conc_fit(:)      ! An array giving particle conc independant of PSD_style [m⁻³]
  REAL(dp), allocatable :: conc_vapour(:)
  INTEGER               :: dmps_sp_min = 0, dmps_sp_max = 0 ! Indices for dmps_special
  REAL(dp), ALLOCATABLE :: losses_fit(:)    ! An array giving particle conc independant of PSD_style [m-3]
  REAL(dp), ALLOCATABLE :: save_measured(:)    ! An array giving particle conc independant of PSD_style [m-3]


  ! MPI send buffer and receive buffer have different lengths; some variables are not updated in the AERO module
  ! and some are not known beforehand (sink, GR). To use less named indices, the send and receive indices are the
  ! same, but the receive buffer has its indices starting later, not from 1
  ! some of the indices can be defined as parameters, others have to be defined later in aerosol_init, to
  ! be more flexible with the number of vapours

  !                     send buffer:
  !           x x x x x x x x x x x x x x x x x x x
  ! buf size: <----------------------------------->
  ! same as
  ! in chem:  o o o (time1, time2, kz)
  !                     recv buffer:
  ! buf size:               <----------------------------------------------------->
  !                         x x x x x x x x x x x x x x x x x x x x x x x x x x x x
  ! common indices:         | | | | | | | | | | | |


  ! only going in
  INTEGER, PARAMETER :: ind_aer_ts    = 4 ! timestep in aero
  INTEGER, PARAMETER :: ind_tempk     = 5 ! temperature at atmosphere level k
  INTEGER, PARAMETER :: ind_pres      = 6 ! pressure at atmosphere level k
  INTEGER, PARAMETER :: ind_rh        = 7 ! RH at atmosphere level k
  INTEGER, PARAMETER :: ind_U10m      = 8 ! RH at atmosphere level k
  INTEGER, PARAMETER :: ind_lsm       = 9 ! RH at atmosphere level k
  INTEGER, PARAMETER :: ind_firstcall = 10  ! is this the first call of the module?
  INTEGER, PARAMETER :: ind_IPR       = 11 ! ion production rate -> ACDC

  ! going in and out, this is where the mpi_recv_aer_buffer gets its first index
  INTEGER, PARAMETER :: ind_kz_tag    = 12          ! last parameter, beginning of return buffer
  INTEGER :: ind1_emiaer,        ind2_emiaer        ! length will come from of emission bins count
  INTEGER :: ind1_seasalt,       ind2_seasalt       ! length will come from n_bins_par
  INTEGER :: ind1_n_conc,        ind2_n_conc        ! length will come from n_bins_par
  INTEGER :: ind1_composition,   ind2_composition   ! length will come from n_bins_par*n_cond_tot
  INTEGER :: ind1_CONS_aer,      ind2_CONS_aer      ! length will come from NSPEC
  INTEGER :: ind1_CONS_clust,    ind2_CONS_clust    ! length will come from sum(G_ACDC(:)%neq_syst)

  ! only coming back
  INTEGER :: ind1_acdc_rates,    ind2_acdc_rates    ! length is 4 (total, neut, pos, neg)
  INTEGER :: ind1_oth_nuc_rates, ind2_oth_nuc_rates ! length varies, for now 1
  INTEGER :: ind1_gr,            ind2_gr            ! length will come from n_bins_par
  INTEGER :: ind1_sink,          ind2_sink          ! length will come from n_cond_tot

  INTEGER :: new_mpi_send_aer_buf_size              ! these are dynamically defined
  INTEGER :: new_mpi_recv_aer_buf_size              ! these are dynamically defined
  REAL(dp), ALLOCATABLE :: new_mpi_send_aer_buf(:)
  REAL(dp), ALLOCATABLE :: new_mpi_recv_aer_buf(:)

  REAL(dp), ALLOCATABLE :: CH_RES_org(:)            ! ??? Currently not used in the chemistry
  ! INTEGER,  ALLOCATABLE :: index_vapor(:)           ! redundant, in VAPOUR_PROP??
  ! INTEGER,  ALLOCATABLE :: sat_vap(:)               ! redundant, in VAPOUR_PROP?
  ! character(len=60), ALLOCATABLE :: vapor_names(:)
  real(dp), ALLOCATABLE :: dp_emi_aer(:,:) ! the upper and lower boundaries of emitted particles, (n_bins_emis,2)

contains

    SUBROUTINE ACDC_INIT(SPC_NAMES_, stage)

        CHARACTER(len=16), INTENT(IN)     :: SPC_NAMES_(:)
        INTEGER :: stage

        if (stage==1) THEN
            ALLOCATE(G_ACDC(nr_of_acdc_modules))
            ALLOCATE(ACDC_SYSTEMS(nr_of_acdc_modules))
            ALLOCATE(ACDC_links(nr_of_acdc_modules))
            ! These are just default values for backwards compatibility and get overrided if the used has defined them.
            ACDC_SYSTEMS        = 0
            ACDC_SYSTEMS(1:2)   = 1
            ACDC_links(1)       = 'A H2SO4 N NH3'
            ACDC_links(2)       = 'A H2SO4 D DMA'
        else
            CALL PARSE_ACDC_SYSTEMS(SPC_NAMES_)
        end if

    END SUBROUTINE ACDC_INIT

    SUBROUTINE AERO_INIT(NSPEC, kz)
        INTEGER :: ii
        INTEGER, intent(IN) :: NSPEC, kz

        CALL PARSE_MULTIMODAL

        ! Initialize the Particle representation
        CALL INITIALIZE_PSD

        ! going in and out, this is where the mpi_recv_aer_buffer gets its first index
        ind1_emiaer = ind_kz_tag + 1
        ind2_emiaer = ind1_emiaer + 9 - 1

        ind1_seasalt = ind2_emiaer + 1
        ind2_seasalt = ind1_seasalt + n_bins_par - 1

        ind1_n_conc = ind2_seasalt + 1
        ind2_n_conc = ind1_n_conc + n_bins_par - 1

        ind1_composition = ind2_n_conc + 1
        ind2_composition = ind1_composition + n_bins_par * n_cond_tot - 1

        ind1_CONS_aer = ind2_composition + 1
        ind2_CONS_aer = ind1_CONS_aer + NSPEC - 1

        ind1_CONS_clust = ind2_CONS_aer + 1
        ind2_CONS_clust = ind1_CONS_clust + (sum(G_ACDC(:)%neq_syst)) - 1

        ! only coming back
        ind1_acdc_rates = ind2_CONS_clust + 1
        ind2_acdc_rates = ind1_acdc_rates + (4*size(G_ACDC)) - 1

        ind1_oth_nuc_rates = ind2_acdc_rates + 1
        ind2_oth_nuc_rates = ind1_oth_nuc_rates

        ind1_gr = ind2_oth_nuc_rates + 1
        ind2_gr = ind1_gr + n_bins_par - 1

        ind1_sink = ind2_gr + 1
        ind2_sink = ind1_sink + n_cond_tot - 1

        ALLOCATE(new_mpi_send_aer_buf(1:ind2_CONS_clust))
        ALLOCATE(new_mpi_recv_aer_buf(ind_kz_tag:ind2_sink))

        new_mpi_send_aer_buf_size = size(new_mpi_send_aer_buf, 1)
        new_mpi_recv_aer_buf_size = size(new_mpi_recv_aer_buf, 1)

        ALLOCATE(VAPOR(kz,n_cond_tot))
        ALLOCATE(SINK(kz,n_cond_tot))
        ALLOCATE(CLUSTERS(kz,sum(G_ACDC(:)%neq_syst)))

        ALLOCATE(Formation_rates(kz,21)) ! 5x4 from ACDC + organic
        Formation_rates = 0d0

        ALLOCATE(MASS_COMPO(kz,n_bins_par,n_cond_tot))
        ALLOCATE(VOL_COMPO(kz,n_bins_par,n_cond_tot))
        ALLOCATE(CH_RES_org(n_cond_tot))
        ! ALLOCATE(index_vapor(n_cond_tot))
        ! ALLOCATE(vapor_names(n_cond_tot))
        ! ALLOCATE(sat_vap(n_cond_tot))

        ! Initialize the nominal diameter vector
        ALLOCATE(nominal_dp(n_bins_par))
        nominal_dp = get_dp()

        ! J_Distr is a vector containing the weighing factor for newly nucleated particles. Majority (~75-90%) of new clustes
        ! are put in the first model bin (whatever that is), and the rest are distributed in the next couple of bins in
        ! diminishing numbers (exp decay), but not further than NPF_DIST (default=1.15) times the smallest diameter.
        ALLOCATE(J_distr(   MINLOC(abs(nominal_dp-nominal_dp(1)*NPF_DIST),1)   ))
        do ii=1,size(J_distr)
            J_distr(ii) = 0.5*exp(-0.7d0*ii)
        end do
        J_distr(1) = J_distr(1) + 1-sum(J_distr)
        print *, 'Distributing nucleated particles over '//TRIM(i2chr(size(J_distr, 1)))//' bins.'

        ! Initialize the dummy for dmps fitting
        ALLOCATE(conc_fit(n_bins_par))
        save_measured = conc_fit
        losses_fit    = conc_fit
        ! Read in background particles
        ! IF (N_MODAL>0) THEN
            call Multimodal(MMODES, get_dp(), conc_fit, N_MODAL)
            conc_fit = conc_fit*1d6
            CAll send_conc(current_PSD%dp_range(1),current_PSD%dp_range(2),conc_fit)
            do ii = 1, n_bins_par
                CALL set_composition(ii,nominal_dp(ii), .false.)
            end do
            ! N_MODAL = -1d0
        ! END IF

        ALLOCATE(dp_emi_aer(9,2))
        dp_emi_aer(:,1) = [3d0,10d0,20d0,30d0,50d0,70d0,100d0,200d0,400d0]*1d-9
        dp_emi_aer(:,2) = [10d0,20d0,30d0,50d0,70d0,100d0,200d0,400d0,1000d0]*1d-9
        ! Allocate the change vectors for integration timestep control
        ALLOCATE(d_dpar(n_bins_par))
        ALLOCATE(d_npar(n_bins_par))
        ALLOCATE(d_npdep(n_bins_par))
        d_dpar = 0d0
        d_npar = 0d0
        d_npdep = 0d0

        if (use_dmps_partial .and. use_dmps) THEN
            ! write(*, *) 'Using dmps_special:'

            ! both limits have to be far enough from the smallest and largest diameters
            if (dmps_lowband_upper_limit > nominal_dp(2)) THEN
                dmps_sp_min = minloc(abs(get_dp()-dmps_lowband_upper_limit),1)
                print *, 'Lower bins replaced from 1 to '//i2chr(dmps_sp_min)
            ELSE IF (dmps_lowband_upper_limit > 0d0) THEN
                STOP '     ---- dmps_lowband_upper_limit is is too small and will not be used ----'
            END IF

            if (dmps_highband_lower_limit < nominal_dp(n_bins_par-1) &
             .and. dmps_highband_lower_limit > nominal_dp(2)) THEN
                dmps_sp_max = minloc(abs(get_dp()-dmps_highband_lower_limit),1)
                print *, 'Upper bins replaced from indices '//TRIM(i2chr(dmps_sp_max))//' to '//TRIM(i2chr(n_bins_par))
            ELSE IF (dmps_highband_lower_limit > 0d0) THEN
                STOP '     ---- dmps_highband_lower_limit is is too large and will not be used ----'
            END IF

        END IF

        ALLOCATE(conc_vapour(n_cond_tot))

        ! Allocate the change vectors for integration timestep control
        ALLOCATE(d_vap(max(VAPOUR_PROP%n_cond_tot,1)))
        d_vap = 0d0

        ! if (Deposition) THEN
        !     ALLOCATE(Depos_composition(max(VAPOUR_PROP%n_cond_tot,1)))
        !     Depos_composition = 0d0
        ! END if
        ! if (Chem_Deposition) THEN
        !     ALLOCATE(c_org_wall(VAPOUR_PROP%n_cond_org))
        !     c_org_wall = 0d0
        !     c_org_wall_old = c_org_wall
        !     ALLOCATE(VAP_DEP_MASS_WALLS(n_cond_org))
        !     VAP_DEP_MASS_WALLS = 0d0
        ! END IF

        do ii=1,n_bins_par
            ! if (Kelvin_taylor) THEN
                pre_Kelvin(ii,:) = 2D0*VAPOUR_PROP%surf_tension * VAPOUR_PROP%molar_mass / ( Rg * 300d0 * VAPOUR_PROP%density * nominal_dp(ii)/2d0)
            ! else
            !     pre_Kelvin(ii,:) = EXP(4D0*VAPOUR_PROP%surf_tension * VAPOUR_PROP%molar_mass / ( Rg * 300d0 * VAPOUR_PROP%density * nominal_dp(ii)) )
            ! end if

        end do

        print *, ''

    END SUBROUTINE AERO_INIT

    SUBROUTINE VAPOUR_INIT(SPC_NAMES_,CHEM_DIR)

        IMPLICIT NONE

        CHARACTER(len=16), INTENT(IN)     :: SPC_NAMES_(:)
        CHARACTER(len=256)                :: buf
        integer, allocatable              :: Natoms(:,:)
        integer                           :: ioi, ii, ind_core
        integer                           :: i, j, k, jj, path_l(2), N_Xtr = 0
        integer                           :: rows, cols
        LOGICAL                           :: elements_missing = .false.
        real(dp)                          :: molar_mass, parameter_A, parameter_B,SURFACE_TENSION_buf=-1d0, dummy
        CHARACTER(len=256)                :: species_name
        CHARACTER(len=*)                  :: CHEM_DIR
        CHARACTER(len=20), ALLOCATABLE    :: atoms_name(:)
        CHARACTER(len=2)    :: noacd

        if (PSD_MODE == 1) write(*,FMT_MSG) 'Using fully stationary PSD scheme with '//TRIM(i2chr(n_bins_par))//' bins.'
        if (PSD_MODE == 2) write(*,FMT_MSG) 'Using fixed grid/moving average PSD scheme with '//TRIM(i2chr(n_bins_par))//' bins.'
        print FMT_LEND,
        write(*,FMT_MSG) 'Reading Vapour name file '// TRIM(Vap_names)
        OPEN(unit=802, File= TRIM(Vap_names) , STATUS='OLD', iostat=ioi)
        call handle_file_io(ioi, Vap_names, &
            'If Condensation is used, "Vapour file" must be defined (in tab "Advanced").')

        H2SO4_ind_in_chemistry = IndexFromName( 'H2SO4', SPC_NAMES_ )

        rows = ROWCOUNT(802)
        cols = COLCOUNT(802)
        VAPOUR_PROP%n_cond_org = 0
        do j = 1,rows
            read(802,*,iostat=ioi) species_name
            if (j<=limit_vapours .or. j==rows) THEN
                if ((TRIM(species_name) == 'HOA')         & ! HOA is sometimes used in PRAM
                    .or.(TRIM(species_name) == 'GENERIC') & ! GENERIC is not active in chemistry
                    ! We want to add H2SO4 to the end of vapours, and it should not be in vapour file anyway
                    .or.(TRIM(species_name) == 'H2SO4') ) &
                    cycle
                k = IndexFromName( TRIM(species_name), SPC_NAMES_ )
                if (k>0) THEN
                    VAPOUR_PROP%n_cond_org = VAPOUR_PROP%n_cond_org + 1
                end if
            end if
        end do
        REWIND(802)

        ! Here we account for the non-volatile pseudocompound GENERIC
        VAPOUR_PROP%n_cond_org = VAPOUR_PROP%n_cond_org + 1
        ! Here we add place for sulfuric acid, non-organic
        VAPOUR_PROP%n_cond_tot = VAPOUR_PROP%n_cond_org + 1

        allocate(VAPOUR_PROP%Vapour_names (VAPOUR_PROP%n_cond_tot) )
        allocate(VAPOUR_PROP%molar_mass   (VAPOUR_PROP%n_cond_tot) )
        allocate(VAPOUR_PROP%psat_a       (VAPOUR_PROP%n_cond_tot) )
        allocate(VAPOUR_PROP%psat_b       (VAPOUR_PROP%n_cond_tot) )
        allocate(VAPOUR_PROP%molec_mass   (VAPOUR_PROP%n_cond_tot) )
        allocate(VAPOUR_PROP%molec_volume (VAPOUR_PROP%n_cond_tot) )
        allocate(VAPOUR_PROP%density      (VAPOUR_PROP%n_cond_tot) )
        allocate(VAPOUR_PROP%surf_tension (VAPOUR_PROP%n_cond_tot) )
        allocate(VAPOUR_PROP%diff         (VAPOUR_PROP%n_cond_tot) )
        allocate(VAPOUR_PROP%c_speed      (VAPOUR_PROP%n_cond_tot) )
        allocate(VAPOUR_PROP%c_sat        (VAPOUR_PROP%n_cond_tot) )
        allocate(VAPOUR_PROP%sink         (VAPOUR_PROP%n_cond_tot) )
        allocate(VAPOUR_PROP%cond_type    (VAPOUR_PROP%n_cond_tot) )
        allocate(VAPOUR_PROP%molec_dia    (VAPOUR_PROP%n_cond_tot) )
        allocate(VAPOUR_PROP%mfractions   (VAPOUR_PROP%n_cond_tot) )
        allocate(VAPOUR_PROP%alpha        (VAPOUR_PROP%n_cond_tot) )
        allocate(VAPOUR_PROP%diff_vol     (VAPOUR_PROP%n_cond_tot) )
        allocate(VAPOUR_PROP%diff_dia     (VAPOUR_PROP%n_cond_tot) )
        allocate(VAPOUR_PROP%ind_ch       (VAPOUR_PROP%n_cond_tot) )
        ! This is the vector that combines chemistry and condensible vapours.
        ! -1 because GENERIC is not in gas phase. Ever. H2SO4 is picked from chemistry manually
        ALLOCATE(index_cond(VAPOUR_PROP%n_cond_org-1))

        index_cond = 0

        if (USE_RH_CORRECTION) THEN
            ! These are only allocated for sulfuric acid, maybe nitric acid in the future.
            ! Uses same index as in VAPOUR_PROP%ind_H2SO4
            allocate(VAPOUR_PROP%wet_dia(VAPOUR_PROP%ind_GENERIC+1:VAPOUR_PROP%n_cond_tot))
            allocate(VAPOUR_PROP%wet_mass(VAPOUR_PROP%ind_GENERIC+1:VAPOUR_PROP%n_cond_tot))
        END IF


        print FMT_SUB, 'Compounds picked from Vapours file: '//TRIM(i2chr(VAPOUR_PROP%n_cond_org))
        print FMT_SUB, 'Total number of condensibles      : '//TRIM(i2chr(VAPOUR_PROP%n_cond_tot))

        ! Reading the vap names and vap vapour_properties
        VAPOUR_PROP%ind_GENERIC = VAPOUR_PROP%n_cond_org
        VAPOUR_PROP%ind_H2SO4   = VAPOUR_PROP%n_cond_tot
        VAPOUR_PROP%Mfractions  = 0.0
        VAPOUR_PROP%Mfractions(VAPOUR_PROP%ind_GENERIC) = 1d0 !
        VAPOUR_PROP%ind_ch = 0

        ! ---------------------------------------------------------------------
        ! ORGANIC VAPOUR PROPERTIES
        ! ---------------------------------------------------------------------
        ii = 1
        do j = 1, rows
            if (cols==4) read(802,*,iostat=ioi)   species_name, molar_mass, parameter_A, parameter_B
            if (cols==5) read(802,*,iostat=ioi)   species_name, molar_mass, parameter_A, parameter_B, SURFACE_TENSION_buf

            if (j<=limit_vapours .or. j==rows) THEN ! j==rows is because "GENERIC" is selected even if limitvapours<rows

                ! Check if the compounds exists in Chemistry and only then add to vapours
                k = IndexFromName( species_name, SPC_NAMES_ )

                if (k>0.and..not.(TRIM(species_name)=='GENERIC')) index_cond(ii) = k
                if (k>0.and..not.(TRIM(species_name)=='GENERIC')) VAPOUR_PROP%ind_ch(ii) = k

                if ((k>0).or.(TRIM(species_name)=='GENERIC')) THEN
                    ! fill the hash table for vapour index -> chemistry index

                    VAPOUR_PROP%molar_mass(ii)   = molar_mass *1D-3 ! kg/mol
                    VAPOUR_PROP%psat_a(ii)       = parameter_A
                    VAPOUR_PROP%psat_b(ii)       = parameter_B
                    VAPOUR_PROP%vapour_names(ii) = TRIM(species_name)
                    VAPOUR_PROP%molec_mass(ii)   = VAPOUR_PROP%molar_mass(ii)/Na  !kg/#
                    VAPOUR_PROP%cond_type(ii)    = 1
                    if (TRIM(species_name)=='GENERIC') THEN
                        VAPOUR_PROP%density(ii)   = HARD_CORE_DENSITY  ! kg/m3
                    ELSE
                        VAPOUR_PROP%density(ii)   = ORGANIC_DENSITY  ! kg/m3
                    END IF
                    VAPOUR_PROP%molec_volume(ii) = VAPOUR_PROP%molec_mass(ii)/VAPOUR_PROP%density(ii)
                    VAPOUR_PROP%diff_vol(ii)     = VAPOUR_PROP%molec_mass(ii)/VAPOUR_PROP%density(ii)

                    if (cols==4) VAPOUR_PROP%surf_tension(ii) = SURFACE_TENSION
                    if (cols==5) THEN
                        if (SURFACE_TENSION_buf<0d0) VAPOUR_PROP%surf_tension(ii) = abs(SURFACE_TENSION_buf)*SURFACE_TENSION
                        if (SURFACE_TENSION_buf>=0d0) VAPOUR_PROP%surf_tension(ii) = SURFACE_TENSION_buf
                    end if
                    SURFACE_TENSION_buf = -1d0

                    VAPOUR_PROP%alpha(ii)         = 1.0
                    ! this is just initial value, always gets updated with T
                    VAPOUR_PROP%c_sat(ii)         = saturation_conc_m3(VAPOUR_PROP%psat_a(ii),VAPOUR_PROP%psat_b(ii), 293.15d0)

                    ii = ii + 1
                END IF
            END IF
        end do
        close(802)

        ! In case GENERIC was not in Vapour file (should not happen if the file was from the GUI), add GENERIC with default values
        if (VAPOUR_PROP%vapour_names(VAPOUR_PROP%n_cond_org) /= 'GENERIC') THEN
            print FMT_WARN0, 'The vapour file did not contain GENERIC, adding it now. You should update your vapour file.'
            ii = VAPOUR_PROP%n_cond_org
            VAPOUR_PROP%vapour_names(ii)  = 'GENERIC'
            VAPOUR_PROP%cond_type(ii)     = 1  ! Generic non-evaporating
            VAPOUR_PROP%molar_mass(ii)    = 437.0 * 1d-3
            VAPOUR_PROP%psat_a(ii)        = 10
            VAPOUR_PROP%psat_b(ii)        = 1d4
            VAPOUR_PROP%density(ii)       = HARD_CORE_DENSITY
            VAPOUR_PROP%molec_mass(ii)    = VAPOUR_PROP%molar_mass(ii)/Na
            VAPOUR_PROP%molec_volume(ii)  = VAPOUR_PROP%molec_mass(ii)/VAPOUR_PROP%density(ii)
            VAPOUR_PROP%diff_vol(ii)      = VAPOUR_PROP%molec_mass(ii)/VAPOUR_PROP%density(ii)
            VAPOUR_PROP%surf_tension(ii)  = SURFACE_TENSION
            VAPOUR_PROP%alpha(ii)         = 1.0
            ! this is just initial value, always gets updated with T
            VAPOUR_PROP%c_sat(ii)         = saturation_conc_m3(VAPOUR_PROP%psat_a(ii),VAPOUR_PROP%psat_b(ii), 293.15d0)
        END IF

        VAPOUR_PROP%Mfractions  = 0.0
        if (TRIM(HARD_CORE) == '') THEN
            VAPOUR_PROP%Mfractions(VAPOUR_PROP%ind_GENERIC) = 1d0 !
        ELSE
            ind_core = IndexFromName(HARD_CORE,VAPOUR_PROP%vapour_names)
            if (ind_core>0 .and. ind_core<=VAPOUR_PROP%n_cond_tot) THEN
                VAPOUR_PROP%Mfractions(ind_core) = 1d0
            ELSE
                VAPOUR_PROP%Mfractions(VAPOUR_PROP%ind_GENERIC) = 1d0
            END IF
        END IF

        if (Use_atoms) THEN
            OPEN(unit=804, File=TRIM(Vap_atoms) , STATUS='OLD', iostat=ioi)
            call handle_file_io(ioi, Vap_atoms, 'Terminating the program.')
            write(*,FMT_MSG) 'Reading the list of elemental composition: '// TRIM(Vap_atoms)

            allocate(Natoms(5,ROWCOUNT(804))) ! C,O,N,H
            allocate(atoms_name(ROWCOUNT(804)))

            Natoms = 0

            DO j=1,ROWCOUNT(804)
                READ(804,*, iostat=ioi) atoms_name(j), molar_mass, Natoms(:,j)
            END DO

            CLOSE(804)

            DO j=1,VAPOUR_PROP%n_cond_org
                jj = IndexFromName(VAPOUR_PROP%vapour_names(j), atoms_name)
                if (jj>0) THEN
                    vapour_prop%diff_vol(j) = (Natoms(1,jj)*15.9D0 + Natoms(2,jj)*6.11D0 &
                                              + Natoms(4,jj)*2.31D0 + Natoms(3,jj)*4.54D0) + Natoms(5,jj) * 22.9D0 ![Å^3]
                ELSE
                    elements_missing = .true.
                END IF
            END DO
            if (elements_missing) print FMT_WARN0, 'Not all organics had atom content, using generic diameter'

            deallocate(Natoms)
            deallocate(atoms_name)

        end if

        ! Now the volumes are updated, the diameter can be calculated
        VAPOUR_PROP%molec_dia = (6D0 * VAPOUR_PROP%molec_volume / pi )**(1D0/3D0)  ! molecular diameter [m]
        if (use_diff_dia_from_diff_vol) THEN
            VAPOUR_PROP%diff_dia = (6D0 * 1d-30 * VAPOUR_PROP%diff_vol / pi )**(1D0/3D0)  ! molecular diameter [m]
        ELSE
            VAPOUR_PROP%diff_dia = VAPOUR_PROP%molec_dia  ! molecular diameter [m]
        END IF
        ! ---------------------------------------------------------------------
        ! Sulfuric acid treated separately
        ! ---------------------------------------------------------------------
        ii = VAPOUR_PROP%n_cond_tot
        VAPOUR_PROP%vapour_names(ii)  = 'H2SO4'
        VAPOUR_PROP%cond_type(ii)     = 2  ! Acid
        VAPOUR_PROP%molar_mass(ii)    = 98.0785 * 1d-3
        VAPOUR_PROP%psat_a(ii)        = 0
        VAPOUR_PROP%psat_b(ii)        = 20000d0
        VAPOUR_PROP%density(ii)       = 1830.5 ! kg/m3
        VAPOUR_PROP%molec_mass(ii)    = VAPOUR_PROP%molar_mass(ii)/Na
        VAPOUR_PROP%molec_volume(ii)  = VAPOUR_PROP%molec_mass(ii)/VAPOUR_PROP%density(ii)
        VAPOUR_PROP%diff_vol(ii)      = 4*6.11D0 + 2*2.31D0 + 22.9D0 ! O=4, H=2, S=1
        VAPOUR_PROP%diff_dia(ii)      = VAPOUR_PROP%molec_dia(ii)
        VAPOUR_PROP%surf_tension(ii)  = 0.07
        VAPOUR_PROP%molec_dia(ii)     = (6D0 * VAPOUR_PROP%molec_volume(ii) / pi )**(1D0/3D0)  ! molecular diameter [m]
        VAPOUR_PROP%alpha(ii)         = 1.0
        VAPOUR_PROP%c_sat(ii)         = 0.0 ! Sulfuric acid stays put
        VAPOUR_PROP%ind_ch(ii)        = IndexFromName('H2SO4', SPC_NAMES_)

        do i = 1, VAPOUR_PROP%n_cond_org-1
            if (TRIM(VAPOUR_PROP%vapour_names(i))/= TRIM(SPC_NAMES_(index_cond(i)))) THEN
                print*, i, index_cond(i), VAPOUR_PROP%vapour_names(i), SPC_NAMES_(index_cond(i))
                stop 'VAPOR INDEXES DO NOT MATCH'
            end if
        end do

        call ORGANIC_NUCL(dummy, [1d0,1d0], SPC_NAMES_, TRIM(CHEM_DIR)//'/INPUT/nucl_homs.txt')

    END SUBROUTINE VAPOUR_INIT


    SUBROUTINE run_arca_psd(CH_GAS, dt_aero, TEMPK, pres, emissions,seasaltparts, J_ORG)

        IMPLICIT NONE

        real(dp), INTENT(INOUT) :: CH_GAS(:)
        real(dp), INTENT(IN) :: TEMPK
        real(dp), INTENT(IN) :: pres
        real(dp), INTENT(IN) :: dt_aero
        real(dp), INTENT(IN) :: emissions(:), seasaltparts(:)
        INTEGER :: i,ii,jj,i_15, i1,i2
        real(dp), INTENT(  OUT) :: J_ORG

        GTEMPK  = TEMPK
        Gpres   = pres

        ! ! Read in background particles
        ! IF (N_MODAL>0) THEN
        !     call Multimodal(MMODES, get_dp(), conc_fit, N_MODAL)
        !     conc_fit = conc_fit*1d6
        !     CAll send_conc(current_PSD%dp_range(1),current_PSD%dp_range(2),conc_fit)
        !     do i = 1, n_bins_par
        !         CALL set_composition(i,nominal_dp(i), .false.)
        !     end do
        !     N_MODAL = -1d0
        ! END IF

        ! Pick the condensibles from chemistry and change units from #/cm^3 to #/m^3
        conc_vapour = 0d0
        conc_vapour(1:VAPOUR_PROP%n_cond_org-1) =  CH_GAS(index_cond)*1D6 ! mol/m3

        ! Poor sulfuric acid always needs special treatment
        conc_vapour(VAPOUR_PROP%ind_H2SO4) = CH_GAS(H2SO4_ind_in_chemistry)*1D6

        ! Update vapour pressures for organics
        VAPOUR_PROP%c_sat(1:VAPOUR_PROP%n_cond_org) =  saturation_conc_m3( &
                            VAPOUR_PROP%psat_a(1:VAPOUR_PROP%n_cond_org),  &
                            VAPOUR_PROP%psat_b(1:VAPOUR_PROP%n_cond_org),  &
                            GTEMPK)

        CALL UPDATE_MOLECULAR_DIFF_AND_CSPEED(VAPOUR_PROP)

        ! IF (CHEM_DEPOSITION) CALL CALCULATE_CHEMICAL_WALL_LOSS(conc_vapour(1:VAPOUR_PROP%n_cond_org),c_org_wall)

        dmass  = 0d0
        d_vap  = 0
        d_dpar = 0

        kelvin_eff = pre_Kelvin ** (300d0/GTEMPK)
        if (Kelvin_taylor) Kelvin_Eff = 1d0 + pre_Kelvin * 300d0 / GTEMPK + (pre_Kelvin * 300d0 / GTEMPK) **2 / 2d0 + (pre_Kelvin * 300d0 / GTEMPK) **3 / 6d0

        if (OPTIONS%NUCLEATION) THEN
            J_TOTAL_M3 = 0d0

            if (OPTIONS%ACDC) THEN
                CALL ACDC_J(CH_GAS, 3d6, dt_aero)

                J_TOTAL_M3 = sum(G_ACDC(:)%J_OUT_M3(1)) ! [particles/s/m^3]

                dmass = 0d0
                dconc_dep_mix = 0d0

                ! Negative mixing ratio makes this nucleation
                mix_ratio = -1d0

                ! New particles are assigned GENERIC composition.
                dmass(1:size(J_distr,1),VAPOUR_PROP%ind_GENERIC) = nominal_dp(1:size(J_distr,1))**3*pi/6d0 * VAPOUR_PROP%density(VAPOUR_PROP%ind_GENERIC)
                ! The new particles are distributed over bins spanning bins that are 15% larger in diameter than the first bin, typically 3 bins
                dconc_dep_mix(1:size(J_distr,1)) = J_TOTAL_M3*dt_aero*J_distr

                CALL Mass_Number_Change('mixing')
                ! Update current_psd
                current_PSD = new_PSD
            END if

            IF (OPTIONS%ORG_NUCL) THEN

                J_ORG = 0d0
                CALL ORGANIC_NUCL(J_ORG, CH_GAS)

                dmass = 0d0
                dconc_dep_mix = 0d0
                ! Negative mixing ratio makes this nucleation
                mix_ratio = -1d0
                i_15 = size(PACK(nominal_dp, nominal_dp-1.55e-9<0),1)
                ! New particles are assigned GENERIC composition.
                dmass(i_15:i_15+1,VAPOUR_PROP%ind_GENERIC) = nominal_dp(i_15:i_15+1)**3*pi/6d0 * VAPOUR_PROP%density(VAPOUR_PROP%ind_GENERIC)
                ! The new particles are equally distributed over 2 bins around 1.5 nm
                dconc_dep_mix(i_15:i_15+1) = J_ORG*dt_aero*0.5d0
                CALL Mass_Number_Change('mixing')
                ! Update current_psd
                current_PSD = new_PSD

                ! Going back to main as cm⁻³
                J_ORG = J_ORG*1d-6

            END IF


        END IF

        IF (OPTIONS%AER_EMISSIONS) THEN

            dmass = 0d0
            dconc_dep_mix = 0d0
            ! Negative mixing ratio makes this emission
            mix_ratio = -1d0

            do ii = 1,size(dp_emi_aer,1)
                i1 = size(PACK(nominal_dp, nominal_dp-dp_emi_aer(ii,1)<0),1)+1
                i2 = size(PACK(nominal_dp, nominal_dp-dp_emi_aer(ii,2)<0),1)
                ! print*, 'dp_emi_aer(ii)',dp_emi_aer(ii,1),dp_emi_aer(ii,2),i1,i2
                ! print*, 'nominal_dp(i1:i2)',nominal_dp(i1:i2)
                dconc_dep_mix(i1:i2) = emissions(ii)*bin_ratio_lg
            end do

            ! Emitted particles are assigned GENERIC composition.
            dmass(:,VAPOUR_PROP%ind_GENERIC) = nominal_dp(:)**3*pi/6d0 * VAPOUR_PROP%density(VAPOUR_PROP%ind_GENERIC)

            dconc_dep_mix(:) = dconc_dep_mix*dt_aero

            CALL Mass_Number_Change('mixing')
            ! Update current_psd
            current_PSD = new_PSD

        END IF

        IF (OPTIONS%SEASALT_EMS) THEN
            ! ----------------------------------------
            ! SEASALT PARTICLES!
            if (sum(seasaltparts)>0d0) THEN
                dmass = 0d0
                dconc_dep_mix = 0d0
                ! Negative mixing ratio makes this emission
                mix_ratio = -1d0

                dconc_dep_mix = seasaltparts*bin_ratio_lg

                ! Emitted particles are assigned GENERIC composition.
                dmass(:,VAPOUR_PROP%ind_GENERIC) = nominal_dp(:)**3*pi/6d0 * VAPOUR_PROP%density(VAPOUR_PROP%ind_GENERIC)

                dconc_dep_mix(:) = dconc_dep_mix*dt_aero

                CALL Mass_Number_Change('mixing')
                ! Update current_psd
                current_PSD = new_PSD
            END IF

        END IF


        if (OPTIONS%CONDENSATION) THEN
            CALL Condensation_apc(VAPOUR_PROP,conc_vapour,dmass, dt_aero,d_dpar,d_vap, kelvin_eff)
            ! Distribute mass
            CALL Mass_Number_Change('condensation')

            ! Update current_psd
            current_PSD = new_PSD

            CH_GAS(index_cond) = conc_vapour(1:VAPOUR_PROP%n_cond_org-1) * 1d-6
        END IF

        if (OPTIONS%COAGULATION) THEN
            Call COAGULATION_ROUTINE(dconc_coag, dt_aero,d_npar, update_K)
            Call Mass_Number_Change('coagulation')

            ! Update current_psd
            current_PSD = new_PSD
        END IF

    END SUBROUTINE run_arca_psd




    PURE CHARACTER(LEN=12) FUNCTION f2chr(number)
        IMPLICIT NONE
        real(dp), INTENT(IN) :: number
        write(f2chr, '(es12.3)') number
        f2chr = ADJUSTL(f2chr)
    END FUNCTION f2chr

    PURE FUNCTION i2chr(number) result(out)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: number
        CHARACTER(len=int(LOG10(MAX(ABS(number)*1d0, 1d0))+2)-min(0,sign(1,number))) :: out
        write(out, '(i0)') number
    END FUNCTION i2chr

    PURE FUNCTION di2chr(number) result(out)
        IMPLICIT NONE
        INTEGER(dint), INTENT(IN) :: number
        CHARACTER(len=int(LOG10(MAX(ABS(number)*1d0, 1d0))+2)-min(0,sign(1_dint,number))) :: out
        write(out, '(i0)') number
    END FUNCTION di2chr

    PURE LOGICAL FUNCTION equal(a,b)
        IMPLICIT NONE
        REAL(dp), INTENT(IN) :: a,b
        equal = ABS(a-b) .lt. 1d-200
    END FUNCTION equal

SUBROUTINE PARSE_MULTIMODAL()
    IMPLICIT NONE
    real(dp) :: buffer(100)
    INTEGER  :: ioi, sze

    buffer = -9999d0
    read(mmodal_input, *, IOSTAT=ioi) buffer
    sze =  SIZE(PACK(buffer,buffer>0))
    if (sze>=3) THEN
        ALLOCATE(MMODES(sze))
        MMODES = PACK(buffer,buffer>0)
    ELSE
        N_MODAL = -1d0
    ENDIF
    IF (mmodal_input_inuse == 0) N_MODAL = -1d0

END SUBROUTINE PARSE_MULTIMODAL

! =================================================================================================
! Parametrisation of organic nucleation as used in Roldin et al., Nat Comm. 2019
! .................................................................................................
SUBROUTINE ORGANIC_NUCL(J_TOTAL_M3, CH_GAS, SPC_NAMES_,nuc_homs_file)
    real(dp) :: J_TOTAL_M3
    real(dp) :: CH_GAS(:)
    real(dp) :: ORGS,kJ,dH,J15
    real(dp), parameter :: dG = -15.1d3   ! [cal/mol]
    real(dp), parameter :: dS = 61.1      ! [cal/mol/K] NOTE this value had a typo in Nat. Comm. and should be positive!
    real(dp), PARAMETER :: Rcal = 1.98720 ! gas constant in calories/mol/K
    integer :: i,jj,n
    integer, save, allocatable  :: inds(:)
    logical, save  :: first_run = .True.
    character(25)  :: name
    character(len=*), optional  :: SPC_NAMES_(:)
    character(len=*), optional  :: nuc_homs_file

    if (first_run) THEN
        OPEN(UNIT=609, FILE=TRIM(nuc_homs_file), STATUS='OLD', ACTION='READ', iostat=jj)
        if (jj/=0) STOP 'Could not find NUCL_HOMS.txt'
        n = rowcount(609)
        print FMT_MSG, 'Using parametrisation for organic nucleation with '//i2chr(n)//' nucleating compounds'
        allocate(inds(n))
        inds = 0
        DO i=1,n
            read(609, *) name
            if (IndexFromName(TRIM(name), SPC_NAMES_)>0) inds(i) = IndexFromName(TRIM(name), SPC_NAMES_)
        END DO
        if (PRODUCT(inds) == 0) THEN
            print FMT_LEND,
            print FMT_FAT0, "List of nucleating organic compounds has vapours which are not in the chemistry."
            print FMT_SUB, "Options: Edit the file 'ModelLib/required/nucl_homs.txt' or turn of organic nucleation."
            print FMT_LEND,
            STOP "  ----------- Bye -------------"
        END IF
        first_run = .False.

    ELSE
        ORGS = sum(CH_GAS(inds))

        dH = dG - dS*GTEMPK
        kJ = 5d-13*EXP(-dH/(Rcal) * (1/GTEMPK - 1/298d0))

        J15 = kJ * ORGS*CH_GAS(VAPOUR_PROP%ind_ch(VAPOUR_PROP%ind_H2SO4))
        ! if (GTIME%printnow) print FMT_MSG, 'Organic formation rate: '//TRIM(ADJUSTL(f2chr(J15)))//' [/s/cm3]'
        J_TOTAL_M3 = J_TOTAL_M3 + J15*1d6
    END IF

END SUBROUTINE ORGANIC_NUCL



! =================================================================================================
! ACDC Nucleation. See that all values here are in SI-UNITS: CUBIC METERS, KELVINS AND PASCALS.
! Written by Tinja Olenius
! .................................................................................................
! Input for get_acdc_J:
! c_acid:            Sulfuric acid concentration [1/m3]
! c_base:            base (ammonia) concentration [1/m3]
! c_org:             Nucleating organic concentration [1/m3]. not in use currently
! cs_H2SO4:          Condensation sink of sulfuric acid [1/s]
! TempK:             Temperature in Kelvins
! IPR:               Ion production rate in ion pairs per second [1/m3/s]. 3d6 is a good guestimate
! dt:                Main time step [s]
! ACDC_solve_ss:     Solve steady state or only to timestep duration (generally makes no difference)
! J_ACDC_NH3_M3:        Particle formation rate due to ammonia [1/s/m3]. Sum of J_NH3_BY_IONS
! acdc_cluster_diam: Outgrowing cluster diameter [m]
! J_NH3_BY_IONS:     Particle formation rate by ions (neutral, negative and positive) [1/s/m3].
!                    The outgrowing cluster typically has 5 to 6 H2SO4 in it
! .................................................................................................
! Input for get_acdc_D:
! c_acid:            Sulfuric acid concentration [1/m3]
! c_dma:             DMA concentration [1/m3]
! c_org:             Nucleating organic concentration [1/m3]. not in use currently
! cs_H2SO4:          Condensation sink of sulfuric acid [1/s]
! TempK:             Temperature in Kelvins
! dt:                Main time step [s]
! time:              model time; used only for outputting information in cluster
! ACDC_solve_ss:     Solve steady state or only to timestep duration (generally makes no difference)
! J_ACDC_DMA_M3:        Particle formation rate due to DMA [1/s/m3]
! acdc_cluster_diam: Outgrowing cluster diameter [m]. The cluster typically has 5 to 6 H2SO4 in it
!
! HOW TO ADD MORE ACDC SUBMODULES?
! - Duplicate any of the src/ACDC/ACDC_* directories and rename it with the next available number.
!   If, for example one adds a sixth module, the new directory is named ACDC_06
! - In the new directory ACDC_06, delete files acdc_equations_0x*.f90 and acdc_system_0x*.f90
! - In the new directory ACDC_06, name all the remaining .f90 files with *_0x#.f90 -> *_0x6.f90
! - In the new system directory, find and replace 0X* with 0X6 (case insensitive search/replace) INSIDE
!   the following files:
!      - acdc_simulation_setup_0x*.f90
!      - driver_acdc_J_0x*.f90 and
!      - get_acdc_J_0x*.f90
! - In this subroutine, update the USE statements and SELECT_PROPER_ACDC_SUBMODULE using the same logic
! - In the GUI->Cluster formation, update the new system settings and press Update ACDC_06 and recompile ARCA
! - The makefile should not need any updating, but in case compiling error, run make from the terminal:
! make cleanmain && make
! =================================================================================================


SUBROUTINE ACDC_J(C, IPR, dt)
    use get_acdc_J_0X1, only : acdc_plugin_0X1=>acdc_plugin
    use get_acdc_J_0X2, only : acdc_plugin_0X2=>acdc_plugin
    use get_acdc_J_0X3, only : acdc_plugin_0X3=>acdc_plugin
    use get_acdc_J_0X4, only : acdc_plugin_0X4=>acdc_plugin
    use get_acdc_J_0X5, only : acdc_plugin_0X5=>acdc_plugin

    implicit none

    REAL(dp), intent(in) :: C(:)
    REAL(dp), intent(in) :: dt
    REAL(dp), intent(in) :: IPR ! Ion pairs / second / m3
    LOGICAL, save        :: first_time = .true., ss_handle
    ! REAL(dp)             :: IPR = 0 ! Ion pairs / second / cm3, converted to Ion pairs / second / m3
    INTEGER :: ii,jj, dd
    REAL(dp) :: CS

    if (first_time) THEN
        ss_handle = ACDC_solve_ss
    end if

    CS = VAPOUR_PROP%SINK(VAPOUR_PROP%ind_H2SO4)
    ! -------- NEW ACDC -------------
    do ii=1,size(G_ACDC)

        if (G_ACDC(ii)%inuse) THEN
            do jj=1,size(G_ACDC(ii)%ACDC_MONOMER_NAMES)
                ! by default, G_ACDC(ii)%ACDC_monConc(jj) = 0d0
                ! If (G_ACDC(ii)%ACDC_LINK_IND(jj)>0) THEN
                    ! if (G_ACDC(ii)%ACDC_LINK_IND(jj) == inm_H2SO4) THEN
                    !     G_ACDC(ii)%ACDC_monConc(jj) = H2SO4*1d6
                    ! ELSE
                        ! G_ACDC(ii)%ACDC_monConc(jj) = C(G_ACDC(ii)%ACDC_LINK_IND(jj))*1d6
                    ! END IF
                ! ELSE if (G_ACDC(ii)%ACDC_LINK_IND(jj)<0) THEN
                    G_ACDC(ii)%ACDC_monConc(jj) = C(-1*G_ACDC(ii)%ACDC_LINK_IND(jj))*1d6
                    ! print*,'CONCE', G_ACDC(ii)%ACDC_MONOMER_NAMES(jj), G_ACDC(ii)%ACDC_monConc(jj)
                ! END IF
            end do
            SELECT_PROPER_ACDC_SUBMODULE: select case (ii)
                case (1)
                    CALL acdc_plugin_0X1(G_ACDC(ii)%ACDC_MONOMER_NAMES,G_ACDC(ii)%ACDC_monConc,CS,GTEMPK,IPR,dt,1d-10,&
                        G_ACDC(ii)%J_OUT_M3,G_ACDC(ii)%Cl_diam,ss_handle, .false., .false., n_clusters=G_ACDC(ii)%ACDC_clusConc)
                case (2)
                    CALL acdc_plugin_0X2(G_ACDC(ii)%ACDC_MONOMER_NAMES,G_ACDC(ii)%ACDC_monConc,CS,GTEMPK,IPR,dt,1d-10,&
                        G_ACDC(ii)%J_OUT_M3,G_ACDC(ii)%Cl_diam,ss_handle, .false., .false., n_clusters=G_ACDC(ii)%ACDC_clusConc)
                case (3)
                    CALL acdc_plugin_0X3(G_ACDC(ii)%ACDC_MONOMER_NAMES,G_ACDC(ii)%ACDC_monConc,CS,GTEMPK,IPR,dt,1d-10,&
                        G_ACDC(ii)%J_OUT_M3,G_ACDC(ii)%Cl_diam,ss_handle, .false., .false., n_clusters=G_ACDC(ii)%ACDC_clusConc)
                case (4)
                    CALL acdc_plugin_0X4(G_ACDC(ii)%ACDC_MONOMER_NAMES,G_ACDC(ii)%ACDC_monConc,CS,GTEMPK,IPR,dt,1d-10,&
                        G_ACDC(ii)%J_OUT_M3,G_ACDC(ii)%Cl_diam,ss_handle, .false., .false., n_clusters=G_ACDC(ii)%ACDC_clusConc)
                case (5)
                    CALL acdc_plugin_0X5(G_ACDC(ii)%ACDC_MONOMER_NAMES,G_ACDC(ii)%ACDC_monConc,CS,GTEMPK,IPR,dt,1d-10,&
                        G_ACDC(ii)%J_OUT_M3,G_ACDC(ii)%Cl_diam,ss_handle, .false., .false., n_clusters=G_ACDC(ii)%ACDC_clusConc)
            end select SELECT_PROPER_ACDC_SUBMODULE

            G_ACDC(ii)%J_OUT_CM3(:) = 1d-6 * G_ACDC(ii)%J_OUT_M3(:)
            ! Missing values are substituted w -1
            WHERE (G_ACDC(ii)%J_OUT_CM3(:)<-9.999d-7 ) G_ACDC(ii)%J_OUT_CM3(:) = -1d0

            if (first_time) print FMT_SUB, 'ACDC system '//TRIM(i2chr(ii))//', Outgrowing cluster size: '//TRIM(f2chr(G_ACDC(ii)%Cl_diam))
            if (first_time) THEN
                if (nominal_dp(1)/G_ACDC(ii)%Cl_diam > 1d0) THEN
                    print FMT_FAT0, 'System nr: '//TRIM(i2chr(ii))
                    print *, 'Cluster diameter ', G_ACDC(ii)%Cl_diam
                    STOP 'Smallest particle bin is larger than outgrowing cluster'
                else if (nominal_dp(1)/G_ACDC(ii)%Cl_diam < 0.9d0) THEN
                    print FMT_WARN0, 'System nr: '//TRIM(i2chr(ii))//', outgrowing cluster diameter:'//TRIM(f2chr(G_ACDC(ii)%Cl_diam))
                    print FMT_WARN0, 'Outgrowing cluster diameter is '//TRIM(f2chr(G_ACDC(ii)%Cl_diam/nominal_dp(1)))//' times larger than the bin where it is put'
                    print FMT_WARN0, 'solution: increase minimum Dp'
                end if
            end if
        END IF
    END DO

    if (first_time) first_time = .false.

END SUBROUTINE ACDC_J  ! END ACDC Nucleation

SUBROUTINE PARSE_ACDC_SYSTEMS(SPC_NAMES_)
    use acdc_simulation_setup_0x1, only : get_system_size_0X1=>get_system_size
    use acdc_system_0x1, only : n_charges_0X1=>n_charges

    use acdc_simulation_setup_0x2, only : get_system_size_0X2=>get_system_size
    use acdc_system_0x2, only : n_charges_0X2=>n_charges

    use acdc_simulation_setup_0x3, only : get_system_size_0X3=>get_system_size
    use acdc_system_0x3, only : n_charges_0X3=>n_charges

    use acdc_simulation_setup_0x4, only : get_system_size_0X4=>get_system_size
    use acdc_system_0x4, only : n_charges_0X4=>n_charges

    use acdc_simulation_setup_0x5, only : get_system_size_0X5=>get_system_size
    use acdc_system_0x5, only : n_charges_0X5=>n_charges

    IMPLICIT NONE
    INTEGER :: ii,jj,kk,ioi=0, sze, counter,mm
    CHARACTER(len=16), INTENT(IN)     :: SPC_NAMES_(:)
    CHARACTER(len=16)   :: name(24)
    CHARACTER(len=256)  :: System,Energies,Dipoles,path,Nickname
    NAMELIST /ACDC_RECORD_NML/ System,Energies,Dipoles,Nickname
    print *, ''
    print FMT_HDR, 'Allocating ACDC systems to the selected input'

    do jj=1,size(G_ACDC)
        ! numbers of clusters and equations, indices of outgoing fluxes,max. diameter in system
        select case (JJ)
        case (1)
            ALLOCATE(G_ACDC(jj)%nout_syst(n_charges_0X1))
            call get_system_size_0X1(G_ACDC(jj)%neq_syst,G_ACDC(jj)%nclust_syst,G_ACDC(jj)%nout_syst,G_ACDC(jj)%Cl_diam)
            print*, 'ACDC',jj,': neq_syst', G_ACDC(jj)%neq_syst
        case (2)
            ALLOCATE(G_ACDC(jj)%nout_syst(n_charges_0X2))
            call get_system_size_0X2(G_ACDC(jj)%neq_syst,G_ACDC(jj)%nclust_syst,G_ACDC(jj)%nout_syst,G_ACDC(jj)%Cl_diam)
            print*, 'ACDC',jj,': neq_syst', G_ACDC(jj)%neq_syst
        case (3)
            ALLOCATE(G_ACDC(jj)%nout_syst(n_charges_0X3))
            call get_system_size_0X3(G_ACDC(jj)%neq_syst,G_ACDC(jj)%nclust_syst,G_ACDC(jj)%nout_syst,G_ACDC(jj)%Cl_diam)
            print*, 'ACDC',jj,': neq_syst', G_ACDC(jj)%neq_syst
        case (4)
            ALLOCATE(G_ACDC(jj)%nout_syst(n_charges_0X4))
            call get_system_size_0X4(G_ACDC(jj)%neq_syst,G_ACDC(jj)%nclust_syst,G_ACDC(jj)%nout_syst,G_ACDC(jj)%Cl_diam)
            print*, 'ACDC',jj,': neq_syst', G_ACDC(jj)%neq_syst
        case (5)
            ALLOCATE(G_ACDC(jj)%nout_syst(n_charges_0X5))
            call get_system_size_0X5(G_ACDC(jj)%neq_syst,G_ACDC(jj)%nclust_syst,G_ACDC(jj)%nout_syst,G_ACDC(jj)%Cl_diam)
            print*, 'ACDC',jj,': neq_syst', G_ACDC(jj)%neq_syst
        end select

        name = '----------------'
        read(ACDC_links(jj),*, iostat=ioi) name(:)

        counter = 0
        do ii=0,11
            if (name((ii*2)+1)/= '----------------') THEN
                counter = counter + 1
            end if
        end do

        if (ACDC_SYSTEMS(jj)==1) THEN
            print FMT_MSG, 'ACDC submodule #'//i2chr(jj)//' is initialized with input definitions.'

            WRITE(path,'(a,i0.2,a)') 'src/ACDC/ACDC_',jj,'/ACDC_RECORD_NML'
            OPEN(UNIT=889, FILE=TRIM(path), STATUS='OLD', ACTION='READ', iostat=ii)
            if (ii==0) THEN
                Nickname = '-not defined-'
                do while (ii==0)
                    READ(889,NML = ACDC_RECORD_NML, IOSTAT=ii)
                end do
                close(889)
                G_ACDC(jj)%SYSTEM_FILE=System
                G_ACDC(jj)%ENERGY_FILE=Energies
                G_ACDC(jj)%DIPOLE_FILE=Dipoles
                G_ACDC(jj)%NICKNAME=Nickname
                print FMT_SUB, 'System name: '//TRIM(G_ACDC(jj)%NICKNAME)
                print FMT_SUB, 'System based on '//TRIM(G_ACDC(jj)%SYSTEM_FILE)
                print FMT_SUB, 'Energies from '//TRIM(G_ACDC(jj)%ENERGY_FILE)
                print FMT_SUB, 'Dipoles from '//TRIM(G_ACDC(jj)%DIPOLE_FILE)
            END IF


            G_ACDC(jj)%inuse = .true.

            ALLOCATE(G_ACDC(jj)%ACDC_LINK_IND(counter))
            ALLOCATE(G_ACDC(jj)%ACDC_monConc(counter))
            G_ACDC(jj)%ACDC_monConc(counter) = 0d0
            ALLOCATE(G_ACDC(jj)%ACDC_MONOMER_NAMES(counter))
            ALLOCATE(G_ACDC(jj)%ACDC_clusConc(G_ACDC(jj)%neq_syst))

            do ii=1,counter
                if (name(((ii-1)*2)+1)/= '----------------') THEN
                    G_ACDC(jj)%ACDC_MONOMER_NAMES(ii)(:) = TRIM(name(((ii-1)*2)+1))
                    ! G_ACDC(jj)%ACDC_LINK_IND(ii) = IndexFromName(TRIM(name(((ii-1)*2)+2)), [(MODS(mm)%NAME, mm=1,size(MODS))] )
                    ! IF (G_ACDC(jj)%ACDC_LINK_IND(ii) == 0) THEN
                        ! print FMT_WARN0, 'In ACDC '//i2chr(jj)//': Compound '//TRIM(name(((ii-1)*2)+2))//' does not exist in input'
                        print FMT_SUB, 'Searching compound '//TRIM(name(((ii-1)*2)+2))//' from chemistry...'
                        G_ACDC(jj)%ACDC_LINK_IND(ii) = -1 * IndexFromName(TRIM(name(((ii-1)*2)+2)), SPC_NAMES_ )
                        if (G_ACDC(jj)%ACDC_LINK_IND(ii)<0) print FMT_SUB, 'Found '//TRIM(name(((ii-1)*2)+2))//' from chemistry. Using as input '//i2chr(ii)//' for ACDC '//i2chr(jj)//'.'
                        if (G_ACDC(jj)%ACDC_LINK_IND(ii)==0) THEN
                            print FMT_FAT0, 'In ACDC 1: Could not find '//TRIM(name(((ii-1)*2)+2))
                            stop
                        END IF
                    ! END IF
                end if
            end do
        ELSE
            print FMT_MSG, 'ACDC submodule #'//i2chr(jj)//' is not used.'
        END IF

    end do

    print FMT_LEND,


END SUBROUTINE PARSE_ACDC_SYSTEMS

end module
