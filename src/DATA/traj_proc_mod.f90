module traj_proc_mod

  use constants_mod
  use tool_mod
  use traj_data_mod

  implicit none

  private

  public :: traj_emission_init, traj_emission_calc


contains


  !============================================================================!
  ! Initialize the variables related to trajectory emission rate calculation
  !============================================================================!
  subroutine traj_emission_init()

  end subroutine traj_emission_init


  !============================================================================!
  ! Calculate the emission rates from the input data
  !============================================================================!
  subroutine traj_emission_calc(infield, ptr, z, time, emis, emi_aer)
    ! NSPEC and chemical indices, ind_xxxx
    use second_Parameters, kpp_dp=>dp

    type(input_field), intent(in   ) :: infield
    integer          , intent(in   ) :: ptr  ! input pointer
    real(dp)         , intent(in   ) :: z(:)
    real(dp)         , intent(in   ) :: time
    real(dp)         , intent(  out) :: emis(:, :)
    real(dp)         , intent(  out) :: emi_aer(:, :)

    integer :: emisbio_levels(3) = [2,3,4]
    integer :: emisant_level1, emisant_level2  ! bottom and top emission levels
    integer :: emisant_nlevel
    integer :: emisid, chemid, sizeid, I, KZ

    real(dp), save :: bio_dz = 0d0 ! height of emission column for biogenics
    real(dp) :: emis_  ! emission rate buffer
    real(dp) :: part_  ! particle emission rate buffer
    real(dp), allocatable :: emis1d_(:)  ! 1d emission rate buffer

    ! integer, allocatable ::

    !------------------------------------------------------!
    ! Obtain biogenic emission rates and set them to the
    ! emission level
    !------------------------------------------------------!

    if (time==0d0) THEN
        do i=1,size(emisbio_levels,1)
            kz = emisbio_levels(i)
            if (kz > 1) then
                bio_dz = bio_dz + (z(kz)-z(kz-1))
            end if
        end do
        print '(a,f0.1,a)', 'biogenic emissions distributed over ',bio_dz, ' meters at layers'
        print*, emisbio_levels
    end if

    !----- CH3CHO -----!
    emisid = infield%id_emisbio_CH3CHO
    chemid = ind_CH3CHO
    call linear_interp_2p( &
      infield%var(infield%id_time)%f1d(ptr), &
      infield%var(infield%id_time)%f1d(ptr+1), &
      infield%var(emisid)%f1d(ptr), &
      infield%var(emisid)%f1d(ptr+1), &
      time, emis_ &
      )
      emis(emisbio_levels, chemid) = emis_/bio_dz

    !----- HCHO -----!
    emisid = infield%id_emisbio_HCHO
    chemid = ind_HCHO
    call linear_interp_2p( &
      infield%var(infield%id_time)%f1d(ptr), &
      infield%var(infield%id_time)%f1d(ptr+1), &
      infield%var(emisid)%f1d(ptr), &
      infield%var(emisid)%f1d(ptr+1), &
      time, emis_ &
      )
      emis(emisbio_levels, chemid) = emis_/bio_dz

    !----- CH3COCH3 -----!
    emisid = infield%id_emisbio_CH3COCH3
    chemid = ind_CH3COCH3
    call linear_interp_2p( &
      infield%var(infield%id_time)%f1d(ptr), &
      infield%var(infield%id_time)%f1d(ptr+1), &
      infield%var(emisid)%f1d(ptr), &
      infield%var(emisid)%f1d(ptr+1), &
      time, emis_ &
      )
      emis(emisbio_levels, chemid) = emis_/bio_dz

    !----- C5H8 -----!
    emisid = infield%id_emisbio_C5H8
    chemid = ind_C5H8
    call linear_interp_2p( &
      infield%var(infield%id_time)%f1d(ptr), &
      infield%var(infield%id_time)%f1d(ptr+1), &
      infield%var(emisid)%f1d(ptr), &
      infield%var(emisid)%f1d(ptr+1), &
      time, emis_ &
      )
      emis(emisbio_levels, chemid) = emis_/bio_dz

    !----- APINENE-----!
    emisid = infield%id_emisbio_APINENE
    chemid = ind_APINENE
    call linear_interp_2p( &
      infield%var(infield%id_time)%f1d(ptr), &
      infield%var(infield%id_time)%f1d(ptr+1), &
      infield%var(emisid)%f1d(ptr), &
      infield%var(emisid)%f1d(ptr+1), &
      time, emis_ &
      )
      emis(emisbio_levels, chemid) = emis_/bio_dz

    !----- BPINENE-----!
    emisid = infield%id_emisbio_BPINENE
    chemid = ind_BPINENE
    call linear_interp_2p( &
      infield%var(infield%id_time)%f1d(ptr), &
      infield%var(infield%id_time)%f1d(ptr+1), &
      infield%var(emisid)%f1d(ptr), &
      infield%var(emisid)%f1d(ptr+1), &
      time, emis_ &
      )
      emis(emisbio_levels, chemid) = emis_/bio_dz

    !----- CH3OH -----!
    emisid = infield%id_emisbio_CH3OH
    chemid = ind_CH3OH
    call linear_interp_2p( &
      infield%var(infield%id_time)%f1d(ptr), &
      infield%var(infield%id_time)%f1d(ptr+1), &
      infield%var(emisid)%f1d(ptr), &
      infield%var(emisid)%f1d(ptr+1), &
      time, emis_ &
      )
      emis(emisbio_levels, chemid) = emis_/bio_dz

    !----- MBO -----!
    emisid = infield%id_emisbio_MBO
    chemid = ind_MBO
    call linear_interp_2p( &
      infield%var(infield%id_time)%f1d(ptr), &
      infield%var(infield%id_time)%f1d(ptr+1), &
      infield%var(emisid)%f1d(ptr), &
      infield%var(emisid)%f1d(ptr+1), &
      time, emis_ &
      )
      emis(emisbio_levels, chemid) = emis_/bio_dz

    !----- OMT -----!
    ! ?? OMT will be distributed to other MT species or directly set to ind_OMT,
    ! this needs to be determined in future
    ! emisid = infield%id_emisbio_OMT
    ! chemid = ind_OMT
    ! call linear_interp_2p( &
    !   infield%var(infield%id_time)%f1d(ptr), &
    !   infield%var(infield%id_time)%f1d(ptr+1), &
    !   infield%var(emisid)%f1d(ptr), &
    !   infield%var(emisid)%f1d(ptr+1), &
    !   time, emis_ &
    !   )

    !----- SQT -----!
    ! ?? SQT will be distributed to other SQT species, and this needs to be
    ! determined in future
    !---------------!
    ! emisid = infield%id_emisbio_SQT
    ! chemid = ind_SQT
    ! call linear_interp_2p( &
    !   infield%var(infield%id_time)%f1d(ptr), &
    !   infield%var(infield%id_time)%f1d(ptr+1), &
    !   infield%var(emisid)%f1d(ptr), &
    !   infield%var(emisid)%f1d(ptr+1), &
    !   time, emis_ &
    !   )

    !----- DMS -----!
    emisid = infield%id_emisbio_DMS
    chemid = ind_DMS
    call linear_interp_2p( &
      infield%var(infield%id_time)%f1d(ptr), &
      infield%var(infield%id_time)%f1d(ptr+1), &
      infield%var(emisid)%f1d(ptr), &
      infield%var(emisid)%f1d(ptr+1), &
      time, emis_ &
      )
    emis(emisbio_levels, chemid) = emis_/bio_dz

    !------------------------------------------------------!
    ! Obtain anthropogenic emission rates and set them to the
    ! emission levels
    !------------------------------------------------------!

    !----- Obtain emission model levels -----!
    ! Here only the levels within the input emission level range will be set
    ! to interpolated emissions.
    ! For example, the input emission levels are like:
    ! x  <-         tlh
    ! o  <-     mlh
    ! x  <- blh     tlh
    ! o  <-     mlh
    ! x  <- blh     tlh
    ! o  <-     mlh
    ! x  <- blh
    ! And the model levels within the lowest blh and highest tlh are considered
    ! as the emission levels, the values are interpolated from the mlh levels,
    ! with extrapolated points using boundary values.
    ! The emission rates need to be added up with biogenic emissions.
    !----------------------------------------!

    emisant_level1 = minloc( z, &
      mask=(z >= infield%var(infield%id_emisant_blh)%f1d(1)), &
      dim=1 &
      )
    emisant_level2 = maxloc( z, &
      mask=(z <= infield%var(infield%id_emisant_tlh)%f1d(infield%nlevel_emisant)), &
      dim=1 &
      )
    emisant_nlevel = emisant_level2 - emisant_level1 + 1
    allocate(emis1d_(emisant_nlevel))

    !----- CO -----!
    emisid = infield%id_emisant_CO
    chemid = ind_CO

     call interp_time2p_level1d( &
      infield%var(infield%id_emisant_mlh)%f1d, infield%var(infield%id_time)%f1d, &
      transpose(infield%var(emisid)%f2d), ptr, &
      z(emisant_level1:emisant_level2), time, &
      emis1d_ &
      )
    emis(emisant_level1:emisant_level2, chemid) = &
      emis(emisant_level1:emisant_level2, chemid) + emis1d_

    !----- NOx -----!
    ! NO and NO2, need to modify how to distribute them in future
    !---------------!
    emisid = infield%id_emisant_NOx
    ! chemid = ind_NO2

     call interp_time2p_level1d( &
      infield%var(infield%id_emisant_mlh)%f1d, infield%var(infield%id_time)%f1d, &
      transpose(infield%var(emisid)%f2d), ptr, &
      z(emisant_level1:emisant_level2), time, &
      emis1d_ &
      )
    emis(emisant_level1:emisant_level2, ind_NO) = &
      emis(emisant_level1:emisant_level2, ind_NO) + emis1d_*0.1d0
    emis(emisant_level1:emisant_level2, ind_NO2) = &
      emis(emisant_level1:emisant_level2, ind_NO2) + emis1d_*0.9d0

    ! ----- NH3 -----!
    emisid = infield%id_emisant_NH3
    chemid = ind_NH3
    call interp_time2p_level1d( &
      infield%var(infield%id_emisant_mlh)%f1d, infield%var(infield%id_time)%f1d, &
      transpose(infield%var(emisid)%f2d), ptr, &
      z(emisant_level1:emisant_level2), time, &
      emis1d_ &
      )
    emis(emisant_level1:emisant_level2, chemid) = &
      emis(emisant_level1:emisant_level2, chemid) + emis1d_

    ! ----- DMA is guessed based on NH3 emissions! -----!
    emis(1:2, ind_DMA) = emis(1:2, ind_NH3) * dma_nh3_fraction

    !----- SO2 -----!
    emisid = infield%id_emisant_SO2
    chemid = ind_SO2

     call interp_time2p_level1d( &
      infield%var(infield%id_emisant_mlh)%f1d, infield%var(infield%id_time)%f1d, &
      transpose(infield%var(emisid)%f2d), ptr, &
      z(emisant_level1:emisant_level2), time, &
      emis1d_ &
      )
    emis(emisant_level1:emisant_level2, chemid) = &
      emis(emisant_level1:emisant_level2, chemid) + emis1d_

    !----- C5H8 -----!
    emisid = infield%id_emisant_C5H8
    chemid = ind_C5H8

    call interp_time2p_level1d( &
      infield%var(infield%id_emisant_mlh)%f1d, infield%var(infield%id_time)%f1d, &
      transpose(infield%var(emisid)%f2d), ptr, &
      z(emisant_level1:emisant_level2), time, &
      emis1d_ &
      )
    emis(emisant_level1:emisant_level2, chemid) = &
      emis(emisant_level1:emisant_level2, chemid) + emis1d_

    ! ----- MT -----!
    emisid = infield%id_emisant_MT

    call interp_time2p_level1d( &
      infield%var(infield%id_emisant_mlh)%f1d, infield%var(infield%id_time)%f1d, &
      transpose(infield%var(emisid)%f2d), ptr, &
      z(emisant_level1:emisant_level2), time, &
      emis1d_ &
      )

    ! Distribution of anth. MT emissions according to TRP1 model species. Data from
    ! Guenther et al (1999). TRP1 is used for biogenicMT emissions, but SAPRC-99 documentation
    ! recommends this over no estimation at all (SOURCE: IMPLEMENTATION OF THE SAPRC-99 CHEMICAL
    ! MECHANISM INTO THE MODELS-3 FRAMEWORK. Report to the United States Environmental Protection Agency
    ! By William P. L. Carter January 29, 2000, page 23)

    emis(emisant_level1:emisant_level2, ind_APINENE) = &
      emis(emisant_level1:emisant_level2, ind_APINENE) + emis1d_ * 0.37719298d0

    emis(emisant_level1:emisant_level2, ind_BPINENE) = &
      emis(emisant_level1:emisant_level2, ind_BPINENE) + emis1d_ * 0.27192982d0

    emis(emisant_level1:emisant_level2, ind_CARENE) = &
      emis(emisant_level1:emisant_level2, ind_CARENE) + emis1d_ * 0.16666667d0

    emis(emisant_level1:emisant_level2, ind_LIMONENE) = &
      emis(emisant_level1:emisant_level2, ind_LIMONENE) + emis1d_ * 0.0877193d0

    emis(emisant_level1:emisant_level2, ind_SABINENE) = &
      emis(emisant_level1:emisant_level2, ind_SABINENE) + emis1d_ * 0.09649123d0

    deallocate(emis1d_)

    !----- PARTICLES -----!
    emisid = infield%id_emisaer_3_10nm
    sizeid = 1

    call interp_time2p_level1d( infield%var(infield%id_emisant_mlh)%f1d, infield%var(infield%id_time)%f1d, &
    transpose(infield%var(emisid)%f2d), ptr, z(emisant_level1:emisant_level2), time, &
    emi_aer(emisant_level1:emisant_level2,sizeid) )

    emisid = infield%id_emisaer_10_20nm
    sizeid = sizeid + 1
    call interp_time2p_level1d( infield%var(infield%id_emisant_mlh)%f1d, infield%var(infield%id_time)%f1d, &
    transpose(infield%var(emisid)%f2d), ptr, z(emisant_level1:emisant_level2), time, &
    emi_aer(emisant_level1:emisant_level2,sizeid) )

    emisid = infield%id_emisaer_20_30nm
    sizeid = sizeid + 1
    call interp_time2p_level1d( infield%var(infield%id_emisant_mlh)%f1d, infield%var(infield%id_time)%f1d, &
    transpose(infield%var(emisid)%f2d), ptr, z(emisant_level1:emisant_level2), time, &
    emi_aer(emisant_level1:emisant_level2,sizeid) )

    emisid = infield%id_emisaer_30_50nm
    sizeid = sizeid + 1
    call interp_time2p_level1d( infield%var(infield%id_emisant_mlh)%f1d, infield%var(infield%id_time)%f1d, &
    transpose(infield%var(emisid)%f2d), ptr, z(emisant_level1:emisant_level2), time, &
    emi_aer(emisant_level1:emisant_level2,sizeid) )

    emisid = infield%id_emisaer_50_70nm
    sizeid = sizeid + 1
    call interp_time2p_level1d( infield%var(infield%id_emisant_mlh)%f1d, infield%var(infield%id_time)%f1d, &
    transpose(infield%var(emisid)%f2d), ptr, z(emisant_level1:emisant_level2), time, &
    emi_aer(emisant_level1:emisant_level2,sizeid) )

    emisid = infield%id_emisaer_70_100nm
    sizeid = sizeid + 1
    call interp_time2p_level1d( infield%var(infield%id_emisant_mlh)%f1d, infield%var(infield%id_time)%f1d, &
    transpose(infield%var(emisid)%f2d), ptr, z(emisant_level1:emisant_level2), time, &
    emi_aer(emisant_level1:emisant_level2,sizeid) )

    emisid = infield%id_emisaer_100_200nm
    sizeid = sizeid + 1
    call interp_time2p_level1d( infield%var(infield%id_emisant_mlh)%f1d, infield%var(infield%id_time)%f1d, &
    transpose(infield%var(emisid)%f2d), ptr, z(emisant_level1:emisant_level2), time, &
    emi_aer(emisant_level1:emisant_level2,sizeid) )

    emisid = infield%id_emisaer_200_400nm
    sizeid = sizeid + 1
    call interp_time2p_level1d( infield%var(infield%id_emisant_mlh)%f1d, infield%var(infield%id_time)%f1d, &
    transpose(infield%var(emisid)%f2d), ptr, z(emisant_level1:emisant_level2), time, &
    emi_aer(emisant_level1:emisant_level2,sizeid) )

    emisid = infield%id_emisaer_400_1000nm
    sizeid = sizeid + 1
    call interp_time2p_level1d( infield%var(infield%id_emisant_mlh)%f1d, infield%var(infield%id_time)%f1d, &
    transpose(infield%var(emisid)%f2d), ptr, z(emisant_level1:emisant_level2), time, &
    emi_aer(emisant_level1:emisant_level2,sizeid) )


  end subroutine traj_emission_calc


  subroutine traj_update_chem_nconc()
  end subroutine traj_update_chem_nconc


  ! subroutine prepare_chemistry_traj()
  ! end subroutine prepare_chemistry_traj

end module traj_proc_mod
