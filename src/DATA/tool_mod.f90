module tool_mod

  use constants_mod

  implicit none

  public

!===== NML_MISC: misc options =====!
! experimental de-normalisation of emissions ... (P.C.)
logical :: weigh_with_srr = .false.
! experimental DMA fraction ... (P.C.)
real(dp) :: dma_nh3_fraction = 0d0

NAMELIST /NML_MISC/ weigh_with_srr, dma_nh3_fraction



contains


!==============================================================================!
!
! Obtain the interpolated value from the time series of data. Usually the data
! are read in from files.
! This is used to substitue this kind of job:
! CH_CONS_ALL(:,ind_O3)  = (CH_gas_hyy(2,nxodrad) + &
!   (time_in_month-(nxodrad-1)*dt_obs) * &
!   (CH_gas_hyy(2,nxodrad+1)-CH_gas_hyy(2,nxodrad))/dt_obs) * air(:) / 1.0d9
!
! Usage:
!   t1: current time, usually is month time
!   ind1: the largest time index before t1
!   dt: time step between two time indices
!   var: data, 1D array
!   linear_interp: the interpolated value
!
!==============================================================================!

function linear_interp(t1, ind1, var, dt)
  REAL(dp) :: t1  ! current time
  INTEGER :: ind1  ! current index
  REAL(dp) :: dt  ! observation time step
  REAL(dp), INTENT(IN) :: var(:)  ! time series of data
  REAL(dp) :: linear_interp  ! interpolated value

  linear_interp = var(ind1) + (t1-(ind1-1)*dt)/dt * (var(ind1+1)-var(ind1))
end function linear_interp


!==============================================================================!
!
! Obtain the interpolated value between two points.
!
! Usage:
!   x0, x1: x at two points
!   y0, y1: y at two points
!   x: current position
!   y: interpolated value at x
!
!==============================================================================!

subroutine linear_interp_2p(x0, x1, y0, y1, x, y)
  REAL(dp), intent(in) :: x0, x1
  REAL(dp), intent(in) :: y0, y1
  REAL(dp), INTENT(in) :: x
  REAL(dp), INTENT(out) :: y

  y = y0 + (x - x0) * (y1 - y0) / (x1 - x0)
end subroutine linear_interp_2p


!==============================================================================!
!
! Linear interpolation for 1D array
!
!==============================================================================!

subroutine interp_1d(nd, xd, yd, ni, xi, yi)
  INTEGER, INTENT(IN) :: nd  ! number of data points
  INTEGER, INTENT(IN) :: ni  ! number of interpolation points

  REAL(dp), INTENT(IN) :: xd(nd), yd(nd)  ! input data points and values, should be monotonic
  REAL(dp), INTENT(IN) :: xi(ni)  ! interpolation points
  REAL(dp), INTENT(OUT) :: yi(ni)  ! interpolation values

  INTEGER :: i, k  ! loop integers
  REAL(dp) :: t
  REAL(dp) :: xd_incr(nd)
  REAL(dp) :: yd_incr(nd)

  ! Check number of data points
  IF (nd < 2) THEN
    WRITE(*,*) 'Interpolation can not be done because there is only one available data point.'
    RETURN
  END IF

  ! Make sure x values are monotonically increasing
  IF (xd(2) > xd(1)) THEN
    xd_incr = xd
    yd_incr = yd
  ELSE
    xd_incr = xd(nd:1:-1)
    yd_incr = yd(nd:1:-1)
  END IF

  ! Default values
  yi(1:ni) = 0.0_dp

  ! When there is only one data point, the interpolation values are equal to that only value
  IF (nd == 1) THEN
    yi(1:ni) = yd_incr(1)
    RETURN
  END IF

  DO i = 1, ni
    ! The point is beyond the left boundary
    IF (xi(i) <= xd_incr(1)) THEN
      ! Extrapolation
      t = (xi(i) - xd_incr(1)) / (xd_incr(2) - xd_incr(1))
      yi(i) = (1.0_dp - t)*yd_incr(1) + t*yd_incr(2)
    ! The point is beyond the right boundary
    ELSE IF (xd_incr(nd) <= xi(i)) THEN
      t = (xi(i) - xd_incr(nd-1)) / (xd_incr(nd) - xd(nd-1))
      yi(i) = (1.0_dp - t)*yd_incr(nd-1) + t*yd_incr(nd)
    ! The point is within the domain
    ELSE
      DO k = 2, nd
        IF (xd_incr(k-1) <= xi(i) .and. xi(i) <= xd_incr(k)) THEN
          t = (xi(i) - xd_incr(k-1)) / (xd_incr(k) - xd_incr(k-1))
          yi(i) = (1.0_dp - t)*yd_incr(k-1) + t*yd_incr(k)
          EXIT
        END IF
      END DO
    END IF
  END DO
END SUBROUTINE interp_1d


!==============================================================================!
!
! Obtain the interpolated value between two time points and then interpolate to
! model levels
!
! Usage:
!   level_in: input 1D level array
!   time_in: input 1D time array
!   data_in: input 2D data with dimension (level, time)
!   ptr_in: current time pointer
!   level_out: levels of outpu data
!   time_out: output time
!   data_out: outpu data, usually with dimension (level)
!
!==============================================================================!

subroutine interp_time2p_level1d( &
  level_in, time_in, data_in, ptr_in, &
  level_out, time_out, data_out &
  )

  real(dp), intent(in ) :: level_in(:)
  real(dp), intent(in ) :: time_in(:)
  real(dp), intent(in ) :: data_in(:, :)
  integer , intent(in ) :: ptr_in

  real(dp), intent(in ) :: level_out(:)
  real(dp), intent(in ) :: time_out
  real(dp), intent(out) :: data_out(:)

  real(dp) :: data_tmp(size(level_in))
  integer :: nlevel_in, nlevel_out
  integer :: i

  nlevel_in  = size(level_in)
  nlevel_out = size(level_out)

  ! Interpolate input data to current time
  do i = 1, nlevel_in
    call linear_interp_2p( &
      time_in(ptr_in), time_in(ptr_in+1), &
      data_in(i, ptr_in), data_in(i, ptr_in+1), &
      time_out, data_tmp(i) &
      )
  end do
  ! Interpolate input to model layer heights
    call interp_1d( &
      nlevel_in, level_in, data_tmp, &
      nlevel_out, level_out, data_out &
      )

end subroutine interp_time2p_level1d

end module tool_mod
