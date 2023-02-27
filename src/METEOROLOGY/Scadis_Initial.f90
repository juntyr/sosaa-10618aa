! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
! Module for initializing meteorology in sosa.                                                  !
! Code for main loop in MT_MainMet                                                              !
! Commented and restructured by Rosa Gierens, University of Helsinki                            !
!                                                                                               !
!                                                                                               !
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

MODULE Scadis_Initial

  ! variables from other modules
  USE SOSA_DATA

  ! subroutines and variables part of the meteorology scheme
  ! USE Scadis_Parameters
  ! USE Scadis_Radiation

  IMPLICIT NONE

  ! will make this privata later!
  !PRIVATE

  !Public subroutines:
  !PUBLIC ::

  !Public variables:
  !PUBLIC ::

  PRIVATE

  PUBLIC :: initiate_meteorology, Beta_PDF


CONTAINS


subroutine initiate_meteorology(flag_model_type)
  integer, intent(in) :: flag_model_type

  ! Use this in future
  select case (flag_model_type)
  case (STATION_MODE)
    call initiate_meteorology_stat()
  case (TRAJECTORY_MODE)
    call initiate_meteorology_traj()
  case default
    write(*,*) 'Wrong value of flag_model_type.'
  end select
end subroutine initiate_meteorology


SUBROUTINE initiate_meteorology_stat()
  integer :: count_nonnan

  ! Roughness at the surface z0 (0.0001-1.0 m)
  z0 = 0.1d0

  ! Cover fraction of low, middle, and high clouds (0 - 10).
  mclob=0.0d0
  mclom=0.0d0
  mclot=0.0d0

  ! Initial data for cloudiness  -> move to Scadis_Parameters
  clob=0.1d0*mclob
  clom=0.1d0*mclom
  clot=0.1d0*mclot

  ! s(k) - vertical leaf area density profile
  dtree=0.0d0
  dtree2=0.0d0
  do k=2,kz
    if(z(k) <= hc) then
      dtree=dtree+0.5*((z(k-1)/hc)**(pa-1)*(1.-z(k-1)/hc)**(3-1)+(z(k)/hc)**(pa-1)*(1.-z(k)/hc)**(3-1))*(z(k)-z(k-1))
      dtree2=dtree2+0.5*((z(k-1)/hc)**(pa2-1)*(1.-z(k-1)/hc)**(3-1)+(z(k)/hc)**(pa2-1)*(1.-z(k)/hc)**(3-1))*(z(k)-z(k-1))
    endif
  enddo

  do k=1,kz
     s1(k)=0.0d0
     s2(k)=0.0d0
     if(k==1) write(*,*) '(message printed only for k1) LAI related', k, LAI_proj, hc, s1(k), dtree, pa, z(k)
     if (z(k) <= hc) s1(k)=LAI_proj*(z(k)/hc)**(pa-1)*(1.-z(k)/hc)**(3-1)/dtree
     if(k==1) write(*,*) '(message printed only for k1) s1(k): ', k, s1(k)
     if (s1(k) < 0.005d0) s1(k)=0.

     if (z(k) <= hc) s2(k)=hc*(z(k)/hc)**(pa2-1)*(1.-z(k)/hc)**(3-1)/dtree2

     if (s2(k) < 0.005d0) s2(k)=0.

     if (hc == 0.0d0) then
        if (z(k) <= 2) kk=k
        nz=min(kk,kz)+1      !  level for downward radiation flux
     else
        if (z(k) <= (hc+2.)) kk=k
        nz=min(kk,kz)+1      !  level for downward radiation flux
     endif

    ! ?? temporary solution
     if (flag_model_type==TRAJECTORY_MODE) then
       nz = 13
     end if

     if (z(k) <= (33.6)) kkk=k
     nsmear=min(kkk,kz)+1 !  the closest level to 33.6 m

     if (z(k) <= (23.3)) kkk5=k
     n233=min(kkk,kz)+1 !  the closest level to Hyytiäla level

  enddo

  s1(1)=0.  ! [m2 m-3], LAD for Meteorology
  s1(2) = 1.32d0  ! [m2 m-3], LAD for understory layer (shallow dwarf shrub layer, h=0.2-0.3 m, projected LAI~0.25 m2 m-2)
                  ! 1.32 = 0.25/((0.38-0.0)/2). So all the understory layer is included in the layer 2.
  s2(1)=0.

  ! Set s1 as ASAM values
  ! s1(2:5) = (/0.077, 0.150, 0.222, 0.150/)


  ! Thus we get s1 and s2 that is now beta-function -
  ! to calculate the integral we should do next to get LAD in % for megan per layer

  do k=2,kz
     if(z(k) <= hc) then
        LAD_P(k) = 0.5 * (s2(k)+s2(k-1)) * (z(k)-z(k-1)) / hc
     endif
  enddo

  ! Now dintegr is close to 1 all time depending on resolution of the model
  ! To calculate real value of s1 and LAD we need to do the next

  do k=1,kz
     s2(k)=s2(k) * LAI_proj / hc   ! now s2 is the real LAD for MEGAN
  enddo

  ! Not used currently, should be deleted in future
  ! alai=0.
  ! do k=2,kz
  !    alai=alai+(s1(k-1)+s1(k))*dz(k)/2.
  ! enddo

  ! Some parameters (?????)
  ! ztop=z(kz)  ! not used
  ! zdown=0.0d0  ! not used
  tau3=300.0d0
  wind=10.0d0
  nua=0.0000141d0
  ! kw=1.0d0  ! not used
  eps=0.001_dp
  profile=9999.0d0
  kva=0.0000224d0
  kta=0.0000214d0

  cc2=0.0436d0
  c833=0.833d0
  c52=0.52d0
  alf=sqrt(cc2)*(c833-c52)/0.1872d0

  face= 2.7_dp
  alph=0.80_dp
  alni=0.80_dp
  roa=1.205
  rh2o=1000.0_dp
  cpa=1009.0_dp
  ch2o=4190.0_dp
  hleaf=0.003
  lv=2256000.0_dp
  lv=2501400.0_dp
  dels=0.99
  delf=0.95
  ! asl=1.5  ! not used
  cd=0.2
  gamma=0.0098_dp

  ! calculation of arrival of total radiation
  proatm=0.747d0   !proatm - average transparency of atmosphere
  lat_rad=lat_deg*deg2rad
  lon_rad=lon_deg*deg2rad
  f1=7.29d-5*2.0d0*sin(lat_rad)
  wsun=0.4145d0*sin(PIx2*(day-79.8d0)/365.24d0)    !  tirgnstrom

  zeit=pi*(time_in_day/c432-1.)
  cos_zenith=sin(lat_rad)*sin(wsun)+cos(lat_rad)*cos(wsun)*cos(zeit)
  optmas=796.*(sqrt(cos_zenith**2+0.0025)-cos_zenith)
  rads=solar_constant*cos_zenith
  tata=tan(lat_rad)*tan(wsun)
  upsn=c432/pi*atan(sqrt(1.-tata**2)/tata)
  downsn=c864*(1.-1./PIx2*atan(sqrt(1.-tata**2)/tata))
  ct=0.4
  if(cos_zenith >= 0.342) ct=0.2
  rads=rads*(1.-0.8*clob-0.5*clom-ct*clot)
  teil=28.5*(asin(cos_zenith))*proatm**4
  rsnt=1.0*rads*teil/(1.+teil)
  rskt=1.0*rads/(1.+teil)
  if(time_in_day <= upsn) rsnt=0.
  if(time_in_day <= upsn) rskt=0.
  if(time_in_day >= downsn) rsnt=0.
  if(time_in_day >= downsn) rsnt=0.
  if(time_in_day <= upsn) cos_zenith=1.
  if(time_in_day <= upsn) cos_zenith=1.
  if(time_in_day >= downsn) cos_zenith=1.
  if(time_in_day >= downsn) cos_zenith=1.


  !*****************************************************
  !boundary conditions and input parameters
  !*****************************************************

  select case (flag_model_type)
  case (STATION_MODE)
    if(abl == 1) then
       tah=border(1,1)
       qah=border(4,1)*roa ! [kg kg-1] --> [kg m-3]
    else
       tah=border_abl(1,1)
       qah=border_abl(4,1)*roa
    endif
  case (TRAJECTORY_MODE)
    tah = sum(infield%var(infield%id_t)%f2d(1, :)) / infield%ntime  ! ??, temporary solution
    qah = sum(infield%var(infield%id_q)%f2d(1, :)) / infield%ntime  ! ??, temporary solution
  end select

  fsoil1=49.0d0
  pp=0.0d0
  temtran=0.000001*0.5
  temdif=0.5

  select case (flag_model_type)
  case (STATION_MODE)
    if(abl == 1) then
       windx=border(2,1)
       windy=border(3,1)
    else
       windx=border_abl(2,1)
       windy=border_abl(3,1)
    endif
  case (TRAJECTORY_MODE)
    windx = sum(infield%var(infield%id_u)%f2d(1, :)) / infield%ntime  ! ??, temporary solution
    windy = sum(infield%var(infield%id_v)%f2d(1, :)) / infield%ntime  ! ??, temporary solution
  end select

  windx=sqrt(windx**2+windy**2)

  windy=0.

  udin=abs(wind)*0.40/log((hh+z0)/z0)

  sk00=0.002
  sks00=0.002

  al1=25.

  wg(1)=0.2
  wg1i=0.2
  wg(2)=0.250
  wg2i=0.250

  select case (flag_model_type)
  case (STATION_MODE)
    if (abl == 1) then
      ta(1)=tah-0.*border(5,1)*z(kz)
    else
      ! ta(1)=tah-border_abl(5,1)*z(kz)
      if (.not. isnan(local_temp(1,1))) then
        ta(1) = local_temp(1,1)
      else
        count_nonnan = count(.not. isnan(local_temp(1, :)))
        write(*,*) '[putian debug] count_nonnan', count_nonnan
        ta(1) = sum( local_temp(1, :), mask=.not. isnan(local_temp(1, :)) ) / count_nonnan
      end if
    endif
  case (TRAJECTORY_MODE)
    ta(1) = sum(infield%var(infield%id_t)%f2d(1, :)) / infield%ntime  ! ??, temporary solution
  end select

  qa(1)=rhovs_buck(ta(1))
  qa(1)=qa(1)*0.7
  bt(1)=udin**2/0.09**0.5

  do k=1,kz
     !
     ts0=ta(1)-10.
     tsoil(k)=ts0

     if(abl == 1) then
        u(k)=windx
        v(k)=windy
     else
        u(k)=windx
        if(z(k) < 300.) u(k)=windx*(z(k)/300.)**0.5
        v(k)=windy
        if(z(k) < 300.) v(k)=windy*(z(k)/300.)

     endif

     w(k)=0.

     qa(k)=qah
     ta(k)=ta(1)+z(k)*(tah-ta(1))/hh
     bt(k)=0.00001
     dbt(k)=0.001/bt(k)

     tsn(k)=ta(k)
     tsd(k)=ta(k)
     l(k)=0.40*(z(k)+z0)/(1+0.40*z(k)/al1)

     dhsn(k)=0.00
     dhsd(k)=0.00
     dvsn(k)=0.00
     dvsd(k)=0.00
     psn(k)=1.
     psk(k)=1.
     firu(k)=0.
     fird(k)=0.
     fniu(k)=0.
     fnid(k)=0.
     fphu(k)=0.
     fphd(k)=0.
     fluxle(k)=0.
     fluxh(k)=0.
     rih(k)=0.
  enddo

  DO k=1,11
     kt(k)=0.40*udin*(z(k)+z0)
  ENDDO
  DO k=12,kz
     kt(k)=kt(11)
  ENDDO

  DO k=1,kz
     ur(k)=0.1
  ENDDO


  nturb=1
  smol=1

  wg1(1)=wg(1)
  wg1(2)=wg(2)

  ! Set new values equal to old values in the initial step
  ! This is done from level 1 to kz
  l1     = l
  tsoil1 = tsoil
  alt    = 1.35
  alt1   = 1.35
  bt1    = bt
  dbt1   = dbt
  u1     = u
  v1     = v
  kt1    = kt
  ta1    = ta
  tsn1   = tsn
  tsd1   = tsd
  qa1    = qa
  w1     = w

  nxod=0
  nxodrad=0


  if(abl == 1) then
     do j=1,1488
        border(4,j)= 0.9*rhovs_buck( border(1,j) )
     enddo

  else

     do j=1,121
        border_abl(4,j)=border_abl(4,j)*roa
     enddo

  endif

  ! About snow
  rou=0.0_dp    !  snow cover is excluded from consideration now
  psnow=300.0_dp

end subroutine initiate_meteorology_stat


subroutine initiate_meteorology_traj()
  call initiate_meteorology_stat()

  ! pressure
  ! pres = sum(input_meteo_1(1, :, 6))/input_meteo_ntime
  ! write(*,*) 'pres: ', pres

end subroutine initiate_meteorology_traj




!************************************************************************************'!
! Beta probability density function for LAD (Eq. 1, Markkanen et al., 2003, BLM)
!************************************************************************************'!
FUNCTION Beta_PDF(x, a, b)
  REAL(dp), INTENT(IN) :: x(:)
  INTEGER, INTENT(IN) :: a, b

  REAL(dp) :: Beta_PDF(SIZE(x))

  INTEGER :: k, nx
  REAL(dp) :: below

  nx = SIZE(x)

  below = 0.0d0
  DO k=2, nx
    below = below + 0.5d0*( x(k-1)**(a-1)*(1.0d0-x(k-1))**(b-1) + x(k)**(a-1)*(1.0d0-x(k))**(b-1) )*(x(k)-x(k-1))
  ENDDO

  Beta_PDF = x**(a-1)*(1.0d0-x)**(b-1) / below
END FUNCTION Beta_PDF

END MODULE Scadis_Initial
