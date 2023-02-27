! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
! Module for meteorology in main loop in sosa, based on SCADIS by Andrej Sogachev.              !
! Initialization done in module MT_InitialMet.                                                  !
! Commented and restructured by Rosa Gierens, University of Helsinki                            !
!                                                                                               !
!                                                                                               !
! Doxygen is used fro ducumenting this code.
!
! LONG EXPLANATION
!  - where is original code to be found
! The original Scadis code that is referred to in the comments, is the code as it was used in Sosa
! 23.08.2014. It can be found from Scadis_Reference.md
!
!
!  - how is the structure   (other modules & subroutines in this one)
!
!
! REFERENCES:                                                                                   !
! - Andrej's papers
! - manual
!
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !


MODULE MT_MainMet

  ! variables from other modules
  USE SOSA_DATA !, ONLY : variablename

  ! subroutines and variables part of the meteorology scheme
  USE Scadis_Radiation
  USE Scadis_Subroutines, ONLY: New2Old,  New2Old_HalfStep
  USE Scadis_Data_Tools, ONLY: MakeNudgingData, make_nudging_data_new
  ! USE Scadis_Parameters
  ! USE Scadis_Tools
  ! USE Scadis_Energy

  IMPLICIT NONE

  ! will make this private later
  !PRIVATE

  !Public subroutines:
  !PUBLIC ::

  !Public variables:
  !PUBLIC ::


CONTAINS


subroutine update_meteorology(dt_obs, flag_model_type)
  real(dp), intent(in) ::  dt_obs ! temporal resolution of surface observations (= input data) [s]
  integer, intent(in) :: flag_model_type

  ! Use this in future
  select case (flag_model_type)
  case (STATION_MODE)
    call update_meteorology_stat(dt_obs)
  case (TRAJECTORY_MODE)
    call update_meteorology_traj(dt_obs)
  case default
    write(*,*) 'Wrong value of flag_model_type.'
  end select
end subroutine update_meteorology


subroutine update_meteorology_traj(dt_obs)
  real(kind = dp), intent(in) ::  dt_obs ! temporal resolution of surface observations (= input data) [s]

  ! call update_meteorology_stat(dt_obs)

  !--------------------------------------------------------------
  ! Get Ktt, Ktb and Kht, Khb so every diffusion term can be calculated more easily
  !--------------------------------------------------------------
  DO k=2,kz-1
    !!!!! Momentum
    Ktt(k) = df(k)*kt(k+1)+(1.0d0-df(k))*kt(k)
    Ktb(k) = (1.0d0-df(k-1))*kt(k-1)+df(k-1)*kt(k)
    !!!!! Heat and scalar
    Kht(k) = df(k)*alt(k+1)*kt(k+1) +(1.0d0-df(k))*alt(k)*kt(k)
    Khb(k) = (1.0d0-df(k-1))*alt(k-1)*kt(k-1) +df(k-1)*alt(k)*kt(k)
  END DO

  !--------------------------------------------------------------
  ! Get the values of ua, va, ta, qa, pres from input data
  !--------------------------------------------------------------

  ! Interpolate input layer heights to current time
  do i = 1, size(infield%var(infield%id_mla)%f2d, 1)
    call linear_interp_2p( &
      infield%var(infield%id_time)%f1d(input_ptr), &
      infield%var(infield%id_time)%f1d(input_ptr+1), &
      infield%var(infield%id_mla)%f2d(i, input_ptr), &
      infield%var(infield%id_mla)%f2d(i, input_ptr+1), &
      time, input_levels_curr(i) &
      )
  end do

  ! Interpolate input data to current time and model levels
  call interp_time2p_level1d( &
    input_levels_curr, infield%var(infield%id_time)%f1d, &
    infield%var(infield%id_u)%f2d, input_ptr, &
    z, time, u1 &
    )

  call interp_time2p_level1d( &
    input_levels_curr, infield%var(infield%id_time)%f1d, &
    infield%var(infield%id_v)%f2d, input_ptr, &
    z, time, v1 &
    )

  call interp_time2p_level1d( &
    input_levels_curr, infield%var(infield%id_time)%f1d, &
    infield%var(infield%id_t)%f2d, input_ptr, &
    z, time, ta1 &
    )

  call interp_time2p_level1d( &
    input_levels_curr, infield%var(infield%id_time)%f1d, &
    infield%var(infield%id_q)%f2d, input_ptr, &
    z, time, qa1 &
    )

  call interp_time2p_level1d( &
    input_levels_curr, infield%var(infield%id_time)%f1d, &
    infield%var(infield%id_lp)%f2d, input_ptr, &
    z, time, pres &
    )

  ! Convert specific humidity to absolute humidity: rhov = qv * rhoa ~ qv * rhod
  qa1 = pres / (Rgas_d * ta1) * qa1


  !--------------------------------------------------------------------
  !
  !      Estimation of Ri number
  !
  !--------------------------------------------------------------------

  !****new part starts ****

     !do k=1,kz
     !   f2(k) = kt(k) * prod1(k)+vtran(k)
     !enddo


     !dtree = 0.
     !do k = 2, kz
     !   if (z(k).le. z(nz)) then
     !
     !      dtree = dtree + 0.5*(z(k) -z(k-1)) *(f2(k-1) +f2(k)) / f2(nz) / hc
     !   endif
     !enddo

     !ddd1 = max(0., (1-dtree)*hc)

  !**** new part ends ****

  do k=2,kz-1    ! loop 9033

    dtz1=diff(ta1,kz,dz,k)+gamma
    dqz1=diff(qa1,kz,dz,k)
    duz1=diff(u1,kz,dz,k)
    dvz1=diff(v1,kz,dz,k)
    dwz1=diff(w1,kz,dz,k)


!        if(k .le. 2) then   ! old
!           duz1=0.4*u1(2)/log(z(2)/z0)/l(k)
!           dvz1=0.4*v1(2)/log(z(2)/z0)/l(k)
!           dtz1=0.4*(ta1(2)-ta1(1)+gamma*dz(2))/log(z(2)/z0)/l(k)
!           dqz1=0.4*(qa1(2)-qa1(1))/log(z(2)/z0)/l(k)
!        endif


    if (k <= 2) then  ! new
      duz1=ur(2)/l(2)
      dvz1=0.
      dtz1=(((ta1(3)-ta1(1))*dz(2)**2   &
           +(ta1(1)-ta1(2))*(dz(2)+dz(3))**2))/   &
           ( (dz(2)+dz(3))*dz(2)**2 -              &
           dz(2)*(dz(2)+dz(3))**2  )+gamma
      dqz1=(((qa1(3)-qa1(1))*dz(2)**2  &
           +(qa1(1)-qa1(2))*(dz(2)+dz(3))**2))/  &
           ( (dz(2)+dz(3))*dz(2)**2 -   &
           dz(2)*(dz(2)+dz(3))**2  )
    endif

    !**** old ****
    !utt=(duz1**2+dvz1**2)
    !if (utt .eq. 0) utt = 1E-15
    !if(utt.ne.0.) rih(k)=(grav*((dtz1  )/ (ta1(k)+gamma*z(k)) + 1.*0.608*dqz1/roa ) ) / utt
    !**** old ****

    !****new part starts ****

    !abb=min(1.1,l(k)/al1)
    abb=abs(alt(k)*rih(k))/(1.+abs(alt(k)*rih(k)))    ! (neutral case =0) <abb < 1, parameter related to convectivity, related to canopy (bounce or                                                          !  mechanic production)

        !abb1=(max(0.,1.-(kt(k)*prod1(k)+1.*vtran(k)) &
        !     /(kt(k)*prod2(k)+kt(k)*prod1(k)+1.*vtran(k)) )   )


    !if(z(k).lt.ddd1) then
    !       abb1=(max(0.,1.-(kt(k)*prod1(k)+0.*vtran(k)) &
    !            /(kt(k)*prod2(k)+kt(k)*prod1(k)+0.*vtran(k)) )   )
        !endif

        a005=2.     ! const value for neutral case
        a0005=0.

        abb3=abb
        rih1(k)=(a0005-(a0005+1./(1.-c52/c833))*abb3)
        !if(z(k).le.ddd1) then
        !   abb3=abb1
        !   a0005=2.
        !   rih1(k)=(a0005-(a0005+1./(1.-c52/c833))*abb3)
        !endif
    if(dtz1.ge.0.) then

           abb3=a005*alt(k)*rih(k)
           rih1(k)=a005

    endif

        utt=max((duz1**2+dvz1**2)  &
             !+bt2(k)/kt(k)  &
             !+2.*det2(k)/kt(k)   &
             +abs(1.*rih1(k)*( alt(k)*grav*((dtz1  )/   &
             (ta1(k)+gamma*z(k))   &
             + 0.608*dqz1/roa)))  &
             ,0.0000001)

!    if(utt.gt.0.) then
           rih(k)=(   grav*((dtz1  )/(ta1(k)+gamma*z(k) )  + 0.608*dqz1/roa ) ) / utt
!    else
!           rih(k)=0.
!    endif


!**** new part ends ****

        rih(1)=0.
        rih(kz)=0.

!        if(rih(k).gt.0.) then
!           alt1(k)=1.35/(1.+1.35*rih(k))

        if (rih(k).ge.0.0d0) then !new gt -> ge
          alt1(k)=1.35d0/(1.0d0+a005*1.35d0*rih(k))        !new 1.35 -> a005
        else
          alt1(k)=1.35d0*(1.0d0- 15.0d0*rih(k))**0.25d0
        endif

    enddo ! 9033

     !--------------------------------------------------------------------
     ! Solution of the turbulent kinetic energy & dissipation rate equation
     !--------------------------------------------------------------------
     do k=2,kz-1    ! loop 933
        !
        dtz1=diff(ta1,kz,dz,k)+gamma
        dqz1=diff(qa1,kz,dz,k)
        duz1=diff(u1,kz,dz,k)
        dvz1=diff(v1,kz,dz,k)
        dwz1=diff(w1,kz,dz,k)

!        if(k.le.2) then  !old
!           duz1=0.4*u1(2)/log(z(2)/z0)/l(k)
!           dvz1=0.4*v1(2)/log(z(2)/z0)/l(k)
!           dtz1=0.4*(ta1(2)-ta1(1)+gamma*dz(2))/log(z(2)/z0)/l(k)
!           dqz1=0.4*(qa1(2)-qa1(1))/log(z(2)/z0)/l(k)
!        endif

        if(k.le.2) then   !new
           duz1=ur(2)/l(2) ! new
           dvz1=0.        ! new

           dtz1=(((ta1(3)-ta1(1))*dz(2)**2   &
                +(ta1(1)-ta1(2))*(dz(2)+dz(3))**2))/   &
                ( (dz(2)+dz(3))*dz(2)**2 -              &
                dz(2)*(dz(2)+dz(3))**2  )+gamma

           dqz1=(((qa1(3)-qa1(1))*dz(2)**2  &
                +(qa1(1)-qa1(2))*(dz(2)+dz(3))**2))/  &
                ( (dz(2)+dz(3))*dz(2)**2 -   &
                dz(2)*(dz(2)+dz(3))**2  )
        endif

        !abb=min(1.1,l(k)/al1)
        prod1(k) = duz1**2. + dvz1**2.
        prod2(k) = -1.*(alt(k)*grav*(dtz1/(ta1(k) +gamma*z(k)) + 0.608*dqz1/roa))

        abb=abs(alt1(k)*rih(k))/(1.+abs(alt1(k)*rih(k)))    ! (neutral case =0) <abb < 1, parameter related to convectivity, related to canopy (bounce or                                                          !  mechanic production)


        !abb1=(max(0.,1.-(kt(k)*prod1(k)+1.*vtran(k)) &
        !     /(kt(k)*prod2(k)+kt(k)*prod1(k)+1.*vtran(k)) )   )

    !if(z(k).lt.ddd1) then
    !       abb1=(max(0.,1.-(kt(k)*prod1(k)+0.*vtran(k)) &
    !            /(kt(k)*prod2(k)+kt(k)*prod1(k)+0.*vtran(k)) )   )
!
 !       endif
        abb3 = abb
        cc22=c52+(c833-c52)*abb3
        cc33=c833
        !a0005 = 0.
        rih1(k)=(a0005-(a0005+1./(1.-c52/c833))*abb3)


        !if(z(k) .le. ddd1) then
        !   abb3 = abb1
        !   a0005 = 2.
        !   rih1(k)=(a0005-(a0005+1./(1.-c52/c833))*abb3)
        !   cc22=c52+(c833-c52)*abb3
        !endif
        beta=(c52-c833)*rih1(k)


        if(rih(k).ge.0.) then
           abb3 = a005*alt(k)*rih(k)
           rih1(k) = a005
           cc22=c52+(c833-c52)*abb3
           cc33 = c833
           beta=(c52-c833)*rih1(k)
        endif

        c15=(12.*cc2**0.5)*cd*s1(k)*( sqrt(u1(k)**2+v1(k)**2)) ! if c16 in use, 10*cc2 is used for c15, if c16 not in use, then 12*c22 is used.

        bt2(k) = (cd*s1(k) * (sqrt(u1(k)**2 + v1(k)**2))**3)*(0. + abb3**0.5)

!**** new part ends ****

        fktt=2.*alf*(df(k)*kt(k+1)+(1.-df(k))*kt(k))
        fktd=2.*alf*((1.-df(k-1))*kt(k-1)+df(k-1)*kt(k))

        if(k.le.2) then    ! new
           fktt=2.*alf*sqrt(kt(k+1)*kt(k)) !new
           ! geometrical mean provides the smallest weight possible, according to Andrey most stable
           fktd=2.*alf*sqrt(kt(k-1)*kt(k))  ! new
        endif    ! new


        !cc22=(c52+(c833-c52)*l(k)/al1)

        !c15=(12.*cc2**0.5)*s(k)*( sqrt(u1(k)**2+v1(k)**2)+sqrt(u1(k)**2+v1(k)**2))/2.
        !cc33=(c833-(c833-c52)*(c15)/dbt(k))

        !beta0=(cc22-c833)


        ! Andrey new addings on 2013.09.17

        vtran(k) = bt(k+1) * fktt/da(k)    - bt(k)*(dz(k)*fktt + dz(k+1)*fktd)/(dz(k+1) + dz(k)) /db(k) + bt(k-1)*fktd/dc(k)

        det2(k) = c15 * bt(k)

        terma(k)= dt_mete*(fktt -w(k)*dz(k))/da(k)
        termb(k)= dt_mete*((dz(k)*fktt+ dz(k+1)*fktd) /(dz(k+1)+dz(k)) +w(k)*(dz(k+1)-dz(k)))/db(k)+ 1.0d0+dt_mete*dbt(k)
        termc(k)= dt_mete*(fktd+ w(k)*dz(k+1))/dc(k)
        !old: termd(k)= bt(k) + dt_mete*kt(k) *(duz1**2+dvz1**2+dwz1**2-alt1(k)*grav*(dtz1/(ta1(k)+gamma*z(k)) + 1.*0.608*dqz1/roa))-0.*(bt(k)-bnud(k))*c_nud !nudging!
        termd(k)= bt(k) + dt_mete*kt(k)*(prod1(k)+prod2(k))

        !old: bg(k)= dt_mete*((dz(k)*fktt+ dz(k+1)*fktd) /(dz(k+1)+dz(k)) +w(k)*(dz(k+1)-dz(k)))/db(k)+ 2.*dt_mete*cc33*dbt(k)+1.0d0
        bg(k)= dt_mete*((dz(k)*fktt+ dz(k+1)*fktd) /(dz(k+1)+dz(k)) +w(k)*(dz(k+1)-dz(k)))/db(k)+ 2.*dt_mete*cc33*dbt(k)+1.0d0 &
             -1.*dt_mete*(c52-cc33)*bt2(k)/bt(k) + 1.*dt_mete*(c52-cc33)*c15

        !old : dg(k)= dbt(k) +dt_mete*cc33*dbt(k)**2 +dt_mete*cc22*cc2 *(duz1**2+dvz1**2+dwz1**2) &
        !     -beta0*dt_mete*cc2*(1.*alt1(k)*grav*(dtz1/(ta1(k)+gamma*z(k))+ 1.*0.608*dqz1/roa ) )
       dg(k)= dbt(k) + dt_mete*cc33*dbt(k)**2 + dt_mete*cc2*(cc22*prod1(k) +beta*prod2(k))

    enddo  ! end of loop 933

     if(abl == 1) then
        bttop=temgrad**2/(cc2**0.5)
        !old CALL gtri(terma,termc,termb,termd,bt1,2,0.d0,3.d0,1,bttop,3.d0,kz,2)
        CALL gtri(terma,termc,termb,termd,bt1,2,0.d0,3.d0,2,0.d0,3.d0,kz,2)
    else

        btbot=ur(1)**2/(cc2**0.5) ! new
        !old: CALL gtri(terma,termc,termb,termd,bt1,2,0.d0,3.d0,2,0.d0,3.d0,kz,2)
        CALL gtri(terma,termc,termb,termd,bt1,1,btbot,3.d0,2,0.d0,3.d0,kz,2)  ! new

     endif

!tke_limit = 1.d-5  ! 1.d-2 in manitou
     do k=1,kz !-1
        !old bt1(k)=max(bt1(k),1.d-5)
        bt1(k)=max(bt1(k),1.d-2) !new
     enddo


     !old CALL gtri(terma, termc, bg, dg, dbt1, 1, (cc2**(0.75)*bt1(2)**0.5/l(1)), 3.d0, 1, 1.*(cc2**(0.75)*bt1(kz)**0.5/l(kz)), 3.d0, kz, 2)
     CALL gtri(terma, termc, bg, dg, dbt1, 1, (cc2**(0.75)*btbot**0.5/l(1)), 3.d0, 2, 0.d0, 3.d0, kz, 2)  !new

     ! limiting disipational rate related to eddy
    do k=1,kz
        dbt1(k) = max(dbt1(k), cc2**0.75*(bt1(k))**0.5/(al1))
     enddo

     !-------------------------------------------------------------
     ! solution for ABL height
     !-------------------------------------------------------------

    call find_boundary_layer(z, rih, nz + 1, pblh, pblh_ind)

!old      if(abl /= 1) al1=max(1.*0.00027*abs(windx)/abs(f1),5.d0)  ! was used for netural cases

     !----------------------------------------------------------
     !  relation for the eddy exchange coefficient of momentum
     !----------------------------------------------------------

!old     kt1(1)=max(cc2*bt1(1)/dbt1(1),1.d-3)
     kt1(1)=max(cc2*bt1(1)/dbt1(1),ur(1)*l(1)) !new
     kt1(kz)=max(cc2*bt1(kz)/dbt1(kz),kt1(1))
 !dublicate    l1(1)=l(1)  !new
     l1(kz)=max(sqrt(bt1(kz))*cc2**0.75/dbt1(kz),l1(1))

     do k=2,kz-1   !  loop 291
        gf=abs(bt1(k))
        gf2=dbt1(k)
        kt1(k)=max(cc2*gf/gf2,kt1(1))
        l1(k)=max(sqrt(gf)*cc2**0.75/gf2,l1(1))
     enddo  ! end of loop 291

     if(abl == 1) then
        do k=2,kz
           if(z(k).gt.hc) l1(k)=l1(k-1)+0.43*dz(k)
        enddo
     endif

     do  k=2,kz-1
        !old   fktt=((df(k)*kt1(k+1) +(1.-df(k))*kt1(k)))
        !old   fktd=(((1.-df(k-1))*kt1(k-1) +df(k-1)*kt1(k)))
        !old   ur(k)=0.5*(sqrt(fktt*abs(sqrt(u1(k+1)**2+v1(k+1)**2) - sqrt(u1(k)**2+v1(k)**2))/dz(k+1))+  &
        !old        sqrt(fktd*abs(sqrt(u1(k)**2+v1(k)**2) - sqrt(u1(k-1)**2+v1(k-1)**2))/dz(k)))



        pa1=(kt1(k+1)+kt1(k))*dz(k)/da(k)/2.  !new
        pa0=(kt1(k)+kt1(k-1))*dz(k+1)/dc(k)/2. !new

        akt_duz=(pa1*(u1(k+1)-u1(k))+pa0*(u1(k)-u1(k-1)))  !new

        akt_dvz=(pa1*(v1(k+1)-v1(k))+pa0*(v1(k)-v1(k-1)))  !new

        ur(k)=(akt_duz**2+akt_dvz**2)**0.25   !new
        ! same thing, but speeds up calculations
     enddo

     !old ur(1)=ur(2)
     ur(1)=(ur(2)*(dz(2)+dz(3))**2-ur(3)*dz(2)**2)/( (dz(2)+dz(3))**2-dz(2)**2) !new cubic parameterization for first level (so that there is no gradient)
     ur(kz)=ur(kz-1)

!**** new part starts ****

     sb=0.
     szb=0.


     do k=2,kz
        szb=(sqrt(abs(bt1(k-1)-0.01))*z(k-1)+ &  ! 0.01 is the tke limit
        sqrt(abs(bt1(k)-0.01))*z(k))*   & ! 0.01 is the tke limit
        dz(k)/2.+szb

        sb=(sqrt(abs(bt1(k-1)-0.01))+sqrt(abs(bt1(k)-0.01)))*  & ! 0.01 is the tke_limit
        dz(k)/2.+sb
     enddo
     !      al0=al1
     if (sb .eq. 0) then
        write(*,*) 'error: divided by zero. sb =', sb
     endif
          al100=0.00027*abs(windx)/abs(f1)  ! blackadar value for netural case
     al1=max(al100,0.075*(szb/sb) )  ! convective case

! szb/sb the relation of the integrals above

!
!**** new part ends ****

     nturb=nturb+1
  !--------------------------------------------------------------
  ! Radiation
  !---------------------------------------------------------------

  ! dhsn, dvsn, dhsd, dvsd

  !--------------------------------------------------------------
  ! Heat fluxes
  !---------------------------------------------------------------

  DO k=1,kz
    if (sl(k) == 0.0d0) then
      tsn1(k )=ta1(k )
      tsd1(k )=ta1(k )
      qsn1(k )=qa1(k )
      qsd1(k )=qa1(k )
    end if
  ENDDO

  DO k=1,nz-1
     fluxle(k)=fluxle(1)
     fluxh(k)=fluxh(1)

     if(sl(k).ne.0.)  &

          tsn1(k)=  &        !  with inertia in leaves i
          (rasn(k) +face*roa*cpa*dhsn(k)*ta1(k)  &
          + 3.*delf*2.*sigma_SB*tsn(k)**4  &
          -face*lv*dvsn(k)*(qsn(k)-qa1(k))  &
          +ch2o*rh2o*hleaf*tsn(k)/dt_mete )/  &
          (face*roa*cpa*dhsn(k)+4.*delf*2.*sigma_SB*tsn(k)**3  &
          +ch2o*rh2o*hleaf/dt_mete   )
     if(sl(k).ne.0.)  &

          tsd1(k)=   &       !  with inertia in leaves i
          (rasd(k) +face*roa*cpa*dhsd(k)*ta1(k)  &
          + 3.*delf*2.*sigma_SB*tsd(k)**4  &
          -face*lv*dvsd(k)*(qsd(k)-qa1(k))  &
          +ch2o*rh2o*hleaf*tsd(k)/dt_mete )/  &
          ( face*roa*cpa*dhsd(k)+4.*delf*2.*sigma_SB*tsd(k)**3  &
          +ch2o*rh2o*hleaf/dt_mete   )  ! same as for surface temperature
  ENDDO

  DO k=2,kz-1
     fktt=(df(k)*alt1(k+1)*kt(k+1)+(1.-df(k))*alt1(k)*kt(k))*dz(k)  !new
     fktd=((1.-df(k-1))*alt1(k-1)*kt(k-1)+df(k-1)*alt1(k)*kt(k))*dz(k+1)  !new

     fluxh3(k)=-cpa*roa*(((ta1(k)-ta1(k-1))/dz(k)+gamma)*fktd+((ta1(k+1)-ta1(k))/dz(k+1)+gamma)*fktt)/(dz(k)+dz(k+1)) !new
     fluxle3(k)=-lv*(((qa1(k)-qa1(k-1))/dz(k))*fktd+((qa1(k+1)-qa1(k))/dz(k+1))*fktt)/(dz(k)+dz(k+1)) !new
  ENDDO

  fluxle3(1)=(fluxle3(2)*(dz(2)+dz(3))**2-fluxle3(3)*dz(2)**2)/( (dz(2)+dz(3))**2-dz(2)**2) ! new
  fluxh3(1)=(fluxh3(2)*(dz(2)+dz(3))**2-fluxh3(3)*dz(2)**2)/( (dz(2)+dz(3))**2-dz(2)**2)  !new

  fluxle3(kz)=fluxle3(kz-1)
  fluxh3(kz)=fluxh3(kz-1)

  tasoil=ta1(1)
  qasoil=qa1(1)

  fluxles =fluxle3(1)  !new
  fluxle(1) =fluxle3(1)  !new
  fluxhs =fluxh3(1)     !new
  fluxh(1) =fluxh3(1)   !new

  sks00= (alt1(1)*kt1(1)+alt1(2)*kt1(2))/2.  !new

  haags=cpa*roa*sks00/dz(2)     !new

  !============ Surface radiation ===============!
  ff1=fnid(1)+fphd(1)+rsnt*psn(1)+ rskt*psk(1)+fird(1) -fniu(1)-fphu(1)
  !tar=ta1(kmix-1)  ! not neeed for anything

  !r soil heat flux
  if ( ISNAN(local_gsoil(nxodrad, 1)) .or. ISNAN(local_gsoil(nxodrad+1, 1)) ) then  ! measurement values missing
      ! based on modelled soil temperatures
      pp=0.5*(3.*tsoil1(10+nmix)-4.*tsoil1(10+nmix-1)+tsoil1(10+nmix-2))/0.10 ! new
  else !- CLEARCUT COMMENT OUT
      ! based on the measurements
      pp = linear_interp(time_in_month, nxodrad, local_gsoil(:, 1), dt_obs)  ! [W m-2], ground heat flux
  endif !- CLEARCUT COMMENT OUT

  kmix = 2  ! ?? temporary solution, should be modified to a more general code
  firu(1)=dels*sigma_SB*ta1(kmix-1)**4+(1.-dels)*fird(1)

  do k=1,kz-1   ! loop 5647
     fluxle(k+1)=fluxle(k)+face*0.5*dz(k+1)* (sl(k+1)*(psn(k+1)*lv*dvsn(k+1)*(qsn(k+1) -qa1(k+1))+  &
          (1.-psn(k+1))*lv*dvsd(k+1)*(qsd(k+1)-qa1(k+1))) +sl(k)*(psn(k)*lv*dvsn(k)*(qsn(k)-qa1(k))+  &
          (1.-psn(k))*lv*dvsd(k)*   (qsd(k)-qa1(k))))   ! latent heat flux

     fluxh(k+1)=fluxh(k)+face*0.5*dz(k+1)* (sl(k+1)*(psn(k+1)*cpa*roa*dhsn(k+1) *(tsn1(k+1)-ta1(k+1))+  &
          (1.-psn(k+1))*cpa*roa*dhsd(k+1) * (tsd1(k+1)-ta1(k+1))) +sl(k)*(psn(k)*cpa*roa*dhsn(k)*(tsn1(k)-ta1(k))+  &
          (1.-psn(k))*cpa*roa*dhsd(k)* (tsd1(k)-ta1(k))))  ! sensible heat flux
  enddo   ! end of 5647

  ! energy balance at the top of the canopy
  balans=fird(nz)+rsnt+rskt -1.*fniu(nz) -firu(nz)-1.*fphu(nz)-pp- 1.*fluxle3(nz)-1.*fluxh3(nz)-0.*tau-0.*tau4   ! C4 & C1 Sog. 2002
  ! energy balance at the surface
  balans1=fnid(1)+fphd(1)+rsnt*psn(1)+ rskt*psk(1)+fird(1) -1.*fniu(1) -firu(1)-1.*fphu(1)-pp- fluxles-fluxhs-1.*tau-tau4   ! C4 & C1 Sog. 2002

  ! pause
  Other_out(:, 1) = ((wee-tniu-rniu)*(1.-phsn)+(wee-tphu-rphu)*phsn)*rsnt*gl/cos_zenith &
                    + ((wee-tniu-rniu)*rskt*(1.-phsk)+ (wee-tphu-rphu)*rskt*phsk)*psk*gd/cos_zenith &
                    + (wee-tniu-rniu)*(fnid+ fniu) &
                    + (wee-tphu-rphu)*(fphd+fphu)  ! Qabs_sn
  Other_out(:, 2) = ((wee-tniu-rniu)*rskt*(1.-phsk) + (wee-tphu-rphu)*rskt*phsk)*psk*gd/cos_zenith &
                    + (wee-tniu-rniu)*(fnid+fniu) &
                    + (wee-tphu-rphu)*(fphd+fphu)  ! Qabs_sd
  Other_out(:, 3) = delf*(fird+firu)  ! Labs
  Other_out(:, 4) = 2.*delf*sigma_SB*tsn1**4  ! Lout_sn
  Other_out(:, 5) = 2.*delf*sigma_SB*tsd1**4  ! Lout_sd
  Other_out(:, 6) = face*cpa*roa*dhsn*(tsn1-ta1)  ! [W (proj. leaf m2)], SH_sn
  Other_out(:, 7) = face*cpa*roa*dhsd*(tsd1-ta1)  ! [W (proj. leaf m2)], SH_sd
  Other_out(:, 8) = face*lv*dvsn*(qsn-qa1)  ! [W (proj. leaf m2)], LH_sn
  Other_out(:, 9) = face*lv*dvsd*(qsd-qa1)  ! [W (proj. leaf m2)], LH_sd
  Other_out(:, 10) = fphu + fniu  ! Su_out
  Other_out(:, 11) = fphd + fnid  ! Sd_out
  Other_out(:, 12) = gl
  Other_out(:, 13) = gd
  Other_out(:, 14) = psn  ! eta
  Other_out(:, 15) = psk  ! mu
  Other_out(:, 16) = rads2  ! radiation on top of the canopy
  Other_out(:, 17) = rsnt  ! direct radiation on top of the canopy
  Other_out(:, 18) = rskt  ! diffuse radiation on top of the canopy
  Other_out(:, 19) = cos_zenith  ! cos(cos_zenithh angle)
  Other_out(:, 20) = shre
  Other_out(:, 21) = re
  Other_out(:, 22) = dhsn
  Other_out(:, 23) = dhsd
  Other_out(:, 24) = dvsn
  Other_out(:, 25) = dvsd
  Other_out(:, 26) = kt1
  Other_out(:, 27) = qsn
  Other_out(:, 28) = qsd
  Other_out(:, 29) = fird
  Other_out(:, 30) = firu
  Other_out(:, 31) = pp
  Other_out(:, 32) = balans
  Other_out(:, 33) = balans1
  ! Other_out(:, 22) = qa1   ! absolute humidity of air

  ! *** AIR molecules ********************************************************************************************
  ! pres(1) = linear_interp(time_in_month, nxodrad, local_pres(:, 1), dt_obs)  ! [Pa]
  ! do k = 2,kz
  !   pres(k)=pres(k-1) * exp(-9.81 * (z(k)-z(k-1))/Rgas_d/(0.5*(ta1(k)+ta1(k-1))))
  ! enddo

  ! Calculation of RH and water vapour pressure used in emission modules
  DO k = 1,kz

     ! Calculate saturation vapor pressure of water
     ! based on Seifield & Pandis, p765, unit in Pascal!!!!!!!!!!!!
     EM_ES(k) = svp(ta1(k))

     ! Water vapour pressure in Pa - different options, Pv = rhovRT/Mv
     EM_EW(k) = (qa1(k)/(18.0153d0/1000.0d0)) * 8.314d0 * ta1(k)
     IF (EM_EW(k)>EM_ES(k)) THEN
       EM_EW(k)=EM_ES(k)
     ENDIF
     ! Water vapour mixing ratio (kg/kg) (based on Megan conversion function), rhov/rhod = PvMv/(PdMd)
     EM_WVM(k) = EM_EW(k) * 18.0153d0/ pres(k) / 28.96d0
     EM_WVM(k) = EM_WVM(k) * 1e3_dp  ! [kg kg-1] --> [g kg-1]
     ! RH(k) = 100.0d0* qa1(k)/ rhovs_buck(ta1(k))
     RH(k) = 100.0d0 * EM_EW(k) / EM_ES(k)

     ! Water vapour pressure in Pa - different options
     !EM_EW(k)  = RH(k) * EM_ES(k) / 100  ! based on RH

     ! Limit EM_EW to maximum value of ES
     !IF (EM_EW(k) .GT. EM_ES(k)) THEN
     !   EM_EW(k) = EM_ES(k)
     !ENDIF

  ENDDO

  ! Replace old values with new values
  CALL New2Old(kz, wg, wg1, tsoil, tsoil1, l, l1,  bt, dbt, bt1, dbt1, &
               u,u1, v,v1, w, w1, kt, kt1, ta, ta1, &
               alt, alt1, tsn, tsn1, tsd, tsd1, qa, qa1)

end subroutine update_meteorology_traj


subroutine update_meteorology_stat(dt_obs)

  REAL(KIND = dp), INTENT(IN   ) ::  dt_obs ! temporal resolution of surface observations (= input data) [s]
  REAL(KIND = dp), DIMENSION(kz) ::  lai_cumulative ! cumulative leaf area index (m^2 / m^2)
  REAL(dp) :: tempv  ! temporary value used for NAN readings
  REAL(dp) :: B_fac,     &
              EoT,   &!Equation of time
              TC,    &!Time Correction factor
              LST,   &!Local Solar Time
              LSN,   &!Local Solar Noon
              UTC

  ! This part seems to be used to calculate something related to snow, but it
  ! is not used elsewhere, so Putian commented it out
  ! s(1:kz) = 0.0d0
  ! do k=1,kz
  !   if (rou.le.z(k)) then! rou = 0 snow cover excluded for now
  !     s(k)=s1(k)*cd+s(k)
  !   else
  !     krelf=k-1
  !     s(k)=profile
  !   end if
  ! end do

  !----------------------------------------------------------------------------!
  ! Calculate the incoming total solar radiation at the top of atmosphere
  ! (rads) to see if it is in the day (rads>0) or night (rads<=0), then read
  ! the incoming radiation from input accordingly.
  !----------------------------------------------------------------------------!
  ! CALL RadiationIntoATM(time_in_day, julian, lat_rad, wsun, cos_zenith, rads)
  rads = solar_constant * cos_zenith

  ! Now new cos_zenith is kind of correct, now I need to find how to deal with
  ! the nighttime ??
  ! declination = get_solar_declination(real(julian_UTC, dp))
  ! write(*,*) 'declination old, new: ', wsun*rad2deg, declination
  ! write(*,*) 'cos_zenith old: ', cos_zenith
  ! cos_zenith = get_cos_zenith( time_in_day_UTC, real(julian_UTC, dp), &
  !   lat_deg, lon_deg)
  ! write(*,*) 'cos_zenith new: ', cos_zenith

  IF (rads > 0.) THEN
    IF ( TRIM(ADJUSTL(STATION)) == 'FCT' ) THEN  ! calculating incoming solar radiation - define here if the case for your station!
      CALL TotalSWRAboveCanopy(cos_zenith, rads, rads2)
      CALL SWRDirectDiffuse(rads, rads2, cos_zenith, rsnt, rskt)
    ELSE IF ( TRIM(ADJUSTL(STATION)) == 'Manitou' ) THEN  ! at the Welgegund station only total global radiation measured
      CALL TotalSWRAboveCanopy(cos_zenith, rads, rads2)
      ! the measured total global radiation above the canopy
      rads2 = linear_interp(time_in_month, nxodrad, dirtop1, dt_obs)

      CALL SWRDirectDiffuse(rads, rads2, cos_zenith, rsnt, rskt)
    ELSE  ! using input values for direct and diffuse global radiation - default option!
      rsnt = linear_interp(time_in_month, nxodrad, dirtop1, dt_obs)
      rskt = linear_interp(time_in_month, nxodrad, difftop1, dt_obs)
    END IF
  ELSE
     rsnt = 0.0_dp
     rskt = 0.0_dp
  END IF

  !Time correction because Hyytiala is not exactly lying on UTC+2
  B_fac=(360.0/365.0*(julian-81.0)*PI)/180.0
  EoT=9.87*sin(2.0*B_fac)-7.53*cos(B_fac)-1.5*sin(B_fac)
  TC=4.0*(EM_Lat*180/PI-UTC*15.0) + EoT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  tata=tan(lat_rad)*tan(wsun)

  ! Old way to calculate sunrise and sunset
  upsn=c432/pi*atan(sqrt(1.-tata**2)/tata)
  downsn=SECONDS_IN_ONE_DAY*(1.-1./PIx2*atan(sqrt(1.-tata**2)/tata))


  ! new way, consider to use this in future
  ! upsn=c432/pi*atan(sqrt(1.-tata**2)/tata)-TC*60.0
  ! downsn=SECONDS_IN_ONE_DAY*(1.-1./PIx2*atan(sqrt(1.-tata**2)/tata))-TC*60.0
  ! write(*,*) 'upsn old', c432/pi*atan(sqrt(1.-tata**2)/tata)
  ! write(*,*) 'upsn new', upsn
  ! write(*,*) 'downsn old', SECONDS_IN_ONE_DAY*(1.-1./PIx2*atan(sqrt(1.-tata**2)/tata))
  ! write(*,*) 'downsn new', downsn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  sumrad=solar_constant*cos_zenith
  clo=min(1.d0,(clob+clom+clot))
  tran=(0.6+0.2*cos_zenith)*(1.-0.4*clob)*(1-0.7*clom)*(1.-0.4*clot)
  clo=0.
  emm2=(1.-clo)*0.0000092*(ta1(nz )**2) +clo
  emm2=max(7.2d-1,emm2)

  if(abl == 1) then

    dtmet=1800.
    nxodmet=int(time_in_month/dtmet)+1
    shift=(time_in_month-(nxodmet-1)*dtmet)/dtmet

    tax=border(1,nxodmet) +shift*(border(1,nxodmet+1)-border(1,nxodmet))


    temgrad=border(5,nxodmet) +shift*(border(5,nxodmet+1)-border(5,nxodmet))

    temgrad=max(temgrad,1.d-2)

    windx=border(2,nxodmet) +shift*(border(2,nxodmet+1)-border(2,nxodmet))

    windy=border(3,nxodmet) +shift *(border(3,nxodmet+1)-border(3,nxodmet))

    qax=qah
    qax=border(4,nxodmet) +shift *(border(4,nxodmet+1)-border(4,nxodmet))

    windy=0.

  else
    ! linearly interpolating for the current time step
    dtmet=10800.0_dp

    CALL interp_input_to_now(time_in_month, dtmet, border_abl(1, :), tax)
    CALL interp_input_to_now(time_in_month, dtmet, border_abl(2, :), windx)
    CALL interp_input_to_now(time_in_month, dtmet, border_abl(3, :), windy)
    CALL interp_input_to_now(time_in_month, dtmet, border_abl(4, :), qax)
    CALL interp_input_to_now(time_in_month, dtmet, border_abl(5, :), temgrad)
    temgrad = -1.0_dp * temgrad

    ! windy=0.0_dp
  endif

  ! LWR data from reanalysis: linearly interpolating for time of the current time step
  dtmet=10800.0_dp  ! time difference btw data points 3h = 10800s

  ! ECMWF data in UTC time, so
  ! time = 0s      -> LWR_in(1)
  ! time = 3600s   -> LWR_in(2)
  ! time = 4*3600s -> LWR_in(3)
  ! etc...
  ! Putian has changed this. Now the input ECMWF data are interpolated to local time in the input files.
  ! nxodmet= floor((time_in_month - 1*3600)/dtmet) + 2 ! "row to read" (index) from the input data
  ! shift=(time_in_month +2*3600 -(nxodmet-1)*dtmet )/dtmet
  ! LWRdown = LWR_in(nxodmet)+shift*(LWR_in(nxodmet+1)-LWR_in(nxodmet))
  CALL interp_input_to_now(time_in_month, dtmet, ECMWF_strd(:), LWRdown)

  temgrad = min(temgrad,1.d-1)

  emm=emm2

  !r soil data: linearly interpolating to current time step
  !r hyy_soil(2,:)  Volumetric water content (m3/m3) in organic layer in -5-0 cm
  surf_soil_hum = linear_interp(time_in_month, nxodrad, stat_sm_1, dt_obs)

  rads2=rskt+rsnt

    if(time_in_day.le.upsn) rsnt=0.
    if(time_in_day.le.upsn) rskt=0.
    if(time_in_day.ge.downsn) rsnt=0.
    if(time_in_day.ge.downsn) rskt=0.
    if(time_in_day.le.upsn) cos_zenith=1.
    if(time_in_day.le.upsn) cos_zenith=1.
    if(time_in_day.ge.downsn) cos_zenith=1.
    if(time_in_day.ge.downsn) cos_zenith=1.
    if(rads2.lt.10.) cos_zenith = 1.

  do k=1,kz    ! loop 223
    qsn(k)=max(qa(k), rhovs_buck(tsn(k)))
    qsd(k)=max(qa(k), rhovs_buck(tsd(k)))
    if(qa(k).ge.qsn(k)) qa(k)=qsn(k)
  enddo

  cor1=0.
  c_nud=0.

  if (abl /= 1) then

    cor1=1.0d0
    c_nud=0.0d0   ! nudged to diurnal cycle
    c_nud1= 0.01d0  ! 0.01, 0.1    !CLEARCUT - IF WE DO NOT NUDGE, THEN THIS SHOULD BE ZERO - else nudging factor is 0.01

    ! calculating the nudgnig variables from observations:  tnud, qnud, unud, vnud
    ! CALL MakeNudgingData(kz, z, time_in_month, nxodrad, STATION, dt_obs, &
    !   tem, hum, uwind, vwind, &
    !   ta, qa, u, v, &
    !   tnud, qnud, unud, vnud, windx)

    ! u, v, ta, qa
    CALL make_nudging_data_new(STATION, 'spline', kz, z, time_in_month, u, 5, loclv_uwind, nxodrad, dt_obs, local_uwind, unud, (/z(1), z(kz)/), (/0.0_dp, windx/))
    CALL make_nudging_data_new(STATION, 'spline', kz, z, time_in_month, v, 5, loclv_vwind, nxodrad, dt_obs, local_vwind, vnud, (/z(1), z(kz)/), (/0.0_dp, windy/))
    CALL make_nudging_data_new(STATION, 'linear', kz, z, time_in_month, ta, 7, loclv_temp, nxodrad, dt_obs, local_temp, tnud)
    CALL make_nudging_data_new(STATION, 'linear', kz, z, time_in_month, qa, 7, loclv_rhov, nxodrad, dt_obs, local_rhov, qnud)
  endif

  !--------------------------------------------------------------
  ! Get Ktt, Ktb and Kht, Khb so every diffusion term can be calculated more easily
  !--------------------------------------------------------------
  DO k=2,kz-1
    !!!!! Momentum
    Ktt(k) = df(k)*kt(k+1)+(1.0d0-df(k))*kt(k)
    Ktb(k) = (1.0d0-df(k-1))*kt(k-1)+df(k-1)*kt(k)
    !!!!! Heat and scalar
    Kht(k) = df(k)*alt(k+1)*kt(k+1) +(1.0d0-df(k))*alt(k)*kt(k)
    Khb(k) = (1.0d0-df(k-1))*alt(k-1)*kt(k-1) +df(k-1)*alt(k)*kt(k)
  END DO

  !--------------------------------------------------------------
  !   solution of U and V - velocity component equations
  !--------------------------------------------------------------

    do k=2,kz-1   ! loop 93
       fktt=2.0d0*Ktt(k)
       fktd=2.0d0*Ktb(k)

       terma(k)= dt_mete*(fktt-w(k)*dz(k))/da(k)

       b1(k)= dt_mete*((dz(k)*fktt+ dz(k+1)*fktd) /(dz(k+1)+dz(k)) +w(k)*(dz(k+1)-dz(k)))/db(k)+ &
            dt_mete*(cd*s1(k)*sqrt(u(k)**2 +v(k)**2+w(k)**2)) +1.0d0

       termc(k)= dt_mete*(fktd+w(k)*dz(k+1))/dc(k)

       !d1(k)= u(k) + dt_mete*(cor1*f1*(v(k)-windy)) -24.*dt_mete*(u(k)-windx)/SECONDS_IN_ONE_DAY       &
       !     *c_nud*(0.5+0.5*sin(0.5*pi*(2.*z(k)/z(kz)-1.)))&
       !     -(u(k)-unud(k))*c_nud1   ! nudging

       d1(k)= u(k) + dt_mete*(cor1*f1*(v(k)-windy)) - dt_mete*(u(k)-windx)/SECONDS_IN_ONE_DAY*c_nud &
       ! *(0.5+0.5*sin(0.5*pi*(2.*z(k)/z(kz)-1.)))&
       -(u(k)-unud(k))*c_nud1     ! nudging


       !d2(k)= v(k)+ dt_mete*(-cor1*f1*(u(k)-u(kz)))-24.*dt_mete*(v(k)-0.)/SECONDS_IN_ONE_DAY   &
       !     *c_nud*(0.5+0.5*sin(0.5*pi*(2.*z(k)/z(kz)-1.)))&
       !     -(v(k)-vnud(k))*c_nud1       ! nudging


       d2(k)= v(k)+ dt_mete*(-cor1*f1*(u(k)-windx)) - dt_mete*(v(k)-0.)/SECONDS_IN_ONE_DAY*c_nud &
       !  *(0.5+0.5*sin(0.5*pi*(2.*z(k)/z(kz)-1.))) &
       -(v(k)-vnud(k))*c_nud1     ! nudging
    enddo !    end of loop 93

    !old      call gtri(terma,termc,b1,d1,u1,1,0.d0,3.d0,2,0.*windx,3.d0,kz,2)
    call gtri(terma,termc,b1,d1,u1,1,0.d0,3.d0,1,1.0d0*windx,3.d0,kz,2)  !new (not much difference in the results)
    call gtri(terma,termc,b1,d2,v1,1,0.d0,3.d0,1,windy,3.d0,kz,2)


     !----------------------------------------------------------
     !  solution of the heat and moisture equations
     !----------------------------------------------------------

     do k=2,kz-1   ! loop 950
        fktt=2.0d0*Kht(k)
        fktd=2.0d0*Khb(k)

        dkz=(fktt-fktd)/(dz(k)+dz(k+1)) !new

        terma(k)= dt_mete*(fktt-w(k)*dz(k))/da(k)
        termc(k)= dt_mete*(fktd+w(k)*dz(k+1))/dc(k)

        b1(k)= dt_mete*((dz(k)*fktt+ dz(k+1)*fktd) /(dz(k+1)+dz(k)) +w(k)*(dz(k+1)-dz(k)))/db(k)+ &
               face*dt_mete*s1(k)*(psn(k)*dhsn(k)+(1.-psn(k))*dhsd(k))+ 1.0d0

!        d1(k)= ta(k) + dt_mete*(face*s1(k)*(psn(k)*dhsn(k)*tsn(k)+ (1.-psn(k))*dhsd(k)*tsd(k)) +gamma*dkz-w(k) *gamma )    &
!             +24*dt_mete*(tax+temgrad*(z(kz)-z(k))-ta(k))/(SECONDS_IN_ONE_DAY)       &
!             *c_nud*(0.5+0.5*sin(0.5*pi*(2.*z(k)/z(kz)-1.)))        &
!             -(ta(k)-tnud(k))*c_nud1 ! nudging
        d1(k)= ta(k) + dt_mete*(face*s1(k)*(psn(k)*dhsn(k)*tsn(k)+ (1.-psn(k))*dhsd(k)*tsd(k)) +gamma*dkz-w(k) *gamma )    &
             + c_nud*dt_mete*(tax+temgrad*(z(kz)-z(k))-ta(k))/(SECONDS_IN_ONE_DAY)       &
                                !       *(0.5+0.5*sin(0.5*pi*(2.*z(k)/z(kz)-1.)))        &
             -(ta(k)-tnud(k))*c_nud1 ! nudging , for temperature

        !  2 - for moisture
        b2(k)= dt_mete*((dz(k)*fktt+ dz(k+1)*fktd) /(dz(k+1)+dz(k)) +w(k)*(dz(k+1)-dz(k)))/db(k)+  &
               face*dt_mete*s1(k)*(psn(k)*dvsn(k)+(1.-psn(k))*dvsd(k))+ 1.0d0
        !
        !  dvsn - dhsn; dvsd - dhsd; qsn - tsn; qsd - qsn
        !
!        d2(k)= qa(k) + dt_mete*(face*s1(k)*(psn(k)*dvsn(k)*qsn(k)+ (1.-psn(k))*dvsd(k)*qsd(k))) &
!             -24.*dt_mete*(qa(k)-qax)/(SECONDS_IN_ONE_DAY)*c_nud*(0.5+0.5*sin(0.5*pi*(2.*z(k)/z(kz)-1.))) &
!             -(qa(k)-qnud(k))*c_nud1                        ! nudging
        d2(k)= qa(k) + dt_mete*(face*s1(k)*(psn(k)*dvsn(k)*qsn(k)+ (1.-psn(k))*dvsd(k)*qsd(k))) &
             - dt_mete*(qa(k)-qax)/(SECONDS_IN_ONE_DAY)*c_nud*(0.5+0.5*sin(0.5*pi*(2.*z(k)/z(kz)-1.))) &
             -(qa(k)-qnud(k))*c_nud1                        ! nudging
             ! -(qa(k)-qnud(k))*0.10                            ! added by Tian due to bad agreement with observation

     enddo    ! end of loop 950

     call gtri(terma,termc,b1,d1,ta1,1,ta(1),3.d0,1,tax,3.d0,kz,2)

     call gtri(terma,termc,b2,d2,qa1,1,qa(1),3.d0,1,qax,3.d0,kz,2)

     !--------------------------------------------------------------------
     !c
     !c      Estimation of Ri number
     !c
     !c--------------------------------------------------------------------

!****new part starts ****

     !do k=1,kz
     !   f2(k) = kt(k) * prod1(k)+vtran(k)
     !enddo


     !dtree = 0.
     !do k = 2, kz
     !   if (z(k).le. z(nz)) then
     !
     !      dtree = dtree + 0.5*(z(k) -z(k-1)) *(f2(k-1) +f2(k)) / f2(nz) / hc
     !   endif
     !enddo

     !ddd1 = max(0., (1-dtree)*hc)

!**** new part ends ****

    do k=2,kz-1    ! loop 9033

        dtz1=diff(ta1,kz,dz,k)+gamma
        dqz1=diff(qa1,kz,dz,k)
        duz1=diff(u1,kz,dz,k)
        dvz1=diff(v1,kz,dz,k)
        dwz1=diff(w1,kz,dz,k)


!        if(k .le. 2) then   ! old
!           duz1=0.4*u1(2)/log(z(2)/z0)/l(k)
!           dvz1=0.4*v1(2)/log(z(2)/z0)/l(k)
!           dtz1=0.4*(ta1(2)-ta1(1)+gamma*dz(2))/log(z(2)/z0)/l(k)
!           dqz1=0.4*(qa1(2)-qa1(1))/log(z(2)/z0)/l(k)
!        endif


        if(k .le. 2) then  ! new
           duz1=ur(2)/l(2)
           dvz1=0.
           dtz1=(((ta1(3)-ta1(1))*dz(2)**2   &
                +(ta1(1)-ta1(2))*(dz(2)+dz(3))**2))/   &
                ( (dz(2)+dz(3))*dz(2)**2 -              &
                dz(2)*(dz(2)+dz(3))**2  )+gamma
           dqz1=(((qa1(3)-qa1(1))*dz(2)**2  &
                +(qa1(1)-qa1(2))*(dz(2)+dz(3))**2))/  &
                ( (dz(2)+dz(3))*dz(2)**2 -   &
                dz(2)*(dz(2)+dz(3))**2  )
        endif

        !**** old ****
        !utt=(duz1**2+dvz1**2)
        !if (utt .eq. 0) utt = 1E-15
        !if(utt.ne.0.) rih(k)=(grav*((dtz1  )/ (ta1(k)+gamma*z(k)) + 1.*0.608*dqz1/roa ) ) / utt
        !**** old ****

!****new part starts ****

        !abb=min(1.1,l(k)/al1)
        abb=abs(alt(k)*rih(k))/(1.+abs(alt(k)*rih(k)))    ! (neutral case =0) <abb < 1, parameter related to convectivity, related to canopy (bounce or                                                          !  mechanic production)

        !abb1=(max(0.,1.-(kt(k)*prod1(k)+1.*vtran(k)) &
        !     /(kt(k)*prod2(k)+kt(k)*prod1(k)+1.*vtran(k)) )   )


    !if(z(k).lt.ddd1) then
    !       abb1=(max(0.,1.-(kt(k)*prod1(k)+0.*vtran(k)) &
    !            /(kt(k)*prod2(k)+kt(k)*prod1(k)+0.*vtran(k)) )   )
        !endif

        a005=2.     ! const value for neutral case
        a0005=0.

        abb3=abb
        rih1(k)=(a0005-(a0005+1./(1.-c52/c833))*abb3)
        !if(z(k).le.ddd1) then
        !   abb3=abb1
        !   a0005=2.
        !   rih1(k)=(a0005-(a0005+1./(1.-c52/c833))*abb3)
        !endif
    if(dtz1.ge.0.) then

           abb3=a005*alt(k)*rih(k)
           rih1(k)=a005

    endif

        utt=max((duz1**2+dvz1**2)  &
             !+bt2(k)/kt(k)  &
             !+2.*det2(k)/kt(k)   &
             +abs(1.*rih1(k)*( alt(k)*grav*((dtz1  )/   &
             (ta1(k)+gamma*z(k))   &
             + 0.608*dqz1/roa)))  &
             ,0.0000001)

!    if(utt.gt.0.) then
           rih(k)=(   grav*((dtz1  )/(ta1(k)+gamma*z(k) )  + 0.608*dqz1/roa ) ) / utt
!    else
!           rih(k)=0.
!    endif


!**** new part ends ****

        rih(1)=0.
        rih(kz)=0.

!        if(rih(k).gt.0.) then
!           alt1(k)=1.35/(1.+1.35*rih(k))

        if (rih(k).ge.0.0d0) then !new gt -> ge
          alt1(k)=1.35d0/(1.0d0+a005*1.35d0*rih(k))        !new 1.35 -> a005
        else
          alt1(k)=1.35d0*(1.0d0- 15.0d0*rih(k))**0.25d0
        endif

    enddo ! 9033

     !--------------------------------------------------------------------
     ! Solution of the turbulent kinetic energy & dissipation rate equation
     !--------------------------------------------------------------------
     do k=2,kz-1    ! loop 933
        !
        dtz1=diff(ta1,kz,dz,k)+gamma
        dqz1=diff(qa1,kz,dz,k)
        duz1=diff(u1,kz,dz,k)
        dvz1=diff(v1,kz,dz,k)
        dwz1=diff(w1,kz,dz,k)

!        if(k.le.2) then  !old
!           duz1=0.4*u1(2)/log(z(2)/z0)/l(k)
!           dvz1=0.4*v1(2)/log(z(2)/z0)/l(k)
!           dtz1=0.4*(ta1(2)-ta1(1)+gamma*dz(2))/log(z(2)/z0)/l(k)
!           dqz1=0.4*(qa1(2)-qa1(1))/log(z(2)/z0)/l(k)
!        endif

        if(k.le.2) then   !new
           duz1=ur(2)/l(2) ! new
           dvz1=0.        ! new

           dtz1=(((ta1(3)-ta1(1))*dz(2)**2   &
                +(ta1(1)-ta1(2))*(dz(2)+dz(3))**2))/   &
                ( (dz(2)+dz(3))*dz(2)**2 -              &
                dz(2)*(dz(2)+dz(3))**2  )+gamma

           dqz1=(((qa1(3)-qa1(1))*dz(2)**2  &
                +(qa1(1)-qa1(2))*(dz(2)+dz(3))**2))/  &
                ( (dz(2)+dz(3))*dz(2)**2 -   &
                dz(2)*(dz(2)+dz(3))**2  )
        endif

        !abb=min(1.1,l(k)/al1)
        prod1(k) = duz1**2. + dvz1**2.
        prod2(k) = -1.*(alt(k)*grav*(dtz1/(ta1(k) +gamma*z(k)) + 0.608*dqz1/roa))

        abb=abs(alt1(k)*rih(k))/(1.+abs(alt1(k)*rih(k)))    ! (neutral case =0) <abb < 1, parameter related to convectivity, related to canopy (bounce or                                                          !  mechanic production)


        !abb1=(max(0.,1.-(kt(k)*prod1(k)+1.*vtran(k)) &
        !     /(kt(k)*prod2(k)+kt(k)*prod1(k)+1.*vtran(k)) )   )

    !if(z(k).lt.ddd1) then
    !       abb1=(max(0.,1.-(kt(k)*prod1(k)+0.*vtran(k)) &
    !            /(kt(k)*prod2(k)+kt(k)*prod1(k)+0.*vtran(k)) )   )
!
 !       endif
        abb3 = abb
        cc22=c52+(c833-c52)*abb3
        cc33=c833
        !a0005 = 0.
        rih1(k)=(a0005-(a0005+1./(1.-c52/c833))*abb3)


        !if(z(k) .le. ddd1) then
        !   abb3 = abb1
        !   a0005 = 2.
        !   rih1(k)=(a0005-(a0005+1./(1.-c52/c833))*abb3)
        !   cc22=c52+(c833-c52)*abb3
        !endif
        beta=(c52-c833)*rih1(k)


        if(rih(k).ge.0.) then
           abb3 = a005*alt(k)*rih(k)
           rih1(k) = a005
           cc22=c52+(c833-c52)*abb3
           cc33 = c833
           beta=(c52-c833)*rih1(k)
        endif

        c15=(12.*cc2**0.5)*cd*s1(k)*( sqrt(u1(k)**2+v1(k)**2)) ! if c16 in use, 10*cc2 is used for c15, if c16 not in use, then 12*c22 is used.

        bt2(k) = (cd*s1(k) * (sqrt(u1(k)**2 + v1(k)**2))**3)*(0. + abb3**0.5)

!**** new part ends ****

        fktt=2.*alf*(df(k)*kt(k+1)+(1.-df(k))*kt(k))
        fktd=2.*alf*((1.-df(k-1))*kt(k-1)+df(k-1)*kt(k))

        if(k.le.2) then    ! new
           fktt=2.*alf*sqrt(kt(k+1)*kt(k)) !new
           ! geometrical mean provides the smallest weight possible, according to Andrey most stable
           fktd=2.*alf*sqrt(kt(k-1)*kt(k))  ! new
        endif    ! new


        !cc22=(c52+(c833-c52)*l(k)/al1)

        !c15=(12.*cc2**0.5)*s(k)*( sqrt(u1(k)**2+v1(k)**2)+sqrt(u1(k)**2+v1(k)**2))/2.
        !cc33=(c833-(c833-c52)*(c15)/dbt(k))

        !beta0=(cc22-c833)


        ! Andrey new addings on 2013.09.17

        vtran(k) = bt(k+1) * fktt/da(k)    - bt(k)*(dz(k)*fktt + dz(k+1)*fktd)/(dz(k+1) + dz(k)) /db(k) + bt(k-1)*fktd/dc(k)

        det2(k) = c15 * bt(k)

        terma(k)= dt_mete*(fktt -w(k)*dz(k))/da(k)
        termb(k)= dt_mete*((dz(k)*fktt+ dz(k+1)*fktd) /(dz(k+1)+dz(k)) +w(k)*(dz(k+1)-dz(k)))/db(k)+ 1.0d0+dt_mete*dbt(k)
        termc(k)= dt_mete*(fktd+ w(k)*dz(k+1))/dc(k)
        !old: termd(k)= bt(k) + dt_mete*kt(k) *(duz1**2+dvz1**2+dwz1**2-alt1(k)*grav*(dtz1/(ta1(k)+gamma*z(k)) + 1.*0.608*dqz1/roa))-0.*(bt(k)-bnud(k))*c_nud !nudging!
        termd(k)= bt(k) + dt_mete*kt(k)*(prod1(k)+prod2(k))

        !old: bg(k)= dt_mete*((dz(k)*fktt+ dz(k+1)*fktd) /(dz(k+1)+dz(k)) +w(k)*(dz(k+1)-dz(k)))/db(k)+ 2.*dt_mete*cc33*dbt(k)+1.0d0
        bg(k)= dt_mete*((dz(k)*fktt+ dz(k+1)*fktd) /(dz(k+1)+dz(k)) +w(k)*(dz(k+1)-dz(k)))/db(k)+ 2.*dt_mete*cc33*dbt(k)+1.0d0 &
             -1.*dt_mete*(c52-cc33)*bt2(k)/bt(k) + 1.*dt_mete*(c52-cc33)*c15

        !old : dg(k)= dbt(k) +dt_mete*cc33*dbt(k)**2 +dt_mete*cc22*cc2 *(duz1**2+dvz1**2+dwz1**2) &
        !     -beta0*dt_mete*cc2*(1.*alt1(k)*grav*(dtz1/(ta1(k)+gamma*z(k))+ 1.*0.608*dqz1/roa ) )
       dg(k)= dbt(k) + dt_mete*cc33*dbt(k)**2 + dt_mete*cc2*(cc22*prod1(k) +beta*prod2(k))

    enddo  ! end of loop 933

     if(abl == 1) then
        bttop=temgrad**2/(cc2**0.5)
        !old CALL gtri(terma,termc,termb,termd,bt1,2,0.d0,3.d0,1,bttop,3.d0,kz,2)
        CALL gtri(terma,termc,termb,termd,bt1,2,0.d0,3.d0,2,0.d0,3.d0,kz,2)
    else

        btbot=ur(1)**2/(cc2**0.5) ! new
        !old: CALL gtri(terma,termc,termb,termd,bt1,2,0.d0,3.d0,2,0.d0,3.d0,kz,2)
        CALL gtri(terma,termc,termb,termd,bt1,1,btbot,3.d0,2,0.d0,3.d0,kz,2)  ! new

     endif

!tke_limit = 1.d-5  ! 1.d-2 in manitou
     do k=1,kz !-1
        !old bt1(k)=max(bt1(k),1.d-5)
        bt1(k)=max(bt1(k),1.d-2) !new
     enddo


     !old CALL gtri(terma, termc, bg, dg, dbt1, 1, (cc2**(0.75)*bt1(2)**0.5/l(1)), 3.d0, 1, 1.*(cc2**(0.75)*bt1(kz)**0.5/l(kz)), 3.d0, kz, 2)
     CALL gtri(terma, termc, bg, dg, dbt1, 1, (cc2**(0.75)*btbot**0.5/l(1)), 3.d0, 2, 0.d0, 3.d0, kz, 2)  !new

     ! limiting disipational rate related to eddy
    do k=1,kz
        dbt1(k) = max(dbt1(k), cc2**0.75*(bt1(k))**0.5/(al1))
     enddo

     !-------------------------------------------------------------
     ! solution for ABL height
     !-------------------------------------------------------------

    call find_boundary_layer(z, rih, nz + 1, pblh, pblh_ind)


!old      if(abl /= 1) al1=max(1.*0.00027*abs(windx)/abs(f1),5.d0)  ! was used for netural cases

     !----------------------------------------------------------
     !  relation for the eddy exchange coefficient of momentum
     !----------------------------------------------------------

!old     kt1(1)=max(cc2*bt1(1)/dbt1(1),1.d-3)
     kt1(1)=max(cc2*bt1(1)/dbt1(1),ur(1)*l(1)) !new
     kt1(kz)=max(cc2*bt1(kz)/dbt1(kz),kt1(1))
 !dublicate    l1(1)=l(1)  !new
     l1(kz)=max(sqrt(bt1(kz))*cc2**0.75/dbt1(kz),l1(1))

     do k=2,kz-1   !  loop 291
        gf=abs(bt1(k))
        gf2=dbt1(k)
        kt1(k)=max(cc2*gf/gf2,kt1(1))
        l1(k)=max(sqrt(gf)*cc2**0.75/gf2,l1(1))
     enddo  ! end of loop 291

     if(abl == 1) then
        do k=2,kz
           if(z(k).gt.hc) l1(k)=l1(k-1)+0.43*dz(k)
        enddo
     endif

     do  k=2,kz-1
        !old   fktt=((df(k)*kt1(k+1) +(1.-df(k))*kt1(k)))
        !old   fktd=(((1.-df(k-1))*kt1(k-1) +df(k-1)*kt1(k)))
        !old   ur(k)=0.5*(sqrt(fktt*abs(sqrt(u1(k+1)**2+v1(k+1)**2) - sqrt(u1(k)**2+v1(k)**2))/dz(k+1))+  &
        !old        sqrt(fktd*abs(sqrt(u1(k)**2+v1(k)**2) - sqrt(u1(k-1)**2+v1(k-1)**2))/dz(k)))



        pa1=(kt1(k+1)+kt1(k))*dz(k)/da(k)/2.  !new
        pa0=(kt1(k)+kt1(k-1))*dz(k+1)/dc(k)/2. !new

        akt_duz=(pa1*(u1(k+1)-u1(k))+pa0*(u1(k)-u1(k-1)))  !new

        akt_dvz=(pa1*(v1(k+1)-v1(k))+pa0*(v1(k)-v1(k-1)))  !new

        ur(k)=(akt_duz**2+akt_dvz**2)**0.25   !new
        ! same thing, but speeds up calculations
     enddo

     !old ur(1)=ur(2)
     ur(1)=(ur(2)*(dz(2)+dz(3))**2-ur(3)*dz(2)**2)/( (dz(2)+dz(3))**2-dz(2)**2) !new cubic parameterization for first level (so that there is no gradient)
     ur(kz)=ur(kz-1)

!**** new part starts ****

     sb=0.
     szb=0.


     do k=2,kz
        szb=(sqrt(abs(bt1(k-1)-0.01))*z(k-1)+ &  ! 0.01 is the tke limit
        sqrt(abs(bt1(k)-0.01))*z(k))*   & ! 0.01 is the tke limit
        dz(k)/2.+szb

        sb=(sqrt(abs(bt1(k-1)-0.01))+sqrt(abs(bt1(k)-0.01)))*  & ! 0.01 is the tke_limit
        dz(k)/2.+sb
     enddo
     !      al0=al1
     if (sb .eq. 0) then
        write(*,*) 'error: divided by zero. sb =', sb
     endif
          al100=0.00027*abs(windx)/abs(f1)  ! blackadar value for netural case
     al1=max(al100,0.075*(szb/sb) )  ! convective case

! szb/sb the relation of the integrals above

!
!**** new part ends ****

     nturb=nturb+1

     !--------------------------------------------------------------
     !   solution of the soil transfer
     !---------------------------------------------------------------

     sand=0.000001*0.50
     sand=temtran  ! above defined temtran=0.000001*0.5, so in this case sand = temtran   ! thermal diffusivity
     if(rou.le.0.0) sand=temtran   ! rou = 0 snow cover excluded for now
     snow=0.000001*0.27
     do  k=1,kz-2
        if(rou.ge.z(k)) kmix=max(1,k+1)  ! rou = 0 snow cover excluded for now, -> kmix = 1
     enddo
     nmix=int(z(kmix-1)/0.05)
     do k=2,kz-1
        ksnow(k)=sand                     ! thermal diffusivity of soil
        if(k.gt.10) ksnow(k)=snow         ! thermal diffusivity
        if(k.gt.(10+nmix-1)) ksnow(k)=0.  ! thermal diffusivity
     enddo
     DO k=2,10-1+nmix
        dksz=0.05*ksnow(k+1)/0.005- 0.05*ksnow(k-1)/0.005  ! depth of each layer is 0.05
        terma(k)= dt_mete*(2.*ksnow(k) +0.05*dksz)/0.005
        termb(k)= dt_mete*2.*ksnow(k) /0.0025+1.
        termc(k)= dt_mete*(2.*ksnow(k) -0.05*dksz)/0.005
        termd(k)=tsoil(k)
     ENDDO
     e(1)=0.   ! constant temperature

     ! taking deep soil temperature from measurements (9  Temperature (K) in 23-60 cm, now it is average of 22-29 cm and 42-58 cm)
     ts0 = linear_interp(time_in_month, nxodrad, local_deep_soiltemp, dt_obs)  ! [K]

     f(1)=ts0  ! condition
     !f(1)=T00 ! condition

     m2=10-1+nmix

     DO m=2,m2
        den=termb(m)-termc(m)*e(m-1)
        f(m)=(termd(m)+termc(m)*f(m-1))/den  ! [K]
        e(m)=terma(m)/den                ! [1]
     ENDDO



     tsoil1(10+nmix)=ta1(kmix-1)
     DO n=1,m2
        m=10-n+nmix
        tsoil1(m)=e(m)*tsoil1(m+1)+f(m)
     ENDDO

     do k=1,kz
        if(k.gt.10+nmix-1) tsoil1(k)=ta1(kmix-1)
     enddo
     !----------------------------------------------------------
     ! radiation blok
     !----------------------------------------------------------

     do k=1,kz-2
        if(rou.ge.z(k)) kmix=max(1,k+1)  ! dublicate, already done above   ! rou = 0 snow cover excluded for now, -> kmix = 1
     enddo

     DO k=1,kz
        sl(k)=s1(k)
     ENDDO

     DO k=2,kz
        gl(k)=1.
        gd(k)=1.
     ENDDO

     DO k=2,kz
        if((sl(k)-sl(k-1)).gt.0.) kf=k
        go to 4015
     ENDDO
4015 ks=max(kf,1)
     lai(1)=0.
     DO k=nz-1,1,-1
        jk=nz-k+1
        lai(jk)=lai(jk-1)+(sl(jk-1)+sl(jk))*dz(jk)/2.
     ENDDO

     lai0=lai(nz)
     do k=nz,kz
        lai(k)=lai0
     enddo
     asl=1.5

     call IntegralFLeafOrientation(kz, nz, rads2, cos_zenith, su, gl)

     DO k=nz,1,-1
!        psn(k)=exp(0.5*(lai(k)-lai0)/cos_zenith)
!        gl(k)=0.5
!        if(rads2.le.10.) gl(k)=cos_zenith
 !       if(rads2.le.5.) gl(k)=0.
 !       if(rads2.le.10.) psn(k)=exp((lai(k)-lai0))
 !       if(su.eq.2) then
 !          psn(k)=exp((lai(k)-lai0))
 !          gl(k)=cos_zenith
 !          if(rads2.le.5.) gl(k)=0.
 !       endif

        lai_cumulative(k) = -(lai(k)-lai0)  ! added by rosa 23.8.2014
     ENDDO

     call PenetrationFDirect(nz, kz, su, rads2, cos_zenith, gl, lai_cumulative, psn)

     ! in theory should be if( su.eq.2 .or. rads2.le.10.), but that produces significanlty different results
     if ( su.eq.2 )  then
         psk=psn
     else
         call PenetrationFDiffuse(kz, nz, gl, lai_cumulative, psk)
     endif

     DO k=nz,1,-1
        if((lai(k)-lai0).ne.0.and.rads2.gt.5.) gd(k)=cos_zenith*log(psk(k))/ (lai(k)-lai0)
        if(su.eq.2) gd(k)=cos_zenith
        if(rads2.le.5.) gd(k)=0.
     ENDDO
     DO k=2,kz-1
        asun(k)=dz(k)*psn(k+1)/da(k)+ (dz(k+1)-dz(k))*psn(k)/db(k)- dz(k+1)*psn(k-1)/dc(k)
        asky(k)=dz(k)*psk(k+1)/da(k)+ (dz(k+1)-dz(k))*psk(k)/db(k)- dz(k+1)*psk(k-1)/dc(k)
     ENDDO


     asun(1)=(psn(2)-psn(1))/dz(2)
     asky(1)=(psk(2)-psk(1))/dz(2)

     alph=0.80_dp
     alni=0.80_dp
     if(tau.gt.0.) then
        alph=alph-0.03*(dt_mete/SECONDS_IN_ONE_DAY)
        alni=alni-0.03*(dt_mete/SECONDS_IN_ONE_DAY)
     endif
     if(rou.gt.0.) then
        alph=max(1.d-1,alph)
        alni=max(1.5d-1,alni)
     endif
     if(rou.le.0.) then
        alph=0.10_dp
        alni=0.15_dp
     endif

     if(ta1(1).ge.T00) then
        alph=0.10_dp
        alni=0.15_dp
     else
        alph=0.50_dp
        alni=0.50_dp
     endif

     !alph = 0.125 !CLEARCUT: include this
     !alni = 0.225 !CLEARCUT: include this

     wee=1.0_dp
     tphu=0.06_dp*(ONE-sqrt( (ONE-0.06_dp)**2-0.09_dp**2))/0.15_dp
     tphd=0.06_dp*(ONE-sqrt( (ONE-0.06_dp)**2-0.09_dp**2))/0.15_dp
     rphu=0.09_dp*(ONE-sqrt( (ONE-0.06_dp)**2-0.09_dp**2))/0.15_dp
     rphd=0.09_dp*(ONE-sqrt( (ONE-0.06_dp)**2-0.09_dp**2))/0.15_dp

     tniu=0.40_dp*(ONE-sqrt( (ONE-0.40_dp)**2-0.35_dp**2))/0.75_dp
     tnid=0.40_dp*(ONE-sqrt( (ONE-0.40_dp)**2-0.35_dp**2))/0.75_dp
     rniu=0.35_dp*(ONE-sqrt( (ONE-0.40_dp)**2-0.35_dp**2))/0.75_dp
     rnid=0.35_dp*(ONE-sqrt( (ONE-0.40_dp)**2-0.35_dp**2))/0.75_dp

     if(su.eq.2) then
        tniu=0.50_dp*(ONE-sqrt( (ONE-0.50_dp)**2-0.35_dp**2))/0.85_dp ! ross
        tnid=0.50_dp*(ONE-sqrt( (ONE-0.50_dp)**2-0.35_dp**2))/0.85_dp ! ross
        rniu=0.35_dp*(ONE-sqrt( (ONE-0.50_dp)**2-0.35_dp**2))/0.85_dp ! ross
        rnid=0.35_dp*(ONE-sqrt( (ONE-0.50_dp)**2-0.35_dp**2))/0.85_dp ! ross
     endif
     phsn=0.4200_dp  ! default: 0.42
     phsk=0.6000_dp  ! 0.6, 0.51 is better in results
     two=1.5
     tu=tphu
     td=tphd
     ru=rphu
     rd=rphd
     al=alph
     pasn=phsn
     pask=phsk
4026 continue
     DO k=1,nz
        fu0(k)=0.001_dp
        fd0(k)=0.001_dp
        fu1(k)=0.001_dp
        fd1(k)=0.001_dp
     ENDDO
4028 continue

     !   ------------------------------------------ haotic --------------

     fu1(1)=al*(fd0(1)+psn(1)*rsnt*pasn +psk(1)*rskt*pask)


     DO k=1,nz-1

        trem1=dz(k+1)*((sl(k+1)+sl(k))/2.0_dp)* (fd0(k+1)*ru-(wee-td)*fu1(k)+  &
             ru*(((psn(k)+psn(k+1))/2.0_dp)*rsnt*pasn*gl(k)/cos_zenith +((psk(k)+psk(k+1))/2.0_dp)*rskt*pask*gd(k)/cos_zenith))

        trem2=dz(k+1)*((sl(k+1)+sl(k))/2.0_dp)* (0.5_dp*(fd0(k+1)+fd0(k))*ru-(wee-td)*(fu1(k)+trem1/2.0_dp)+  &
             ru*(((psn(k)+psn(k+1))/2.0_dp)*rsnt*pasn*gl(k)/cos_zenith +((psk(k)+psk(k+1))/2.0_dp)*rskt*pask*gd(k)/cos_zenith))

        trem3=dz(k+1)*((sl(k+1)+sl(k))/2.0_dp)* (0.5_dp*(fd0(k+1)+fd0(k))*ru-(wee-td)*(fu1(k)+trem2/2.0_dp)+  &
             ru*(((psn(k)+psn(k+1))/2.0_dp)*rsnt*pasn*gl(k)/cos_zenith +((psk(k)+psk(k+1))/2.0_dp)*rskt*pask*gd(k)/cos_zenith))

        trem4=dz(k+1)*((sl(k+1)+sl(k))/2.0_dp)* (fd0(k+1)*ru-(wee-td)*(fu1(k)+trem3)+  &
             ru*(((psn(k)+psn(k+1))/2.0_dp)*rsnt*pasn*gl(k)/cos_zenith +((psk(k)+psk(k+1))/2.0_dp)*rskt*pask*gd(k)/cos_zenith))

        fu1(k+1)=fu1(k)+(trem1+2.0_dp*trem2+2.0_dp*trem3+trem4)/6.0_dp
     ENDDO

     DO k=nz,2,-1
        trek1=dz(k)*((sl(k)+sl(k-1))/2.0_dp)* (fu0(k-1)*rd-(wee-tu)*fd1(k)+  &
             td*(((psn(k-1)+psn(k))/2.0_dp)*rsnt*pasn *gl(k)/cos_zenith +((psk(k-1)+psk(k))/2.0_dp)*rskt*pask*gd(k)/cos_zenith))

        trek2=dz(k)*((sl(k)+sl(k-1))/2.0_dp)* (0.5_dp*(fu0(k-1)+fu0(k))*rd-(wee-tu)*(fd1(k)+trek1/2.0_dp)+ &
             td*(((psn(k-1)+psn(k))/2.0_dp)*rsnt*pasn*gl(k)/cos_zenith +((psk(k-1)+psk(k))/2.0_dp)*rskt*pask*gd(k)/cos_zenith))

        trek3=dz(k)*((sl(k)+sl(k-1))/2.0_dp)* (0.5_dp*(fu0(k-1)+fu0(k))*rd-(wee-tu)*(fd1(k)+trek2/2.0_dp)+  &
             td*(((psn(k-1)+psn(k))/2.0_dp)*rsnt*pasn*gl(k)/cos_zenith +((psk(k-1)+psk(k))/2.0_dp)*rskt*pask*gd(k)/cos_zenith))

        trek4=dz(k)*((sl(k)+sl(k-1))/2.0_dp)* (fu0(k-1)*rd-(wee-tu)*(fd1(k)+trek3)+  &
             td*(((psn(k-1)+psn(k))/2.0_dp)*rsnt*pasn*gl(k)/cos_zenith +((psk(k-1)+psk(k))/2.0_dp)*rskt*pask*gd(k)/cos_zenith))

        fd1(k-1)=fd1(k)+(trek1+2.0_dp*trek2+2.0_dp*trek3+trek4)/6.0_dp

     ENDDO

     DO k=1,nz
        IF(ABS(fu0(k)/fu1(k)-1.0_dp).GT.eps) go to 4033
     ENDDO
     DO k=1,nz-1
        IF(ABS(fd0(k)/fd1(k)-1.0_dp).GT.eps) go to 4033
     ENDDO
     !-----------------------------------------haotic------------------------
     go to 4035
4033 continue
     DO k=1,nz
        fu0(k)=fu1(k)
        fd0(k)=fd1(k)
     ENDDO
     go to 4028
4035 if(two.lt.2.) go to 4036
     go to 4038
4036 continue
     DO k=1,nz
        fphu(k)=fu1(k)
        fphd(k)=fd1(k)
     ENDDO
     tu=tniu
     td=tnid
     ru=rniu
     rd=rnid
     al=alni
     pasn=1.0_dp-phsn
     pask=1.0_dp-phsk
     two=two+1.0_dp
     go to 4026

4038 continue
     DO k=1,nz
        fniu(k)=fu1(k)
        fnid(k)=fd1(k)
     ENDDO

     iru0=(dels*sigma_SB*ta1(kmix-1)**4+1.*(1.-dels)*fird(1))/pi ! net LWR at surface
     iruh=sigma_SB*(ta1(nz  )**4)*emm/pi

     !r using reanalysis value instead
     iruh = LWRdown/pi   !CLEARCUT - comment this part out

     firu(1)=iru0*pi
     fird(nz)=iruh*pi
     DO k=1,nz
        fu0(k)=0.001
        fd0(k)=0.001
        fu1(k)=0.001
        fd1(k)=0.001
     ENDDO

40278 continue

     fu1(1)=iru0*pi
     fd1(nz)=iruh*pi

     DO k=1,nz-1

        trom1=dz(k+1)*((1.*sl(k)+1.*sl(k+1))/2.)*delf*(sigma_SB*(psn(k)*(tsn(k)**4)+ (1.-psn(k))*(tsd(k)**4))-fu1(k))

        trom2=dz(k+1)*((1.*sl(k)+1.*sl(k+1))/2.)*delf*(sigma_SB*(0.5*(psn(k+1)*(tsn(k+1)**4) +psn(k)*(tsn(k)**4))+  &
             0.5*((1.-psn(k+1))*(tsd(k+1)**4) +(1.-psn(k))*(tsd(k)**4)))-(fu1(k)+trom1/2.))

        trom3=dz(k+1)*((1.*sl(k)+1.*sl(k+1))/2.)*delf*(sigma_SB*(0.5*(psn(k+1)*(tsn(k+1)**4) +psn(k)*(tsn(k)**4))+  &
             0.5*((1.-psn(k+1))*(tsd(k+1)**4) +(1.-psn(k))*(tsd(k)**4)))-(fu1(k)+trom2/2.)  )

        trom4=dz(k+1)*((1.*sl(k)+1.*sl(k+1))/2.)*delf*(sigma_SB*(psn(k+1)*(tsn(k+1)**4)+ (1.-psn(k+1))*(tsd(k+1)**4))-(fu1(k)+trom3) )
        !
        fu1(k+1)=fu1(k)+(trom1+2.*trom2+2.*trom3+trom4)/6.

     ENDDO
     !
     DO k=nz,2,-1

        tram1=dz(k)*((1.*sl(k)+1.*sl(k-1))/2.)*delf*(sigma_SB*(psn(k)*(tsn(k)**4)+(1.-psn(k))*(tsd(k)**4))-fd1(k))

        tram2=dz(k)*((1.*sl(k)+1.*sl(k-1))/2.)*delf*(sigma_SB*(0.5*(psn(k-1)*(tsn(k-1)**4) +psn(k)*(tsn(k)**4))+  &
             0.5*((1.-psn(k-1))*(tsd(k-1)**4) +(1.-psn(k))*(tsd(k)**4)))-(fd1(k)+tram1/2.)  )

        tram3=dz(k)*((1.*sl(k)+1.*sl(k-1))/2.)*delf*(sigma_SB*(0.5*(psn(k-1)*(tsn(k-1)**4) +psn(k)*(tsn(k)**4))+  &
             0.5*((1.-psn(k-1))*(tsd(k-1)**4) +(1.-psn(k))*(tsd(k)**4)))-(fd1(k)+tram2/2.)  )

        tram4=dz(k)*((1.*sl(k)+1.*sl(k-1))/2.)*delf*(sigma_SB*(psn(k-1)*(tsn(k-1)**4)+ (1.-psn(k-1))*(tsd(k-1)**4))-(fd1(k)+tram3) )

        fd1(k-1)=fd1(k)+(tram1+2.*tram2+2.*tram3+tram4)/6.

     ENDDO

     DO k=1,nz
        IF(ABS(fu0(k)/fu1(k)-1.).GT.eps) THEN
          go to 40373
        END IF

     ENDDO

     DO k=1,nz-1
        IF(ABS(fd0(k)/fd1(k)-1.).GT.eps) go to 40373
     ENDDO

     !-----------------------------------------haotic------------------------
     go to 40376

40373 continue

     DO k=1,nz
        fu0(k)=fu1(k)
        fd0(k)=fd1(k)
     ENDDO
     go to 40278
40376 continue

     DO k=1,nz
        firu(k)=fu1(k)
        fird(k)=fd1(k)
     ENDDO

     fird(nz)=iruh*pi
     firu(1)=iru0*pi

     nnna=1

     nnna=nnna+1
     pgr=0.6305_dp
     pre=0.6056_dp

     nua=0.0000141_dp
     DO k=1,nz
        rlf=0.03_dp
        if(su.eq.2) then
           rlf=0.05_dp
        endif
        grsn(k)=grav*(pgr*rlf)**3*abs(tsn(k)-ta1(k))/(T00*nua**2)
        grsd(k)=grav*(pgr*rlf)**3*abs(tsd(k)-ta1(k))/(T00*nua**2)
        re(k)=sqrt(u1(k)**2+v1(k)**2)*pre*rlf/nua
     ENDDO

     DO k=1,nz
        if(re(k).ge.0..and.re(k).le.h(1)) nure(k)=h(2)*re(k)**h(3)
        if(re(k).gt.h(1).and.re(k).le.h(4)) nure(k)=h(5)*re(k)**h(6)
        if(re(k).gt.h(4)) nure(k)=h(7)*re(k)**h(8)
     ENDDO
     DO k=1,nz
        if(grsn(k).ge.0..and.grsn(k).le.h(10)) nusn(k)=h(3)
        if(grsn(k).gt.h(10).and.grsn(k).le.h(11)) nusn(k)=h(12)*grsn(k)**h(13)
        if(grsn(k).gt.h(11).and.grsn(k).le.h(14)) nusn(k)=h(15)*grsn(k)**h(16)
        if(grsn(k).gt.h(14)) nusn(k)=h(17)*grsn(k)**h(18)
     ENDDO
     DO k=1,nz
        if(grsd(k).ge.0..and.grsd(k).le.h(10)) nusd(k)=h(3)
        if(grsd(k).gt.h(10).and.grsd(k).le.h(11)) nusd(k)=h(12)*grsd(k)**h(13)
        if(grsd(k).gt.h(11).and.grsd(k).le.h(14)) nusd(k)=h(15)*grsd(k)**h(16)
        if(grsd(k).gt.h(14)) nusd(k)=h(17)*grsd(k)**h(18)
     ENDDO

     DO k=1,nz
        dhsn(k)=1.*kta*nusn(k)/(pgr*rlf)

        if((re(k)*re(k)).ge.grsn(k)) dhsn(k)=1.*kta*nure(k)/(pre*rlf)
     ENDDO

     DO k=1,nz
        dhsd(k)=1.*kta*nusd(k)/(pgr*rlf)
        if((re(k)*re(k)).ge.grsd(k)) dhsd(k)=1.*kta*nure(k)/(pre*rlf)
     ENDDO

     DO k=1,nz
        if(re(k).ge.0..and.re(k).le.h(1)) shre(k)=h(19)*re(k)**h(3)
        if(re(k).gt.h(1).and.re(k).le.h(4)) shre(k)=h(20)*re(k)**h(6)
        if(re(k).gt.h(4)) shre(k)=h(21)*re(k)**h(8)
     ENDDO

     DO k=1,nz
        if(grsn(k).ge.0..and.grsn(k).le.h(23)) shsn(k)=h(3)
        if(grsn(k).gt.h(23).and.grsn(k).le.h(24)) shsn(k)=h(25)*grsn(k)**h(13)
        if(grsn(k).gt.h(24).and.grsn(k).le.h(26)) shsn(k)=h(27)*grsn(k)**h(16)
        if(grsn(k).gt.h(26)) shsn(k)=h(28)*grsn(k)**h(18)
     ENDDO

     DO k=1,nz
        if(grsd(k).ge.0..and.grsd(k).le.h(23)) shsd(k)=h(3)
        if(grsd(k).gt.h(23).and.grsd(k).le.h(24)) shsd(k)=h(25)*grsd(k)**h(13)
        if(grsd(k).gt.h(24).and.grsd(k).le.h(26)) shsd(k)=h(27)*grsd(k)**h(16)
        if(grsd(k).gt.h(26)) shsd(k)=h(28)*grsd(k)**h(18)
     ENDDO


     DO k=1,nz
        shsn(k)=kva*shsn(k)/(pgr*rlf)
        if((re(k)*re(k)).ge.grsn(k)) then
          shsn(k)=kva*shre(k)/(pre*rlf)
        end if
     ENDDO


     DO k=1,nz
        shsd(k)=kva*shsd(k)/(pgr*rlf)
        if((re(k)*re(k)).ge.grsd(k)) then
          shsd(k)=kva*shre(k)/(pre*rlf)
        end if
     ENDDO

     fniu(nz)=fniu(nz-1)
     fphu(nz)=fphu(nz-1)

     DO k=2,nz-1
        asun(k)=dz(k)*psn(k+1)/da(k)+ (dz(k+1)-dz(k))*psn(k)/db(k)- dz(k+1)*psn(k-1)/dc(k)
        asky(k)=dz(k)*psk(k+1)/da(k)+ (dz(k+1)-dz(k))*psk(k)/db(k)- dz(k+1)*psk(k-1)/dc(k)
        aphu(k)=dz(k)*fphu(k+1)/da(k)+ (dz(k+1)-dz(k))*fphu(k)/db(k)- dz(k+1)*fphu(k-1)/dc(k)
        aphd(k)=dz(k)*fphd(k+1)/da(k)+ (dz(k+1)-dz(k))*fphd(k)/db(k)- dz(k+1)*fphd(k-1)/dc(k)
        aniu(k)=dz(k)*fniu(k+1)/da(k)+ (dz(k+1)-dz(k))*fniu(k)/db(k)- dz(k+1)*fniu(k-1)/dc(k)
        anid(k)=dz(k)*fnid(k+1)/da(k)+ (dz(k+1)-dz(k))*fnid(k)/db(k)- dz(k+1)*fnid(k-1)/dc(k)
        airu(k)=dz(k)*firu(k+1)/da(k)+ (dz(k+1)-dz(k))*firu(k)/db(k)- dz(k+1)*firu(k-1)/dc(k)
        aird(k)=dz(k)*fird(k+1)/da(k)+ (dz(k+1)-dz(k))*fird(k)/db(k)- dz(k+1)*fird(k-1)/dc(k)
     ENDDO

     asun(1)=(psn(2)-psn(1))/dz(2)
     asky(1)=(psk(2)-psk(1))/dz(2)
     aphu(1)=(fphu(2)-fphu(1))/dz(2)
     aphd(1)=(fphd(2)-fphd(1))/dz(2)
     aniu(1)=(fniu(2)-fniu(1))/dz(2)
     anid(1)=(fnid(2)-fnid(1))/dz(2)
     airu(1)=(firu(2)-firu(1))/dz(2)
     aird(1)=(fird(2)-fird(1))/dz(2)

     face=2.7
     if(su.eq.1) then
        dmax=0.0026
        dmin=0.0000723
        topt=17.
        tmin=0.
        tmax=50.
        par50=20.
        humg=0.05
     endif

     if(su.eq.2) then
        dmax=0.0056
        dmin=0.0000723
        topt=18.
        tmin=0.
        tmax=50.
        par50=20.
        humg=0.05
        face=2.
     endif


     if(hc.eq.0.) then
        dmax=0.0056
        dmin=0.0000723
        topt=18.
        tmin=0.
        tmax=70.
        par50=50.
        humg=0.05
     endif



     DO k=1,nz-1  ! loop 4063
        if(sl(k).eq.0.0) go to 4072
        rosn(k)=(wee-tphu-rphu)*( (phsn*rsnt*gl(k)/cos_zenith) +(phsk*rskt*psk(k)*gd(k)/cos_zenith+fphu(k)+fphd(k)))/face
        go to 4073
4072    rosn(k)=0.
4073    if(sl(k).eq.0.) go to 4074

        rosd(k)=rosn(k)

        rosd(k)=(wee-tphu-rphu)* (phsk*rskt*psk(k)*gd(k)/cos_zenith+fphu(k)+fphd(k))/face

        go to 4075
4074    rosd(k)=0.
4075    continue
        if(rosn(k).le.0.) rosn(k)=0.
        if(rosd(k).le.0.) rosd(k)=0.
        tc=(tmax-topt)/(topt-tmin)

        rhg= qa1(1) / rhovs_buck(ta1(1))

 ! by Sampo, to avoid division by zero
        if( (wgf-wgwilt) == 0 ) then
           defsoil = 1
        else
           defsoil=min(1.d0, ((1.0*wg(2)+0.*wg(1))-wgwilt)/(wgf-wgwilt))
        end if

        if(wg(2).le.wgwilt) defsoil=0.

        if(rou.le.0.) goto 99

        if(ta1(k).lt.T00.or.tsd(k).le.T00) then
           defsoil=0.
        else
           defsoil = min(4.d-1,(ta1(k)-T00)**2.5/T00)

        endif
99      continue
        temlim=max(0.0_dp,(1.0_dp-0.0016_dp*(T00+18.0_dp-ta1(k))**2))
        dvsn(k)= dmin+  &
             dmax*(1.0_dp-exp(-0.01_dp*rosn(k)))  &
             * temlim
        dvsd(k)= dmin +  &
             dmax*(1.0_dp-exp(-0.03_dp*rosd(k)))  &
             * temlim
     ENDDO  ! end of loop 4063
     gstm_h2o_sn = dvsn
     gstm_h2o_sd = dvsd

     ! Integral coefficients of water vapour exchange between the surface of phytoelements and air between leaves:
     do k=2,nz-1
        dvsn(k)=(dvsn(k)*shsn(k)/(dvsn(k)+shsn(k)))
        dvsd(k)=(dvsd(k)*shsd(k)/(dvsd(k)+shsd(k)))
     enddo
     ! **************************************************
     DO k=1,nz  ! loop 4076
        rasn(k)=0.
        rasd(k)=0.
        goto 5252
        psn(k)=max(psn(k),1.d-2)
        !
        !   ************************+++  method 2 +++*********************


    if(sl(k).ne.0.0)   &
             rasn(k)=        &
             (rsnt*asun(k)   &

             /psn(k)         &

             + (rskt*asky(k)-1.*aphu(k)+1.*aphd(k)-   &
             1.*aniu(k)+1.*anid(k)-1.*airu(k)+1.*aird(k)))  &
             /(sl(k))
        if(sl(k).ne.0.0)   &
             rasd(k)=        &
             (rsnt*asun(k)*0. &
             +(rskt*asky(k)-1.*aphu(k)+1.*aphd(k)-    &
             1.*aniu(k)+1.*anid(k)-1.*airu(k)+1.*aird(k)))   &
             /sl(k)
        !   ************************+++  method 2 +++*********************

        !  *************************+++  method 1  +++++++*************************
5252    continue
        if(sl(k).ne.0.)  &
             rasn(k)=  &
             ( ((wee-tniu-rniu)*(1.-phsn)+(wee-tphu-rphu)*phsn)  &

             *rsnt*gl(k)/cos_zenith   &

             +( ((wee-tniu-rniu)*rskt*(1.-phsk)+  &
             (wee-tphu-rphu)*rskt*phsk)*psk(k)*gd(k)/cos_zenith  &
             +(wee-tniu-rniu)*(fnid(k)+ fniu(k))  &
             +(wee-tphu-rphu)*(fphd(k)+fphu(k))  &
             +delf*(fird(k)+firu(k))             ))
        if(sl(k).ne.0.)  &

             rasd(k)=  &
             ((((wee-tniu-rniu)*rskt*(1.-phsk)  &
             +(wee-tphu-rphu)*rskt*phsk)*psk(k)*gd(k)/cos_zenith  &

             + (wee-tniu-rniu)*(fnid(k)+fniu(k))  &
             +(wee-tphu-rphu)*(fphd(k)+fphu(k))  &
             +delf*(fird(k)+firu(k))             ))

     ENDDO  ! end of loop 4076

     DO k=1,kz
        if(sl(k).eq.0.) tsn1(k )=ta1(k )
        if(sl(k).eq.0.) tsd1(k )=ta1(k )
     ENDDO

     DO k=1,nz-1
        fluxle(k)=fluxle(1)
        fluxh(k)=fluxh(1)

        if(sl(k).ne.0.)  &

             tsn1(k)=  &        !  with inertia in leaves i
             (rasn(k) +face*roa*cpa*dhsn(k)*ta1(k)  &
             + 3.*delf*2.*sigma_SB*tsn(k)**4  &
             -face*lv*dvsn(k)*(qsn(k)-qa1(k))  &
             +ch2o*rh2o*hleaf*tsn(k)/dt_mete )/  &
             (face*roa*cpa*dhsn(k)+4.*delf*2.*sigma_SB*tsn(k)**3  &
             +ch2o*rh2o*hleaf/dt_mete   )
        if(sl(k).ne.0.)  &

             tsd1(k)=   &       !  with inertia in leaves i
             (rasd(k) +face*roa*cpa*dhsd(k)*ta1(k)  &
             + 3.*delf*2.*sigma_SB*tsd(k)**4  &
             -face*lv*dvsd(k)*(qsd(k)-qa1(k))  &
             +ch2o*rh2o*hleaf*tsd(k)/dt_mete )/  &
             ( face*roa*cpa*dhsd(k)+4.*delf*2.*sigma_SB*tsd(k)**3  &
             +ch2o*rh2o*hleaf/dt_mete   )  ! same as for surface temperature
     ENDDO

     DO k=2,kz-1
        fktt=(df(k)*alt1(k+1)*kt(k+1)+(1.-df(k))*alt1(k)*kt(k))*dz(k)  !new
        fktd=((1.-df(k-1))*alt1(k-1)*kt(k-1)+df(k-1)*alt1(k)*kt(k))*dz(k+1)  !new

        fluxh3(k)=-cpa*roa*(((ta1(k)-ta1(k-1))/dz(k)+gamma)*fktd+((ta1(k+1)-ta1(k))/dz(k+1)+gamma)*fktt)/(dz(k)+dz(k+1)) !new
        fluxle3(k)=-lv*(((qa1(k)-qa1(k-1))/dz(k))*fktd+((qa1(k+1)-qa1(k))/dz(k+1))*fktt)/(dz(k)+dz(k+1)) !new
     ENDDO

     fluxle3(1)=(fluxle3(2)*(dz(2)+dz(3))**2-fluxle3(3)*dz(2)**2)/( (dz(2)+dz(3))**2-dz(2)**2) ! new
     fluxh3(1)=(fluxh3(2)*(dz(2)+dz(3))**2-fluxh3(3)*dz(2)**2)/( (dz(2)+dz(3))**2-dz(2)**2)  !new

     fluxle3(kz)=fluxle3(kz-1)
     fluxh3(kz)=fluxh3(kz-1)

     gw=0.2
     g0=0.32
     tasoil=ta1(1)
     qasoil=qa1(1)

     fluxles =fluxle3(1)  !new
     fluxle(1) =fluxle3(1)  !new
     fluxhs =fluxh3(1)     !new
     fluxh(1) =fluxh3(1)   !new

     sks00= (alt1(1)*kt1(1)+alt1(2)*kt1(2))/2.  !new

     haags=cpa*roa*sks00/dz(2)     !new

     !============ Surface radiation ===============!
     ff1=fnid(1)+fphd(1)+rsnt*psn(1)+ rskt*psk(1)+fird(1) -fniu(1)-fphu(1)
     !tar=ta1(kmix-1)  ! not neeed for anything


     if(rou.gt.0.05) then
        g0=0.45
        if((ff1-firu(1)).lt.0.) g0=0.1

     else
        g0=0.32
        if((ff1-firu(1)).lt.0.) g0=0.19
     endif
     ! old  pp=g0*(ff1-firu(1))

     !pp=0.5*(3.*tsoil1(10+nmix)-4.*tsoil1(10+nmix-1)+tsoil1(10+nmix-2))/0.10 ! new

     !r soil heat flux
     if ( ISNAN(local_gsoil(nxodrad, 1)) .or. ISNAN(local_gsoil(nxodrad+1, 1)) ) then  ! measurement values missing
         ! based on modelled soil temperatures
         pp=0.5*(3.*tsoil1(10+nmix)-4.*tsoil1(10+nmix-1)+tsoil1(10+nmix-2))/0.10 ! new
     else !- CLEARCUT COMMENT OUT
         ! based on the measurements
         pp = linear_interp(time_in_month, nxodrad, local_gsoil(:, 1), dt_obs)  ! [W m-2], ground heat flux
     endif !- CLEARCUT COMMENT OUT


     !fluxc=(ff1-firu(1))-pp-fluxle(1)  ! not needed for anything

     !old ta1(kmix-1)=(ff1 -pp-fluxles-(1.-dels)*fird(1)+ haags*((ta1(kmix)+ta1(kmix))/2.+gamma*(dz(kmix)))  &
     !old     +3.*sigma_SB*dels*(ta(kmix-1)**4))/ (4.*sigma_SB*dels*(ta(kmix-1)**3))/ (1.+haags/(4.*sigma_SB*dels*(ta(kmix-1)**3)))

     ta1(kmix-1)=(ff1 -pp-fluxles-0.*(1.-dels)*fird(1)+ haags*(ta1(kmix)+gamma*(dz(kmix))) &  !new
          +3.*sigma_SB*dels*(ta(kmix-1)**4))/ (4.*sigma_SB*dels*(ta(kmix-1)**3))/(1.+haags/(4.*sigma_SB*dels*(ta(kmix-1)**3))) !new

     ! fluxhs = -roa*cpa*sks00*( ((ta1(kmix)+ta1(kmix))/2. -ta1(kmix-1))/(dz(kmix))+gamma) !don't needed
     firu(1)=dels*sigma_SB*ta1(kmix-1)**4+(1.-dels)*fird(1)
     tau=0.

     ! guess: snow melting?
     if(ta1(kmix-1).ge.T00.and.rou.gt.0.) then  ! if air temperature > 0C and snow cover present
        ta1(kmix-1)=T00
        fluxhs = -roa*cpa*sks00* ( ((ta1(kmix)+ta1(kmix))/2. -ta1(kmix-1))/(dz(kmix))+gamma )
        firu(1)=dels*sigma_SB*T00**4+(1.-dels)*fird(1)
        tau=max(0.d0,(((ff1-firu(1))-pp-fluxles-fluxhs)+0.*tau ))  ! net radiation on surface - soil heat flux - turbulent heat fluxes  ! this is the offset from balance
        tau2=1.*dt_mete*(sks00*qmgs-tau*rh2o/333500./psnow)/rh2o       ! rh2o = 1000  ! psnow = 300
        rou=rou+tau2
        rou=max(rou,0.d0)
     endif
     tau4=0.

     !     soil parameters for loamy sand
     wgmax=0.410
     wgf=0.150  !  the field capacity of the soil which corresponds to an hydraulic conductivity of 0.1 mm d−1
     wgwilt=0.075
     bsoil=4.38
     cgsat=3.057*0.000001
     psoil=4.
     asoil=0.404
     cw2ref=3.7
     cw1sat=0.098
     wps=-9.  ! moisture potential when the soil is saturated
     tau3=tau/333500.   ! 333500 could be sulamislämpö, and I'm guessing that tau is the energy used for melting snow
                        ! -> water from tau3 melting snow
     wg1(2)=wg(2)+0.*dt_mete*(tau3 -fluxle(nz)/lv)/(rh2o*1.)+0.*dt_mete*(wg2i-wg(2))/SECONDS_IN_ONE_DAY  ! wg1 volumetric water content of the soil  ! eq. (9) in Sog. 2002
     !wg1(2) is constant!
     ! volumetric water content of the soil, time step
     ! d(wg)/dt_mete = f_m -> wg1 = wg + f_hum * dt_mete
     wg1(2)=min(wg1(2),wgmax)
     wgeq=wg1(2)-wgmax*asoil*(wg1(2)/wgmax)**psoil* (1.-wg1(2)/wgmax)**(8.*psoil)
     wg(1)=max(wg(1),0.d0)

     if (wg(1) /= 0) then
        cw1=cw1sat*(wgmax/wg(1))**(0.5*bsoil+1.)
        cw1=min(cw1,4.d0)
     else
        cw1=4
     end if

     if(fluxles.le.0.) cw1=0.*cw1sat
     cw2=cw2ref*wg1(2)/(wgmax+wgf-wg1(2))
     if(rou.gt.0.0) cw1=0.098
     if(rou.gt.0.0) cw2=0.
     wg1(1)=wg(1)+1.*dt_mete*tau3/(rh2o*0.010) -dt_mete* cw1*fluxles/lv/(rh2o*0.010) -dt_mete*cw2*(wg(1)-wgeq)/SECONDS_IN_ONE_DAY+0.*dt_mete*(wg1i-wg(1))/SECONDS_IN_ONE_DAY
     ! volumetric water content of the soil, time step
     ! d(wg)/dt_mete = f_m -> wg1 = wg + f_hum * dt_mete
     ! equation (8) in Sog. 2002
     !   term 1: tau3/(rh2o*0.010)     ! rh2o = density of liquid water  ! 0.010 thickness of the soil layer
                                       ! tau3 could be water from melting snow
     !   term 2: - cw1*fluxles/lv/(rh2o*0.010) ! fluxles evaporation  ! lv = 2501400 ! cw1 hydraulic coefficient C_1
     !   term 3: - cw2*(wg(1)-wgeq)/SECONDS_IN_ONE_DAY   ! wgeq = the equilibrium value  ! cw2 hydraulic coefficient C_2

     !   term 4: 0.*(wg1i-wg(1))/SECONDS_IN_ONE_DAY
     !           wg1i = 0.2

     wg1(1)=max(wg1(1),0.d0)
     if(rou.gt.1.d-2) wg1(1)=min(wg1(1),3.d-2)

     !r using measured value instead
  if (.not. isnan(surf_soil_hum)) then
     wg1(1) = surf_soil_hum
  endif
     cgw=cgsat*exp(0.5*bsoil*log(wgmax/wg1(2))/log(10.))
     temdif= max(418.*exp(-log(abs(wps*(wgmax/wg1(2))**bsoil))-2.7),1.72d-1)     ! thermal conductivity of soil, eq. below C5 in Sog. 2002
     ! wps*(wgmax/wg1(2))**bsoil) = moisture potential of soil tension
     ! wps =  moisture potential when the soil is saturated
     ! wgmax = maximum volumetric water content that a given soil type can hold (saturation water con)
     ! wg1(2) = volumetric water content of the soil, layer 1 is 0.01 m thick and layer 2 is 1 m thick
     ! bsoil = the slope of the retention curve  (?)
     ! see Appendix C in Sog. 2002
     temtran=temdif/((1.-wg1(2))*1270000.+wg1(2)*4180000.)     ! thermal diffusivity
     ! 1270000 = density of the solid soil * specific heat of the solid soil

     if(rou.gt.0.0) then

        if(wg1(1).le.0.03) then
           vla=1.
        else
           vla=1.
        endif

     else


        ! equation below C5 in Sog. 2002
        if(wg1(1).le.wgf) then
           vla=(0.5*(1.-cos(pi*wg1(1)/wgf))+0.*vla)/1.  ! relative humidity
        else
           vla=1.   ! fraction of humidity
        endif

     endif


 !     vla=0.2  !r works better this way

     if(ta1(kmix-1).lt.T00) vla=0.



     ! border condition for humidity
     qa1(kmix-1)=vla* rhovs_buck(ta1(kmix-1))+1.*(1.-vla)*qa1(kmix)

     qa1(kmix-1)=min(qa1(kmix-1), rhovs_buck(ta1(kmix-1)))

     do k=1,kz-1   ! loop 5647
        fluxle(k+1)=fluxle(k)+face*0.5*dz(k+1)* (sl(k+1)*(psn(k+1)*lv*dvsn(k+1)*(qsn(k+1) -qa1(k+1))+  &
             (1.-psn(k+1))*lv*dvsd(k+1)*(qsd(k+1)-qa1(k+1))) +sl(k)*(psn(k)*lv*dvsn(k)*(qsn(k)-qa1(k))+  &
             (1.-psn(k))*lv*dvsd(k)*   (qsd(k)-qa1(k))))   ! latent heat flux

        fluxh(k+1)=fluxh(k)+face*0.5*dz(k+1)* (sl(k+1)*(psn(k+1)*cpa*roa*dhsn(k+1) *(tsn1(k+1)-ta1(k+1))+  &
             (1.-psn(k+1))*cpa*roa*dhsd(k+1) * (tsd1(k+1)-ta1(k+1))) +sl(k)*(psn(k)*cpa*roa*dhsn(k)*(tsn1(k)-ta1(k))+  &
             (1.-psn(k))*cpa*roa*dhsd(k)* (tsd1(k)-ta1(k))))  ! sensible heat flux
     enddo   ! end of 5647

     ! energy balance at the top of the canopy
     balans=fird(nz)+rsnt+rskt -1.*fniu(nz) -firu(nz)-1.*fphu(nz)-pp- 1.*fluxle3(nz)-1.*fluxh3(nz)-0.*tau-0.*tau4   ! C4 & C1 Sog. 2002
     ! energy balance at the surface
     balans1=fnid(1)+fphd(1)+rsnt*psn(1)+ rskt*psk(1)+fird(1) -1.*fniu(1) -firu(1)-1.*fphu(1)-pp- fluxles-fluxhs-1.*tau-tau4   ! C4 & C1 Sog. 2002
     ! tau = 0., tau4 = 0., tau might be energy used for melting snow

    ! pause
    Other_out(:, 1) = ((wee-tniu-rniu)*(1.-phsn)+(wee-tphu-rphu)*phsn)*rsnt*gl/cos_zenith &
                      + ((wee-tniu-rniu)*rskt*(1.-phsk)+ (wee-tphu-rphu)*rskt*phsk)*psk*gd/cos_zenith &
                      + (wee-tniu-rniu)*(fnid+ fniu) &
                      + (wee-tphu-rphu)*(fphd+fphu)  ! Qabs_sn
    Other_out(:, 2) = ((wee-tniu-rniu)*rskt*(1.-phsk) + (wee-tphu-rphu)*rskt*phsk)*psk*gd/cos_zenith &
                      + (wee-tniu-rniu)*(fnid+fniu) &
                      + (wee-tphu-rphu)*(fphd+fphu)  ! Qabs_sd
    Other_out(:, 3) = delf*(fird+firu)  ! Labs
    Other_out(:, 4) = 2.*delf*sigma_SB*tsn1**4  ! Lout_sn
    Other_out(:, 5) = 2.*delf*sigma_SB*tsd1**4  ! Lout_sd
    Other_out(:, 6) = face*cpa*roa*dhsn*(tsn1-ta1)  ! [W (proj. leaf m2)], SH_sn
    Other_out(:, 7) = face*cpa*roa*dhsd*(tsd1-ta1)  ! [W (proj. leaf m2)], SH_sd
    Other_out(:, 8) = face*lv*dvsn*(qsn-qa1)  ! [W (proj. leaf m2)], LH_sn
    Other_out(:, 9) = face*lv*dvsd*(qsd-qa1)  ! [W (proj. leaf m2)], LH_sd
    Other_out(:, 10) = fphu + fniu  ! Su_out
    Other_out(:, 11) = fphd + fnid  ! Sd_out
    Other_out(:, 12) = gl
    Other_out(:, 13) = gd
    Other_out(:, 14) = psn  ! eta
    Other_out(:, 15) = psk  ! mu
    Other_out(:, 16) = rads2  ! radiation on top of the canopy
    Other_out(:, 17) = rsnt  ! direct radiation on top of the canopy
    Other_out(:, 18) = rskt  ! diffuse radiation on top of the canopy
    Other_out(:, 19) = cos_zenith  ! cos(zenith angle)
    Other_out(:, 20) = shre
    Other_out(:, 21) = re
    Other_out(:, 22) = dhsn
    Other_out(:, 23) = dhsd
    Other_out(:, 24) = dvsn
    Other_out(:, 25) = dvsd
    Other_out(:, 26) = kt1
    Other_out(:, 27) = qsn
    Other_out(:, 28) = qsd
    Other_out(:, 29) = fird
    Other_out(:, 30) = firu
    Other_out(:, 31) = pp
    Other_out(:, 32) = balans
    Other_out(:, 33) = balans1
    ! Other_out(:, 22) = qa1   ! absolute humidity of air

    !===== AIR molecules =====!
    pres(1) = linear_interp(time_in_month, nxodrad, local_pres(:, 1), dt_obs)  ! [Pa]
    do k = 2,kz
      pres(k)=pres(k-1) * exp(-9.81 * (z(k)-z(k-1))/Rgas_d/(0.5*(ta1(k)+ta1(k-1))))
    enddo

     ! Calculation of RH and water vapour pressure used in emission modules
     DO k = 1,kz

        ! Calculate saturation vapor pressure of water
        ! based on Seifield & Pandis, p765, unit in Pascal!!!!!!!!!!!!
        EM_ES(k) = svp(ta1(k))

        ! Water vapour pressure in Pa - different options, Pv = rhovRT/Mv
        EM_EW(k) = (qa1(k)/(18.0153d0/1000.0d0)) * 8.314d0 * ta1(k)
        IF (EM_EW(k)>EM_ES(k)) THEN
          EM_EW(k)=EM_ES(k)
        ENDIF
        ! Water vapour mixing ratio (kg/kg) (based on Megan conversion function), rhov/rhod = PvMv/(PdMd)
        EM_WVM(k) = EM_EW(k) * 18.0153d0/ pres(k) / 28.96d0
        EM_WVM(k) = EM_WVM(k) * 1e3_dp  ! [kg kg-1] --> [g kg-1]
        ! RH(k) = 100.0d0* qa1(k)/ rhovs_buck(ta1(k))
        RH(k) = 100.0d0 * EM_EW(k) / EM_ES(k)

        ! Water vapour pressure in Pa - different options
        !EM_EW(k)  = RH(k) * EM_ES(k) / 100  ! based on RH

        ! Limit EM_EW to maximum value of ES
        !IF (EM_EW(k) .GT. EM_ES(k)) THEN
        !   EM_EW(k) = EM_ES(k)
        !ENDIF

     ENDDO

  ! Replace old values with new values
  CALL New2Old(kz, wg, wg1, tsoil, tsoil1, l, l1,  bt, dbt, bt1, dbt1, u,u1, v,v1,  &
               w, w1, kt, kt1, ta, ta1,alt, alt1, tsn, tsn1, tsd, tsd1, qa, qa1)

end subroutine update_meteorology_stat


FUNCTION diff(p,kz,dz,k)
  ! some kind of difference derivative(?), uses only indices k-1, k and k+1
  INTEGER :: kz,k
  REAL(kind=dp) :: diff,p(kz),dz(kz),da1,db1,dc1
  da1=dz(k+1)*(dz(k)+dz(k+1))
  db1=dz(k)*dz(k+1)
  dc1=dz(k)*(dz(k)+dz(k+1))
  diff=dz(k)*p(k+1)/da1+(dz(k+1)-dz(k))*p(k)/db1-dz(k+1)*p(k-1)/dc1
  RETURN
END FUNCTION diff

FUNCTION diffall(p)
  REAL(KIND=dp) :: diffall(kz), p(kz)

  diffall(1) = (p(2)-p(1))/dz(2)
  diffall(2:kz-1) = dz(2:kz-1)*p(3:kz)/da(2:kz-1) + (dz(3:kz)-dz(2:kz-1))*p(2:kz-1)/db(2:kz-1) - dz(3:kz)*p(1:kz-2)/dc(2:kz-1)
  diffall(kz) = (p(kz)-p(kz-1))/dz(kz)
END FUNCTION diffall

!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine to solve the vertical exchange of gases - Thomas algorithm
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   a	=	array of matrix elements left of the diagonal
!   c	=	array of matrix diagonal elements
!   b	=	array of matrix elements right of the diagonal
!   d	=	input right-hand-side vector
!   p	=	output vector
!   l1	=	lower boundary condition option indicator
!   a1	=	lower boundary value
!   q1	=	lower boundary value
!   lm	=	upper boundary condition option indicator
!   am	=	upper boundary value
!   qm	=	upper boundary value
!   kz	=	upper level number
!   l2	=	lower level start number

SUBROUTINE gtri(a, c, b, d, p, l1, a1, q1, lm, am, qm, kz, l2)
  IMPLICIT NONE
  ! interface:
  INTEGER, INTENT(in) :: l1,l2,lm,kz
  REAL(kind=dp), INTENT(in) :: a1,q1,am,qm
  REAL(kind=dp), INTENT(in) :: a(kz),b(kz),c(kz),d(kz)
  REAL(kind=dp), INTENT(inout) :: p(kz)

  ! local:
  INTEGER :: m2,m,n
  REAL(kind=dp) :: den
  REAL(kind=dp) :: e(kz),f(kz)

  IF(l1.EQ.1) e(l2-1)=0.0d0
  IF(l1.EQ.1) f(l2-1)=a1
  IF(l1.EQ.2) e(l2-1)=1.0d0
  IF(l1.EQ.2) f(l2-1)=-a1
  IF(l1.EQ.3) e(l2-1)=a1/(a1-1.0d0)
  IF(l1.EQ.3) f(l2-1)=q1/(1.0d0-a1)

  m2=kz-1

  DO m=l2,m2
    den=b(m)-c(m)*e(m-1)
    f(m)=(d(m)+c(m)*f(m-1))/den
    e(m)=a(m)/den
  ENDDO

  IF(lm.EQ.1) p(kz)=am
  IF(lm.EQ.2) p(kz)=(f(m2)+am)/(1.0d0-e(m2))
  IF(lm.EQ.3) p(kz)=(f(m2)+qm/am)/((1.0d0+am)/am-e(m2))

  DO n=1,m2-l2+2
    m=kz-n
    p(m)=e(m)*p(m+1)+f(m)
  ENDDO
END SUBROUTINE gtri


!-------------------------------------------------------------------------------
! Get the data at current time from the input data array.
! Time of month for data_in: 0, dt, dt*2, dt*3, ...
! time_in_month: current time
! ind1: left index of the current time period, ind1 = INT(time_in_month/dt) + 1
! ind2: right index of the current time period, ind2 = ind1 + 1
!-------------------------------------------------------------------------------
SUBROUTINE interp_input_to_now(time_in_month, dt, data_in, data_out)
  REAL(dp), INTENT(IN) :: time_in_month
  REAL(dp), INTENT(IN) :: dt
  REAL(dp), INTENT(IN) :: data_in(:)
  REAL(dp), INTENT(OUT) :: data_out

  INTEGER :: ind1, ind2

  ind1 = INT(time_in_month/dt)+1
  ind2 = ind1 + 1

  data_out = data_in(ind1) + ( time_in_month - (ind1-1)*dt ) / dt * ( data_in(ind2) - data_in(ind1) )
END SUBROUTINE interp_input_to_now


subroutine find_boundary_layer(z, Ri, bottom_layer, pblh, pblh_ind)
  ! Layer heights
  real(dp), intent(in) :: z(:)

  ! Richardson number
  real(dp), intent(in) :: Ri(:)

  ! The pblh should be at this layer or above
  ! station: usually set to top layer of canopy
  ! trajectory: usually set to 2 which is just above the ground
  integer , intent(in) :: bottom_layer

  real(dp), intent(out) :: pblh
  integer , intent(out) :: pblh_ind

  ! Critical Richardson number which is used to identify the boundary layer
  real(dp), parameter :: Ri_c = 0.25d0

  ! Number of layers, equals to size(z) and size(Ri)
  integer :: nlayer

  integer :: k

  nlayer = size(z)
  pblh_ind = bottom_layer
  pblh     = z(bottom_layer)

  ! Search from the top to find a layer starting turbulence
  do k=nlayer-1, bottom_layer, -1
    if ( Ri(k) <= Ri_c ) then
      pblh_ind = k
      pblh     = z(k)
      exit
    end if
  end do
end subroutine find_boundary_layer


END MODULE MT_MainMet
