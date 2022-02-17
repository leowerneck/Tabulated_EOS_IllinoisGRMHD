subroutine calc_taus(eos_params, eos_tempmin, grhydro_hot_atmo_temp, &
     grhydro_y_e_max, grhydro_y_e_min, rho_gf, kappa_gf, &
     rho,temp,ye,oldtauruff,tauruff,chiross, &
     heatflux,heaterms,heateave,nzones,rad,ds,compos,xentropy)

  !-----------------------------------
  ! ADDED FOR STANDALONE COMPILATION
  use NRPyEOS
  implicit none

  interface
     subroutine nrpyeos_full(eos_params,rho,ye,temperature,P,eps, &
          S,cs2,depsdT,dPdeps,dPdrho,X_a,X_h,X_n,X_p,Abar,Zbar, &
          mu_e,mu_n,mu_p,muhat) bind(C)
       use, intrinsic :: iso_c_binding, only: c_double
       use NRPyEOS
       implicit none
       type(NRPyEOS_params), intent(in) :: eos_params
       real(kind=c_double),  intent(in) :: rho,ye,temperature
       real(kind=c_double), intent(out) :: P, eps, S, cs2
       real(kind=c_double), intent(out) :: depsdT, dPdeps, dPdrho
       real(kind=c_double), intent(out) :: X_a, X_h, X_n, X_p, Abar, Zbar
       real(kind=c_double), intent(out) :: mu_e, mu_n, mu_p, muhat
     end subroutine nrpyeos_full
  end interface

  type(NRPyEOS_params), intent(in) :: eos_params
  real*8, intent(in) :: eos_tempmin
  real*8, intent(in) :: grhydro_hot_atmo_temp
  real*8, intent(in) :: grhydro_y_e_max
  real*8, intent(in) :: grhydro_y_e_min
  real*8, intent(in) :: rho_gf
  real*8, intent(in) :: kappa_gf
  !-----------------------------------

  real*8, intent(inout) :: rho(nzones) ! density in g/cm^3
  real*8, intent(inout) :: temp(nzones) ! temperature in MeV
  real*8, intent(inout) :: ye(nzones) ! ye, dimensionless

  integer, intent(in) :: nzones ! number of radial zones
  real*8, intent(in) :: rad(nzones) !radial points, cm
  real*8, intent(in) :: ds(nzones) !line element, sqrt(g_{rr}) * dr

  real*8, intent(inout) :: oldtauruff(nzones,3) !tau used for leakage from last iteration, one for each nu
  real*8, intent(out) :: tauruff(nzones,3) !new tau used for leakage from current iteration, one for each nu

  real*8, intent(out) :: compos(nzones,4)
  real*8, intent(out) :: chiross(nzones,3) !chi calculated from Rosswog scheme, one for each nu, to be interpolated throughout 3D

  real*8, intent(out) :: heatflux(nzones,3) !flux of neutrinos used in heating
  real*8, intent(out) :: heaterms(3) !rms energy of neutrinos at neutrinosphere, one for each nu, to be interpolated throughout 3D
  ! real*8 :: lum_local(nzones,3) ! erg/s/cm^3 energy output or input
  real*8, intent(out) :: heateave(3) !average energy of neutrinos at neutrinosphere, one for each nu, to be interpolated throughout 3D
  real*8, intent(out) :: xentropy(nzones)
  ! real*8, intent(in) :: rho_min !atmosphere density
  logical :: have_old_tau !whether we have rufftau or need to calculate from stratch

  !EOS & local variables
  integer keytemp, keyerr
  ! real*8 :: precision = 1.0d-10
  real*8 :: matter_rho,matter_temperature,matter_ye
  real*8 :: matter_enr,matter_prs,matter_ent
  real*8 :: matter_cs2,matter_dedt,matter_dpdrhoe
  real*8 :: matter_dpderho,matter_xa,matter_xh
  real*8 :: matter_xn,matter_xp,matter_abar
  real*8 :: matter_zbar,matter_mue,matter_mun
  real*8 :: matter_mup,matter_muhat

  real*8 :: eta_e(nzones),eta_p(nzones),eta_n(nzones)
  real*8 :: eta_nue(nzones),eta_nua(nzones),eta_nux(nzones),eta_hat(nzones)
  real*8 :: eta_pn(nzones),eta_np(nzones)
  real*8 :: massfracs_xa(nzones),massfracs_xh(nzones)
  real*8 :: massfracs_xp(nzones),massfracs_xn(nzones)
  real*8 :: massfracs_abar(nzones),massfracs_zbar(nzones)

  integer :: i !counter
  ! integer :: rl

  !heating variables
  real*8 :: radial_luminosity(2),lumrad(nzones,3)
  integer :: ns_location(3)
  ! real*8 :: lum(3)
  ! real*8 :: leak_dummy1,leak_dummy2,leak_dummy3,leak_dummy4(3),leak_dummy5(3)
  ! real*8 :: leak_dummy6(3)
  ! real*8 :: lepton_blocking(nzones,2)

  !constants & parameters
  real*8, parameter :: Qnp = 1.293333d0 !m_n - m_p, MeV
  real*8, parameter :: Cv = 0.5d0 + 2.0d0*0.23d0 !vector coupling
  real*8, parameter :: Ca = 0.5d0 !axial coupling
  real*8, parameter :: alpha = 1.23d0 !gA
  real*8, parameter :: me_mev = 0.510998910d0 !mass of electron in MeV
  real*8, parameter :: sigma_0 = 1.76d-44 ! in units of cm^2
  real*8, parameter :: avo = 6.0221367d23 !Avogadro's number
  real*8, parameter :: pi = 4d0*atan2(1d0,1d0) !3.1415926535897932384d0
  real*8, parameter :: twothirds = 2.0d0/3.0d0

  !Ruffert tau stuff
  real*8 :: kappa_tot(nzones,3) !total kappa, 1/cm
  real*8 :: kappa_tot_p(nzones,3) !previous kappa, for comparing, 1/cm
  real*8 :: kappa_scat_n(nzones,3) !total scattering kappa on neutrons, 1/cm
  real*8 :: kappa_scat_p(nzones,3) !total scattering kappa on protons, 1/cm
  real*8 :: kappa_abs_n(nzones) !total absorption kappa on neutrons, 1/cm
  real*8 :: kappa_abs_p(nzones) !total absorption kappa on protons, 1/cm

  real*8 :: local_eta_nue(nzones) !interpolated eta's, nue
  real*8 :: local_eta_nua(nzones) !interpolated eta's, nua
  real*8 :: local_eta_nux(nzones) !interpolated eta's, nux

  real*8 :: csn_0, csp_0 !crosssection preambles
  real*8 :: kappa_const_scat_p(nzones)
  real*8 :: kappa_const_scat_n(nzones)
  real*8 :: kappa_const_abs(nzones)

  real*8 :: xerr
  real*8 :: xerr_out = 1.0d-10
  real*8 :: xlye,xyn,xyp,xynp,xypn,t1,t2

  integer :: icount
  integer, parameter :: icount_max = 200
  integer, parameter :: eoskey = 4
  ! integer :: anyerr
  integer, parameter :: npoints = 1

  !function declarations
  real*8 :: get_fermi_integral_leak

  !Rosswog chi variables
  real*8 :: scattering_kappa
  real*8 :: kappa_tilde_nu_scat(nzones,3,3)
  real*8 :: kappa_tilde_nu_abs(nzones,3,3)
  real*8 :: block_factor
  real*8 :: abs_kappa

  real*8 :: zeta(nzones,3) !zeta
  real*8 :: chi(nzones,3) !chi

  ! cactus related stuff
  ! character(len=512) :: warnline

!#############################

  if(oldtauruff(1,1).gt.0.0d0) then
     have_old_tau = .true.
  else
     have_old_tau = .false.
  endif

  !first get all the EOS variables
  keytemp = 1
  keyerr = 0
  do i=1,nzones

     ! fix potential undershoots near shock
     if(rho(i).lt.rho(nzones)) rho(i) = rho(nzones)
     if(rho(i).lt.1.0d3.and.i.lt.nzones-1.and.i.gt.1) then
        if(rho(i+1).lt.rho(nzones)) rho(i+1) = rho(nzones)
        rho(i) = 0.5d0*(rho(i-1)+rho(i+1))
     endif

     ! fix potential undershoots near shock
     ! need to make sure that the temperature in particular does
     ! nothing crazy near the shock
     ! 1) make sure it is not within a factor of 2 of the table
     ! bound
     if(temp(i).lt.2.0d0*eos_tempmin.and.i.lt.nzones-1.and.i.gt.1) then
        temp(i) = 0.5d0*(temp(i-1)+temp(i+1))
     endif
     ! 2) make sure it does not drop by more than 50% from one zone to the
     !    next
     if(i.lt.nzones .and. i.gt.1) then
        if( temp(i) .lt. 0.5d0*temp(i-1) ) then
           temp(i) = 0.5d0*(temp(i-1)+temp(i+1))
        endif
     endif

     ! fix potential undershoots near shock
     if(ye(i).lt.0.0d0.or.ye(i).gt.0.53d0.and.(i.lt.nzones-1.and.i.gt.1)) then
        ye(i) = 0.5d0*(ye(i-1)+ye(i+1))
     endif

     matter_rho = rho(i)
     matter_temperature = max(temp(i),GRHydro_hot_atmo_temp)
     matter_ye = max(GRHydro_Y_e_min,min(GRHydro_Y_e_max,ye(i))) ! ye(i)

     ! call eos_omni_full(eoskey,keytemp,precision,npoints,&
     !      matter_rho*rho_gf,matter_enr,matter_temperature,matter_ye, &
     !      matter_prs,matter_ent,matter_cs2,matter_dedt,matter_dpderho, &
     !      matter_dpdrhoe,matter_xa,matter_xh,matter_xn,matter_xp,matter_abar, &
     !      matter_zbar,matter_mue,matter_mun,matter_mup,matter_muhat, &
     !      keyerr,anyerr)
     call nrpyeos_full(eos_params,&
          matter_rho*rho_gf,matter_ye,matter_temperature,matter_prs, &
          matter_enr,matter_ent,matter_cs2,matter_dedt,matter_dpderho,matter_dpdrhoe, &
          matter_xa,matter_xh,matter_xn,matter_xp,matter_abar,matter_zbar,&
          matter_mue,matter_mun,matter_mup,matter_muhat)

     write(*,"(A20,1P10e22.15)") "(ZelmaniLeak) rho = ", matter_rho*rho_gf
     write(*,"(A20,1P10e22.15)") "(ZelmaniLeak) Y_e = ", matter_ye
     write(*,"(A20,1P10e22.15)") "(ZelmaniLeak)  T  = ", matter_temperature
     write(*,"(A20,1P10e22.15)") "(ZelmaniLeak)  P  = ", matter_prs
     write(*,"(A20,1P10e22.15)") "(ZelmaniLeak) eps = ", matter_enr
     write(*,"(A20,1P10e22.15)") "(ZelmaniLeak)  S  = ", matter_ent
     write(*,"(A20,1P10e22.15)") "(ZelmaniLeak) mue = ", matter_mue
     write(*,"(A20,1P10e22.15)") "(ZelmaniLeak) mup = ", matter_mup
     write(*,"(A20,1P10e22.15)") "(ZelmaniLeak) mun = ", matter_mun
     write(*,"(A20,1P10e22.15)") "(ZelmaniLeak) muh = ", matter_muhat
     write(*,"(A20,1P10e22.15)") "(ZelmaniLeak) X_a = ", matter_xa
     write(*,"(A20,1P10e22.15)") "(ZelmaniLeak) X_h = ", matter_xh
     write(*,"(A20,1P10e22.15)") "(ZelmaniLeak) X_n = ", matter_xn
     write(*,"(A20,1P10e22.15)") "(ZelmaniLeak) X_p = ", matter_xp
     

     if (keyerr.ne.0) then
        !$OMP CRITICAL
        write(*,"(A15,1P10E15.6)") "rho: ", matter_rho
        ! call CCTK_WARN(1,warnline)
        write(*,"(A15,1P10E15.6)") "rho(i): ", rho(i)
        ! call CCTK_WARN(1,warnline)
        write(*,"(A15,1P10E15.6)") "rho(i-1): ", rho(i-1)
        ! call CCTK_WARN(1,warnline)
        write(*,"(A15,1P10E15.6)") "rho(i+1): ", rho(i+1)
        ! call CCTK_WARN(1,warnline)
        write(*,"(A15,1P10E15.6)") "temperature: ", matter_temperature
        ! call CCTK_WARN(1,warnline)
        write(*,"(A15,1P10E15.6)") "ye: ", matter_ye
        ! call CCTK_WARN(1,warnline)
        write(*,"(A15,i10)") "eos error", keyerr
        ! call CCTK_WARN(1,warnline)
        write(*,"(A15,i10,1P10E15.6)") "location: ", i, rad(i)
        ! call CCTK_WARN(1,warnline)
        ! call CCTK_WARN(0,"set_eos_variables: EOS error in leakage calc_tau")
        !$OMP END CRITICAL
     endif
     xentropy(i) = matter_ent
     compos(i,1) = matter_xh
     compos(i,2) = matter_xn
     compos(i,3) = matter_xp
     compos(i,4) = matter_xa

     ! in our EOS the rest mass difference is in the chemical potentials of the neucleons
     eta_e(i) = matter_mue/matter_temperature
     eta_p(i) = matter_mup/matter_temperature
     eta_n(i) = matter_mun/matter_temperature

     massfracs_xa(i) = matter_xa
     massfracs_xh(i) = matter_xh
     massfracs_xp(i) = matter_xp
     massfracs_xn(i) = matter_xn
     massfracs_abar(i) = matter_abar
     massfracs_zbar(i) = matter_zbar

     eta_hat(i) = eta_n(i)-eta_p(i) - Qnp/matter_temperature
     eta_nue(i) = eta_e(i) - eta_n(i) + eta_p(i) !fully includes effects of rest masses
     eta_nua(i) = -eta_nue(i)
     eta_nux(i) = 0.0d0

  enddo

!#############################

  !now calculate Ruffert tau

  !use previous tauruff
  if (have_old_tau) tauruff = oldtauruff

  !initialize
  kappa_tot(:,:)   = 1.0d0 ! 1/cm
  kappa_tot_p(:,:) = 1.0d0 ! 1/cm
  kappa_scat_n(:,:) = 1.0d-5 ! 1/cm
  kappa_scat_p(:,:) = 1.0d-5 ! 1/cm
  kappa_abs_n(:)  = 1.0d-5 ! 1/cm
  kappa_abs_p(:)  = 1.0d-5 ! 1/cm

  local_eta_nux(:) = 0.0d0
  local_eta_nue(:) = 0.0d0
  local_eta_nua(:) = 0.0d0

  !cross section coeffs
  csn_0 = (1.0d0 + 5.0d0*alpha**2) / 24.0d0
  csp_0 = (4.0d0*(Cv-1.0d0)**2 + 5.0d0*alpha**2) / 24.0d0

  do i=1,nzones
     ! constant parts of kappa (A6)
     t1 = sigma_0 * avo * rho(i) * (temp(i)/me_mev)**2
     kappa_const_scat_n(i) = csn_0 * t1
     kappa_const_scat_p(i) = csp_0 * t1
     ! (A11) constant part
     kappa_const_abs(i) = (1.0d0+3.0d0*alpha**2)/4.0d0 * t1
  enddo

  ! Loop to get converged result for tau.
  ! This is discussed in the text between equations (A5) and
  ! (A6). Note that for the initial iteration tau is set to
  ! 1 and the etas are set to 1.0d-5

  icount = 1
  xerr = 1.0d0

  do while(xerr.gt.xerr_out .and. icount.lt.icount_max)
     ! copy over current into previous kappa
     kappa_tot_p = kappa_tot

     ! set up new kappa based on individual contributions
     ! nu_e; (A17)
     kappa_tot(:,1) = &
          kappa_scat_p(:,1) &
          + kappa_scat_n(:,1) &
          + kappa_abs_n(:)
     ! antis; (A18)
     kappa_tot(:,2) = &
              kappa_scat_p(:,2) &
            + kappa_scat_n(:,2) &
            + kappa_abs_p(:)

     ! nu_xs (A19)
     kappa_tot(:,3) = &
          + kappa_scat_p(:,3) &
          + kappa_scat_n(:,3) !&

     ! Integrate optical depths: Equation (A20)
     ! Note that this is not done for particle and energy transport
     if(icount.gt.2.or..not.have_old_tau) then
        tauruff(:,:) = 0.0d0
        do i=nzones-1,1,-1
           tauruff(i,1:3) = tauruff(i+1,1:3) + kappa_tot(i,1:3)*ds(i)
           compos(i,4) = kappa_tot(i,1)
        enddo
     endif

     do i=1,nzones
        ! (A5) equilibrium eta, we have rest masses in our chemical potentials
        ! no need to include mass difference in eta
        local_eta_nux(i) = 0.0d0   ! (A2)
        local_eta_nue(i) = eta_nue(i) * (1.0d0-exp(-tauruff(i,1))) ! (A3); note that the ^0 etas are set to 0.0d0
        local_eta_nua(i) = eta_nua(i) * (1.0d0-exp(-tauruff(i,2))) ! (A4)


        !assuming completely dissociated
        xlye = ye(i)
        xyn = (1.0d0-xlye) / (1.0d0 + 2.0d0/3.0d0 * max(eta_n(i),0.0d0)) ! (A8)
        xyp = xlye / (1.0d0 + 2.0d0/3.0d0*max(eta_p(i),0.0d0))

        if(massfracs_xh(i).lt.0.5d0) then
           t1 = exp(-eta_hat(i))
           xynp = max((2.0d0*xlye-1.0d0)/ (t1-1.0d0),0.0d0) ! (A13)
           xypn = max(xynp * t1,0.0d0) ! (A14)
        else
           xypn = massfracs_xp(i)
           xynp = massfracs_xn(i)
        endif


        ! electron neutrinos
        t1 = get_fermi_integral_leak(5,local_eta_nue(i)) / &
             get_fermi_integral_leak(3,local_eta_nue(i))
        t2 = 1.0d0 + exp(eta_e(i)-get_fermi_integral_leak(5,local_eta_nue(i)) / &
             get_fermi_integral_leak(4,local_eta_nue(i))) ! (A15)

        kappa_scat_n(i,1) = kappa_const_scat_n(i) * xyn  * t1
        kappa_scat_p(i,1) = kappa_const_scat_p(i) * xyp  * t1
        kappa_abs_n(i) = kappa_const_abs(i) * xynp * t1 / t2 ! (A11)

        ! anti-electron neutrinos
        t1 = get_fermi_integral_leak(5,local_eta_nua(i)) / &
             get_fermi_integral_leak(3,local_eta_nua(i))
        t2 = 1.0d0 + exp(-eta_e(i)-get_fermi_integral_leak(5,local_eta_nua(i)) / &
             get_fermi_integral_leak(4,local_eta_nua(i))) ! (A16)

        kappa_scat_n(i,2) = kappa_const_scat_n(i) * xyn  * t1 ! (A6)
        kappa_scat_p(i,2) = kappa_const_scat_p(i) * xyp  * t1
        kappa_abs_p(i) = kappa_const_abs(i) * xypn * t1 / t2 ! (A12)

        ! nux neutrinos
        t1 = get_fermi_integral_leak(5,local_eta_nux(i)) / &
             get_fermi_integral_leak(3,local_eta_nux(i))
        kappa_scat_n(i,3) = kappa_const_scat_n(i) * xyn * t1 ! (A6)
        kappa_scat_p(i,3) = kappa_const_scat_p(i) * xyp * t1

     enddo

     ! compute relative change xerr
     xerr = 0.0d0
     do i=1,nzones
        xerr = max(xerr,abs(kappa_tot(i,1)/kappa_tot_p(i,1)-1.0d0))
        xerr = max(xerr,abs(kappa_tot(i,2)/kappa_tot_p(i,2)-1.0d0))
        xerr = max(xerr,abs(kappa_tot(i,3)/kappa_tot_p(i,3)-1.0d0))
     enddo

     icount = icount + 1

  enddo

  if(icount.ge.icount_max) then
     write(6,"(i5,1P10E15.6)") icount,xerr,xerr_out
     stop "icount > icount_max in leakage; leak.F90"
  endif

  ! Recompute tau based on the most recent kappa_tot
  tauruff = 0.0d0

  do i=nzones-1,1,-1
     tauruff(i,1:3) = tauruff(i+1,1:3) + kappa_tot(i,1:3)*ds(i)
  enddo

  have_old_tau = .true.

  !set degeneracy factors to interpolated values
  eta_nue(:) = local_eta_nue(:)
  eta_nua(:) = local_eta_nua(:)
  eta_nux(:) = local_eta_nux(:)

  write(*,*) "kappa_n_nue  = ",kappa_scat_n(1,1)
  write(*,*) "kappa_p_nue  = ",kappa_scat_p(1,1)
  write(*,*) "kappa_a_nue  = ",kappa_abs_p(1)

  write(*,*) "kappa_n_anue = ",kappa_scat_n(1,2)
  write(*,*) "kappa_p_anue = ",kappa_scat_p(1,2)
  write(*,*) "kappa_a_anue = ",kappa_abs_n(1)

  write(*,*) "kappa_n_nux  = ",kappa_scat_n(1,3)
  write(*,*) "kappa_p_nux  = ",kappa_scat_p(1,3)

!#######################################

!Rosswog chi

  zeta(:,:) = 0.0d0
  chi(:,:) = 0.0d0

  do i=1,nzones

     !scattering
     scattering_kappa = rho(i)*avo*0.25d0*sigma_0/me_mev**2
     kappa_tilde_nu_scat(i,1,1) = massfracs_xn(i)*scattering_kappa
     kappa_tilde_nu_scat(i,1,2) = massfracs_xp(i)*scattering_kappa
     kappa_tilde_nu_scat(i,2,1) = massfracs_xn(i)*scattering_kappa
     kappa_tilde_nu_scat(i,2,2) = massfracs_xp(i)*scattering_kappa
     kappa_tilde_nu_scat(i,3,1) = massfracs_xn(i)*scattering_kappa
     kappa_tilde_nu_scat(i,3,2) = massfracs_xp(i)*scattering_kappa


     scattering_kappa = rho(i)*avo*0.0625d0*sigma_0/me_mev**2* &
          massfracs_abar(i)*(1.0d0-massfracs_zbar(i)/massfracs_abar(i))**2 ! only have 1 factor of A because kappa multiples the number fraction, not mass fractions
     kappa_tilde_nu_scat(i,1,3) = massfracs_xh(i)*scattering_kappa
     kappa_tilde_nu_scat(i,2,3) = massfracs_xh(i)*scattering_kappa
     kappa_tilde_nu_scat(i,3,3) = massfracs_xh(i)*scattering_kappa

     eta_pn(i) = avo*rho(i)*(massfracs_xn(i)-massfracs_xp(i))/(exp(eta_hat(i))-1.0d0)
     eta_pn(i) = max(eta_pn(i),0.0d0)
     eta_np(i) = avo*rho(i)*(massfracs_xp(i)-massfracs_xn(i))/(exp(-eta_hat(i))-1.0d0)
     eta_np(i) = max(eta_np(i),0.0d0)

     if (rho(i).lt.1.0d11) then
        !non degenerate here, use mass fractions as chemical potentials fail at low densities
        eta_pn(i) = avo*rho(i)*massfracs_xp(i)
        eta_np(i) = avo*rho(i)*massfracs_xn(i)
     endif

     !absorption
     abs_kappa = (1.0d0+3.0d0*alpha**2)*0.25d0*sigma_0/me_mev**2
     block_factor = 1.0d0 + exp(eta_e(i)-get_fermi_integral_leak(5,eta_nue(i))/ &
          get_fermi_integral_leak(4,eta_nue(i)))
     kappa_tilde_nu_abs(i,1,1) = eta_np(i)*abs_kappa/block_factor
     kappa_tilde_nu_abs(i,2,1) = 0.0d0 !no absorption of a-type on neutrons
     kappa_tilde_nu_abs(i,3,1) = 0.0d0 !no absorption of x-type neutrinos
     kappa_tilde_nu_abs(i,1,2) = 0.0d0 !no absorption of e-type on protons
     block_factor = 1.0d0 + exp(-eta_e(i)-get_fermi_integral_leak(5,eta_nua(i))/ &
          get_fermi_integral_leak(4,eta_nua(i)))
     kappa_tilde_nu_abs(i,2,2) = eta_pn(i)*abs_kappa/block_factor
     kappa_tilde_nu_abs(i,3,2) = 0.0d0 !no absorption of x-type neutrinos
     kappa_tilde_nu_abs(i,1,3) = 0.0d0 !no absorption on nuclei
     kappa_tilde_nu_abs(i,2,3) = 0.0d0 !no absorption on nuclei
     kappa_tilde_nu_abs(i,3,3) = 0.0d0 !no absorption on nuclei

     !sum up opacities to get zeta (again, factoring out energy dependence)
     zeta(i,1) = kappa_tilde_nu_scat(i,1,1) + kappa_tilde_nu_scat(i,1,2) + &
          kappa_tilde_nu_scat(i,1,3) + kappa_tilde_nu_abs(i,1,1) + &
          kappa_tilde_nu_abs(i,1,2) + kappa_tilde_nu_abs(i,1,3)

     zeta(i,2) = kappa_tilde_nu_scat(i,2,1) + kappa_tilde_nu_scat(i,2,2) + &
          kappa_tilde_nu_scat(i,2,3) + kappa_tilde_nu_abs(i,2,1) + &
          kappa_tilde_nu_abs(i,2,2) + kappa_tilde_nu_abs(i,2,3)

     zeta(i,3) = kappa_tilde_nu_scat(i,3,1) + kappa_tilde_nu_scat(i,3,2) + &
          kappa_tilde_nu_scat(i,3,3) + kappa_tilde_nu_abs(i,3,1) + &
          kappa_tilde_nu_abs(i,3,2) + kappa_tilde_nu_abs(i,3,3)

  enddo

  do i=nzones-1,1,-1
     !integrate zeta to get chi, tau with energy dependence factored out
     chi(i,1) = chi(i+1,1) + zeta(i,1)*ds(i)
     chi(i,2) = chi(i+1,2) + zeta(i,2)*ds(i)
     chi(i,3) = chi(i+1,3) + zeta(i,3)*ds(i)
  enddo

  chi(nzones,:) = chi(nzones-1,:)

  chiross = chi

  !neutrinosphere located at tau = 2/3
  ns_location(:) = 1
  do i=1,nzones
     if (tauruff(i,1).gt.twothirds) then
        ns_location(1) = i
     endif
     if (tauruff(i,2).gt.twothirds) then
        ns_location(2) = i
     endif
     if (tauruff(i,3).gt.twothirds) then
        ns_location(3) = i
     endif
  enddo

  !rms energy at neutrino sphere
  heaterms(1) = temp(ns_location(1))* &
       sqrt(get_fermi_integral_leak(5,eta_nue(ns_location(1)))/ &
       get_fermi_integral_leak(3,eta_nue(ns_location(1))))
  heaterms(2) = temp(ns_location(2))* &
       sqrt(get_fermi_integral_leak(5,eta_nua(ns_location(2)))/ &
       get_fermi_integral_leak(3,eta_nua(ns_location(2))))
  heaterms(3) = temp(ns_location(3))* &
       sqrt(get_fermi_integral_leak(5,eta_nux(ns_location(3)))/ &
       get_fermi_integral_leak(3,eta_nux(ns_location(3))))

  !mean energy at neutrino sphere
  heateave(1) = temp(ns_location(1))* &
       get_fermi_integral_leak(5,eta_nue(ns_location(1)))/ &
       get_fermi_integral_leak(4,eta_nue(ns_location(1)))
  heateave(2) = temp(ns_location(2))* &
       get_fermi_integral_leak(5,eta_nua(ns_location(2)))/ &
       get_fermi_integral_leak(4,eta_nua(ns_location(2)))
  heateave(3) = temp(ns_location(3))* &
       get_fermi_integral_leak(5,eta_nux(ns_location(3)))/ &
       get_fermi_integral_leak(4,eta_nux(ns_location(3)))

  !now we leak along this lines to determine luminosity, will use this for heating!

  radial_luminosity(:) = 0.0d0
  lumrad(:,:) = 0.0d0
  heatflux = 0.0d0

  ! if(do_heat.ne.0) then
  !    do i =1,nzones
  !       lepton_blocking(i,1) = 1.0d0/(1.0d0 + exp(eta_e(i) - &
  !            get_fermi_integral_leak(5,eta_nue(ns_location(1)))/ &
  !            get_fermi_integral_leak(4,eta_nue(ns_location(1)))))
  !       lepton_blocking(i,2) = 1.0d0/(1.0d0 + exp(-eta_e(i) - &
  !            get_fermi_integral_leak(5,eta_nua(ns_location(2)))/ &
  !            get_fermi_integral_leak(4,eta_nua(ns_location(2)))))
  !    enddo

  !    rl = -1
  !    do i=1,nzones-1

  !       leak_dummy1 = 0.0d0
  !       leak_dummy2 = 0.0d0
  !       leak_dummy3 = 0.0d0 !this is ldt, the actual leakage will take care of overflow of dyedt
  !       leak_dummy4 = 0.0d0
  !       leak_dummy5 = 0.0d0
  !       leak_dummy6 = 0.0d0

  !       !call leak with heating turned on to get the change in
  !       !luminosity, then we'll interpolate that and use for heating
  !       call calc_leak(rho(i),temp(i),ye(i),chi(i,:),tauruff(i,:),heatflux(i,:), &
  !            heaterms,heateave,leak_dummy1,leak_dummy2,leak_dummy3,&
  !            lum,leak_dummy4,leak_dummy5,leak_dummy6,rho_min,rl,rad(i))

  !       !lum coming from calc_leak is not actually luminosity, rather,
  !       !must multiply by volume (of spherical shell)
  !       !then add to integrated luminosity
  !       radial_luminosity(1) = radial_luminosity(1) + lum(1)*4.0d0*pi*((rad(i)+rad(i+1))/2.0d0)**2*ds(i)
  !       radial_luminosity(2) = radial_luminosity(2) + lum(2)*4.0d0*pi*((rad(i)+rad(i+1))/2.0d0)**2*ds(i)
  !       lum_local(i,1:3) = lum(1:3)
  !       lumrad(i,1) = leak_dummy5(1)
  !       lumrad(i,2) = leak_dummy5(2)

  !       !this is what gets interpolated and put into the leak routine.
  !       heatflux(i+1,1) = radial_luminosity(1)/(4.0d0*pi*((rad(i)+rad(i+1))/2.0d0)**2)&
  !            * lepton_blocking(i,1) !luminosity/4pir^2
  !       heatflux(i+1,2) = radial_luminosity(2)/(4.0d0*pi*((rad(i)+rad(i+1))/2.0d0)**2)&
  !            * lepton_blocking(i,2) !luminosity/4pir^2

  !    enddo

  ! else
  !    heaterms = 0.0d0
  !    heateave = 0.0d0
  ! endif


end subroutine calc_taus


function get_fermi_integral_leak(ifermi,eta)

  implicit none
  integer ifermi
  real*8 get_fermi_integral_leak
  real*8 eta
  real*8 fermi_integral_analytical

  fermi_integral_analytical = 0.0d0

  ! Expressions for Fermi integrals given in Takahashi et al. 1978
  if (eta.gt.1.D-3) then
     select case (ifermi)
     case (0)
        fermi_integral_analytical = &
             log10(1.0d0+exp(eta))
     case (1)
        fermi_integral_analytical = &
             (eta**2/2.0D0 + 1.6449d0)/(1.0D0+EXP(-1.6855d0*eta))
     case (2)
        fermi_integral_analytical = &
             (eta**3/3.0D0 + 3.2899d0*eta)/(1.0D0-EXP(-1.8246d0*eta))
     case (3)
        fermi_integral_analytical = &
             (eta**4/4.0D0 + 4.9348d0*eta**2+11.3644d0) / &
             (1.0D0+EXP(-1.9039d0*eta))
     case (4)
        fermi_integral_analytical = &
             (eta**5/5.0D0 + 6.5797d0*eta**3+45.4576d0*eta) / &
             (1.0D0-EXP(-1.9484d0*eta))
     case (5)
        fermi_integral_analytical = &
             (eta**6/6.0D0 + 8.2247d0*eta**4 + 113.6439d0*eta**2 + &
             236.5323d0)/(1.0D0+EXP(-1.9727d0*eta))
     end select

  else
     select case (ifermi)
     case (0)
        fermi_integral_analytical = &
             log10(1.0d0+exp(eta))
     case (1)
        fermi_integral_analytical = &
             EXP(eta)/(1.0D0+0.2159d0*EXP(0.8857d0*eta))
     case (2)
        fermi_integral_analytical = &
             2.0D0*EXP(eta)/(1.0D0+0.1092d0*EXP(0.8908d0*eta))
     case (3)
        fermi_integral_analytical = &
             6.0D0*EXP(eta)/(1.0D0+0.0559d0*EXP(0.9069d0*eta))
     case (4)
        fermi_integral_analytical = &
             24.0D0*EXP(eta)/(1.0D0+0.0287d0*EXP(0.9257d0*eta))
     case (5)
        fermi_integral_analytical = &
             120.0D0*EXP(eta) / (1.0D0 + 0.0147d0*EXP(0.9431d0*eta))
     end select

  endif
  get_fermi_integral_leak = fermi_integral_analytical

  return
end function get_fermi_integral_leak
