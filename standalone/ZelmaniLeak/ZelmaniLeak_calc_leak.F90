subroutine calc_leak(eos_params,do_heat,rho_gf,inv_eps_gf,&
     GRHydro_hot_atmo_temp,GRHydro_Y_e_min, GRHydro_Y_e_max,&
     rho,temp,ye,chi,tau,heatflux,heaterms,heateave,&
     depsdt,dyedt,ldt,lum,eave,heatout,netheatout,reflevel,rad)
! WARNING: Be careful when changing the arguments to this function; it is also
! called from calc_tau to get the luminosity available for heating
! along the rays.

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

  type(NRPyEOS_params), intent(in):: eos_params
  integer, intent(in) :: do_heat
  real*8, intent(in)  :: rho_gf, inv_eps_gf
  real*8, intent(in)  :: GRHydro_hot_atmo_temp, GRHydro_Y_e_min, GRHydro_Y_e_max
  real*8, parameter   :: f_heat = 1.0d0
  !-----------------------------------

  real*8, intent(in) :: rho ! density in g/cm^3
  real*8, intent(in) :: temp ! temperature in MeV
  real*8, intent(in) :: ye ! ye, dimensionless

  real*8, intent(in) :: chi(3) !chi calculated from Rosswog scheme, one for each nu, to be interpolated throughout 3D
  real*8, intent(in) :: tau(3) !tau calculated from Rosswog scheme, one for each nu, to be interpolated throughout 3D
  real*8, intent(in) :: heatflux(3) !flux used for heating, one for each nu, interpolated throughout 3D
  real*8, intent(in) :: heaterms(3) !rms neutrino energy at neutrinosphere, one for each nu, interpolated throughtout 3D
  real*8, intent(in) :: heateave(3) !average neutrino energy at neutrinosphere, one for each nu, interpolated throughtout 3D
  real*8, intent(out) :: lum(3) ! local dE/dt (luminosities)
  real*8, intent(out) :: eave(3) ! local <E> (average energy)
  real*8, intent(out) :: heatout(3) ! local Q+ ergs/s/cm^3 (heating)
  real*8, intent(out) :: netheatout(3) ! local net heating ergs/s/cm^3 (heating)
  real*8, intent(out) :: depsdt !change in the internal energy, ergs/cm^3/s
  real*8, intent(out) :: dyedt !change in electron fraction
  real*8, intent(in) :: rad ! radius of the point we are dealing with
  integer, intent(in) :: reflevel

  !EOS & local variables
  ! real*8 :: precision = 1.0d-10
  real*8 :: matter_rho,matter_temperature,matter_ye
  real*8 :: matter_enr,matter_prs,matter_ent
  real*8 :: matter_cs2,matter_dedt,matter_dpdrhoe
  real*8 :: matter_dpderho,matter_xa,matter_xh
  real*8 :: matter_xn,matter_xp,matter_abar
  real*8 :: matter_zbar,matter_mue,matter_mun
  real*8 :: matter_mup,matter_muhat

  real*8 :: eta_e,eta_p,eta_n
  real*8 :: eta_hat,eta_nue,eta_nua,eta_nux
  real*8 :: eta_pn, eta_np
  real*8 :: fermi_fac

  !for heating
  real*8 :: heat_const !preamble
  ! real*8 :: f_heat !factor for increasing heating -- set via Cactus parameter
  real*8 :: F(2) !Inverse flux factor
  real*8 :: heat_eff(2) !effective heating

  !constants & parameters
  real*8, parameter :: Qnp = 1.293333d0 !m_n - m_p, MeV
  real*8, parameter :: Cv = 0.5d0 + 2.0d0*0.23d0 !vector coupling
  real*8, parameter :: Ca = 0.5d0 !axial coupling
  real*8, parameter :: alpha = 1.23d0 !gA
  real*8, parameter :: me_mev = 0.510998910d0 !mass of electron in MeV
  real*8, parameter :: sigma_0 = 1.76d-44 ! in units of cm^2
  real*8, parameter :: avo = 6.0221367d23 !Avogadro's number
  real*8, parameter :: pi = 4d0*atan2(1d0,1d0) !3.141592653589793238462d0 !pi
  real*8, parameter :: clite = 29979245800.0d0 !speed of light
  real*8, parameter :: hc_mevcm = 1.23984172d-10 !hc in units of MeV*cm
  real*8, parameter :: gamma_0 = 5.565d-2 ! dimensionless
  real*8, parameter :: fsc = 1.0d0/137.036d0 ! fine structure constant, dimensionless
  real*8, parameter :: mev_to_erg = 1.60217733d-6 !conversion
  real*8, parameter :: massn_cgs = 1.674927211d-24 !neutron mass in grams
  real*8, parameter :: verysmall = 1.0d-20
  integer, parameter :: eoskey = 4
  integer, parameter :: npoints = 1

  !diffusion
  real*8 :: scattering_kappa,abs_kappa
  real*8 :: kappa_tilde_nu_scat(3,3)
  real*8 :: kappa_tilde_nu_abs(3,3)
  real*8 :: block_factor
  real*8 :: zeta(3)
  real*8 :: rate_const
  real*8 :: R_diff(3),Q_diff(3)

  !free emission
  real*8 :: beta
  real*8 :: gamma_const,gamma, R_gamma
  real*8 :: block_factor_e,block_factor_a,block_factor_x
  real*8 :: R_pair, Q_pair, pair_const
  real*8 :: enr_m,enr_p,enr_tilde_m,enr_tilde_p
  real*8 :: R_loc(3),Q_loc(3)

  !leakage
  real*8 :: R_eff(3),Q_eff(3)

  !function
  real*8 :: get_fermi_integral_leak

  real*8 :: ldt

  !cactus relates stuff
  ! character(len=512) :: warnline

  matter_rho = rho
  matter_temperature = max(temp,GRHydro_hot_atmo_temp)
  matter_ye = min(max(ye,GRHydro_Y_e_min),GRHydro_Y_e_max)

  ! call EOS_Omni_full(eoskey,keytemp,GRHydro_eos_rf_prec,npoints,&
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
  matter_enr = matter_enr * inv_eps_gf
  ! if (anyerr.ne.0) then
  !    !$OMP CRITICAL
  !    write(*,"(A15,1P10E15.6)") "rho: ", matter_rho
  !    ! call CCTK_WARN(1,warnline)
  !    write(*,"(A15,1P10E15.6)") "temperature: ", matter_temperature
  !    ! call CCTK_WARN(1,warnline)
  !    write(*,"(A15,1P10E15.6)") "ye: ", matter_ye
  !    ! call CCTK_WARN(1,warnline)
  !    ! write(*,"(A15,i10)") "eos error", keyerr
  !    ! call CCTK_WARN(1,warnline)
  !    ! call CCTK_WARN(0,"set_eos_variables: EOS error in leakage calc_leak")
  !    !$OMP END CRITICAL
  ! endif

  ! don't do anything outside the shock
  if(matter_xh.gt.0.5.and.matter_rho.lt.1.0d13) then
     lum(1:3) = 0.0d0
     depsdt = 0.0d0
     dyedt = 0.0d0
     eave = 0.0d0
     return
  endif

  !Don't do anything in the atmosphere
  ! if(matter_rho*rho_gf.lt. &
  !      (1.0d0+GRHydro_atmo_tolerance)*rho_min) then
  !    lum(1:3) = 0.0d0
  !    depsdt = 0.0d0
  !    dyedt = 0.0d0
  !    eave = 0.0d0
  !    return
  ! endif

  eta_e = matter_mue/temp
  eta_p = matter_mup/temp
  eta_n = matter_mun/temp


  eta_hat = eta_n - eta_p - Qnp/temp
  eta_nue = eta_e - eta_n + eta_p !fully includes effects of rest masses
  eta_nua = -eta_nue
  eta_nux = 0.0d0

  !interpolate etas like we do in GR1D
  eta_nue = eta_nue*(1.0d0-exp(-tau(1)))
  eta_nua = eta_nua*(1.0d0-exp(-tau(2)))

  scattering_kappa = rho*avo*0.25d0*sigma_0/me_mev**2
  kappa_tilde_nu_scat(1,1) = matter_xn*scattering_kappa
  kappa_tilde_nu_scat(1,2) = matter_xp*scattering_kappa
  kappa_tilde_nu_scat(2,1) = matter_xn*scattering_kappa
  kappa_tilde_nu_scat(2,2) = matter_xp*scattering_kappa
  kappa_tilde_nu_scat(3,1) = matter_xn*scattering_kappa
  kappa_tilde_nu_scat(3,2) = matter_xp*scattering_kappa

  scattering_kappa = rho*avo*0.0625d0*sigma_0/me_mev**2* &
       matter_abar*(1.0d0-matter_zbar/matter_abar)**2

  kappa_tilde_nu_scat(1,3) = matter_xh*scattering_kappa
  kappa_tilde_nu_scat(2,3) = matter_xh*scattering_kappa
  kappa_tilde_nu_scat(3,3) = matter_xh*scattering_kappa

  eta_pn = avo*rho*(matter_xn-matter_xp)/(exp(eta_hat)-1.0d0)
  eta_pn = max(eta_pn,0.0d0)
  eta_np = avo*rho*(matter_xp-matter_xn)/(exp(-eta_hat)-1.0d0)
  eta_np = max(eta_np,0.0d0)

  if (rho.lt.1.0d11) then
     !non degenerate here, use mass fractions as chemical potentials fail at low densities
     eta_pn = avo*rho*matter_xp
     eta_np = avo*rho*matter_xn
  endif

  !absorption
  abs_kappa = (1.0d0+3.0d0*alpha**2)*0.25d0*sigma_0/me_mev**2
  block_factor = 1.0d0 + exp(eta_e-get_fermi_integral_leak(5,eta_nue)/ &
       (get_fermi_integral_leak(4,eta_nue)+verysmall))
  kappa_tilde_nu_abs(1,1) = eta_np*abs_kappa/block_factor
  kappa_tilde_nu_abs(2,1) = 0.0d0 !no absorption of a-type on neutrons
  kappa_tilde_nu_abs(3,1) = 0.0d0 !no absorption of x-type neutrinos
  kappa_tilde_nu_abs(1,2) = 0.0d0 !no absorption of e-type on protons
  block_factor = 1.0d0 + exp(-eta_e-get_fermi_integral_leak(5,eta_nua)/ &
       (get_fermi_integral_leak(4,eta_nua)+verysmall))
  kappa_tilde_nu_abs(2,2) = eta_pn*abs_kappa/block_factor
  kappa_tilde_nu_abs(3,2) = 0.0d0 !no absorption of x-type neutrinos
  kappa_tilde_nu_abs(1,3) = 0.0d0 !no absorption on nuclei
  kappa_tilde_nu_abs(2,3) = 0.0d0 !no absorption on nuclei
  kappa_tilde_nu_abs(3,3) = 0.0d0 !no absorption on nuclei

  !sum up opacities to get zeta (again, factoring out energy dependence)
  zeta(1) = kappa_tilde_nu_scat(1,1) + kappa_tilde_nu_scat(1,2) + &
       kappa_tilde_nu_scat(1,3) + kappa_tilde_nu_abs(1,1) + &
       kappa_tilde_nu_abs(1,2) + kappa_tilde_nu_abs(1,3)

  zeta(2) = kappa_tilde_nu_scat(2,1) + kappa_tilde_nu_scat(2,2) + &
       kappa_tilde_nu_scat(2,3) + kappa_tilde_nu_abs(2,1) + &
       kappa_tilde_nu_abs(2,2) + kappa_tilde_nu_abs(2,3)

  zeta(3) = kappa_tilde_nu_scat(3,1) + kappa_tilde_nu_scat(3,2) + &
       kappa_tilde_nu_scat(3,3) + kappa_tilde_nu_abs(3,1) + &
       kappa_tilde_nu_abs(3,2) + kappa_tilde_nu_abs(3,3)

  rate_const = 4.0d0*pi*clite*zeta(1)/(hc_mevcm**3*6.0d0*chi(1)**2)
  R_diff(1) = rate_const*temp*get_fermi_integral_leak(0,eta_nue)
  Q_diff(1) = rate_const*temp**2*get_fermi_integral_leak(1,eta_nue)

  rate_const = 4.0d0*pi*clite*zeta(2)/(hc_mevcm**3*6.0d0*chi(2)**2)
  R_diff(2) = rate_const*temp*get_fermi_integral_leak(0,eta_nua)
  Q_diff(2) = rate_const*temp**2*get_fermi_integral_leak(1,eta_nua)

  rate_const = 16.0d0*pi*clite*zeta(3)/(hc_mevcm**3*6.0d0*chi(3)**2)
  R_diff(3) = rate_const*temp*get_fermi_integral_leak(0,eta_nux)
  Q_diff(3) = rate_const*temp**2*get_fermi_integral_leak(1,eta_nux)

  !now for the free emission
  beta = pi*clite*(1.0d0+3.0d0*alpha**2)*sigma_0/(hc_mevcm**3*me_mev**2)

  R_loc = 0.0d0
  Q_loc = 0.0d0

  R_loc(1) = beta*eta_pn*temp**5*get_fermi_integral_leak(4,eta_e)
  Q_loc(1) = beta*eta_pn*temp**6*get_fermi_integral_leak(5,eta_e)
  R_loc(2) = beta*eta_np*temp**5*get_fermi_integral_leak(4,-eta_e)
  Q_loc(2) = beta*eta_np*temp**6*get_fermi_integral_leak(5,-eta_e)

  !e-e+ pair processes from Ruffert et al.
  fermi_fac = (get_fermi_integral_leak(4,eta_e)/ &
       (get_fermi_integral_leak(3,eta_e)+verysmall) + &
       get_fermi_integral_leak(4,-eta_e)/&
       (get_fermi_integral_leak(3,-eta_e)+verysmall))
  block_factor_e = 1.0d0+exp(eta_nue-0.5d0*fermi_fac)
  block_factor_a = 1.0d0+exp(eta_nua-0.5d0*fermi_fac)
  block_factor_x = 1.0d0+exp(eta_nux-0.5d0*fermi_fac)

  enr_m = 8.0d0*pi/hc_mevcm**3*temp**4*get_fermi_integral_leak(3,eta_e)
  enr_p = 8.0d0*pi/hc_mevcm**3*temp**4*get_fermi_integral_leak(3,-eta_e)

  enr_tilde_m = 8.0d0*pi/hc_mevcm**3*temp**5*get_fermi_integral_leak(4,eta_e)
  enr_tilde_p = 8.0d0*pi/hc_mevcm**3*temp**5*get_fermi_integral_leak(4,-eta_e)

  pair_const = sigma_0*clite/me_mev**2*enr_m*enr_p

  R_pair =  pair_const/(36.0d0*block_factor_e*block_factor_a)* &
       ((Cv-Ca)**2+(Cv+Ca)**2)
  R_loc(1) = R_loc(1) + R_pair
  Q_loc(1) = Q_loc(1) + R_pair*0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)
  R_loc(2) = R_loc(2) + R_pair
  Q_loc(2) = Q_loc(2) + R_pair*0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)

  R_pair =  pair_const/(9.0d0*block_factor_x**2)*((Cv-Ca)**2+(Cv+Ca-2.0d0)**2)
  R_loc(3) = R_loc(3) + R_pair
  Q_loc(3) = Q_loc(3) + R_pair*0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)

  !plasmon decay from Ruffert et al.
  gamma = gamma_0*sqrt((pi**2+3.0d0*eta_e**2)/3.0d0)
  block_factor_e = 1.0d0 + exp(eta_nue-(1.0d0+0.5d0*gamma**2/(1.0d0+gamma)))
  block_factor_a = 1.0d0 + exp(eta_nua-(1.0d0+0.5d0*gamma**2/(1.0d0+gamma)))
  block_factor_x = 1.0d0 + exp(eta_nux-(1.0d0+0.5d0*gamma**2/(1.0d0+gamma)))

  gamma_const = pi**3*sigma_0*clite*temp**8/(me_mev**2*3.0d0*fsc*hc_mevcm**6)* &
       gamma**6*exp(-gamma)*(1.0d0+gamma)

  R_gamma = Cv**2*gamma_const/(block_factor_e*block_factor_a)
  R_loc(1) = R_loc(1) + R_gamma
  Q_loc(1) = Q_loc(1) + R_gamma*0.5d0*temp*(2.0d0+gamma**2/(1.0d0+gamma))
  R_loc(2) = R_loc(2) + R_gamma
  Q_loc(2) = Q_loc(2) + R_gamma*0.5d0*temp*(2.0d0+gamma**2/(1.0d0+gamma))

  R_gamma = (Cv-1.0d0)**2*4.0d0*gamma_const/block_factor_x**2
  R_loc(3) = R_loc(3) + R_gamma
  Q_loc(3) = Q_loc(3) + R_gamma*0.5d0*temp*(2.0d0+gamma**2/(1.0d0+gamma))


  R_pair = 0.231d0*(2.0778d2/mev_to_erg)*0.5d0* &
       (matter_xn**2+matter_xp**2+28.0d0/3.0d0*matter_xn*matter_xp)* &
       rho**2*temp**(4.5d0)
  Q_pair = R_pair*temp/0.231d0*0.504d0

  R_loc(1) = R_loc(1) + R_pair
  Q_loc(1) = Q_loc(1) + Q_pair

  R_loc(2) = R_loc(2) + R_pair
  Q_loc(2) = Q_loc(2) + Q_pair
  R_loc(3) = R_loc(3) + 4.0d0*R_pair
  Q_loc(3) = Q_loc(3) + 4.0d0*Q_pair

  !now calculate leakage!!!!!!!
  R_eff(:) = R_loc(:)/(1.0d0+R_loc(:)/R_diff(:))
  Q_eff(:) = Q_loc(:)/(1.0d0+Q_loc(:)/Q_diff(:))

  eave(1:3) = Q_eff(1:3)/R_eff(1:3)


  if (do_heat.ne.0) then

     heat_const = f_heat*(1.0d0+3.0d0*alpha**2) * sigma_0 / (me_mev**2 * massn_cgs) * 0.25d0 !cm^2/MeV^2/g

     F(1:2) = (4.275d0*tau(1:2)+1.15d0)*exp(-2.0d0*tau(1:2))

     heat_eff(1) = heat_const * rho * matter_xn * heatflux(1) * &
          heaterms(1)**2 * F(1) / mev_to_erg ! cm^2/MeV^2/g * g/cm^3 erg/s/cm^2 MeV^2 MeV/erg = MeV/cm^3/s
     heat_eff(2) = heat_const * rho * matter_xp * heatflux(2) * &
          heaterms(2)**2 * F(2) / mev_to_erg ! cm^2/MeV^2/g * g/cm^3 erg/s/cm^2 MeV^2 MeV/erg = MeV/cm^3/s

     lum(1:2) = (Q_eff(1:2)-heat_eff(1:2))*mev_to_erg !ergs/cm^3/s
! debug:
!     lum(1:2) = (Q_eff(1:2))*mev_to_erg !ergs/cm^3/s
     lum(3) = Q_eff(3)*mev_to_erg !ergs/cm^3/s

     heatout(1:2) = heat_eff(1:2)*mev_to_erg !ergs/cm^3/s
     heatout(3) = 0.0d0

     netheatout(1) = abs(min(0.0d0,lum(1)))
     netheatout(2) = abs(min(0.0d0,lum(2)))
     netheatout(3) = 0.0d0

     depsdt = -sum(lum(1:3))/rho
     dyedt = -(R_eff(1)-R_eff(2)+heat_eff(2)/heateave(2)-heat_eff(1)/heateave(1))*massn_cgs/rho
! debug:
!     dyedt = -(R_eff(1)-R_eff(2))*massn_cgs/rho

  else
     lum(1:3) = Q_eff(1:3)*mev_to_erg
     depsdt = -sum(Q_eff(1:3))*mev_to_erg/rho
     dyedt = -(R_eff(1)-R_eff(2))*massn_cgs/rho
  endif

  if (ye+dyedt*ldt.le.eos_params%eos_yemin*1.05d0) then
     dyedt = 0.0d0
     depsdt = 0.0d0
  endif

  if (ye+dyedt*ldt.ge.eos_params%eos_yemax*0.95d0) then
     dyedt = 0.0d0
     depsdt = 0.0d0
  endif


  if (matter_enr+depsdt*ldt.lt.-eos_params%energy_shift*inv_eps_gf) then !.and.&
!       (reflevel.eq.-1.or.reflevel.ge.grhydro_c2p_warn_from_reflevel)) then
     !$OMP CRITICAL
     ! call CCTK_WARN(1,"Problem in leakage; energy change too large")
     write(*,"(A15,i10)") "reflevel: ", reflevel
     ! call CCTK_WARN(1,warnline)
     write(*,"(A15,1P10E15.6)") "rho: ", rho
     ! call CCTK_WARN(1,warnline)
     write(*,"(A15,1P10E15.6)") "temp: ", temp
     ! call CCTK_WARN(1,warnline)
     write(*,"(A15,1P10E15.6)") "Y_e: ", ye
     ! call CCTK_WARN(1,warnline)
     write(*,"(A15,1P10E15.6)") "dyedt*ldt ", dyedt*ldt
     ! call CCTK_WARN(1,warnline)
     write(*,"(A15,1P10E15.6)") "X_n,X_p,X_a,X_h: ", matter_xn, matter_xp, &
          matter_xa, matter_xh
     ! call CCTK_WARN(1,warnline)
     write(*,"(A15,1P10E15.6)") "eps: ", matter_enr
     ! call CCTK_WARN(1,warnline)
     write(*,"(A15,1P10E15.6)") "depsdt*ldt:", depsdt*ldt
     ! call CCTK_WARN(1,warnline)
     write(*,"(A30,1P10E15.6)") "matter_enr+depsdt*ldt:", matter_enr+depsdt*ldt
     ! call CCTK_WARN(1,warnline)
     write(*,"(A30,1P10E15.6)") "radius:", rad
     ! call CCTK_WARN(1,warnline)
     ! call CCTK_WARN(0,"aborting")
     !$OMP END CRITICAL
  endif


end subroutine calc_leak
