#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

! eoskey:
! 1 --- polytropic EOS
! 2 --- gamma-law EOS
! 3 --- hybrid EOS
! 4 --- finite-T microphysical NSE EOS

subroutine EOS_Omni_EOS_Press(eoskey,keytemp,rf_precision,npoints,&
                              rho,eps,temp,ye,press,keyerr,anyerr)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS
  

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
  CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
  CCTK_REAL, intent(out)   :: press(npoints)

  ! local vars
  integer          :: i,p
  character(256)   :: warnstring
  real*8           :: hybrideos_local_gamma, hybrideos_local_k
  real*8           :: hybrideos_local_eps, hybrideos_p_poly, hybrideos_p_th
  real*8,parameter :: zero = 0.0d0
  ! temporary vars for nuc_eos
  real*8, dimension(1) :: xrho,xye,xtemp,xenr,xent
  real*8, dimension(1) :: xprs,xmunu,xcs2
  real*8, dimension(1) :: xdedt,xdpderho,xdpdrhoe
  ! temporary vars for coldeos + gamma law
  integer :: ir
  real*8 :: press_cold, press_th
  real*8 :: eps_cold, eps_th
  real*8 :: gamma

  anyerr    = 0
  keyerr(:) = 0

  select case (eoskey)
     case (1)
        ! polytropic EOS
        if(keytemp.eq.1) then
           do i=1,npoints
              eps(i) = press_gf * poly_k_cgs * &
                   (rho(i)*inv_rho_gf)**(poly_gamma) / &
                   (poly_gamma - 1.0d0) / rho(i)
           enddo
        endif
        do i=1,npoints
           press(i) = press_gf * poly_k_cgs * &
                (rho(i)*inv_rho_gf)**poly_gamma
        enddo
     case (2)
        ! gamma-law EOS
        if(keytemp.eq.1) then
           do i=1,npoints
              eps(i) = press_gf * gl_k_cgs * &
                   (rho(i)*inv_rho_gf)**(gl_gamma) / &
                   (gl_gamma - 1.0d0) / rho(i)
           enddo
        endif
        do i=1,npoints
           press(i) = (gl_gamma - 1.0d0) * rho(i) * eps(i)

        enddo
     case (3)
        ! hybrid EOS
        do i=1,npoints
           hybrideos_local_gamma = hybrideos_gamma(1)
           hybrideos_local_k = hybrideos_k(1)
           hybrideos_local_eps = hybrideos_eps(1)
           if (n_pieces .gt. 1) then
              do p = 1,n_pieces-1
                 if (rho(i) .gt. hybrideos_rho(p)) then
                    hybrideos_local_gamma = hybrideos_gamma(p+1)
                    hybrideos_local_k = hybrideos_k(p+1)
                    hybrideos_local_eps = hybrideos_eps(p+1)
                 end if
              end do
           end if
           if (keytemp .eq. 1) then
              !Set epsilon to the cold part only
              eps(i)=hybrideos_local_k/(hybrideos_local_gamma - 1.0)*rho(i)**(hybrideos_local_gamma - 1.0) + hybrideos_local_eps
           end if
           hybrideos_p_poly = hybrideos_local_k * (rho(i))**hybrideos_local_gamma
           hybrideos_p_th = (hybrid_gamma_th-1.0)* (rho(i)*eps(i) - &
                hybrideos_p_poly/(hybrideos_local_gamma-1) - &
                hybrideos_local_eps * rho(i))
           hybrideos_p_th = max(zero, hybrideos_p_th)
           press(i) = hybrideos_p_poly + hybrideos_p_th
        enddo

     case (4)
        ! nuc eos
        if(keytemp.eq.1) then
           call nuc_eos_m_kt1_press_eps(npoints,rho,temp,ye,&
                eps,press,keyerr,anyerr)
        else if(keytemp.eq.0) then
           call nuc_eos_m_kt0_press(npoints,rho,temp,ye,eps,press,&
                rf_precision,keyerr,anyerr)
        else
           call CCTK_ERROR("This keytemp is not suppported!")
           STOP
        endif
     case (5)
        ! cold tabular EOS with gamma law
        do i=1,npoints
           if(rho(i).lt.coldeos_rhomin) then
              press(i) = coldeos_low_kappa * rho(i)**coldeos_low_gamma
              eps(i) = press(i) / (coldeos_low_gamma - 1.0) / rho(i) 
              cycle
           else if(rho(i).gt.coldeos_rhomax) then
              keyerr(i) = 103
              anyerr = 1
           else
              xrho(1) = log10(rho(i))
              ir = 2 + INT( (xrho(1) - coldeos_logrho(1) - 1.0d-10) * coldeos_dlrhoi )
           endif
           ir = max(2, min(ir,coldeos_nrho))

           gamma = (coldeos_gamma(ir) - coldeos_gamma(ir-1)) / &
                (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                (xrho(1) - coldeos_logrho(ir-1)) + coldeos_gamma(ir-1)

           eps_cold = (coldeos_eps(ir) - coldeos_eps(ir-1)) / &
                    (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                    (xrho(1) - coldeos_logrho(ir-1)) + coldeos_eps(ir-1)

           if(keytemp.eq.0) then
              eps_th = max(0.0d0,eps(i) - eps_cold)
           else if (keytemp.eq.1) then
              eps_th = 0.0d0
              eps(i) = eps_cold
           endif

           press_cold = coldeos_kappa * rho(i)**gamma
           press_th = coldeos_thfac*(coldeos_gammath - 1.0d0)*rho(i)*eps_th
           press(i) = press_cold + press_th
        enddo

     case (6)
        ! barotropic tabular EOS with gamma law
        do i=1,npoints
           if(rho(i).lt.barotropiceos_rhomin) then
              press(i) = barotropiceos_low_kappa * rho(i)**barotropiceos_low_gamma
              eps(i) = press(i) / (barotropiceos_low_gamma - 1.0) / rho(i) 
              cycle
           else if(rho(i).gt.barotropiceos_rhomax) then
              keyerr(i) = 103
              anyerr = 1
           else
              xrho(1) = log10(rho(i))
              ir = 2 + &
                   INT( (xrho(1) - barotropiceos_logrho(1) - 1.0d-10) &
                   * barotropiceos_dlrhoi )
           endif
           ir = max(2, min(ir,barotropiceos_nrho))

           xprs(1) = (barotropiceos_logpress(ir) - barotropiceos_logpress(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho(1) - barotropiceos_logrho(ir-1)) + barotropiceos_logpress(ir-1)

           xenr(1) = (barotropiceos_logeps(ir) - barotropiceos_logeps(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho(1) - barotropiceos_logrho(ir-1)) + barotropiceos_logeps(ir-1)

           xtemp(1) = (barotropiceos_temp(ir) - barotropiceos_temp(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho(1) - barotropiceos_logrho(ir-1)) + barotropiceos_temp(ir-1)

           xye(1) = (barotropiceos_ye(ir) - barotropiceos_ye(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho(1) - barotropiceos_logrho(ir-1)) + barotropiceos_ye(ir-1)

           press_cold = 10.0**xprs(1)
           eps_cold = 10.0**xenr(1) - barotropiceos_energyshift

           if(keytemp.eq.0) then
              eps_th = max(0.0d0,eps(i) - eps_cold)
           else if (keytemp.eq.1) then
              eps_th = 0.0d0
              eps(i) = eps_cold
           endif
           press_th = coldeos_thfac*(barotropiceos_gammath - 1.0d0)*rho(i)*eps_th
           press(i) = press_cold + press_th

        enddo

     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
     end select

end subroutine EOS_Omni_EOS_Press

subroutine EOS_Omni_EOS_PressOMP(eoskey,keytemp,rf_precision,npoints,&
                              rho,eps,temp,ye,press,keyerr,anyerr)

  use EOS_Omni_Module
#ifdef _OPENMP
  use OMP_LIB
#endif

  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
  CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
  CCTK_REAL, intent(out)   :: press(npoints)

  ! local vars
  integer          :: i,p
  character(256)   :: warnstring
  real*8           :: hybrideos_local_gamma,  hybrideos_local_k_cgs, hybrideos_local_k, hybrideos_local_eps, &
                      hybrideos_p_poly, hybrideos_p_th
  real*8,parameter :: zero = 0.0d0
  ! temporary vars for nuc_eos
  real*8, dimension(1) :: xrho,xye,xtemp,xenr,xent
  real*8, dimension(1) :: xprs,xmunu,xcs2
  real*8, dimension(1) :: xdedt,xdpderho,xdpdrhoe
  ! temporary vars for coldeos + gamma law
  integer :: ir
  real*8 :: press_cold, press_th
  real*8 :: eps_cold, eps_th
  real*8 :: gamma
  integer :: num_threads, my_thread_num, slice_len, slice_min, slice_max
  
  CCTK_INT :: my_anyerr

  anyerr    = 0
  keyerr(:) = 0

  select case (eoskey)
     case (1)
        ! polytropic EOS
        if(keytemp.eq.1) then
           !$OMP PARALLEL DO
           do i=1,npoints
              eps(i) = press_gf * poly_k_cgs * &
                   (rho(i)*inv_rho_gf)**(poly_gamma) / &
                   (poly_gamma - 1.0d0) / rho(i)
           enddo
           !$OMP END PARALLEL DO
        endif
        !$OMP PARALLEL DO
        do i=1,npoints
           press(i) = press_gf * poly_k_cgs * &
                (rho(i)*inv_rho_gf)**poly_gamma
        enddo
        !$OMP END PARALLEL DO
     case (2)
        ! gamma-law EOS
        if(keytemp.eq.1) then
           !$OMP PARALLEL DO
           do i=1,npoints
              eps(i) = press_gf * gl_k_cgs * &
                   (rho(i)*inv_rho_gf)**(gl_gamma) / &
                   (gl_gamma - 1.0d0) / rho(i)
           enddo
           !$OMP END PARALLEL DO
        endif
        !$OMP PARALLEL DO
        do i=1,npoints
           press(i) = (gl_gamma - 1.0d0) * rho(i) * eps(i)
        enddo
        !$OMP END PARALLEL DO

     case (3)
        ! hybrid EOS
        !$OMP PARALLEL DO PRIVATE(hybrideos_local_gamma,hybrideos_local_k_cgs, hybrideos_local_k, &
        !$OMP                     hybrideos_local_eps, hybrideos_p_poly, hybrideos_p_th)
        do i=1,npoints
           hybrideos_local_gamma = hybrideos_gamma(1)
           hybrideos_local_k = hybrideos_k(1)
           hybrideos_local_eps = hybrideos_eps(1)
           if (n_pieces .gt. 1) then
              do p = 1,n_pieces-1
                 if (rho(i) .gt. hybrideos_rho(p)) then
                    hybrideos_local_gamma = hybrideos_gamma(p+1)
                    hybrideos_local_k = hybrideos_k(p+1)
                    hybrideos_local_eps = hybrideos_eps(p+1)
                 end if
              end do
           end if
           if (keytemp .eq. 1) then
              !Set epsilon to the cold part only
              eps(i)=hybrideos_local_k/(hybrideos_local_gamma - 1.0)* &
                   rho(i)**(hybrideos_local_gamma - 1.0) + hybrideos_local_eps
           end if
           hybrideos_p_poly = hybrideos_local_k * (rho(i))**hybrideos_local_gamma
           hybrideos_p_th = (hybrid_gamma_th-1)* (rho(i)*eps(i) - &
                hybrideos_p_poly/(hybrideos_local_gamma-1) - &
                hybrideos_local_eps * rho(i))
           hybrideos_p_th = max(zero, hybrideos_p_th)
           press(i) = hybrideos_p_poly + hybrideos_p_th
        enddo
        !$OMP END PARALLEL DO

     case (4)
        ! nuc eos
        my_anyerr = 0
        !$OMP PARALLEL REDUCTION(+: my_anyerr)
#ifdef _OPENMP
        num_threads = omp_get_num_threads()
        my_thread_num = omp_get_thread_num()
#else
        num_threads = 1
        my_thread_num = 0
#endif
        slice_len = (npoints + num_threads-1)/num_threads
        slice_min = 1 + my_thread_num * slice_len
        slice_max = min(slice_min + slice_len - 1, npoints)
        if(keytemp.eq.1) then
           call nuc_eos_m_kt1_press_eps(slice_max-slice_min+1,&
                rho(slice_min:slice_max),temp(slice_min:slice_max),&
                ye(slice_min:slice_max),eps(slice_min:slice_max),&
                press(slice_min:slice_max),keyerr(slice_min:slice_max),&
                my_anyerr)
        else
           call nuc_eos_m_kt0_press(slice_max-slice_min+1,&
                rho(slice_min:slice_max),temp(slice_min:slice_max),&
                ye(slice_min:slice_max),eps(slice_min:slice_max),&
                press(slice_min:slice_max),&
                rf_precision,keyerr(slice_min:slice_max),anyerr)
        endif
        !$OMP END PARALLEL
        ! return 0/1 for false/true rather than zero/nonzero just in case a
        ! caller relied on this
        if(my_anyerr.ne.0) then
           anyerr = 1
        end if
     case (5)
        ! cold tabular EOS with gamma law
        !$OMP PARALLEL DO PRIVATE(xrho,ir,gamma,eps_cold,eps_th, &
        !$OMP                     press_cold,press_th)
        do i=1,npoints
           if(rho(i).lt.coldeos_rhomin) then
              press(i) = coldeos_low_kappa * rho(i)**coldeos_low_gamma
              eps(i) = press(i) / (coldeos_low_gamma - 1.0) / rho(i) 
              cycle
           else if(rho(i).gt.coldeos_rhomax) then
              keyerr(i) = 103
              !$OMP CRITICAL
              anyerr = 1
              !$OMP END CRITICAL
           else
              xrho(1) = log10(rho(i))
              ir = 2 + INT( (xrho(1) - coldeos_logrho(1) - 1.0d-10) * coldeos_dlrhoi )
           endif
           ir = max(2, min(ir,coldeos_nrho))

           gamma = (coldeos_gamma(ir) - coldeos_gamma(ir-1)) / &
                (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                (xrho(1) - coldeos_logrho(ir-1)) + coldeos_gamma(ir-1)

           eps_cold = (coldeos_eps(ir) - coldeos_eps(ir-1)) / &
                    (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                    (xrho(1) - coldeos_logrho(ir-1)) + coldeos_eps(ir-1)

           if(keytemp.eq.0) then
              eps_th = max(0.0d0,eps(i) - eps_cold)
           else if (keytemp.eq.1) then
              eps_th = 0.0d0
              eps(i) = eps_cold
           endif

           press_cold = coldeos_kappa * rho(i)**gamma
           press_th = coldeos_thfac*(coldeos_gammath - 1.0d0)*rho(i)*eps_th
           press(i) = press_cold + press_th
        enddo
        !$OMP END PARALLEL DO

     case (6)
        ! barotropic tabular EOS with gamma law
        !$OMP PARALLEL DO PRIVATE(xrho,xprs,xenr,xye,ir, &
        !$OMP press_cold,eps_cold,eps_th,press_th)
        do i=1,npoints
           if(rho(i).lt.barotropiceos_rhomin) then
              press(i) = barotropiceos_low_kappa * rho(i)**barotropiceos_low_gamma
              eps(i) = press(i) / (barotropiceos_low_gamma - 1.0) / rho(i) 
              cycle
           else if(rho(i).gt.barotropiceos_rhomax) then
              keyerr(i) = 103
              anyerr = 1
           else
              xrho(1) = log10(rho(i))
              ir = 2 + &
                   INT( (xrho(1) - barotropiceos_logrho(1) - 1.0d-10) &
                   * barotropiceos_dlrhoi )
           endif
           ir = max(2, min(ir,barotropiceos_nrho))

           xprs(1) = (barotropiceos_logpress(ir) - barotropiceos_logpress(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho(1) - barotropiceos_logrho(ir-1)) + barotropiceos_logpress(ir-1)

           xenr(1) = (barotropiceos_logeps(ir) - barotropiceos_logeps(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho(1) - barotropiceos_logrho(ir-1)) + barotropiceos_logeps(ir-1)

           xtemp(1) = (barotropiceos_temp(ir) - barotropiceos_temp(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho(1) - barotropiceos_logrho(ir-1)) + barotropiceos_temp(ir-1)

           xye(1) = (barotropiceos_ye(ir) - barotropiceos_ye(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho(1) - barotropiceos_logrho(ir-1)) + barotropiceos_ye(ir-1)

           press_cold = 10.0**xprs(1)
           eps_cold = 10.0**xenr(1) - barotropiceos_energyshift

           if(keytemp.eq.0) then
              eps_th = max(0.0d0,eps(i) - eps_cold)
           else if (keytemp.eq.1) then
              eps_th = 0.0d0
              eps(i) = eps_cold
           endif
           press_th = coldeos_thfac*(barotropiceos_gammath - 1.0d0)*rho(i)*eps_th
           press(i) = press_cold + press_th

        enddo
        !$OMP END PARALLEL DO


     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
     end select

end subroutine EOS_Omni_EOS_PressOMP


subroutine EOS_Omni_EOS_DPressByDRho(eoskey,keytemp,rf_precision,npoints,&
                              rho,eps,temp,ye,dpdrhoe,keyerr,anyerr)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
  CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
  CCTK_REAL, intent(out)   :: dpdrhoe(npoints)

  ! local vars
  integer          :: i,p
  character(256)   :: warnstring
  real*8           :: hybrideos_local_gamma, hybrideos_local_k_cgs, hybrideos_local_k, hybrideos_local_eps, &
                      hybrideos_dp_poly, hybrideos_dp_th1, hybrideos_dp_th2, hybrideos_dp_th
  real*8,parameter :: zero = 0.0d0
  ! temporary vars for nuc_eos
  real*8, dimension(1) :: xrho,xye,xtemp,xenr,xent
  real*8, dimension(1) :: xprs,xmunu,xcs2
  real*8, dimension(1) :: xdedt,xdpderho,xdpdrhoe
  ! temporary vars for cold tabulated EOS + gamma law
  integer :: ir
  real*8 :: press_cold, press_th
  real*8 :: eps_cold, eps_th
  real*8 :: gamma, cs2, h

  anyerr    = 0
  keyerr(:) = 0

  select case (eoskey)
     case (1)
        ! polytropic EOS
        if(keytemp.eq.1) then
           do i=1,npoints
              eps(i) = press_gf * poly_k_cgs * &
                   (rho(i)*inv_rho_gf)**(poly_gamma) / &
                   (poly_gamma - 1.0d0) / rho(i)
           enddo
        endif
        do i=1,npoints
           dpdrhoe(i) = press_gf * poly_k_cgs *  &
                poly_gamma * inv_rho_gf *        & 
                (rho(i)*inv_rho_gf) ** (poly_gamma - 1.d0)
        enddo
     case (2)
        ! gamma-law EOS
        if(keytemp.eq.1) then
           do i=1,npoints
              eps(i) = press_gf * gl_k_cgs * &
                   (rho(i)*inv_rho_gf)**(gl_gamma) / &
                   (gl_gamma - 1.0d0) / rho(i)
           enddo
        endif
        do i=1,npoints
           dpdrhoe(i) = (gl_gamma-1.0d0) * &
                eps(i)
        enddo
     case (3)
        ! hybrid EOS
        do i=1,npoints
           hybrideos_local_gamma = hybrideos_gamma(1)
           hybrideos_local_k = hybrideos_k(1)
           hybrideos_local_eps = hybrideos_eps(1)
           if (n_pieces .gt. 1) then
              do p = 1,n_pieces-1
                 if (rho(i) .gt. hybrideos_rho(p)) then
                    hybrideos_local_gamma = hybrideos_gamma(p+1)
                    hybrideos_local_k = hybrideos_k(p+1)
                    hybrideos_local_eps = hybrideos_eps(p+1)
                 end if
              end do
           end if
           if (keytemp .eq. 1) then
              !Set epsilon to the cold part only
              eps(i)=hybrideos_local_k/(hybrideos_local_gamma - 1.0)* &
                   rho(i)**(hybrideos_local_gamma - 1.0) + hybrideos_local_eps
           end if
           hybrideos_dp_poly = hybrideos_local_gamma * hybrideos_local_k * &
                rho(i)**(hybrideos_local_gamma - 1.0d0)
           
           hybrideos_dp_th = (hybrid_gamma_th - 1.0d0) * (eps(i) - &
                hybrideos_dp_poly/(hybrideos_local_gamma - 1.0d0) - &
                hybrideos_local_eps)

           dpdrhoe(i) = hybrideos_dp_poly + max(0.0d0,hybrideos_dp_th)

        enddo

     case (4)
        if(keytemp.eq.1) then
           call nuc_eos_m_kt1_short(npoints,rho,temp,ye,&
                eps,xprs,xent,xcs2,xdedt,xdpderho,dpdrhoe,xmunu,&
                keyerr,anyerr)
        else if(keytemp.eq.0) then
           call nuc_eos_m_kt0_short(npoints,rho,temp,ye,&
                eps,xprs,xent,xcs2,xdedt,xdpderho,dpdrhoe,xmunu,rf_precision,&
                keyerr,anyerr)
        else
           call CCTK_ERROR("This keytemp is not supported!")
           STOP
        endif
        
     case (5)
        ! with the cold eos we have to assume P = P(rho), so
        ! by definition dPdrho is at constant internal energy
        ! and entropy (the latter, because T=0)
        do i=1,npoints
           if(rho(i).lt.coldeos_rhomin) then
              dpdrhoe(i) = coldeos_low_kappa * coldeos_low_gamma * &
                   rho(i)**(coldeos_low_gamma-1.0d0)
              cycle
           else if(rho(i).gt.coldeos_rhomax) then
              keyerr(i) = 103
              anyerr = 1
           else
              xrho(1) = log10(rho(i))
              ir = 2 + INT( (xrho(1) - coldeos_logrho(1) - 1.0d-10) * coldeos_dlrhoi )
           endif
           ir = max(2, min(ir,coldeos_nrho))

           gamma = (coldeos_gamma(ir) - coldeos_gamma(ir-1)) / &
                (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                (xrho(1) - coldeos_logrho(ir-1)) + coldeos_gamma(ir-1)
           
           ! this is the cold speed of sound squared
           cs2 = (coldeos_cs2(ir) - coldeos_cs2(ir-1)) / &
                (coldeos_cs2(ir) - coldeos_cs2(ir-1)) * &
                (xrho(1) - coldeos_logrho(ir-1)) + coldeos_cs2(ir-1)

           eps_cold = (coldeos_eps(ir) - coldeos_eps(ir-1)) / &
                    (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                    (xrho(1) - coldeos_logrho(ir-1)) + coldeos_eps(ir-1)

           if(keytemp.eq.0) then
              eps_th = max(0.0d0,eps(i) - eps_cold)
           else if (keytemp.eq.1) then
              eps_th = 0.0d0
              eps(i) = eps_cold
           endif

           press_cold = coldeos_kappa * rho(i)**gamma
           press_th = coldeos_thfac*(coldeos_gammath - 1.0d0)*rho(i)*eps_th
           
           ! this is the cold enthalpy, because it belongs to the
           ! cold speed of sound squared
           h = 1.0d0 + eps_cold + press_cold / rho(i)
           dpdrhoe(i) = cs2*h + press_th / rho(i)

        enddo

     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
   end select
     
end subroutine EOS_Omni_EOS_DPressByDRho

subroutine EOS_Omni_EOS_DPressByDEps(eoskey,keytemp,rf_precision,npoints,&
                              rho,eps,temp,ye,dpdepsrho,keyerr,anyerr)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
  CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
  CCTK_REAL, intent(out)   :: dpdepsrho(npoints)

  ! local vars
  integer          :: i,p
  character(256)   :: warnstring
  ! temporary vars for nuc_eos
  real*8, dimension(1) :: xrho,xye,xtemp,xenr,xent
  real*8, dimension(1) :: xprs,xmunu,xcs2
  real*8, dimension(1) :: xdedt,xdpderho,xdpdrhoe
  ! temporary vars for coldeos + gamma law                                                              
  integer :: ir
  real*8 :: eps_cold, eps_th

  anyerr    = 0
  keyerr(:) = 0

  select case (eoskey)
     case (1)
        ! polytropic EOS
        if(keytemp.eq.1) then
           do i=1,npoints
              eps(i) = press_gf * poly_k_cgs * &
                   (rho(i)*inv_rho_gf)**(poly_gamma) / &
                   (poly_gamma - 1.0d0) / rho(i)
           enddo
        endif
        do i=1,npoints
           dpdepsrho(i) = 0.0d0
        enddo
     case (2)
        ! gamma-law EOS
        if(keytemp.eq.1) then
           do i=1,npoints
              eps(i) = press_gf * gl_k_cgs * &
                   (rho(i)*inv_rho_gf)**(gl_gamma) / &
                   (gl_gamma - 1.0d0) / rho(i)
           enddo
        endif
        do i=1,npoints
           dpdepsrho(i) = (gl_gamma - 1.0d0) * &
                rho(i)
        enddo
     case (3)
        ! hybrid EOS
        do i=1,npoints
           dpdepsrho(i) = (hybrid_gamma_th - 1.0d0) * rho(i)
        enddo
     case (4)
        if(keytemp.eq.1) then
           call nuc_eos_m_kt1_short(npoints,rho,temp,ye,&
                eps,xprs,xent,xcs2,xdedt,dpdepsrho,xdpdrhoe,xmunu,&
                keyerr,anyerr)
        else if(keytemp.eq.0) then
           call nuc_eos_m_kt0_short(npoints,rho,temp,ye,&
                eps,xprs,xent,xcs2,xdedt,dpdepsrho,xdpdrhoe,xmunu,rf_precision,&
                keyerr,anyerr)
        else 
           call CCTK_ERROR("This keytemp is not supported!")
           STOP
        endif
     case (5)
        ! with the cold eos we have to assume P = P(rho), so
        ! only the gamma law has non-zero dPdeps
        do i=1,npoints
           if(rho(i).lt.coldeos_rhomin) then
              dpdepsrho(i) = 0.0d0
              cycle
           else if(rho(i).gt.coldeos_rhomax) then
              keyerr(i) = 103
              anyerr = 1
           else
              xrho(1) = log10(rho(i))
              ir = 2 + INT( (xrho(1) - coldeos_logrho(1) - 1.0d-10) * coldeos_dlrhoi )
           endif
           ir = max(2, min(ir,coldeos_nrho))

           eps_cold = (coldeos_eps(ir) - coldeos_eps(ir-1)) / &
                    (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                    (xrho(1) - coldeos_logrho(ir-1)) + coldeos_eps(ir-1)

           if(keytemp.eq.0) then
              eps_th = max(0.0d0,eps(i) - eps_cold)
           else if (keytemp.eq.1) then
              eps_th = 0.0d0
              eps(i) = eps_cold
           endif

           dpdepsrho(i) = coldeos_thfac * (coldeos_gammath - 1.0d0)*rho(i)

        enddo
     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
   end select
     
end subroutine EOS_Omni_EOS_DPressByDEps

subroutine EOS_Omni_EOS_cs2(eoskey,keytemp,rf_precision,npoints,&
                            rho,eps,temp,ye,cs2,keyerr,anyerr)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints)
  CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
  CCTK_REAL, intent(out)   :: cs2(npoints)

  ! local vars
  integer          :: i,p
  character(256)   :: warnstring
  real*8, dimension(1) :: xpress,xdpdrhoe,xdpderho
  real*8           :: hybrideos_local_gamma, hybrideos_local_k, hybrideos_local_eps, &
                      hybrideos_p_poly, hybrideos_p_th
  real*8,parameter :: zero = 0.0d0
  ! temporary vars for nuc_eos
  real*8, dimension(1) :: xrho,xye,xtemp,xenr,xent
  real*8, dimension(1) :: xprs,xmunu,xcs2
  real*8, dimension(1) :: xdedt
  ! temporary vars for cold tabulated EOS + gamma law                                                   
  integer :: ir
  real*8 :: press_cold, press_th
  real*8 :: eps_cold, eps_th
  real*8 :: gamma, cs2_cold, cs2_th, h, h_cold

  anyerr    = 0
  keyerr(:) = 0

  select case (eoskey)
     case (1)
        ! polytropic EOS
        if(keytemp.eq.1) then
           do i=1,npoints
              eps(i) = press_gf * poly_k_cgs * &
                   (rho(i)*inv_rho_gf)**(poly_gamma) / &
                   (poly_gamma - 1.0d0) / rho(i)
           enddo
        endif
        do i=1,npoints
           xpress = press_gf*poly_k_cgs * &
                (rho(i)*inv_rho_gf)**(poly_gamma)
           cs2(i) = poly_gamma * xpress(1) / rho(i) / &
                (1 + eps(i) + xpress(1)/rho(i))
        enddo
     case (2)
        ! gamma-law EOS
        if(keytemp.eq.1) then
           do i=1,npoints
              eps(i) = press_gf * gl_k_cgs * &
                   (rho(i)*inv_rho_gf)**(gl_gamma) / &
                   (gl_gamma - 1.0d0) / rho(i)
           enddo
        endif
        do i=1,npoints
           xpress(1) = (gl_gamma-1.0d0)*rho(i)*eps(i)
           xdpdrhoe(1) = (gl_gamma-1.0d0)*eps(i)
           xdpderho(1) = (gl_gamma-1.0d0)*rho(i)
           cs2(i) = (xdpdrhoe(1) + xpress(1) * xdpderho(1) / (rho(i)**2)) / &
                (1.0d0 + eps(i) + xpress(1)/rho(i))
        enddo
     case(3)
        ! hybrid EOS
        do i=1,npoints
           hybrideos_local_gamma = hybrideos_gamma(1)
           hybrideos_local_k = hybrideos_k(1)
           hybrideos_local_eps = hybrideos_eps(1)
           if (n_pieces .gt. 1) then
              do p = 1,n_pieces-1
                 if (rho(i) .gt. hybrideos_rho(p)) then
                    hybrideos_local_gamma = hybrideos_gamma(p+1)
                    hybrideos_local_k = hybrideos_k(p+1)
                    hybrideos_local_eps = hybrideos_eps(p+1)
                 end if
              end do
           end if
           if (keytemp .eq. 1) then
              !Set epsilon to the cold part only
              eps(i)=hybrideos_local_k/(hybrideos_local_gamma - 1.0)* &
                   rho(i)**(hybrideos_local_gamma - 1.0) + hybrideos_local_eps
           end if
           hybrideos_p_poly = hybrideos_local_k * (rho(i))**hybrideos_local_gamma
           hybrideos_p_th = (hybrid_gamma_th-1)* (rho(i)*eps(i) - &
                hybrideos_p_poly/(hybrideos_local_gamma-1) - &
                hybrideos_local_eps * rho(i))
           hybrideos_p_th = max(zero, hybrideos_p_th)
           xpress(1) = hybrideos_p_poly + hybrideos_p_th
           cs2(i) = (hybrideos_local_gamma * hybrideos_p_poly + &
                     hybrid_gamma_th * hybrideos_p_th) / rho(i) / &
                    (1.0d0 + eps(i) + xpress(1)/rho(i))
        enddo
     case(4)
        ! nuc_eos
        if(keytemp.eq.1) then
           call nuc_eos_m_kt1_press_eps_cs2(npoints,rho,temp,ye,&
                eps,xprs,cs2,keyerr,anyerr)
        else if(keytemp.eq.0) then
           call nuc_eos_m_kt0_press_cs2(npoints,rho,temp,ye,&
                eps,xprs,cs2,rf_precision,keyerr,anyerr)
        else
           call CCTK_ERROR("This keytemp is not supported!")
           STOP
        endif
     case (5)
        ! with the cold eos we have to assume P = P(rho), so
        ! by definition dPdrho is at constant internal energy
        ! and entropy (the latter, because T=0)
        do i=1,npoints
           if(rho(i).lt.coldeos_rhomin) then
              xprs(1) = coldeos_low_kappa * rho(i)**coldeos_low_gamma
              eps(i) = xprs(1) / (coldeos_low_gamma - 1.0d0) / rho(i)
              cs2(i) = coldeos_low_kappa * coldeos_low_gamma * &
                   rho(i)**(coldeos_low_gamma-1.0d0) / &
                   (1.0d0 + eps(i) + xprs(1))
              cycle
           else if(rho(i).gt.coldeos_rhomax) then
              keyerr(i) = 103
              anyerr = 1
                 else
              xrho(1) = log10(rho(i))
              ir = 2 + INT( (xrho(1) - coldeos_logrho(1) - 1.0d-10) * coldeos_dlrhoi )
           endif
           ir = max(2, min(ir,coldeos_nrho))

           gamma = (coldeos_gamma(ir) - coldeos_gamma(ir-1)) / &
                (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                (xrho(1) - coldeos_logrho(ir-1)) + coldeos_gamma(ir-1)

           ! this is the cold speed of sound squared                                                    
           cs2_cold = (coldeos_cs2(ir) - coldeos_cs2(ir-1)) / &
                (coldeos_cs2(ir) - coldeos_cs2(ir-1)) * &
                (xrho(1) - coldeos_logrho(ir-1)) + coldeos_cs2(ir-1)

           eps_cold = (coldeos_eps(ir) - coldeos_eps(ir-1)) / &
                    (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                    (xrho(1) - coldeos_logrho(ir-1)) + coldeos_eps(ir-1)

           if(keytemp.eq.0) then
              eps_th = max(0.0d0,eps(i) - eps_cold)
           else if (keytemp.eq.1) then
              eps_th = 0.0d0
              eps(i) = eps_cold
           endif

           press_cold = coldeos_kappa * rho(i)**gamma
           press_th = coldeos_thfac*(coldeos_gammath - 1.0d0)*rho(i)*eps_th

           xdpdrhoe(1) = coldeos_thfac*(coldeos_gammath - 1.0d0)*eps_th
           xdpderho(1) = coldeos_thfac*(coldeos_gammath - 1.0d0)*rho(i)
           cs2_th = (xdpdrhoe(1) + press_th * xdpderho(1) / (rho(i)**2))

           h = 1.0d0 + eps(i) + (press_cold+press_th) / rho(i)
           h_cold = 1.0d0 + eps_cold + press_cold / rho(i)
           
           cs2(i) = (cs2_cold * h_cold + cs2_th) / h
        enddo
           
     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
   end select
     
end subroutine EOS_Omni_EOS_cs2


subroutine EOS_Omni_EOS_eps_from_press(eoskey,keytemp,rf_precision,npoints,&
                              rho,eps,temp,ye,press,xeps,keyerr,anyerr)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: rho(npoints),ye(npoints),press(npoints)
  CCTK_REAL, intent(inout) :: eps(npoints), temp(npoints)
  CCTK_REAL, intent(out)   :: xeps(npoints)

  ! local vars
  integer        :: i,p
  character(256) :: warnstring
  ! temporary vars for nuc_eos
  real*8, dimension(1) :: xrho,xye,xtemp,xenr,xent
  real*8, dimension(1) :: xprs,xmunu,xcs2
  real*8, dimension(1) :: xdedt,xdpderho,xdpdrhoe

  if(keytemp.eq.1) then
     anyerr = 1
     keyerr(:) = -1
  else
     anyerr    = 0
     keyerr(:) = 0
  endif

  select case (eoskey)
     case (1)
        ! polytropic EOS
        do i=1,npoints
           xeps(i) = press(i) / (poly_gamma - 1.0d0) / rho(i)
        enddo
     case (2)
        ! gamma-law EOS
        do i=1,npoints
           xeps(i) = press(i) / (gl_gamma - 1.0d0) / rho(i)
        enddo
     case (3)
        ! hybrid EOS
        write(warnstring,*) "EOS_Omni_EpsFromPress call not supported for hybrid EOS"
        call CCTK_ERROR(warnstring)
        STOP
     case (4)
        ! nuc EOS
        write(warnstring,*) "EOS_Omni_EpsFromPress call not supported for nuc_eos yet"
        call CCTK_ERROR(warnstring)
        STOP
     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
     end select

end subroutine EOS_Omni_EOS_eps_from_press


subroutine EOS_Omni_EOS_RhoFromPressEpsTempEnt(eoskey,&
                           ikeytemp,rf_precision,&
                           npoints,rho,eps,temp,entropy,ye,press,keyerr,anyerr)

  ! This routine returns the density and spec. internal energy based on inputs of
  ! pressure, temperature, and Y_e.
  ! The current implementation is robust, but very slow; it should be used only
  ! for initial data where speed does not matter.

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(in)     :: eoskey,ikeytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: ye(npoints),press(npoints)
  CCTK_REAL, intent(inout) :: rho(npoints),temp(npoints),eps(npoints)
  CCTK_REAL, intent(inout) :: entropy(npoints)


  ! local vars
  integer        :: i,p
  character(256) :: warnstring
  integer :: keytemp
  ! temporary vars for nuc_eos
  real*8, dimension(1) :: xrho,xye,xtemp,xenr,xent
  real*8, dimension(1) :: xprs,xmunu,xcs2
  real*8, dimension(1) :: xdedt,xdpderho,xdpdrhoe
  ! temporary vars for cold tabulated EOS
  integer :: ir
  real*8  :: gamma
  
  ! helper vars
  logical :: cont
  integer :: counter
  real*8, dimension(1) :: rho_guess, rho_guess2
  real*8, dimension(1) :: press_guess,press_guess2, mydpdrho
  real*8  :: fac
  keytemp = ikeytemp

  anyerr    = 0
  keyerr(:) = 0

  select case (eoskey)
     case (1)
        ! polytropic EOS
        do i=1,npoints
           rho(i) = (press(i) / poly_k)**(1.0/poly_gamma)
           eps(i) = press(i) / (poly_gamma - 1.0d0) / rho(i)
        enddo
     case (2)
        ! gamma-law EOS
        do i=1,npoints
           rho(i) = press(i) / ((gl_gamma - 1.0d0)*eps(i))
        enddo
     case (3)
        ! hybrid EOS
        write(warnstring,*) "EOS_Omni_RestMassDensityFromPressEpsTemp not supported for hybrid EOS"
        call CCTK_ERROR(warnstring)
        STOP
     case (4)
        ! nuc EOS
        if(ikeytemp.eq.2) then
           call CCTK_ERROR("This function does not work yet when coming in with entropy")
           STOP
        else if(ikeytemp.eq.1) then
           keytemp = 1
        else
           call CCTK_ERROR("This function does not work yet when coming in with this keytemp")
           STOP
        endif
        
        call CCTK_ERROR("This routine does not work with the new nuc_eos_cxx [yet]")
        STOP

        keytemp = 1
        do i=1,npoints
           fac = 1.0d0
           counter = 0
           cont = .true.
           xprs(1) = press(i)
           press_guess(1) = xprs(1)
           rho_guess(1) = rho(i)
           xprs(1) = press(i)
           xtemp(1) = temp(i)
           xye(1) = ye(i)
           do while(cont)
              counter = counter + 1
              rho_guess2 = rho_guess * 1.0001d0
              call nuc_eos_m_kt1_short(1,rho_guess2,xtemp,xye,xenr,&
                   press_guess2,xent,xcs2,xdedt,xdpderho,xdpdrhoe,xmunu,&
                   keyerr(i),anyerr)
              call nuc_eos_m_kt1_short(1,rho_guess2,xtemp,xye,xenr,&
                   press_guess2,xent,xcs2,xdedt,xdpderho,xdpdrhoe,xmunu,&
                   keyerr(i),anyerr)
              mydpdrho(1) = (press_guess2(1)-press_guess(1))/(rho_guess2(1)-rho_guess(1))
              if (mydpdrho(1).lt.0.0d0) then
                 !$OMP CRITICAL
                 write(warnstring,"(A25,1P10E15.6)") "Issue with table, dpdrho.lt.0",&
                      rho_guess,xtemp,xye
                 call CCTK_ERROR(warnstring)
                 !$OMP END CRITICAL
                 STOP
              endif
              
              if (abs(1.0d0-press_guess(1)/xprs(1)).lt.rf_precision) then
                 cont = .false.
                 rho(i) = rho_guess(1)
                 eps(i) = xenr(1)
              else
                 if (fac*(xprs(1)-press_guess(1))/mydpdrho(1)/rho_guess(1).gt.0.1d0) then
                    rho_guess = 0.99d0*rho_guess
                 else
                    rho_guess(1) = rho_guess(1) + fac*(xprs(1)-press_guess(1))/mydpdrho(1)
                 endif
              endif

              if (counter.gt.100) then
                 fac = 0.01d0
              endif

              if (rho_guess(1).lt.1.0d3.or.counter.gt.100000) then
                 keyerr(i) = 473
                 anyerr = 1
                 return
              endif
           enddo
        enddo
        
     case (5)
        ! cold tabulated EOS
        ! deal only with case in which thermal pressure is zero
        if(keytemp.ne.1) then
           call CCTK_ERROR("finding rho(press) for tabulated cold EOS only possible if keytemp=1")
           STOP
        endif

        do i=1,npoints
           fac = 1.0d0
           counter = 0
           cont = .true.
           xprs(1) = press(i)
           press_guess(1) = xprs(1)
           rho_guess(1) = rho(i)
           xprs(1) = press(i)
           do while(cont)
              counter = counter + 1
              rho_guess2(1) = rho_guess(1) * 1.0001d0

              ! press 2
              if(rho_guess2(1).lt.coldeos_rhomin) then
                 keyerr(i) = 104
                 anyerr = 1
              else if(rho_guess2(1).gt.coldeos_rhomax) then
                 keyerr(i) = 103
                 anyerr = 1
              else
                 xrho(1) = log10(rho_guess2(1))
                 ir = 2 + INT( (xrho(1) - coldeos_logrho(1) - 1.0d-10) * coldeos_dlrhoi )
              endif
              ir = max(2, min(ir,coldeos_nrho))
           
              gamma = (coldeos_gamma(ir) - coldeos_gamma(ir-1)) / &
                (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                (xrho(1) - coldeos_logrho(ir-1)) + coldeos_gamma(ir-1)

              press_guess2(1) = coldeos_kappa * rho_guess2(1)**gamma
              ! end press 2
              ! press 1
              if(rho_guess(1).lt.coldeos_rhomin) then
                 keyerr(i) = 104
                 anyerr = 1
              else if(rho_guess(1).gt.coldeos_rhomax) then
                 keyerr(i) = 103
                 anyerr = 1
              else
                 xrho(1) = log10(rho_guess(1))
                 ir = 2 + INT( (xrho(1) - coldeos_logrho(1) - 1.0d-10) * coldeos_dlrhoi )
              endif
              ir = max(2, min(ir,coldeos_nrho))
           
              gamma = (coldeos_gamma(ir) - coldeos_gamma(ir-1)) / &
                (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                (xrho(1) - coldeos_logrho(ir-1)) + coldeos_gamma(ir-1)

              press_guess(1) = coldeos_kappa * rho_guess(1)**gamma
              ! end press 1

              ! derivative
              mydpdrho(1) = (press_guess2(1)-press_guess(1))/(rho_guess2(1)-rho_guess(1))
              if (mydpdrho(1).lt.0.0d0) then
                 !$OMP CRITICAL
                 write(warnstring,"(A25,1P10E15.6)") "Issue with table, dpdrho.lt.0",&
                      rho_guess(1),rho_guess2(1)
                 call CCTK_ERROR(warnstring)
                 !$OMP END CRITICAL
                 STOP
              endif
              
              if (abs(1.0d0-press_guess(1)/xprs(1)).lt.rf_precision) then
                 cont = .false.
                 rho(i) = rho_guess(1)
                 eps(i) = (coldeos_eps(ir) - coldeos_eps(ir-1)) / &
                    (coldeos_logrho(ir) - coldeos_logrho(ir-1)) * &
                    (xrho(1) - coldeos_logrho(ir-1)) + coldeos_eps(ir-1)
              else
                 if (fac*(xprs(1)-press_guess(1))/mydpdrho(1)/rho_guess(1).gt.0.1d0) then
                    rho_guess(1) = 0.99d0*rho_guess(1)
                 else
                    rho_guess(1) = rho_guess(1) + fac*(xprs(1)-press_guess(1))/mydpdrho(1)
                 endif
              endif

              if (counter.gt.100) then
                 fac = 0.01d0
              endif

              if (rho_guess(1).lt.1.0d-15.or.counter.gt.100000) then
                 keyerr(i) = 473
                 anyerr = 1
                 return
              endif
           enddo

        enddo

     case (6)
        ! barotropic tabulated EOS
        ! deal only with case in which thermal pressure is zero
        if(keytemp.ne.1) then
           call CCTK_ERROR("finding rho(press) for tabulated barotropic EOS only possible if keytemp=1")
           STOP
        endif

        do i=1,npoints
           fac = 1.0d0
           counter = 0
           cont = .true.
           xprs(1) = press(i)
           press_guess(1) = xprs(1)
           rho_guess(1) = rho(i)
           do while(cont)
              counter = counter + 1
              rho_guess2(1) = rho_guess(1) * 1.0001d0

              ! press 2
              if(rho_guess2(1).lt.barotropiceos_rhomin) then
                 keyerr(i) = 104
                 anyerr = 1
              else if(rho_guess2(1).gt.barotropiceos_rhomax) then
                 keyerr(i) = 103
                 anyerr = 1
              else
                 xrho(1) = log10(rho_guess2(1))
                 ir = 2 + INT( (xrho(1) - barotropiceos_logrho(1) - 1.0d-10) * barotropiceos_dlrhoi )
              endif
              ir = max(2, min(ir,barotropiceos_nrho))
           
              xprs(1) = (barotropiceos_logpress(ir) - barotropiceos_logpress(ir-1)) / &
                   (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                   (xrho(1) - barotropiceos_logrho(ir-1)) + barotropiceos_logpress(ir-1)
              
              press_guess2(1) = 10.0d0**xprs(1)
              ! end press 2
              ! press 1
              if(rho_guess(1).lt.barotropiceos_rhomin) then
                 keyerr(i) = 104
                 anyerr = 1
              else if(rho_guess(1).gt.barotropiceos_rhomax) then
                 keyerr(i) = 103
                 anyerr = 1
              else
                 xrho(1) = log10(rho_guess(1))
                 ir = 2 + INT( (xrho(1) - barotropiceos_logrho(1) - 1.0d-10) * barotropiceos_dlrhoi )
              endif
              ir = max(2, min(ir,barotropiceos_nrho))
           
              xprs(1) = (barotropiceos_logpress(ir) - barotropiceos_logpress(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho(1) - barotropiceos_logrho(ir-1)) + barotropiceos_logpress(ir-1)

              press_guess(1) = 10.0d0**xprs(1)
              ! end press 1

 !             write(6,"(i7,1P10E15.6)") counter, rho_guess, press_guess, press_guess2, press(i)

              ! derivative
              mydpdrho(1) = (press_guess2(1)-press_guess(1))/(rho_guess2(1)-rho_guess(1))
              if (mydpdrho(1).lt.0.0d0) then
                 !$OMP CRITICAL
                 write(warnstring,"(A25,1P10E15.6)") "Issue with table, dpdrho.lt.0",&
                      rho_guess(1),rho_guess2(1)
                 call CCTK_ERROR(warnstring)
                 !$OMP END CRITICAL
                 STOP
              endif
              
              if (abs(1.0d0-press_guess(1)/press(i)).lt.rf_precision) then
                 cont = .false.
                 rho(i) = rho_guess(1)
                 xenr(1) = (barotropiceos_logeps(ir) - barotropiceos_logeps(ir-1)) / &
                      (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                      (xrho(1) - barotropiceos_logrho(ir-1)) + barotropiceos_logeps(ir-1)
                 eps(i) = 10.0d0**xenr(1) - barotropiceos_energyshift
              else
                 if (fac*(press(i)-press_guess(1))/mydpdrho(1)/rho_guess(1).gt.0.1d0) then
                    rho_guess(1) = 0.99d0*rho_guess(1)
                 else
                    rho_guess(1) = rho_guess(1) + fac*(press(i)-press_guess(1))/mydpdrho(1)
                 endif
              endif

              if (counter.gt.100) then
                 fac = 0.01d0
              endif

              if (rho_guess(1).lt.1.0d-15.or.counter.gt.100000) then
!                 !$OMP CRITICAL
!                 write(warnstring,"(A25,1P10E15.6)") "rho(p) issue", rho_guess,press(i),press_guess
!                 call CCTK_ERROR(warnstring)
!                 STOP
!                 !$OMP END CRITICAL
                 keyerr(i) = 473
                 anyerr = 1
                 return
              endif
           enddo

        enddo

     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented!"
        call CCTK_ERROR(warnstring)
        STOP
     end select

end subroutine EOS_Omni_EOS_RhoFromPressEpsTempEnt


subroutine EOS_Omni_EOS_PressEpsTempYe_from_Rho(eoskey,keytemp,rf_precision,npoints,&
                              rho,eps,temp,ye,press,keyerr,anyerr)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, intent(in)     :: eoskey,keytemp,npoints
  CCTK_INT, intent(out)    :: keyerr(npoints)
  CCTK_INT, intent(out)    :: anyerr
  CCTK_REAL, intent(in)    :: rf_precision
  CCTK_REAL, intent(in)    :: rho(npoints)
  CCTK_REAL, intent(out)   :: eps(npoints), ye(npoints) 
  CCTK_REAL, intent(out)   :: temp(npoints),press(npoints)

  ! local vars
  integer          :: i,p
  character(256)   :: warnstring
  ! temporary vars for nuc_eos
  real*8, dimension(1):: xrho,xye,xtemp,xenr
  real*8, dimension(1):: xprs
  integer :: ir
  real*8 :: press_cold, press_th
  real*8 :: eps_cold, eps_th

  anyerr    = 0
  keyerr(:) = 0

  select case (eoskey)
     case (6)
        ! barotropic tabular EOS with gamma law
        do i=1,npoints
           if(rho(i).lt.barotropiceos_rhomin) then
              press(i) = barotropiceos_low_kappa * rho(i)**barotropiceos_low_gamma
              eps(i) = press(i) / (barotropiceos_low_gamma - 1.0) / rho(i) 
              cycle
           else if(rho(i).gt.barotropiceos_rhomax) then
              keyerr(i) = 103
              anyerr = 1
           else
              xrho(1) = log10(rho(i))
              ir = 2 + &
                   INT( (xrho(1) - barotropiceos_logrho(1) - 1.0d-10) &
                   * barotropiceos_dlrhoi )
           endif
           ir = max(2, min(ir,barotropiceos_nrho))

           xprs(1) = (barotropiceos_logpress(ir) - barotropiceos_logpress(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho(1) - barotropiceos_logrho(ir-1)) + barotropiceos_logpress(ir-1)

           xenr(1) = (barotropiceos_logeps(ir) - barotropiceos_logeps(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho(1) - barotropiceos_logrho(ir-1)) + barotropiceos_logeps(ir-1)

           xtemp(1) = (barotropiceos_temp(ir) - barotropiceos_temp(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho(1) - barotropiceos_logrho(ir-1)) + barotropiceos_temp(ir-1)

           xye(1) = (barotropiceos_ye(ir) - barotropiceos_ye(ir-1)) / &
                (barotropiceos_logrho(ir) - barotropiceos_logrho(ir-1)) * &
                (xrho(1) - barotropiceos_logrho(ir-1)) + barotropiceos_ye(ir-1)

           press_cold = 10.0**xprs(1)
           eps_cold = 10.0**xenr(1) - barotropiceos_energyshift

           if(keytemp.eq.0) then
              eps_th = max(0.0d0,eps(i) - eps_cold)
           else if (keytemp.eq.1) then
              eps_th = 0.0d0
              eps(i) = eps_cold
           endif
           press_th = coldeos_thfac*(barotropiceos_gammath - 1.0d0)*rho(i)*eps_th
           press(i) = press_cold + press_th

           temp(i) = xtemp(1)
           ye(i) = xye(1)

        enddo

     case DEFAULT
        write(warnstring,*) "eoskey ",eoskey," not implemented for EOS_Omni_EOS_PressEpsTempYe_from_Rho!"
        call CCTK_ERROR(warnstring)
        STOP
     end select

   end subroutine EOS_Omni_EOS_PressEpsTempYe_from_Rho
