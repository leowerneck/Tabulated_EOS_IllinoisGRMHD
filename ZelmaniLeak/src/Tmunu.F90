#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "SpaceMask.h"
#define PRESS_NU_CONSTANT 3.52127727d24
#define PI 3.14159265358979d0
#define PI4 97.4090910340024d0
#define PI2 9.86960440108936d0
#define MEV_TO_ERG 1.60217733d-6
#define SPATIAL_DETERMINANT(gxx_,gxy_,gxz_,gyy_,gyz_,gzz_) \
  (-(gxz_)**2*(gyy_) + 2*(gxy_)*(gxz_)*(gyz_) - (gxx_)*(gyz_)**2 - (gxy_)**2*(gzz_) \
   + (gxx_)*(gyy_)*(gzz_))


! include neutrino pressure in stress energy tensor
subroutine ZelmaniLeak_Tmunu(CCTK_ARGUMENTS)

  use EOS_Omni_module
#ifdef HAVE_CAPABILITY_Fortran
  use cctk
#endif
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i,j,k
  CCTK_INT :: nx,ny,nz

  ! EOS 
  character(len=512) :: warnline
  integer   :: eoskey,keytemp,n,anyerr,keyerr(1)
  CCTK_REAL :: xdummy(1),xeps(1)

  ! pnu vars:
  CCTK_REAL :: eta
  CCTK_REAL :: pnufac
  CCTK_REAL :: pnuconst, F3const, F3

  ! Tmunu helpers
  CCTK_REAL :: betaxlow,betaylow,betazlow,beta2
  CCTK_REAL :: velxlow,velylow,velzlow,tmunutemp
  CCTK_REAL :: utlow,uxlow,uylow,uzlow
  CCTK_REAL :: enufac

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  if(global_rho_max*inv_rho_gf < pnu_rho_start / 10.0d0) return;
  
  pnuconst = PRESS_NU_CONSTANT * press_gf
  F3const = 7.0d0*PI4/60.0d0

  ! EOS Omni setup:
  eoskey=4;keytemp=1;n=1;anyerr=0;keyerr=0
  xdummy=0.0d0

  if(include_enu_in_tmunu.ne.0) then
     enufac = 1.0d0
  else
     enufac = 0.0d0
  endif
     
  
  !compute neutrino pressure
  !$OMP PARALLEL DO PRIVATE(i, j, k, xdummy, eta, pnufac, F3, xeps,&
  !$OMP                     keyerr,anyerr,velxlow,velylow,velzlow,beta2,&
  !$OMP                     betaxlow,betaylow,betazlow,tmunutemp,&
  !$OMP                     utlow,uxlow,uylow,uzlow)
  do k=1,nz
     do j=1,ny
        do i=1,nx
            keyerr = 0
            anyerr = 0
            xeps = 0.0d0
            call EOS_Omni_short(eoskey,keytemp,igm_eos_root_finding_precision,&
                 n,rho(i,j,k),xeps,temperature(i,j,k),y_e(i,j,k),&
                 xdummy,xdummy,xdummy,xdummy,&
                 xdummy,xdummy,munu(i,j,k),keyerr,anyerr)

            if(anyerr.ne.0) then
               write(warnline,"(3i5,1P10E15.6)") i,j,k,&
                    x(i,j,k),y(i,j,k),z(i,j,k)
               write(warnline,"(3i5,i6)") i,j,k,&
                    keyerr
               write(warnline,"(3i5,1P10E15.6)") i,j,k,&
                    rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k)
               call CCTK_WARN(1,warnline)
               call CCTK_WARN(0,"Problem in Pnu in Tmunu EOS call 1!")
            endif

            ! let's compute the neutrino pressure
            eta = munu(i,j,k)/temperature(i,j,k)
            F3 = F3const + 0.5d0*eta*eta*(PI2 + 0.5d0*eta*eta)
! calculate pnu, attentuate contribution at lower densities
            pnu(i,j,k) = F3 * pnuconst * temperature(i,j,k)**4 * & 
                 exp(- pnu_rho_start*rho_gf / rho(i,j,k))

            velxlow = gxx(i,j,k)*vel(i,j,k,1) + gxy(i,j,k)*vel(i,j,k,2) + gxz(i,j,k)*vel(i,j,k,3)
            velylow = gxy(i,j,k)*vel(i,j,k,1) + gyy(i,j,k)*vel(i,j,k,2) + gyz(i,j,k)*vel(i,j,k,3)
            velzlow = gxz(i,j,k)*vel(i,j,k,1) + gyz(i,j,k)*vel(i,j,k,2) + gzz(i,j,k)*vel(i,j,k,3)
            
            betaxlow = gxx(i,j,k)*betax(i,j,k) + gxy(i,j,k)*betay(i,j,k) + gxz(i,j,k)*betaz(i,j,k)
            betaylow = gxy(i,j,k)*betax(i,j,k) + gyy(i,j,k)*betay(i,j,k) + gyz(i,j,k)*betaz(i,j,k)
            betazlow = gxz(i,j,k)*betax(i,j,k) + gyz(i,j,k)*betay(i,j,k) + gzz(i,j,k)*betaz(i,j,k)

            beta2 = betax(i,j,k)*betaxlow + betay(i,j,k)*betaylow + betaz(i,j,k)*betazlow

!            tmunutemp = w_lorentz(i,j,k)**2 * (pnu(i,j,k))

            tmunutemp = w_lorentz(i,j,k)**2 * (3.0D0*pnu(i,j,k)*enufac + pnu(i,j,k))
#if 0
! old, from 2006; don't want to do this in the future
            tmunutemp = w_lorentz(i,j,k)**2 * (3.0D0*pnu(i,j,k) + pnu(i,j,k))
#endif
            utlow = (-alp(i,j,k) + vel(i,j,k,1)*betaxlow + vel(i,j,k,2)*betaylow + vel(i,j,k,3)*betazlow)
            uxlow = velxlow
            uylow = velylow
            uzlow = velzlow
! Calculate Tmunu (the lower components!)
            eTtt(i,j,k) = eTtt(i,j,k) + tmunutemp*utlow**2 + pnu(i,j,k)*(beta2 - alp(i,j,k)**2)
            
            eTtx(i,j,k) = eTtx(i,j,k) + tmunutemp*utlow*uxlow + pnu(i,j,k)*betaxlow
            eTty(i,j,k) = eTty(i,j,k) + tmunutemp*utlow*uylow + pnu(i,j,k)*betaylow
            eTtz(i,j,k) = eTtz(i,j,k) + tmunutemp*utlow*uzlow + pnu(i,j,k)*betazlow
            
            eTxx(i,j,k) = eTxx(i,j,k) + tmunutemp*uxlow**2 + pnu(i,j,k)*gxx(i,j,k)
            eTyy(i,j,k) = eTyy(i,j,k) + tmunutemp*uylow**2 + pnu(i,j,k)*gyy(i,j,k)
            eTzz(i,j,k) = eTzz(i,j,k) + tmunutemp*uzlow**2 + pnu(i,j,k)*gzz(i,j,k)
            
            eTxy(i,j,k) = eTxy(i,j,k) + tmunutemp*uxlow*uylow + pnu(i,j,k)*gxy(i,j,k)
            eTxz(i,j,k) = eTxz(i,j,k) + tmunutemp*uxlow*uzlow + pnu(i,j,k)*gxz(i,j,k)
            eTyz(i,j,k) = eTyz(i,j,k) + tmunutemp*uylow*uzlow + pnu(i,j,k)*gyz(i,j,k)
            
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine ZelmaniLeak_Tmunu
