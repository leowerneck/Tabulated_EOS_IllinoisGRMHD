#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#define SPATIAL_DETERMINANT(gxx_,gxy_,gxz_,gyy_,gyz_,gzz_) \
  (-(gxz_)**2*(gyy_) + 2*(gxy_)*(gxz_)*(gyz_) - (gxx_)*(gyz_)**2 - (gxy_)**2*(gzz_) \
   + (gxx_)*(gyy_)*(gzz_))



subroutine ZelmaniLeak_ye_of_rho(CCTK_ARGUMENTS)

  use EOS_Omni_module
#ifdef HAVE_CAPABILITY_Fortran
  use cctk
#endif
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: i,j,k
  CCTK_INT :: nx,ny,nz
  CCTK_REAL :: yeofrho_logrho1 
  CCTK_REAL :: yeofrho_logrho2
  CCTK_REAL :: yeofrho_high_logrho

  ! EOS parameters
  integer   :: eoskey 
  integer   :: keytemp 
  integer   :: n 
  integer   :: anyerr 
  CCTK_INT  :: keyerr(1)
  CCTK_REAL :: xdummy,xtemp
  CCTK_REAL :: xrho
  CCTK_REAL :: xye_new, xye, delta_ye, delta_s
  CCTK_REAL :: rdetg, w, tenthalpy, vlowx, vlowy, vlowz
  CCTK_REAL :: sold,sup,epsold
  character(len=512) :: warnline
  integer :: rl

  n = 1
  anyerr = 0
  eoskey = 4

  yeofrho_logrho1 = log10(yeofrho_rho1)
  yeofrho_logrho2 = log10(yeofrho_rho2)
  yeofrho_high_logrho = log10(yeofrho_high_correction_rho)

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  rl = nint(log10(dble(cctk_levfac(1)))/log10(2.0d0))


  !$OMP PARALLEL DO PRIVATE(k, i, j, xdummy, keyerr, anyerr, warnline, keytemp)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           keytemp = 0
           call EOS_Omni_short(eoskey,keytemp,igm_eos_root_finding_precision,&
                n,rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                xdummy,entropy(i,j,k),xdummy,xdummy,&
                xdummy,xdummy,munu(i,j,k),keyerr,anyerr)
           if(anyerr.ne.0.and.rl.ge.-1) then
              call CCTK_WARN(1,"EOS error in ZelmaniLeak::ye_of_rho 1")
              write(warnline,"(4i6)") i,j,k,rl
              call CCTK_WARN(1,warnline)
              write(warnline,"(1P10E15.6)") x(i,j,k),y(i,j,k),z(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,"(1P10E15.6)") rho(i,j,k),eps(i,j,k),y_e(i,j,k),&
                   temperature(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,"(A8,i8)") "keyerr:", keyerr(1)
              call CCTK_WARN(1,warnline)
              call CCTK_WARN(0,"aborting!")
           endif
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  if(bounce.ne.1.and.do_ye_of_rho.ne.0) then
     if(do_ye_of_rho_from_profile.eq.0) then
        !$OMP PARALLEL DO PRIVATE(i, j, k, xrho, xye, xye_new, delta_ye, sold, sup,&
        !$OMP                     keyerr, anyerr, keytemp, delta_s, xdummy, epsold,&
        !$OMP                     w, rdetg, warnline, vlowx, vlowy, vlowz, tenthalpy)
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 
                 ! going to cgs for density
                 xrho = rho(i,j,k) * inv_rho_gf
                 xye = y_e(i,j,k)
                 
                 call fit_ye(xrho,xye_new,yeofrho_logrho2,yeofrho_logrho1,&
                      yeofrho_do_high_correction,&
                      yeofrho_high_correction_ye,yeofrho_high_logrho)
                 delta_ye = min(0.0d0, xye_new-xye)

!                 if(delta_ye < 0.0d0.and.entropy(i,j,k).gt.0.1d0) then
                 if(delta_ye < 0.0d0.and.entropy(i,j,k).gt.0.1d0) then
                    
                    if(munu(i,j,k).gt.10.0d0.and.xrho.lt.2.0d12) then
                       delta_s = -delta_ye*(munu(i,j,k)-10.0d0) / &
                            temperature(i,j,k)
                    else
                       delta_s = 0.0d0
                    endif
                    sold = entropy(i,j,k)
                    entropy(i,j,k) = entropy(i,j,k) + delta_s
                    epsold = eps(i,j,k)
                    sup = entropy(i,j,k)
                    y_e(i,j,k) = y_e(i,j,k) + delta_ye
                    
                    ! now use the entropy to update things
                    keytemp = 2
                    call EOS_Omni_short(eoskey,keytemp,igm_eos_root_finding_precision,&
                         n,rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                         press(i,j,k),entropy(i,j,k),xdummy,xdummy,&
                         xdummy,xdummy,xdummy,keyerr,anyerr)
!                      xdummy,xdummy,munu(i,j,k),keyerr,anyerr)

!# bookkeeping -- we turn this of for performance reasons
!#                 energycheck(i,j,k) = energycheck(i,j,k) + (eps(i,j,k)-epsold)*rho(i,j,k)

                    if(anyerr.ne.0) then
                       write(warnline,"(3i5,1P10E15.6)") i,j,k,&
                            x(i,j,k),y(i,j,k),z(i,j,k)
                       call CCTK_WARN(1,warnline)
                       write(warnline,"(3i5,1P10E15.6)") i,j,k,&
                            entropy(i,j,k),y_e(i,j,k)
                       call CCTK_WARN(1,warnline)
                       call CCTK_WARN(0,"Problem in ye_of_rho EOS call 2!")
                    endif
                 
                 endif

              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
     else
        !$OMP PARALLEL DO PRIVATE(i, j, k, xrho, xye, xye_new, delta_ye, sold, sup,&
        !$OMP                     keyerr, anyerr, keytemp, delta_s, xdummy, epsold,&
        !$OMP                     w, rdetg, warnline, vlowx, vlowy, vlowz, tenthalpy)
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 
                 ! going to cgs for density
                 xrho = rho(i,j,k) * inv_rho_gf
                 xye = y_e(i,j,k)
                 
                 call yeofrho(xrho,xye_new)
                 delta_ye = min(0.0d0, xye_new-xye)

                 if(delta_ye < 0.0d0) then

                    if(munu(i,j,k).gt.10.0d0.and.xrho.lt.2.0d12) then
                       delta_s = -delta_ye*(munu(i,j,k)-10.0d0) / &
                            temperature(i,j,k)
                    else
                       delta_s = 0.0d0
                    endif
                    sold = entropy(i,j,k)
                    entropy(i,j,k) = entropy(i,j,k) + delta_s
                    epsold = eps(i,j,k)
                    sup = entropy(i,j,k)
                    y_e(i,j,k) = y_e(i,j,k) + delta_ye
                    
                    ! now use the entropy to update things
                    keytemp = 2
                    call EOS_Omni_short(eoskey,keytemp,igm_eos_root_finding_precision,&
                         n,rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                         press(i,j,k),entropy(i,j,k),xdummy,xdummy,&
                         xdummy,xdummy,xdummy,keyerr,anyerr)
!                      xdummy,xdummy,munu(i,j,k),keyerr,anyerr)

!# bookkeeping -- we turn this of for performance reasons
!#                 energycheck(i,j,k) = energycheck(i,j,k) + (eps(i,j,k)-epsold)*rho(i,j,k)

                    if(anyerr.ne.0) then
                       write(warnline,"(3i5,1P10E15.6)") i,j,k,&
                            x(i,j,k),y(i,j,k),z(i,j,k)
                       call CCTK_WARN(1,warnline)
                       write(warnline,"(3i5,1P10E15.6)") i,j,k,&
                            entropy(i,j,k),y_e(i,j,k)
                       call CCTK_WARN(1,warnline)
                       call CCTK_WARN(0,"Problem in ye_of_rho EOS call 2!")
                    endif
                 
                 endif

              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
     endif
  endif


     
end subroutine ZelmaniLeak_ye_of_rho

subroutine fit_ye(incomingrho,outgoingye,yeofrho_logrho2,yeofrho_logrho1, &
     do_high_correction,high_correction_ye,high_correction_logrho)

  ! in CGS
#ifdef HAVE_CAPABILITY_Fortran
  use cctk
#endif
  implicit none
  DECLARE_CCTK_PARAMETERS
    
  integer do_high_correction
  real*8 yeofrho_logrho2, yeofrho_logrho1
  real*8 high_correction_ye,high_correction_logrho
  real*8 incomingrho,outgoingye
  real*8 logincomingrho
    
  real*8 x,absx,slope

  if (incomingrho.lt.1.03d0) then
     stop "incoming rho is < 1.03d0 g/ccm "
  end if
  
  logincomingrho = log10(incomingrho)

  ! Leo says: Part I of Eq. (1) in https://arxiv.org/pdf/astro-ph/0504072.pdf
  x = max(-1.0d0,min(1.0d0,(2.0d0*logincomingrho - yeofrho_logrho2 - yeofrho_logrho1) &
       /(yeofrho_logrho2-yeofrho_logrho1)))
  
  absx = abs(x)

  ! Leo says: Part II of Eq. (1) in https://arxiv.org/pdf/astro-ph/0504072.pdf
  outgoingye = 0.5d0*(yeofrho_ye2+yeofrho_ye1) + x/2.0d0*(yeofrho_ye2-yeofrho_ye1) &
       + yeofrho_yec*(1.0d0-absx+4.0d0*absx*(absx-0.5d0)*(absx-1.0d0))
  
  if(do_high_correction.ne.0 & 
       .and.logincomingrho.gt.yeofrho_logrho2) then
     slope = (high_correction_ye-yeofrho_ye2) /  &
          (high_correction_logrho-yeofrho_logrho2)
     outgoingye = yeofrho_ye2+slope*(logincomingrho-yeofrho_logrho2)
  endif

end subroutine fit_ye

subroutine ZelmaniLeak_calc_entropy(CCTK_ARGUMENTS)

! we don't use this  routine any more, but keep it around
! in case we need the funcitonality somewhere else

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
  CCTK_REAL :: yeofrho_logrho1 
  CCTK_REAL :: yeofrho_logrho2

  ! EOS parameters
  integer   :: eoskey 
  integer   :: keytemp 
  integer   :: n 
  integer   :: anyerr 
  CCTK_INT  :: keyerr(1)
  CCTK_REAL :: xdummy(1),xtemp
  CCTK_REAL :: xrho
  CCTK_REAL :: xye_new, xye, delta_ye, delta_s
  CCTK_REAL :: rdetg, w, tenthalpy, vlowx, vlowy, vlowz
  CCTK_REAL :: sold,sup,epsold
  character(len=512) :: warnline
  integer :: rl

  n = 1
  anyerr = 0
  eoskey = 4

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)


  !$OMP PARALLEL DO PRIVATE(k, i, j, xdummy, keyerr, anyerr, warnline, keytemp)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           keytemp = 0
           call EOS_Omni_short(eoskey,keytemp,igm_eos_root_finding_precision,&
                n,rho(i,j,k),eps(i,j,k),temperature(i,j,k),y_e(i,j,k),&
                xdummy,entropy(i,j,k),xdummy,xdummy,&
                xdummy,xdummy,munu(i,j,k),keyerr,anyerr)
           if(anyerr.ne.0) then
              call CCTK_WARN(1,"EOS error in ZelmaniLeak::ye_of_rho 1")
              write(warnline,"(4i5)") i,j,k,aint(log10(dble(cctk_levfac(1)))/log10(2.0d0))
              call CCTK_WARN(1,warnline)
              write(warnline,"(1P10E15.6)") x(i,j,k),y(i,j,k),z(i,j,k)
              call CCTK_WARN(1,warnline)
              write(warnline,"(1P10E15.6)") rho(i,j,k),eps(i,j,k),y_e(i,j,k),&
                   temperature(i,j,k)
              call CCTK_WARN(1,warnline)
              call CCTK_WARN(0,"aborting!")
           endif
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO
  
end subroutine ZelmaniLeak_calc_entropy
