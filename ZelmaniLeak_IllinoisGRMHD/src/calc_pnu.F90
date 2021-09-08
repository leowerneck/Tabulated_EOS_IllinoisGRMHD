#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
! 4/3*pi*(1MeV)^4/(hc)^3
#define PRESS_NU_CONSTANT 3.52127727d24
#define PI 3.14159265358979d0
#define PI4 97.4090910340024d0
#define PI2 9.86960440108936d0
#define MEV_TO_ERG 1.60217733d-6
#define SPATIAL_DETERMINANT(gxx_,gxy_,gxz_,gyy_,gyz_,gzz_) \
  (-(gxz_)**2*(gyy_) + 2*(gxy_)*(gxz_)*(gyz_) - (gxx_)*(gyz_)**2 - (gxy_)**2*(gzz_) \
   + (gxx_)*(gyy_)*(gzz_))

subroutine ZelmaniLeak_CalcPnu(taurhs,sxrhs,syrhs,szrhs,&
                                   rho,eps,temperature, &
                                   entropy,y_e,munu,pnu,&
                                   wlorentz, vel,&
                                   alp, gxx,gxy,gxz,gyy,gyz,gzz,&
                                   x,y,z,dx,dy,dz,&
                                   nx,ny,nz,dtime)

  use EOS_Omni_module
#ifdef HAVE_CAPABILITY_Fortran
  use cctk
#endif
  implicit none
  

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer   :: nx,ny,nz
  CCTK_REAL :: taurhs(nx,ny,nz)
  CCTK_REAL,target :: sxrhs(nx,ny,nz)
  CCTK_REAL,target :: syrhs(nx,ny,nz)
  CCTK_REAL,target :: szrhs(nx,ny,nz)

  CCTK_REAL :: rho(nx,ny,nz),eps(nx,ny,nz)
  CCTK_REAL :: temperature(nx,ny,nz), y_e(nx,ny,nz)
  CCTK_REAL :: entropy(nx,ny,nz),munu(nx,ny,nz)
  CCTK_REAL :: wlorentz(nx,ny,nz),pnu(nx,ny,nz)
  CCTK_REAL :: vel(nx,ny,nz,3)
  CCTK_REAL :: alp(nx,ny,nz),gxx(nx,ny,nz),gxy(nx,ny,nz),gxz(nx,ny,nz)
  CCTK_REAL :: gyy(nx,ny,nz),gyz(nx,ny,nz),gzz(nx,ny,nz)
  CCTK_REAL :: x(nx,ny,nz),y(nx,ny,nz),z(nx,ny,nz),dx,dy,dz
  CCTK_REAL :: dtime

  ! EOS and other shit
  integer   :: eoskey,keytemp,n,anyerr,keyerr(1)
  CCTK_REAL :: xdummy(1),xeps(1)

  ! local vars
  integer :: i,j,k
  character(len=512) :: warnline
  integer :: iflux,ixoffset,iyoffset,izoffset
  CCTK_REAL :: invdxyz
 
  ! pnu vars:
  CCTK_REAL :: eta
  CCTK_REAL :: pnufac
  CCTK_REAL :: pnuconst, F3const, F3
  CCTK_REAL,allocatable :: dpnudr(:,:,:)
  CCTK_REAL, pointer :: srhs(:,:,:)
  CCTK_REAL :: detg, stress
!  CCTK_REAL :: attfactor


  pnuconst = PRESS_NU_CONSTANT * press_gf
  F3const = 7.0d0*PI4/60.0d0

  ! EOS Omni setup:
  eoskey=4;keytemp=1;n=1;anyerr=0;keyerr=0
  xdummy=0.0d0

  allocate(dpnudr(nx,ny,nz))

  ! compute neutrino pressure
  ! do so only for densities about 1/10 of the density
  ! at which we want the neutrino pressure to be present
  !$OMP PARALLEL DO PRIVATE(i, j, k, xdummy, eta, pnufac, F3, xeps,&
  !$OMP                     keyerr,anyerr)
  do k=1,nz
     do j=1,ny
        do i=1,nx

           if(rho(i,j,k).gt.pnu_rho_start*rho_gf/10.0d0) then           
              keyerr = 0
              anyerr = 0
              xeps = 0.0d0
              call EOS_Omni_short(eoskey,keytemp,GRHydro_eos_rf_prec,&
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
                 call CCTK_WARN(0,"Problem in Pnu EOS call 1!")
              endif
              
              ! let's compute the neutrino pressure
              eta = munu(i,j,k)/temperature(i,j,k)
              F3 = F3const + 0.5d0*eta*eta*(PI2 + 0.5d0*eta*eta)
              pnu(i,j,k) = F3 * pnuconst * temperature(i,j,k)**4
           else
              pnu(i,j,k) = 0.0d0
           endif
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  ! compute the neutrino pressure gradient
  ! and apply it to scon and tau
  ! do 3 swipes for x, y, z
  do iflux=3,1,-1
     dpnudr(:,:,:) = 0.0d0

     select case (iflux)

     case (1)
        ixoffset = 1
        iyoffset = 0
        izoffset = 0
        invdxyz = 1.0d0/(2.0d0*dx);
        srhs => sxrhs

     case (2)
        ixoffset = 0
        iyoffset = 1
        izoffset = 0
        invdxyz = 1.0d0/(2.0d0*dy);
        srhs => syrhs

     case (3)
        ixoffset = 0
        iyoffset = 0
        izoffset = 1  
        invdxyz = 1.0d0/(2.0d0*dz);
        srhs => szrhs

     case DEFAULT
        stop "We should never get here: calc_pnu.F90"

     end select

     ! apply neutrino stress; smoothly taper it on
     ! via a tanh at the transition density
     !
     !$OMP PARALLEL DO PRIVATE(i, j, k, detg, stress)
     do k=GRHydro_stencil-1,nz-GRHydro_stencil+1
        do j=GRHydro_stencil-1,ny-GRHydro_stencil+1
           do i=GRHydro_stencil-1,nx-GRHydro_stencil+1

                 dpnudr(i,j,k) =  &
                      (  pnu(i+ixoffset,j+iyoffset,k+izoffset) -  &
                      pnu(i-ixoffset,j-iyoffset,k-izoffset) )  &
                      * invdxyz 
              
                 detg = SPATIAL_DETERMINANT(gxx(i,j,k),gxy(i,j,k),
                 gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))
                 stress = dpnudr(i,j,k)*alp(i,j,k)*sqrt(detg) * &
                      (tanh( (rho(i,j,k) - pnu_rho_start*rho_gf) / &
                      (pnu_rho_start*rho_gf/10.0d0)) + 1.0d0)/2.0d0

                 srhs(i,j,k) = srhs(i,j,k) - stress
                 taurhs(i,j,k) = taurhs(i,j,k) - vel(i,j,k,iflux)*stress

           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  enddo

  deallocate(dpnudr)

end subroutine ZelmaniLeak_CalcPnu
