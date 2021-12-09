#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
! 4/3*pi*(1MeV)^4/(hc)^3
! Leo says: convert everything to cgs first!
! Python snippet:
! import numpy as np
! import astropy.constants as ct
! import astropy.units as units
! c          = ct.c.cgs.value
! h          = ct.h.cgs.value
! MeV_to_erg = units.MeV.to(units.erg)
! const      = (4/3)*np.pi/(h*c)**3 * MeV_to_erg**4
! print("Computed value:",const)
! print("In the code   : 3.52127727d24")
! print("Relative error:",np.abs((const-3.52127727e24)/const))
! Output:
! Computed value: 3.521275355092826e+24
! In the code   : 3.52127727d24
! Relative error: 5.438106881631993e-07
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
                               nx,ny,nz,ngx,ngy,ngz,dtime)

  use EOS_Omni_module
#ifdef HAVE_CAPABILITY_Fortran
  use cctk
#endif
  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer   :: nx,ny,nz
  integer   :: ngx,ngy,ngz
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

  ! EOS parameters
  integer   :: eoskey,keytemp,n,anyerr,keyerr(1)
  CCTK_REAL :: xdummy(1),xeps(1)

  ! local vars
  integer :: i,j,k
  character(len=512) :: warnline
  integer :: iflux,ixoffset,iyoffset,izoffset
  CCTK_REAL :: invdxyz

  ! pnu vars
  CCTK_REAL :: eta
  CCTK_REAL :: pnufac
  CCTK_REAL :: pnuconst, F3const, F3
  CCTK_REAL,allocatable :: dpnudr(:,:,:)
  CCTK_REAL, pointer :: srhs(:,:,:)
  CCTK_REAL :: detg, stress

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
              call EOS_Omni_short(eoskey,keytemp,ZL_eos_root_finding_precision,&
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
              ! Leo says:
              ! The expression implemented below is obtained using
              ! Eq. 3.19 of http://adsabs.harvard.edu/pdf/1977ApJ...212..859B.
              ! Note that:
              ! -> F3const  = (7/60)*pi^{4} = 11.364... = (5.682) * 2
              ! -> pi^{2}/2 = 4.935... = (2.467) * 2
              ! -> 1/4 = (1/8) * 2
              ! The factors in parenthesis above are the ones that appear
              ! in Eq. 3.19 of the paper, and they are multiplied by a factor
              ! of 2 because the neutrino pressure contains
              !
              ! P_{nu} ~ F_{3}(eta) + F_{3}(-eta),
              !
              ! which eliminates all the terms with odd powers of eta and
              ! multiplies the factors multiplying the even powers of eta
              ! by two.
              F3 = F3const + 0.5d0*eta*eta*(PI2 + 0.5d0*eta*eta)

              ! Leo says: Eq. 36 of https://arxiv.org/pdf/0912.2393.pdf:
              !
              ! P_{nu} = (4pi/3(hc)^3) T^{4} [ F_{3}(eta) + F_{3}(-eta) ]
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
     do k=ngz+1,nz-ngz+1
        do j=ngy+1,ny-ngy+1
           do i=ngx+1,nx-ngx+1

              ! Leo says: second-order centered finite differences
              dpnudr(i,j,k) =  &
                   (  pnu(i+ixoffset,j+iyoffset,k+izoffset) -  &
                   pnu(i-ixoffset,j-iyoffset,k-izoffset) )  &
                   * invdxyz

              detg = SPATIAL_DETERMINANT(gxx(i,j,k),gxy(i,j,k),
              gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k))

              ! Leo says: The source terms for tau and S^{i} are given
              !           by Eqs. (37) in https://arxiv.org/pdf/0912.2393.pdf.
              !
              ! We note, however, that the implementation above is not exactly
              ! the same as the one appearing in the paper, since the paper
              ! focuses on the case of 1D GRHD. For GRMHD, one would start with
              !
              ! nabla_{\mu}T^{\mu\nu} = <Source terms>,
              !
              ! which explains the factor sqrt(-g) = alpha*sqrt(gamma) appearing
              ! below. The base source term is just partial_{i}P_{nu}.
              !
              ! As mentioned above, the function
              !
              ! f(rho) = 0.5 * [ tanh( (rho_{b} - rho_{start}) / (rho_{start}/10) ) + 1 ]
              !
              ! is used to transition smoothly from the regime where the source terms
              ! do not contribute (rho_{b} < rho_{start}) and the regime where they
              ! do contribute (rho_{b} > rho_{start}).
              stress = dpnudr(i,j,k)*alp(i,j,k)*sqrt(detg) * &
                   (tanh( (rho(i,j,k) - pnu_rho_start*rho_gf) / &
                   (pnu_rho_start*rho_gf/10.0d0)) + 1.0d0)/2.0d0

              ! Leo says: First of Eqs. (37) in https://arxiv.org/pdf/0912.2393.pdf.
              srhs(i,j,k) = srhs(i,j,k) - stress

              ! Leo says: Second of Eqs. (37) in https://arxiv.org/pdf/0912.2393.pdf.
              taurhs(i,j,k) = taurhs(i,j,k) - vel(i,j,k,iflux)*stress

           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  enddo

  deallocate(dpnudr)

end subroutine ZelmaniLeak_CalcPnu
