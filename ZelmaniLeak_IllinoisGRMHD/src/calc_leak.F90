#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine ZelmaniLeak_CalcLeak(CCTK_ARGUMENTS)

  use EOS_Omni_Module
#ifdef HAVE_CAPABILITY_Fortran
  use cctk
#endif
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8,parameter :: inv_length_gf3 = 3.220809592832d15
  real*8,parameter :: pi = 4d0*atan2(1d0,1d0) !pi
  real*8,parameter :: pi2 = 2*4d0*atan2(1d0,1d0) ! 2pi
  real*8,parameter :: pi3o2 = 3d0/2d0*4d0*atan2(1d0,1d0) ! 3pi/2
  real*8,parameter :: pio2 = 0.5d0*4d0*atan2(1d0,1d0) ! pi/2
  real*8,parameter :: Binv = 1.0d-51

  integer :: nx,ny,nz
  integer :: i,j,k
  integer :: ii,jj,kk,ll,jj_offset, jj_origin
  integer :: rl
  real*8  :: xchi(3),xtau(3)
  real*8  :: xr,xt,xp,xtemp
  real*8  :: xheatflux(3), xlum(3), xeave(3), xheat(3), xnetheat(3)
  real*8  :: xheaterms(3),xheateave(3)
  real*8  :: depsdt
  real*8  :: dyedt
  real*8  :: ldt
  real*8  :: lumfac
  ! CCTK_REAL :: t,p
  CCTK_REAL, parameter :: tiny  = 1.0d-10
  character(len=512) :: warnline

  ! for 2D (r,theta) interpolation
  integer, parameter :: nfs = 9
  integer, parameter :: nfs2 = 6
  real*8 :: fint(4,nfs),fint3D(8,nfs),fint_out(nfs)
  real*8 :: rr(2),tt(2),pp(2)
  integer isym

  ! Leo says: Spritz MOD - nu variables
  real*8 :: det,radial_velocity
  det = 0.0d0
  radial_velocity = 0.0d0
  lum_inf_nue = 0.0d0
  lum_inf_nua = 0.0d0
  lum_inf_nux = 0.0d0

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  if(do_tau.eq.0) return
  if(leak_in_premerger.ne.0) then
     ! if(in_premerger.eq.0 .and. in_merger.eq.0) return
  else
     if(in_merger.eq.0 .or. (have_interp_data.ne.1) ) then
        !$OMP PARALLEL DO PRIVATE(i,j,k)
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 lum_local(i,j,k,1:3) = 0.0d0
                 heat_local(i,j,k,1:3) = 0.0d0
                 net_heat_local(i,j,k,1:3) = 0.0d0
                 eave_local(i,j,k,1:3) = 0.0d0
                 ! Leo says: Spritz MOD by LSADD
                 lum_nue(i,j,k) = 0.0d0
                 lum_nua(i,j,k) = 0.0d0
                 lum_nux(i,j,k) = 0.0d0
                 ! Leo says: end Spritz MOD by LSADD
              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
        return
     endif
  endif


  lumfac = cctk_delta_space(1)*cctk_delta_space(2)*cctk_delta_space(3)*&
       inv_length_gf3


  if(CCTK_EQUALS(symm,"octant")) then
     lumfac = lumfac * 8.0d0
     isym = 3
  else if (CCTK_EQUALS(symm,"bitant")) then
     lumfac = lumfac * 2.0d0
     isym = 2
  else if (CCTK_EQUALS(symm,"full")) then
     isym = 1
  else
     call CCTK_WARN(0,"Symmetry not known")
  endif

  ! setting up...
  ldt = CCTK_DELTA_TIME*inv_time_gf
  rl = aint(log10(dble(cctk_levfac(1)))/log10(2.0d0))

  ! remember symmetries and ghost zones when implementing this!!!
  if(ntheta.eq.1 .and. nphi.eq.1) then
     ! simplest case; spherically symmetric leakage
     !$OMP PARALLEL DO PRIVATE(i,j,k,xr,xt,ii,jj,ll,xchi,xtau,depsdt,dyedt,xlum, &
     !$OMP                     fint,fint_out,rr,tt,pp,xheatflux,xeave,xheat,xheateave,&
     !$OMP                     xnetheat,xheaterms,xtemp,radial_velocity,det)
     do k=1,nz
        do j=1,ny
           do i=1,nx

              lum_local(i,j,k,1:3) = 0.0d0
              heat_local(i,j,k,1:3) = 0.0d0
              net_heat_local(i,j,k,1:3) = 0.0d0

              ! no leakage in the atmosphere (if there is one)
              if(rho(i,j,k) .le. ZL_rho_threshold) cycle

              if(rl.ge.-1) then

                 ! only leak inside max leak radius
                 if(r(i,j,k).lt.rad_max-drad*1.1d0) then
                    xr = r(i,j,k)
                    ii = floor(xr/drad)+1
                    ii = max(1,min(nrad-1,ii))

                    xlum(1:3) = 0.0d0  
                    xeave(1:3) = 0.0d0
                    xheat(1:3) = 0.0d0
                    xnetheat(1:3) = 0.0d0

                    ! this means that our data point is between ii and ii+1
                    call linterp(rad(ii),rad(ii+1), &
                         zi_xiross(ii,1,1,1),zi_xiross(ii+1,1,1,1),xr,xchi(1))
                    call linterp(rad(ii),rad(ii+1), &
                         zi_xiross(ii,1,1,2),zi_xiross(ii+1,1,1,2),xr,xchi(2))
                    call linterp(rad(ii),rad(ii+1), &
                         zi_xiross(ii,1,1,3),zi_xiross(ii+1,1,1,3),xr,xchi(3))

                    call linterp(rad(ii),rad(ii+1), &
                         zi_tauruff(ii,1,1,1),zi_tauruff(ii+1,1,1,1),xr,xtau(1))
                    call linterp(rad(ii),rad(ii+1), &
                         zi_tauruff(ii,1,1,2),zi_tauruff(ii+1,1,1,2),xr,xtau(2))
                    call linterp(rad(ii),rad(ii+1), &
                         zi_tauruff(ii,1,1,3),zi_tauruff(ii+1,1,1,3),xr,xtau(3))

                    call linterp(rad(ii),rad(ii+1), &
                         zi_heatflux(ii,1,1,1),zi_heatflux(ii+1,1,1,1),xr,xheatflux(1))
                    call linterp(rad(ii),rad(ii+1), &
                         zi_heatflux(ii,1,1,2),zi_heatflux(ii+1,1,1,2),xr,xheatflux(2))
                    call linterp(rad(ii),rad(ii+1), &
                         zi_heatflux(ii,1,1,3),zi_heatflux(ii+1,1,1,3),xr,xheatflux(3))
                    tau3D(i,j,k,1:3) = xtau(1:3)

                    ! if the temperature is below our minimum, floor it
                    xtemp = max(temperature(i,j,k),ZL_T_atmosphere)
                    call calc_leak(rho(i,j,k)*inv_rho_gf,xtemp,&
                         y_e(i,j,k),xchi,xtau,xheatflux,zi_heaterms(1,1,:),&
                         zi_heateave(1,1,:),depsdt,dyedt,ldt,xlum,xeave,xheat,&
                         xnetheat,ZL_rho_atmosphere, rl,r(i,j,k))

                    if(isnan(xlum(1)) .or. &
                         isnan(xlum(2)) .or. &
                         isnan(xlum(3)) .or. &
                         isnan(xnetheat(1)) .or. &
                         isnan(xnetheat(2))) then
                       call CCTK_WARN(1,"LEAKAGE ERROR!")
                       write(warnline,"(A5,i6,1P10E15.6)") "rl: ",rl,x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k)
                       call CCTK_WARN(1,warnline)
                       write(warnline,"(1P10E15.6)") rho(i,j,k),eps(i,j,k),y_e(i,j,k)
                       call CCTK_WARN(1,warnline)
                       write(warnline,"(1P10E15.6)") xtau(1),xtau(2),xtau(3)
                       call CCTK_WARN(1,warnline)
                       write(warnline,"(1P10E15.6)") xlum(1),xlum(2),xlum(3)
                       call CCTK_WARN(1,warnline)
                       write(warnline,"(1P10E15.6)") xnetheat(1),xnetheat(2),xnetheat(3)
                       call CCTK_WARN(1,warnline)
                       write(warnline,"(1P10E15.6)") depsdt*ldt*eps_gf,dyedt*ldt
                       call CCTK_WARN(1,warnline)
                       call CCTK_WARN(0,"Aborting!")
                    endif


                    y_e(i,j,k) = y_e(i,j,k) + dyedt * ldt
                    eps(i,j,k) = eps(i,j,k) + depsdt * ldt * eps_gf

                    lum_local(i,j,k,1:3) = xlum(1:3) * lumfac * Binv
                    lum_int_local(i,j,k,1:3) = lum_int_local(i,j,k,1:3) + lum_local(i,j,k,1:3)*ldt

                    ! Leo says: Spritz MOD by LSADD - lum_nue,a,x
                    lum_nue(i,j,k) = xlum(1) * Binv
                    lum_nua(i,j,k) = xlum(2) * Binv
                    lum_nux(i,j,k) = xlum(3) * Binv
                    ! Leo says: end of Spritz MOD by LSADD

                    net_heat_local(i,j,k,1:3) = xnetheat(1:3) * lumfac * Binv
                    heat_local(i,j,k,1:3) = xheat(1:3) * lumfac * Binv

                    ! the real average energy is obtained by dividing the sum
                    ! of the below over the entire domain by the total luminosity
                    eave_local(i,j,k,1:3) = xeave(1:3) * lum_local(i,j,k,1:3)

                    ! Leo says: original code
                    ! call SpatialDeterminant(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det) ! This subroutine is in Utils.F90
                    ! This is a copy & past of the SpatialDeterminant found in Utils.F90
                    det = -gxz(i,j,k)**2*gyy(i,j,k) + 2*gxy(i,j,k)*gxz(i,j,k)*gyz(i,j,k) - gxx(i,j,k)*gyz(i,j,k)**2 - gxy(i,j,k)**2*gzz(i,j,k) + gxx(i,j,k)*gyy(i,j,k)*gzz(i,j,k)

                    ! Leo says: spherically symmetric leakage
                    radial_velocity = (x(i, j, k) * vel(i, j, k, 1) + y(i, j, k) * vel(i, j, k, 2) + z(i, j, k) * vel(i, j, k, 3)) / ((x(i, j, k)**2 + y(i, j, k)**2 + z(i, j, k)**2) ** 0.5d0 + 1.0d-8)
                    if (isnan(radial_velocity)) then
                       call CCTK_WARN(1,"Radial velocity is NaN!")
                       !call CCTK_WARN(1,"x, y, z") x(i, j, k), y(i, j, k), z(i, j, k)
                    endif
                    lum_inf_nue(i, j, k) = xlum(1) * alp(i, j, k)**2 * w_lorentz(i, j, k) * (1.0d0 + radial_velocity) * sqrt(det) * lumfac * Binv
                    lum_inf_nua(i, j, k) = xlum(2) * alp(i, j, k)**2 * w_lorentz(i, j, k) * (1.0d0 + radial_velocity) * sqrt(det) * lumfac * Binv
                    lum_inf_nux(i, j, k) = xlum(3) * alp(i, j, k)**2 * w_lorentz(i, j, k) * (1.0d0 + radial_velocity) * sqrt(det) * lumfac * Binv
                 endif
              endif
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  else if(ntheta.gt.1.and.nphi.eq.1) then
     ! axisymmetric case, 2D interpolation
     !$OMP PARALLEL DO PRIVATE(i,j,k,xr,xt,ii,jj,ll,xchi,xtau,depsdt,dyedt,xlum, &
     !$OMP                     fint,fint_out,rr,tt,xheatflux,xeave,xheat,xheateave,&
     !$OMP                     xnetheat,xheaterms,xtemp)
     do k=1,nz
        do j=1,ny
           do i=1,nx

              lum_local(i,j,k,1:3) = 0.0d0
              heat_local(i,j,k,1:3) = 0.0d0
              net_heat_local(i,j,k,1:3) = 0.0d0

              ! no leakage in the atmosphere (if there is one)
              if(rho(i,j,k) .le. ZL_rho_threshold) cycle


              ! only leak inside max leak radius
              if(r(i,j,k).lt.rad_max-drad*1.1d0) then
                 xr = r(i,j,k)
                 ii = floor(xr/drad)+1
                 ii = max(1,min(nrad-1,ii))

                 ! theta with symmetry handling
#ifdef SYMMETRIC_OPERATORS
                 xt = asin(z(i,j,k)/(xr+tiny))
                 if(isym.eq.2.or.isym.eq.3) then
                    if(xt.lt.0.) then
                       xt = -xt
                    endif
                 endif
#else
                 xt = acos(z(i,j,k)/(xr+tiny))
                 if(isym.eq.2.or.isym.eq.3) then
                    if(xt.gt.pio2) then
                       xt = pi - xt
                    endif
                 endif
#endif
                 jj = floor(xt/dtheta)+1
                 jj = max(1,min(ntheta-1,jj))

                 ! interpolation setup
                 rr(1) = rad(ii)
                 rr(2) = rad(ii+1)
                 tt(1) = theta(jj)
                 tt(2) = theta(jj+1)

                 fint(1,1:3) = zi_xiross(ii,jj,1,1:3)
                 fint(1,4:6) = zi_tauruff(ii,jj,1,1:3)
                 fint(1,7:9) = zi_heatflux(ii,jj,1,1:3)

                 fint(2,1:3) = zi_xiross(ii+1,jj,1,1:3)
                 fint(2,4:6) = zi_tauruff(ii+1,jj,1,1:3)
                 fint(2,7:9) = zi_heatflux(ii+1,jj,1,1:3)

                 fint(3,1:3) = zi_xiross(ii,jj+1,1,1:3)
                 fint(3,4:6) = zi_tauruff(ii,jj+1,1,1:3)
                 fint(3,7:9) = zi_heatflux(ii,jj+1,1,1:3)

                 fint(4,1:3) = zi_xiross(ii+1,jj+1,1,1:3)
                 fint(4,4:6) = zi_tauruff(ii+1,jj+1,1,1:3)
                 fint(4,7:9) = zi_heatflux(ii+1,jj+1,1,1:3)

                 call linterp2Dn(rr,tt,fint,nfs,xr,xt,fint_out)
                 xchi = fint_out(1:3)
                 xtau = fint_out(4:6)
                 xheatflux = fint_out(7:9)
                 tau3D(i,j,k,1:3) = xtau(1:3)


                 xlum(1:3)  = 0.0d0
                 xeave(1:3) = 0.0d0
                 xheat(1:3)  = 0.0d0
                 xnetheat(1:3)  = 0.0d0

                 call linterp3(theta(jj),theta(jj+1), &
                      zi_heateave(jj,1,1:3),zi_heateave(jj+1,1,1:3),xt,xheateave(1:3))

                 call linterp3(theta(jj),theta(jj+1), &
                      zi_heaterms(jj,1,1:3),zi_heaterms(jj+1,1,1:3),xt,xheaterms(1:3))

                 ! if the temperature is below our minimum, floor it
                 xtemp = max(temperature(i,j,k),ZL_T_atmosphere)
                 call calc_leak(rho(i,j,k)*inv_rho_gf,xtemp,&
                      y_e(i,j,k),xchi,xtau,xheatflux,xheaterms,&
                      xheateave,depsdt,dyedt,ldt,xlum,xeave,xheat,xnetheat,&
                      ZL_rho_atmosphere,rl,r(i,j,k))

                 if(isnan(xlum(1)) .or. &
                      isnan(xlum(2)) .or. &
                      isnan(xlum(3)) .or. isnan(xnetheat(1)) .or. &
                      isnan(xnetheat(2))) then
                    !$OMP CRITICAL
                    call CCTK_WARN(1,"LEAKAGE ERROR!")
                    write(warnline,"(A5,i6,1P10E15.6)") "rl: ",rl,x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k)
                    call CCTK_WARN(1,warnline)
                    write(warnline,"(3i6,1P2E15.6,2i6)") i,j,k,xr,xt,ii,jj
                    call CCTK_WARN(1,warnline)
                    write(warnline,"(1P10E15.6)") rho(i,j,k),eps(i,j,k),y_e(i,j,k)
                    call CCTK_WARN(1,warnline)
                    write(warnline,"(1P10E15.6)") xtau(1),xtau(2),xtau(3)
                    call CCTK_WARN(1,warnline)
                    write(warnline,"(1P10E15.6)") xlum(1),xlum(2),xlum(3)
                    call CCTK_WARN(1,warnline)
                    write(warnline,"(1P10E15.6)") xnetheat(1),xnetheat(2),xnetheat(3)
                    call CCTK_WARN(1,warnline)
                    write(warnline,"(1P10E15.6)") depsdt*ldt*eps_gf,dyedt*ldt
                    call CCTK_WARN(1,warnline)
                    call CCTK_WARN(0,"Aborting!")
                    !$OMP END CRITICAL
                 endif

#if 0
                 if(CCTK_MyProc(cctkGH).eq.0) then
                    !$OMP CRITICAL
                    open(666,file="heat.dat")
                    do ll=1,nrad
                       write(666,"(1P10E15.6)") rad(ll),zi_heatflux(ll,1,1,1),zi_heatflux(ll,1,1,2),&
                            zi_heatflux(ll,1,1,3),zi_heateave(ll,1,1:3)
                    enddo
                    close(666)
                    !                    call CCTK_WARN(0,"end debug")                                                      
                    !$OMP END CRITICAL
                 endif
#endif


                 y_e(i,j,k) = y_e(i,j,k) + dyedt * ldt
                 eps(i,j,k) = eps(i,j,k) + depsdt * ldt * eps_gf

                 lum_local(i,j,k,1:3) = xlum(1:3) * lumfac * Binv
                 lum_int_local(i,j,k,1:3) = lum_int_local(i,j,k,1:3) + lum_local(i,j,k,1:3)*ldt

                 ! Leo says: Spritz MOD by LSADD - lum_nue,a,x
                 lum_nue(i,j,k) = xlum(1) * Binv
                 lum_nua(i,j,k) = xlum(2) * Binv
                 lum_nux(i,j,k) = xlum(3) * Binv
                 ! Leo says: end of Spritz MOD by LSADD

                 net_heat_local(i,j,k,1:3) = xnetheat(1:3) * lumfac * Binv
                 heat_local(i,j,k,1:3) = xheat(1:3) * lumfac * Binv

                 ! the real average energy is obtained by dividing the sum
                 ! of the below over the entire domain by the total luminosity
                 eave_local(i,j,k,1:3) = xeave(1:3) * lum_local(i,j,k,1:3)

                 ! Leo says: original code
                 ! call SpatialDeterminant(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det) ! This subroutine is in Utils.F90
                 ! This is a copy & past of the SpatialDeterminant found in Utils.F90
                 det = -gxz(i,j,k)**2*gyy(i,j,k) + 2*gxy(i,j,k)*gxz(i,j,k)*gyz(i,j,k) - gxx(i,j,k)*gyz(i,j,k)**2 - gxy(i,j,k)**2*gzz(i,j,k) + gxx(i,j,k)*gyy(i,j,k)*gzz(i,j,k)

                 ! Leo says: axisymmetric leakage
                 radial_velocity = (x(i, j, k) * vel(i, j, k, 1) + y(i, j, k) * vel(i, j, k, 2) + z(i, j, k) * vel(i, j, k, 3)) / ((x(i, j, k)**2 + y(i, j, k)**2 + z(i, j, k)**2) ** 0.5d0 + 1.0d-8)
                 if (isnan(radial_velocity)) then
                    call CCTK_WARN(1,"Radial velocity is NaN!")
                    !call CCTK_WARN(1,"x, y, z") x(i, j, k), y(i, j, k), z(i, j, k)
                 endif
                 lum_inf_nue(i, j, k) = xlum(1) * alp(i, j, k)**2 * w_lorentz(i, j, k) * (1.0d0 + radial_velocity) * sqrt(det) * lumfac * Binv
                 lum_inf_nua(i, j, k) = xlum(2) * alp(i, j, k)**2 * w_lorentz(i, j, k) * (1.0d0 + radial_velocity) * sqrt(det) * lumfac * Binv
                 lum_inf_nux(i, j, k) = xlum(3) * alp(i, j, k)**2 * w_lorentz(i, j, k) * (1.0d0 + radial_velocity) * sqrt(det) * lumfac * Binv
              endif
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  else 
     ! full 3D case
     !$OMP PARALLEL DO PRIVATE(i,j,k,xr,xt,xp,ii,jj,kk,ll,jj_offset,jj_origin, &
     !$OMP                     xchi,xtau,depsdt,dyedt,xlum, &
     !$OMP                     fint,fint3D,fint_out,rr,tt,pp,xheatflux,xeave,xheat,xheateave,&
     !$OMP                     xnetheat,xheaterms,xtemp)
     do k=1,nz
        do j=1,ny
           do i=1,nx

              lum_local(i,j,k,1:3) = 0.0d0
              heat_local(i,j,k,1:3) = 0.0d0
              net_heat_local(i,j,k,1:3) = 0.0d0
              eave_local(i,j,k,1:3) = 0.0d0

              ! no leakage in the atmosphere (if there is one)
              if(rho(i,j,k) .le. ZL_rho_threshold) cycle

              ! only leak inside max leak radius
              if(r(i,j,k).lt.rad_max-drad*1.1d0) then
                 xr = r(i,j,k)
                 ii = floor(xr/drad)+1
                 ii = max(1,min(nrad-1,ii))

                 ! theta with symmetry handling
#ifdef SYMMETRIC_OPERATORS
                 if(xr .ne. 0d0) then
                    xt = asin(z(i,j,k)/xr)
                 else
                    xt = 0d0
                 end if
                 if(isym.eq.2.or.isym.eq.3) then
                    if(xt.lt.0.) then
                       xt = -xt
                    endif
                    jj_origin = 1
                 else
                    jj_origin = ntheta/2
                 endif
                 ! in tau.cc the theta coordinates of the nodes are computed
                 ! such that for i = 0,...,ntheta-1 and for a full run
                 !   theta(ntheta/2  -i) = -(2*i+1)*(dtheta/2)
                 !   theta(ntheta/2+1+i) = +(2*i+1)*(dtheta/2)
                 ! and for bitant/octant runs we only have a single point theta<0
                 !   theta(i) = (2*i-1)*(dtheta/2)
                 ! There is a spaceing of dtheta and there is no point
                 ! with theta .eq. 0.  We do have points at the poles theta =
                 ! PI and theta = -PI.  Here we need to invert the relation to
                 ! arrive at i given theta such that the node picked is beyond
                 ! theta as seen from the origin (more positive when theta > 0
                 ! and more negative when theta < 0).
                 ! We want to do this
                 ! * symmetric around theta = 0
                 ! * without floating point additions/subtractions since they
                 !   are inaccurate
                 jj = ceiling(abs(xt)/(0.5d0*dtheta)) ! measured in halfs of dtheta away from xt=0
                 if(mod(jj,2) .eq. 0) then ! round up to next odd integer
                    jj = jj + 1
                 end if
                 ! at this point jj is (2*i+1) in the formulae above. We now
                 ! have to convert to the actual index into the array.
                 if(xt.ge.0d0) then ! we always interpolate from xt=0 outwards
                    ! since jj is beyond theta, pick an offset so that
                    ! jj+jj_offset is between theta and 0
                    jj_offset = -1
                    ! 2*i+1 = jj => i = (jj-1)/2 and we need ntheta/2+1 + i
                    jj = jj_origin+1 + (jj-1)/2
                    jj = max(2,min(ntheta,jj))                ! limit to valid range
                 else
                    ! since jj is beyond theta, pick an offset so that
                    ! jj+jj_offset is between theta and 0
                    jj_offset = +1
                    ! 2*i+1 = jj => i = (jj-1)/2 and we need ntheta/2 - i
                    jj = jj_origin - (jj-1)/2
                    jj = max(1,min(ntheta-1,jj))                ! limit to valid range
                 end if

                 ! phi measured either from +x or -x axis but always [0:pi]
                 if(y(i,j,k) .gt. 0d0 .or. y(i,j,k) .eq. 0d0 .and. x(i,j,k) .ge. 0d0) then
                    xp = atan2(y(i,j,k),x(i,j,k))
                 else
                    xp = atan2(-y(i,j,k),-x(i,j,k))
                 end if
                 if(isym.eq.3) then
                    ! phi repeats every pi/2, and phi is already limited to [0:pi]
                    if(xp.ge.pio2) then
                       xp = xp - pio2
                    end if
                 end if
                 kk = int(xp/dphi) + 1
                 if(isym.eq.3) then
                    kk = max(1,min(nphi-1,kk))
                 else
                    kk = max(1,min(nphi/2-1,kk))
                    if(.not. (y(i,j,k) .gt. 0d0 .or. &
                         y(i,j,k) .eq. 0d0 .and. x(i,j,k) .ge. 0d0)) then
                       kk = nphi/2 + kk
                       xp = -xp
                    end if
                 end if
#else
                 xt = acos(z(i,j,k)/(xr+tiny))
                 if(isym.eq.2.or.isym.eq.3) then
                    if(xt.gt.pio2) then
                       xt = pi - xt
                    endif
                 endif
                 jj = floor(xt/dtheta)+1
                 jj = max(1,min(ntheta-1,jj))

                 ! phi with symmetry handling
                 xp = atan2(y(i,j,k),x(i,j,k)+tiny)
                 if(xp.lt.0) then
                    ! make mapping unique
                    xp = xp + pi2 
                 endif
                 if(isym.eq.3) then
                    if(xp.gt.pi3o2) then
                       xp = xp - pi3o2
                    else if(xp.gt.pi) then
                       xp = xp - pi
                    else if(xp.gt.pio2) then
                       xp = xp - pio2
                    endif
                 endif
                 kk = floor(xp/dphi)+1
                 kk = max(1,min(nphi-1,kk))

                 jj_offset = 1
                 ! jj_origin us not used
#endif

                 ! interpolation setup
                 rr(1) = rad(ii)
                 rr(2) = rad(ii+1)
                 tt(1) = theta(jj)
                 tt(2) = theta(jj+jj_offset)
                 pp(1) = phi(kk)
                 pp(2) = phi(kk+1)

                 fint3D(1,1:3) = zi_xiross(ii,jj,kk,1:3)
                 fint3D(1,4:6) = zi_tauruff(ii,jj,kk,1:3)
                 fint3D(1,7:9) = zi_heatflux(ii,jj,kk,1:3)
                 !
                 fint3D(2,1:3) = zi_xiross(ii,jj+jj_offset,kk,1:3)
                 fint3D(2,4:6) = zi_tauruff(ii,jj+jj_offset,kk,1:3)
                 fint3D(2,7:9) = zi_heatflux(ii,jj+jj_offset,kk,1:3)
                 !
                 fint3D(3,1:3) = zi_xiross(ii,jj,kk+1,1:3)
                 fint3D(3,4:6) = zi_tauruff(ii,jj,kk+1,1:3)
                 fint3D(3,7:9) = zi_heatflux(ii,jj,kk+1,1:3)
                 !
                 fint3D(4,1:3) = zi_xiross(ii,jj+jj_offset,kk+1,1:3)
                 fint3D(4,4:6) = zi_tauruff(ii,jj+jj_offset,kk+1,1:3)
                 fint3D(4,7:9) = zi_heatflux(ii,jj+jj_offset,kk+1,1:3)
                 !
                 fint3D(5,1:3) = zi_xiross(ii+1,jj,kk,1:3)
                 fint3D(5,4:6) = zi_tauruff(ii+1,jj,kk,1:3)
                 fint3D(5,7:9) = zi_heatflux(ii+1,jj,kk,1:3)
                 !
                 fint3D(6,1:3) = zi_xiross(ii+1,jj+jj_offset,kk,1:3)
                 fint3D(6,4:6) = zi_tauruff(ii+1,jj+jj_offset,kk,1:3)
                 fint3D(6,7:9) = zi_heatflux(ii+1,jj+jj_offset,kk,1:3)
                 !
                 fint3D(7,1:3) = zi_xiross(ii+1,jj,kk+1,1:3)
                 fint3D(7,4:6) = zi_tauruff(ii+1,jj,kk+1,1:3)
                 fint3D(7,7:9) = zi_heatflux(ii+1,jj,kk+1,1:3)
                 !
                 fint3D(8,1:3) = zi_xiross(ii+1,jj+jj_offset,kk+1,1:3)
                 fint3D(8,4:6) = zi_tauruff(ii+1,jj+jj_offset,kk+1,1:3)
                 fint3D(8,7:9) = zi_heatflux(ii+1,jj+jj_offset,kk+1,1:3)

                 if((rr(1) - xr)*(rr(2) - xr) .gt. 0d0 .or. &
                      (tt(1) - xt)*(tt(2) - xt) .gt. 0d0 .or. &
                      (pp(1) - xp)*(pp(2) - xp) .gt. 0d0) then
                    write(warnline,"('suspicious interpolation interval (',E15.6,',',E15.6,')x(',E15.6,',',E15.6,')x(',E15.6,',',E15.6,') for point (',3E15.6,')')") rr(1),rr(2),tt(1),tt(2),pp(1),pp(2),xr,xt,xp
                    call CCTK_WARN(1,warnline)
                 end if
                 call linterp3Dn(tt,pp,rr,fint3D,nfs,xt,xp,xr,fint_out)
                 xchi = fint_out(1:3)
                 xtau = fint_out(4:6)
                 xheatflux = fint_out(7:9)
                 tau3D(i,j,k,1:3) = xtau(1:3)

                 if(isnan(xtau(1))) then
                    if(CCTK_MyProc(cctkGH).eq.0) then
                       !$OMP CRITICAL
                       write(warnline,"(3i6,1P10E15.6)") i,j,k,x(i,j,k),y(i,j,k),z(i,j,k)
                       call CCTK_WARN(1,warnline)
                       write(warnline,"(i5,1P10E15.6)") ii,rad(ii),xr,rad(ii+1)
                       call CCTK_WARN(1,warnline)
                       write(warnline,"(i5,1P10E15.6)") jj,theta(jj),xt,theta(jj+jj_offset)
                       call CCTK_WARN(0,warnline)
                       write(warnline,"(i5,1P10E15.6)") kk,phi(kk),xp,phi(kk+1)
                       call CCTK_WARN(0,warnline)
                       write(warnline,"(1P10E15.6)") xtau(1:3)
                       call CCTK_WARN(0,warnline)
                       call CCTK_WARN(0,"aborting!!!")
                       !$OMP END CRITICAL
                    endif
                 endif

                 xlum(1:3)  = 0.0d0  
                 xeave(1:3) = 0.0d0
                 xheat(1:3)  = 0.0d0
                 xnetheat(1:3)  = 0.0d0

                 fint(:,:) = 0.0d0

                 fint(1,1:3) = zi_heateave(jj,kk,1:3)
                 fint(1,4:6) = zi_heaterms(jj,kk,1:3)

                 fint(2,1:3) = zi_heateave(jj+jj_offset,kk,1:3)
                 fint(2,4:6) = zi_heaterms(jj+jj_offset,kk,1:3)

                 fint(3,1:3) = zi_heateave(jj,kk+1,1:3)
                 fint(3,4:6) = zi_heaterms(jj,kk+1,1:3)

                 fint(4,1:3) = zi_heateave(jj+jj_offset,kk+1,1:3)
                 fint(4,4:6) = zi_heaterms(jj+jj_offset,kk+1,1:3)

                 ! Leo says: double check that these quantities are
                 !           the ones given by Eqs. (A3) and (A4) of
                 !           https://arxiv.org/abs/astro-ph/0302301
                 call linterp2Dn(tt,pp,fint,nfs2,xt,xp,fint_out)
                 xheateave(1:3) = fint_out(1:3)
                 xheaterms(1:3) = fint_out(4:6)

                 ! if the temperature is below our minimum, floor it
                 xtemp = max(temperature(i,j,k),ZL_T_atmosphere)
                 ! Leo says: Spritz MOD by LSADD
                 if(isnan(xchi(1)).or.isnan(xchi(2)).or.isnan(xchi(3))) then
                    write(*,*) "xchi is NaN. ",rr,tt,pp
                 endif
                 ! Leo says: end of Spritz MOD by LSADD
                 call calc_leak(rho(i,j,k)*inv_rho_gf,xtemp,&
                      y_e(i,j,k),xchi,xtau,xheatflux,xheaterms,&
                      xheateave,depsdt,dyedt,ldt,xlum,xeave,xheat,xnetheat,&
                      ZL_rho_atmosphere,rl,r(i,j,k))

                 if(isnan(xlum(1)) .or. &
                      isnan(xlum(2)) .or. &
                      isnan(xlum(3)) .or. isnan(xnetheat(1)) .or. &
                      isnan(xnetheat(2))) then
                    !$OMP CRITICAL
                    call CCTK_WARN(1,"LEAKAGE ERROR!")
                    write(warnline,"(A5,i6,1P10E15.6)") "rl: ",rl,x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k)
                    call CCTK_WARN(1,warnline)
                    write(warnline,"(3i6,1P3E15.6,3i6)") i,j,k,xr,xt,xp,ii,jj,kk
                    call CCTK_WARN(1,warnline)
                    write(warnline,"(1P10E15.6)") rho(i,j,k),eps(i,j,k),y_e(i,j,k)
                    call CCTK_WARN(1,warnline)
                    write(warnline,"(1P10E15.6)") xtau(1),xtau(2),xtau(3)
                    call CCTK_WARN(1,warnline)
                    write(warnline,"(1P10E15.6)") xlum(1),xlum(2),xlum(3)
                    call CCTK_WARN(1,warnline)
                    write(warnline,"(1P10E15.6)") xnetheat(1:3)
                    call CCTK_WARN(1,warnline)
                    write(warnline,"(1P10E15.6)") depsdt*ldt*eps_gf,dyedt*ldt
                    call CCTK_WARN(1,warnline)
                    call CCTK_WARN(0,"Aborting!")
                    !$OMP END CRITICAL
                 endif

                 y_e(i,j,k) = y_e(i,j,k) + dyedt * ldt
                 eps(i,j,k) = eps(i,j,k) + depsdt * ldt * eps_gf

                 lum_local(i,j,k,1:3) = xlum(1:3) * lumfac * Binv
                 lum_int_local(i,j,k,1:3) = lum_int_local(i,j,k,1:3) + lum_local(i,j,k,1:3)*ldt

                 !LSADD lum_nue,a,x
                 lum_nue(i,j,k) = xlum(1) * Binv
                 lum_nua(i,j,k) = xlum(2) * Binv
                 lum_nux(i,j,k) = xlum(3) * Binv
                 !end LSADD

                 net_heat_local(i,j,k,1:3) = xnetheat(1:3) * lumfac * Binv
                 heat_local(i,j,k,1:3) = xheat(1:3) * lumfac * Binv

                 ! the real average energy is obtained by dividing the sum
                 ! of the below over the entire domain by the total luminosity
                 eave_local(i,j,k,1:3) = xeave(1:3) * lum_local(i,j,k,1:3)

                 ! Leo says: original code
                 ! call SpatialDeterminant(gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k),det) ! This subroutine is in Utils.F90
                 ! This is a copy & past of the SpatialDeterminant found in Utils.F90
                 det = -gxz(i,j,k)**2*gyy(i,j,k) + 2*gxy(i,j,k)*gxz(i,j,k)*gyz(i,j,k) - gxx(i,j,k)*gyz(i,j,k)**2 - gxy(i,j,k)**2*gzz(i,j,k) + gxx(i,j,k)*gyy(i,j,k)*gzz(i,j,k)

                 ! Leo says: Full 3D
                 radial_velocity = (x(i, j, k) * vel(i, j, k, 1) + y(i, j, k) * vel(i, j, k, 2) + z(i, j, k) * vel(i, j, k, 3)) / ((x(i, j, k)**2 + y(i, j, k)**2 + z(i, j, k)**2) ** 0.5d0 + 1.0d-8)
                 if (isnan(radial_velocity)) then
                    call CCTK_WARN(1,"Radial velocity is NaN!")
                    !call CCTK_WARN(1,"x, y, z") x(i, j, k), y(i, j, k), z(i, j, k)
                 endif
                 lum_inf_nue(i, j, k) = xlum(1) * alp(i, j, k)**2 * w_lorentz(i, j, k) * (1.0d0 + radial_velocity) * sqrt(det) * lumfac * Binv
                 lum_inf_nua(i, j, k) = xlum(2) * alp(i, j, k)**2 * w_lorentz(i, j, k) * (1.0d0 + radial_velocity) * sqrt(det) * lumfac * Binv
                 lum_inf_nux(i, j, k) = xlum(3) * alp(i, j, k)**2 * w_lorentz(i, j, k) * (1.0d0 + radial_velocity) * sqrt(det) * lumfac * Binv
              endif
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
  endif

contains 
  subroutine linterp2Dn(xin,yin,fin,n,x,y,fout)

#ifdef HAVE_CAPABILITY_Fortran
    use cctk
#endif
    implicit none
    real*8  :: xin(2),yin(2)
    integer :: n
    real*8  :: fin(4,n)
    real*8  :: x,y,fout(n)
    real*8  :: invdxdy,dx,dy

    !      y2 f3           f4
    !
    !
    !      y1 f1           f2
    !         x1           x2
    ! Leo says: original code
    ! invdxdy = 1.0d0/( (xin(2)-xin(1))*(yin(2)-yin(1)) )
    ! Leo says: Spritz MOD by LSADD
    dy = yin(2)-yin(1)
    dx = xin(2)-xin(1)
    if(dy.eq.0.0d0) dy = 1.0e-8
    if(dx.eq.0.0d0) dx = 1.0e-8
    invdxdy = 1.0d0/( dx*dy )
    ! Leo says: end of Spritz MOD by LSADD

#ifdef SYMMETRIC_OPERATORS
    fout(1:n) = invdxdy * (  (fin(1,1:n) * (xin(2)-x) &
         + fin(2,1:n) * (x-xin(1)))*(yin(2)-y) &
         + (fin(3,1:n) * (xin(2)-x) &
         + fin(4,1:n) * (x-xin(1)))*(y-yin(1)) )
#else
    fout(1:n) = invdxdy * (  fin(1,1:n) * (xin(2)-x)*(yin(2)-y) &
         + fin(2,1:n) * (x-xin(1))*(yin(2)-y) &
         + fin(3,1:n) * (xin(2)-x)*(y-yin(1)) &
         + fin(4,1:n) * (x-xin(1))*(y-yin(1)) )
#endif


  end subroutine linterp2Dn

  subroutine linterp3Dn(xin,yin,zin,fin,n,x,y,z,fout)

#ifdef HAVE_CAPABILITY_Fortran
    use cctk
#endif
    implicit none
    real*8  :: xin(2),yin(2),zin(2)
    integer :: n
    real*8  :: fin(8,n),fint2D(4,n)
    real*8  :: x,y,z,fout(n)
    real*8  :: invdxdydz,dx,dy,dz

    !        y2 f7           f8
    !
    !  z2
    !        y1 f5           f6
    !           x1           x2

    !        y2 f3           f4
    !
    !  z1
    !        y1 f1           f2
    !           x1           x2
    ! Leo says: original code
    ! invdxdydz = 1.0d0/( (xin(2)-xin(1))*(yin(2)-yin(1))*(zin(2)-zin(1)) )
    ! Leo says: Spritz MOD by LSADD
    dz = zin(2)-zin(1)
    dy = yin(2)-yin(1)
    dx = xin(2)-xin(1)
    if(dz.eq.0.0d0) dz = 1.0e-8
    if(dy.eq.0.0d0) dy = 1.0e-8
    if(dx.eq.0.0d0) dx = 1.0e-8
    invdxdydz = 1.0d0/( dx*dy*dz )
    ! Leo says: end of Spritz MOD by LSADD

    ! auto-generated via mathematica, terms slightly reordered:
#ifdef SYMMETRIC_OPERATORS
    fout(1:n) = invdxdydz * (                                     &
         ((fin(8,1:n)*(x - xin(1)) +                              &
         fin(7,1:n)*(xin(2) - x))*(y - yin(1)) +                &
         (fin(6,1:n)*(x - xin(1)) +                              &
         fin(5,1:n)*(xin(2) - x))*(yin(2) - y))*(z - zin(1)) +  &
         ((fin(4,1:n)*(x - xin(1)) +                              &
         fin(3,1:n)*(xin(2) - x))*(y - yin(1)) +                &
         (fin(2,1:n)*(x - xin(1)) +                              &
         fin(1,1:n)*(xin(2) - x))*(yin(2) - y))*(zin(2) - z)    &
         )
#else
    fout(1:n) = invdxdydz * (                                 &
         fin(8,1:n)*(x - xin(1))*(y - yin(1))*(z - zin(1)) +  &
         fin(7,1:n)*(xin(2) - x)*(y - yin(1))*(z - zin(1)) +  &
         fin(5,1:n)*(x - xin(2))*(y - yin(2))*(z - zin(1)) +  &
         fin(6,1:n)*(x - xin(1))*(yin(2) - y)*(z - zin(1)) +  &
         fin(4,1:n)*(x - xin(1))*(y - yin(1))*(zin(2) - z) +  &
         fin(3,1:n)*(xin(2) - x)*(y - yin(1))*(zin(2) - z) +  &
         fin(2,1:n)*(x - xin(1))*(yin(2) - y)*(zin(2) - z) +  &
         fin(1,1:n)*(xin(2) - x)*(yin(2) - y)*(zin(2) - z) ) 
#endif
  end subroutine linterp3Dn

  subroutine linterp(x1,x2,y1,y2,x,y)

#ifdef HAVE_CAPABILITY_Fortran
    use cctk
#endif
    implicit none

    CCTK_REAL slope,x1,x2,y1,y2,x,y

    if (x2.lt.x1) then
       call CCTK_WARN(0,"Error in linterp!")
    endif

    slope = (y2 - y1) / (x2 - x1)
    y = slope*(x-x1) + y1

  end subroutine  linterp

  subroutine linterp3(x1,x2,y1,y2,x,y)

#ifdef HAVE_CAPABILITY_Fortran
    use cctk
#endif
    implicit none

    CCTK_REAL slope(3),x1(3),x2(3),y1(3),y2(3),x,y(3)

    slope = (y2 - y1) / (x2 - x1)
    y = slope*(x-x1) + y1

  end subroutine  linterp3

end subroutine ZelmaniLeak_CalcLeak

subroutine calc_leak(rho,temp,ye,chi,tau,heatflux,heaterms,heateave,&
     depsdt,dyedt,ldt,lum,eave,heatout,netheatout,rho_min,reflevel,rad)
  ! WARNING: Be careful when changing the arguments to this function; it is also
  ! called from calc_tau to get the luminosity available for heating
  ! along the rays.


  !  use eosmodule, only : energy_shift, eos_yemin, eos_yemax
  use EOS_Omni_module, only: energy_shift, eos_yemin, eos_yemax, &
       rho_gf, inv_rho_gf, inv_eps_gf

#ifdef HAVE_CAPABILITY_Fortran
  use cctk
#endif
  implicit none
  DECLARE_CCTK_PARAMETERS

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
  real*8, intent(in) :: rho_min !atmosphere density
  real*8, intent(in) :: rad ! radius of the point we are dealing with
  integer, intent(in) :: reflevel

  !EOS & local variables
  integer keytemp, keyerr
  real*8 :: precision = 1.0d-10
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
  integer :: anyerr

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
  character(len=512) :: warnline

  matter_rho = rho
  matter_temperature = max(temp,ZL_T_atmosphere)
  matter_ye = min(max(ye,ZL_yemin),ZL_yemax)

  keytemp = 1
  keyerr = 0
  anyerr = 0
  call EOS_Omni_full(eoskey,keytemp,ZL_eos_root_finding_precision,npoints,&
       matter_rho*rho_gf,matter_enr,matter_temperature,matter_ye, &
       matter_prs,matter_ent,matter_cs2,matter_dedt,matter_dpderho, &
       matter_dpdrhoe,matter_xa,matter_xh,matter_xn,matter_xp,matter_abar, &
       matter_zbar,matter_mue,matter_mun,matter_mup,matter_muhat, &
       keyerr,anyerr)
  matter_enr = matter_enr * inv_eps_gf
  if (anyerr.ne.0) then
     !$OMP CRITICAL
     write(warnline,"(A15,1P10E15.6)") "rho: ", matter_rho
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "temperature: ", matter_temperature
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "ye: ", matter_ye
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,i10)") "eos error", keyerr
     call CCTK_WARN(1,warnline)
     call CCTK_WARN(0,"set_eos_variables: EOS error in leakage calc_leak")
     !$OMP END CRITICAL
  endif

  ! don't do anything outside the shock
  if(matter_xh.gt.0.5.and.matter_rho.lt.1.0d13) then 
     lum(1:3) = 0.0d0
     depsdt = 0.0d0
     dyedt = 0.0d0
     eave = 0.0d0
     return
  endif

  ! Don't do anything in the atmosphere
  if(matter_rho*rho_gf.lt.ZL_rho_threshold) then 
     lum(1:3) = 0.0d0
     depsdt = 0.0d0
     dyedt = 0.0d0
     eave = 0.0d0
     return
  endif

  eta_e = matter_mue/temp
  eta_p = matter_mup/temp
  eta_n = matter_mun/temp

  eta_hat = eta_n - eta_p - Qnp/temp
  eta_nue = eta_e - eta_n + eta_p !fully includes effects of rest masses
  eta_nua = -eta_nue
  eta_nux = 0.0d0

  !interpolate etas like we do in GR1D
  ! Leo says: Eq. (24) of https://arxiv.org/abs/0912.2393
  eta_nue = eta_nue*(1.0d0-exp(-tau(1)))
  eta_nua = eta_nua*(1.0d0-exp(-tau(2)))

  ! Leo says: this is sigma_{nui,nuc} * rho * N_{avo},
  !           where sigma_{nui,nuc} is given just below
  !           Eq. (A17) of https://arxiv.org/abs/astro-ph/0302301
  scattering_kappa = rho*avo*0.25d0*sigma_0/me_mev**2
  kappa_tilde_nu_scat(1,1) = matter_xn*scattering_kappa ! Leo says: electron neutrino scattering off neutrons
  kappa_tilde_nu_scat(1,2) = matter_xp*scattering_kappa ! Leo says: electron neutrino scattering off protons
  kappa_tilde_nu_scat(2,1) = matter_xn*scattering_kappa ! Leo says: electron antineutrino scattering off neutrons
  kappa_tilde_nu_scat(2,2) = matter_xp*scattering_kappa ! Leo says: electron antineutrino scattering off protons
  kappa_tilde_nu_scat(3,1) = matter_xn*scattering_kappa ! Leo says: heavy-lepton neutrinos scattering off neutrons
  kappa_tilde_nu_scat(3,2) = matter_xp*scattering_kappa ! Leo says: heavy-lepton neutrinos scattering off protons

  ! Leo says: this is sigma_{nui,A} * rho * N_{avo},
  !           where sigma_{nui,A} is given just below
  !           Eq. (A18) of https://arxiv.org/abs/astro-ph/0302301
  !           One factor of A is missing because kappa multiples the
  !           number fraction, not mass fractions.
  scattering_kappa = rho*avo*0.0625d0*sigma_0/me_mev**2* &
       matter_abar*(1.0d0-matter_zbar/matter_abar)**2

  ! Leo says: Spritz comment:
  ! LSMOD sometimes this generated errors, due to matter_abar = 0
  if(isnan(scattering_kappa)) then
     scattering_kappa = 0.0d0
  endif
  ! Leo says: end of spritz comment

  kappa_tilde_nu_scat(1,3) = matter_xh*scattering_kappa ! Leo says: electron neutrino scattering off heavy nuclei
  kappa_tilde_nu_scat(2,3) = matter_xh*scattering_kappa ! Leo says: electron antineutrino scattering off heavy nuclei
  kappa_tilde_nu_scat(3,3) = matter_xh*scattering_kappa ! Leo says: heavy-lepton neutrinos scattering off heavy nuclei

  ! Leo says: Eq. (A8) in https://arxiv.org/abs/1306.4953
  !           There is an extra multiplicative factor of
  !           N_{avo} * rho
  eta_pn = avo*rho*(matter_xn-matter_xp)/(exp(eta_hat)-1.0d0)
  eta_pn = max(eta_pn,0.0d0)

  ! Leo says: Eq. (A8) in https://arxiv.org/abs/1306.4953
  !           with n <-> p. Note again the extra multiplicative
  !           factor of N_{avo} * rho
  eta_np = avo*rho*(matter_xp-matter_xn)/(exp(-eta_hat)-1.0d0)
  eta_np = max(eta_np,0.0d0)

  if (rho.lt.1.0d11) then
     !non degenerate here, use mass fractions as chemical potentials fail at low densities
     eta_pn = avo*rho*matter_xp
     eta_np = avo*rho*matter_xn
  endif

  !absorption
  ! Leo says: Constants in Eq. (A11) of Ruffert et al. (1996)
  !           https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  abs_kappa = (1.0d0+3.0d0*alpha**2)*0.25d0*sigma_0/me_mev**2

  ! Leo says: Eq. (A15) of Ruffert et al. (1996)
  !           https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  block_factor = 1.0d0 + exp(eta_e-get_fermi_integral_leak(5,eta_nue)/ &
       (get_fermi_integral_leak(4,eta_nue)+verysmall))

  ! Leo says: Eq. (A11) of Ruffert et al. (1996) without Fermi-Dirac integrals
  !           https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  kappa_tilde_nu_abs(1,1) = eta_np*abs_kappa/block_factor
  kappa_tilde_nu_abs(2,1) = 0.0d0 !no absorption of a-type on neutrons
  kappa_tilde_nu_abs(3,1) = 0.0d0 !no absorption of x-type neutrinos
  kappa_tilde_nu_abs(1,2) = 0.0d0 !no absorption of e-type on protons

  ! Leo says: Eq. (A16) of Ruffert et al. (1996)
  !           https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  block_factor = 1.0d0 + exp(-eta_e-get_fermi_integral_leak(5,eta_nua)/ &
       (get_fermi_integral_leak(4,eta_nua)+verysmall))

  ! Leo says: Eq. (A16) of Ruffert et al. (1996) without Fermi-Dirac integrals
  !           https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  kappa_tilde_nu_abs(2,2) = eta_pn*abs_kappa/block_factor
  kappa_tilde_nu_abs(3,2) = 0.0d0 !no absorption of x-type neutrinos
  kappa_tilde_nu_abs(1,3) = 0.0d0 !no absorption on nuclei
  kappa_tilde_nu_abs(2,3) = 0.0d0 !no absorption on nuclei
  kappa_tilde_nu_abs(3,3) = 0.0d0 !no absorption on nuclei

  ! Leo says: zeta = lambda^{-1}E^{-2}, as defined in Eqs. (A32) and
  !           (A33) of https://arxiv.org/abs/1306.4953
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

  ! Leo says: useful Spritz comment:
  ! Federico: Q_diff components below should be compared with the following references:
  ! - Rosswog, Liebendorfer, mnras (2003) - Eq (A35) - https://arxiv.org/abs/astro-ph/0302301
  ! - Galeazzi, Kastaun, Rezzolla, Font, PRD (2013) - Eq (A35) - https://arxiv.org/abs/1306.4953
  ! NOTE1: the factor 6 at the denominator, seems to be different from the
  !        references, where it is set to 3.
  ! NOTE2: one possible explanation for NOTE1 can be found in reference:
  !        - O'Connor, Ott, (2010) - point (ii) before Eq. (25)
!!! TO BE CLARIFIED !!!
  rate_const = 4.0d0*pi*clite*zeta(1)/(hc_mevcm**3*6.0d0*chi(1)**2)
  R_diff(1) = rate_const*temp*get_fermi_integral_leak(0,eta_nue)
  Q_diff(1) = rate_const*temp**2*get_fermi_integral_leak(1,eta_nue)

  rate_const = 4.0d0*pi*clite*zeta(2)/(hc_mevcm**3*6.0d0*chi(2)**2)
  R_diff(2) = rate_const*temp*get_fermi_integral_leak(0,eta_nua)
  Q_diff(2) = rate_const*temp**2*get_fermi_integral_leak(1,eta_nua)

  rate_const = 16.0d0*pi*clite*zeta(3)/(hc_mevcm**3*6.0d0*chi(3)**2)
  R_diff(3) = rate_const*temp*get_fermi_integral_leak(0,eta_nux)
  Q_diff(3) = rate_const*temp**2*get_fermi_integral_leak(1,eta_nux)

  ! Leo says: Spritz MOD by LSADD
  if((isnan(R_diff(1)).or.isnan(R_diff(2)).or.isnan(R_diff(3))).and.reflevel.ge.0) then
     write(*,*) "Problem in diffusion rates, R_diff(nue),R_diff(nua),zeta(1),zeta(2),zeta(3): ",R_diff(1),R_diff(2),R_diff(3),rate_const
  endif
  ! Leo says: end of Spritz MOD by LSADD
  
  !now for the free emission
  ! Leo says: Eq. (A6) in https://arxiv.org/abs/1306.4953
  beta = pi*clite*(1.0d0+3.0d0*alpha**2)*sigma_0/(hc_mevcm**3*me_mev**2)

  R_loc = 0.0d0
  Q_loc = 0.0d0

  ! Leo says: useful Spritz comment:
  ! electron capture rates
  ! Federico: Following equations for Q_loc(i), with i = 1,2 should be compared to:
  ! - Galeazzi, Kastaun, Rezzolla, Font, PRD (2013) - i=1 --> Eq (A4)
  !                                                   i=2 --> Eq (A2)
  !           Alternatively, check:
  ! - Rosswog, Liebendorfer, mnras (2003) - i=1 --> Eq (A10)
  !                                         i=2 --> Eq (A13)

  ! Leo says: Neutrino emission rate due to electron capture
  ! Eq. (A6) in https://arxiv.org/abs/astro-ph/0302301
  R_loc(1) = beta*eta_pn*temp**5*get_fermi_integral_leak(4,eta_e)

  ! Leo says: Neutrino cooling rate due to electron capture
  ! Eq. (A10) in https://arxiv.org/abs/astro-ph/0302301
  Q_loc(1) = beta*eta_pn*temp**6*get_fermi_integral_leak(5,eta_e)

  ! Leo says: Neutrino emission rate due to positron capture
  ! Eq. (A12) in https://arxiv.org/abs/astro-ph/0302301
  R_loc(2) = beta*eta_np*temp**5*get_fermi_integral_leak(4,-eta_e)

  ! Leo says: Neutrino cooling rate due to positron capture
  ! Eq. (A13) in https://arxiv.org/abs/astro-ph/0302301
  Q_loc(2) = beta*eta_np*temp**6*get_fermi_integral_leak(5,-eta_e)

  ! Leo says: Spritz MOD by LSADD
  if(temp < 0) then
     R_loc(1) = 0.0d0
     R_loc(2) = 0.0d0
     Q_loc(1) = 0.0d0
     Q_loc(2) = 0.0d0
  endif
  ! Leo says: end of Spritz MOD by LSADD

  ! Leo says: useful Spritz comment:
  !e-e+ pair processes from Ruffert et al.
  ! Federico: The following equations for Q_loc(i) should be compared with:
  ! - Ruffert, Janka, Schafer (1995) - Eq. (B16)
  !           Also check:
  ! - Galeazzi, Kastaun, Rezzolla, Font, PRD (2013) - Eq. (A13), (A16) and (A18)
  ! NOTE: In this last Eq. (A18) there is a typo because a factor 1/2 is missing.

  ! Leo says: Fermi-Dirac integrals in
  ! Eq. (A15) in https://arxiv.org/abs/astro-ph/0302301
  fermi_fac = (get_fermi_integral_leak(4,eta_e)/ &
       (get_fermi_integral_leak(3,eta_e)+verysmall) + &
       get_fermi_integral_leak(4,-eta_e)/&
       (get_fermi_integral_leak(3,-eta_e)+verysmall))

  ! Eq. (B9) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  block_factor_e = 1.0d0+exp(eta_nue-0.5d0*fermi_fac)
  block_factor_a = 1.0d0+exp(eta_nua-0.5d0*fermi_fac)
  block_factor_x = 1.0d0+exp(eta_nux-0.5d0*fermi_fac)

  ! Leo says: Electron energy moment
  ! Eq. (B5) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  enr_m = 8.0d0*pi/hc_mevcm**3*temp**4*get_fermi_integral_leak(3,eta_e)

  ! Leo says: Positron energy moment
  ! Eq. (B5) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  enr_p = 8.0d0*pi/hc_mevcm**3*temp**4*get_fermi_integral_leak(3,-eta_e)

  ! Leo says: "Tilded" electron energy moment
  ! Eq. (B6) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  enr_tilde_m = 8.0d0*pi/hc_mevcm**3*temp**5*get_fermi_integral_leak(4,eta_e)

  ! Leo says: "Tilded" positron energy moment
  ! Eq. (B6) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  enr_tilde_p = 8.0d0*pi/hc_mevcm**3*temp**5*get_fermi_integral_leak(4,-eta_e)

  ! Leo says: constants appearing in
  ! Eq. (B8) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  pair_const = sigma_0*clite/me_mev**2*enr_m*enr_p

  ! Leo says: Electron neutrino/antineutrino emission rate
  !           due to electron-positron pair annihilation
  ! Eq. (B8) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  R_pair =  pair_const/(36.0d0*block_factor_e*block_factor_a)* &
       ((Cv-Ca)**2+(Cv+Ca)**2) ! Leo says: (C_{1}+C_{2})

  ! Leo says: total electron neutrino emission rate is the sum of all emission rates
  R_loc(1) = R_loc(1) + R_pair

  ! Leo says: total electron antineutrino emission rate is the sum of all emission rates
  R_loc(2) = R_loc(2) + R_pair

  ! Leo says: total electron neutrino cooling rate is the sum of all cooling rates
  ! This uses Eq. (B16) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  Q_loc(1) = Q_loc(1) + R_pair*0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)

  ! Leo says: total electron antineutrino cooling rate is the sum of all cooling rates
  ! This uses Eq. (B16) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  Q_loc(2) = Q_loc(2) + R_pair*0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)

  ! Leo says: Spritz MOD by LSADD
  if(isnan(R_loc(1)).or.isnan(R_loc(2))) then
     write(*,*) "Problem in e-e+ processes, R_pair: ",R_pair
  endif
  ! Leo says: end of Spritz MOD by LSADD

  ! Leo says: Heavy lepton neutrino/antineutrino emission rate
  !           due to electron-positron pair annihilation
  ! Eq. (B10) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  R_pair =  pair_const/(9.0d0*block_factor_x**2)*((Cv-Ca)**2+(Cv+Ca-2.0d0)**2)

  ! Leo says: total heavy lepton neutrino emission rate is the sum of all emission rates
  R_loc(3) = R_loc(3) + R_pair

  ! Leo says: total heavy lepton neutrino cooling rate is the sum of all cooling rates
  ! This uses Eq. (B16) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  Q_loc(3) = Q_loc(3) + R_pair*0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)

  ! Leo says: Spritz MOD by LSADD
  if(isnan(R_loc(3))) then
     write(*,*) "Problem in e-e+ processes for nux, R_pair: ",R_pair
  endif
  ! Leo says: end of Spritz MOD by LSADD

  ! Leo says: useful Spritz comment
  ! plasmon decay from Ruffert et al.
  ! Federico: The equations for the following Q_loc(i) should be compared with:
  ! - Ruffert, Janka, Schafer (1995) - Eq. (B17)
  ! - Galeazzi, Kastaun, Rezzolla, Font, PRD (2013) - Eq. (A25)

  ! Leo says: Eq. (A23) in https://arxiv.org/abs/1306.4953, implemented
  ! using the notation of https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  ! (see paragraph below Eq. B12).
  gamma = gamma_0*sqrt((pi**2+3.0d0*eta_e**2)/3.0d0)

  ! Leo says: Eq. (B13) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  block_factor_e = 1.0d0 + exp(eta_nue-(1.0d0+0.5d0*gamma**2/(1.0d0+gamma)))
  block_factor_a = 1.0d0 + exp(eta_nua-(1.0d0+0.5d0*gamma**2/(1.0d0+gamma)))
  block_factor_x = 1.0d0 + exp(eta_nux-(1.0d0+0.5d0*gamma**2/(1.0d0+gamma)))

  ! Leo says: Constants appearing in Eq. (B12) in
  ! https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  gamma_const = pi**3*sigma_0*clite*temp**8/(me_mev**2*3.0d0*fsc*hc_mevcm**6)* &
       gamma**6*exp(-gamma)*(1.0d0+gamma)

  ! Leo says: Electron neutrino/antineutrino cooling rate due to plasmon decay
  ! Eq. (B11) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  R_gamma = Cv**2*gamma_const/(block_factor_e*block_factor_a)

  ! Leo says: total electron neutrino emission rate is the sum of all emission rates
  R_loc(1) = R_loc(1) + R_gamma

  ! Leo says: total electron antineutrino emission rate is the sum of all emission rates
  R_loc(2) = R_loc(2) + R_gamma

  ! Leo says: total electron neutrino cooling rate is the sum of all cooling rates
  ! This uses Eq. (B17) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  Q_loc(1) = Q_loc(1) + R_gamma*0.5d0*temp*(2.0d0+gamma**2/(1.0d0+gamma))

  ! Leo says: total electron antineutrino cooling rate is the sum of all cooling rates
  ! This uses Eq. (B17) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  Q_loc(2) = Q_loc(2) + R_gamma*0.5d0*temp*(2.0d0+gamma**2/(1.0d0+gamma))

  ! Leo says: Heavy nuclei neutrino/antineutrino cooling rate due to plasmon decay
  ! Eq. (B12) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  R_gamma = (Cv-1.0d0)**2*4.0d0*gamma_const/block_factor_x**2

  ! Leo says: total heavy lepton neutrino emission rate is the sum of all emission rates
  R_loc(3) = R_loc(3) + R_gamma

  ! Leo says: total heavy lepton neutrino cooling rate is the sum of all cooling rates
  ! This uses Eq. (B17) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  Q_loc(3) = Q_loc(3) + R_gamma*0.5d0*temp*(2.0d0+gamma**2/(1.0d0+gamma))

  ! Leo says: useful Spritz comment
  ! Nucleon-nucleon bremsstrahlung, Burrows et al.
  ! Federico: The following equations for Q_loc(i) should be compared with:
  ! - Burrows, Reddy, Thomphson, Nucl. Phys. A (2006) - Eq. (142) and the following text.

  ! Leo says: Eq. (2.1) in https://pos.sissa.it/100/154/pdf
  ! TODO: why 0.231?
  R_pair = 0.231d0*(2.0778d2/mev_to_erg)*0.5d0* &
       (matter_xn**2+matter_xp**2+28.0d0/3.0d0*matter_xn*matter_xp)* &
       rho**2*temp**(4.5d0)
  ! Leo says: Spritz MOD by LSADD
  if(temp < 0) then
     R_pair = 0
  endif
  ! Leo says: end of Spritz MOD by LSADD

  ! Leo says: useful Spritz comment:
  ! LSADD: probably from <omega> = 4.364 T, Burrows et al.
  ! Leo says: Eq. (2.1) in https://pos.sissa.it/100/154/pdf
  ! TODO: why 0.504?
  Q_pair = R_pair*temp/0.231d0*0.504d0

  ! Leo says: total electron neutrino emission/cooling rates
  !           are the sum of all emission/cooling rates
  R_loc(1) = R_loc(1) + R_pair
  Q_loc(1) = Q_loc(1) + Q_pair

  ! Leo says: total electron antineutrino emission/cooling rates
  !           are the sum of all emission/cooling rates
  R_loc(2) = R_loc(2) + R_pair
  Q_loc(2) = Q_loc(2) + Q_pair

  ! Leo says: total heavy nuclei neutrino emission/cooling rates
  !           are the sum of all emission/cooling rates
  !           The factor of 4 appears to be related to g_{nu}, as described below
  !           Eq. (B19) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  R_loc(3) = R_loc(3) + 4.0d0*R_pair
  Q_loc(3) = Q_loc(3) + 4.0d0*Q_pair

  ! Leo says: Spritz MOD by LSADD
  if((isnan(R_loc(1)).or.isnan(R_loc(2)).or.isnan(R_loc(3))).and.reflevel.ge.0) then
     write(*,*) "Problem in plasmon decay processes, gamma_const,gamma,R_gamma,R_pair: ",gamma_const,gamma,R_gamma,R_pair
  endif
  ! Leo says: end of Spritz MOD by LSADD

  ! Leo says: useful Spritz comment
  ! now calculate leakage!!!!!!!
  ! Federico: Check the following references for Q_eff:
  ! - Rosswog, Liebendorfer, mnras (2003) - Eq. (A2)
  ! - Galeazzi, Kastaun, Rezzolla, Font, PRD (2013) - Eq. (62)
  ! - O'Connor, Ott (2010) - Eq (25).

  ! Leo says: Eq. (A1) in https://arxiv.org/abs/astro-ph/0302301
  R_eff(:) = R_loc(:)/(1.0d0+R_loc(:)/R_diff(:))

  ! Leo says: Eq. (A2) in https://arxiv.org/abs/astro-ph/0302301
  Q_eff(:) = Q_loc(:)/(1.0d0+Q_loc(:)/Q_diff(:))

  ! Leo says: Eq. (A11) in https://arxiv.org/abs/astro-ph/0302301
  eave(1:3) = Q_eff(1:3)/R_eff(1:3)

  if (do_heat.ne.0) then
     ! Leo says: useful Spritz comment
     ! LSMOD: heatflux includes luminosity from inner radius, shell volume and lepton blocking factor.

     ! Leo says: Constants in Eq. (31) of https://arxiv.org/abs/0912.2393
     heat_const = f_heat*(1.0d0+3.0d0*alpha**2) * sigma_0 / (me_mev**2 * massn_cgs) * 0.25d0 !cm^2/MeV^2/g

     ! Leo says: See paragraph below Eq. (31) of https://arxiv.org/abs/0912.2393
     F(1:2) = (4.275d0*tau(1:2)+1.15d0)*exp(-2.0d0*tau(1:2))

     ! Leo says: Eq. (31) of https://arxiv.org/abs/0912.2393 (i=n)
     heat_eff(1) = heat_const * rho * matter_xn * heatflux(1) * &
          heaterms(1)**2 * F(1) / mev_to_erg ! cm^2/MeV^2/g * g/cm^3 erg/s/cm^2 MeV^2 MeV/erg = MeV/cm^3/s

     ! Leo says: Eq. (31) of https://arxiv.org/abs/0912.2393 (i=p)
     heat_eff(2) = heat_const * rho * matter_xp * heatflux(2) * &
          heaterms(2)**2 * F(2) / mev_to_erg ! cm^2/MeV^2/g * g/cm^3 erg/s/cm^2 MeV^2 MeV/erg = MeV/cm^3/s

     ! Leo says: Electron neutrino/antineutrino luminosity:
     !           leakage related integrand of Eq. (32)
     !           of https://arxiv.org/abs/0912.2393
     lum(1:2) = (Q_eff(1:2)-heat_eff(1:2))*mev_to_erg !ergs/cm^3/s

     ! Leo says: Heavy lepton neutrino luminosity:
     !           leakage related integrand of Eq. (32)
     !           of https://arxiv.org/abs/0912.2393
     lum(3) = Q_eff(3)*mev_to_erg !ergs/cm^3/s

     ! Leo says: Heating due to electron neutrinos/antineutrinos
     heatout(1:2) = heat_eff(1:2)*mev_to_erg !ergs/cm^3/s

     ! Leo says: Heating due to heavy lepton neutrinos
     heatout(3) = 0.0d0

     netheatout(1) = abs(min(0.0d0,lum(1)))
     netheatout(2) = abs(min(0.0d0,lum(2)))
     netheatout(3) = 0.0d0

     ! Leo says: Energy source term. See e.g.,
     ! Eqs. (21) and (22) of https://arxiv.org/abs/2012.10174
     depsdt = -sum(lum(1:3))/rho

     ! Leo says: Electron fraction source term. See e.g.,
     ! Eqs. (19) and (20) of https://arxiv.org/abs/2012.10174
     dyedt = -(R_eff(1)-R_eff(2)+heat_eff(2)/heateave(2)-heat_eff(1)/heateave(1))*massn_cgs/rho
     ! debug:
     !     dyedt = -(R_eff(1)-R_eff(2))*massn_cgs/rho

  else
     lum(1:3) = Q_eff(1:3)*mev_to_erg

     ! Leo says: Energy source term without heating
     depsdt = -sum(Q_eff(1:3))*mev_to_erg/rho

     ! Leo says: Electron fraction source term without heating
     dyedt = -(R_eff(1)-R_eff(2))*massn_cgs/rho
  endif

  ! Leo says: Spritz addition
  !if (lum(1).ge.1.0d52.and.reflevel.ge.Spritz_c2p_warn_from_reflevel) then
  ! Federico: we are deciding the maximum luminosity to be computed via the parameter ZelmaniLeak::maximum_luminosity
  !           which is by default a huge value (1.0e60 [erg/s])
  if (lum(1).ge.maximum_luminosity.and.reflevel.ge.0) then
     call CCTK_WARN(1,"In my opinion, luminosity is too high. Setting to zero.")
     write(warnline,"(A15,i10)") "reflevel: ", reflevel
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "rho: ", rho
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "temp: ", temp
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "Y_e: ", ye
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "dyedt*ldt ", dyedt*ldt
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "X_n,X_p,X_a,X_h: ", matter_xn, matter_xp, &
          matter_xa, matter_xh
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "eps: ", matter_enr
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "depsdt*ldt:", depsdt*ldt
     call CCTK_WARN(1,warnline)
     write(warnline,"(A30,1P10E15.6)") "matter_enr+depsdt*ldt:", matter_enr+depsdt*ldt
     call CCTK_WARN(1,warnline)
     write(warnline,"(A30,1P10E15.6)") "radius:", rad
     call CCTK_WARN(1,warnline)
     write(warnline,"(A30,1P10E15.6)") "R_pair,R_gamma,pair_const,zeta(1):", R_pair,R_gamma,pair_const,zeta(1)
     call CCTK_WARN(1,warnline)
     write(warnline,"(A30,1P10E15.6)") "Q_loc(1),Q_diff(1):", Q_loc(1),Q_diff(1)
     call CCTK_WARN(1,warnline)
     write(warnline,"(A30,1P10E15.6)") "rate_const,chi(1):", rate_const,chi(1)
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "lum(1:3): ", lum(1), lum(2), lum(3)
     call CCTK_WARN(1,warnline)
     lum(1) = 0.0d0
     lum(2) = 0.0d0
     lum(3) = 0.0d0
     depsdt = 0.0d0
     dyedt  = 0.0d0
     eave(1) = 0.0d0
     eave(2) = 0.0d0
     eave(3) = 0.0d0
     heatout(1) = 0.0d0
     heatout(2) = 0.0d0
     heatout(3) = 0.0d0
     netheatout(1) = 0.0d0
     netheatout(2) = 0.0d0
     netheatout(3) = 0.0d0
  endif

  if (ye+dyedt*ldt.le.eos_yemin*1.05d0) then
     dyedt = 0.0d0
     depsdt = 0.0d0
  endif

  if (ye+dyedt*ldt.ge.eos_yemax*0.95d0) then
     dyedt = 0.0d0
     depsdt = 0.0d0
  endif


  if (matter_enr+depsdt*ldt.lt.-energy_shift*inv_eps_gf.and.&
       (reflevel.eq.-1.or.reflevel.ge.0)) then
     !$OMP CRITICAL
     call CCTK_WARN(1,"Problem in leakage; energy change too large")
     write(warnline,"(A15,i10)") "reflevel: ", reflevel
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "rho: ", rho
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "temp: ", temp
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "Y_e: ", ye
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "dyedt*ldt ", dyedt*ldt
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "X_n,X_p,X_a,X_h: ", matter_xn, matter_xp, &
          matter_xa, matter_xh
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "eps: ", matter_enr
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "depsdt*ldt:", depsdt*ldt
     call CCTK_WARN(1,warnline)
     write(warnline,"(A30,1P10E15.6)") "matter_enr+depsdt*ldt:", matter_enr+depsdt*ldt
     call CCTK_WARN(1,warnline)
     write(warnline,"(A30,1P10E15.6)") "radius:", rad
     call CCTK_WARN(1,warnline)
     write(warnline,"(A30,1P10E15.6)") "R_pair,R_gamma,pair_const,zeta(1):", R_pair,R_gamma,pair_const,zeta(1)
     call CCTK_WARN(1,warnline)
     call CCTK_WARN(0,"aborting")
     !$OMP END CRITICAL
  endif

  if ((isnan(lum(1)).or.isnan(lum(2)).or.isnan(lum(3))).and.reflevel.ge.0) then
     !$OMP CRITICAL
     call CCTK_WARN(1,"Problem in leakage; NaN in luminosity")
     write(warnline,"(A15,i10)") "reflevel: ", reflevel
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "rho: ", rho
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "temp: ", temp
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "Y_e: ", ye
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "dyedt*ldt ", dyedt*ldt
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "X_n,X_p,X_a,X_h: ", matter_xn, matter_xp, &
          matter_xa, matter_xh
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "eps: ", matter_enr
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "depsdt*ldt:", depsdt*ldt
     call CCTK_WARN(1,warnline)
     write(warnline,"(A30,1P10E15.6)") "matter_enr+depsdt*ldt:", matter_enr+depsdt*ldt
     call CCTK_WARN(1,warnline)
     write(warnline,"(A30,1P10E15.6)") "radius:", rad
     call CCTK_WARN(1,warnline)
     write(warnline,"(A30,1P10E15.6)") "R_pair,R_gamma,pair_const,zeta(1):", R_pair,R_gamma,pair_const,zeta(1)
     call CCTK_WARN(1,warnline)
     write(warnline,"(A30,1P10E15.6)") "Q_loc(1),Q_diff(1):", Q_loc(1),Q_diff(1)
     call CCTK_WARN(1,warnline)
     write(warnline,"(A30,1P10E15.6)") "rate_const,chi(1):", rate_const,chi(1)
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "lum(1:3): ", lum(1), lum(2), lum(3)
     call CCTK_WARN(1,warnline)
     call CCTK_WARN(0,"aborting")
     !$OMP END CRITICAL
  endif

  if (lum(1).ge.1.0d52.and.reflevel.ge.0) then
     !$OMP CRITICAL
     call CCTK_WARN(1,"In my opinion, luminosity is too high.")
     write(warnline,"(A15,i10)") "reflevel: ", reflevel
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "rho: ", rho
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "temp: ", temp
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "Y_e: ", ye
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "dyedt*ldt ", dyedt*ldt
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "X_n,X_p,X_a,X_h: ", matter_xn, matter_xp, &
          matter_xa, matter_xh
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "eps: ", matter_enr
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "depsdt*ldt:", depsdt*ldt
     call CCTK_WARN(1,warnline)
     write(warnline,"(A30,1P10E15.6)") "matter_enr+depsdt*ldt:", matter_enr+depsdt*ldt
     call CCTK_WARN(1,warnline)
     write(warnline,"(A30,1P10E15.6)") "radius:", rad
     call CCTK_WARN(1,warnline)
     write(warnline,"(A30,1P10E15.6)") "R_pair,R_gamma,pair_const,zeta(1):", R_pair,R_gamma,pair_const,zeta(1)
     call CCTK_WARN(1,warnline)
     write(warnline,"(A30,1P10E15.6)") "Q_loc(1),Q_diff(1):", Q_loc(1),Q_diff(1)
     call CCTK_WARN(1,warnline)
     write(warnline,"(A30,1P10E15.6)") "rate_const,chi(1):", rate_const,chi(1)
     call CCTK_WARN(1,warnline)
     write(warnline,"(A15,1P10E15.6)") "lum(1:3): ", lum(1), lum(2), lum(3)
     call CCTK_WARN(1,warnline)
     call CCTK_WARN(0,"aborting")
     !$OMP END CRITICAL
  endif

end subroutine calc_leak
