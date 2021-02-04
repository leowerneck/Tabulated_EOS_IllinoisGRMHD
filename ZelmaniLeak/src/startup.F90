#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


subroutine ZelmaniLeak_startup(CCTK_ARGUMENTS)

  use EOS_Omni_module
#ifdef HAVE_CAPABILITY_Fortran
  use cctk
#endif
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS


!# bookkeeping turned off for performance
!#  energycheck = 0.0d0

  if(do_pnu.ne.0) then
     pnu = 0.0d0
  endif


  
end subroutine ZelmaniLeak_startup

subroutine ZelmaniLeak_startup_global(CCTK_ARGUMENTS)

  use yeofrhomod
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  character(len=1024) profilefilename
  integer length
  length = 1024

  if(do_ye_of_rho_from_profile.ne.0) then
     ! don't try doing this more than once
     if(allocated(rhotable)) return

     call CCTK_INFO("Reading Y_e(rho) profile!")
     ! setting module variables                                                                                     
      density_factor = rho_gf
      nzones = Zones
      rho_min = 5.0d5
      rho_max = 5.0d14
      ! allocating memory for profile file data   
      allocate(prhoa(Profile_Zones))
      allocate(pyea(Profile_zones))
      allocate(yetable(nzones))
      allocate(rhotable(nzones))
      ! Let's figure out the profile filename                                                                        
      call CCTK_FortranString(length,Profile_File,profilefilename)
      call readprofile_ye(profilefilename,Profile_Zones,prhoa,pyea)
      call setuprho_ye(nzones,rho_min,rho_max,rhotable)
      call setupye_ye(nzones,Profile_Zones,rhotable,yetable,prhoa,pyea)

  endif

  if(rho_abs_min .lt. 0.0) then
     call CCTK_WARN(0,"Must set GRHydro::rho_abs_min!")
  endif



end subroutine ZelmaniLeak_startup_global
