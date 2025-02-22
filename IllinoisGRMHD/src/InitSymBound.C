/*
  Set the symmetries for the IllinoisGRMHD variables
*/

#include "cctk.h"
#include <cstdio>
#include <cstdlib>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"
#include "IllinoisGRMHD_headers.h"

extern "C" void IllinoisGRMHD_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if( ( CCTK_EQUALS(Matter_BC,"frozen") && !CCTK_EQUALS(EM_BC,"frozen") ) ||
      ( !CCTK_EQUALS(Matter_BC,"frozen") && CCTK_EQUALS(EM_BC,"frozen") ) )
    CCTK_VError(VERR_DEF_PARAMS,"If Matter_BC or EM_BC is set to FROZEN, BOTH must be set to frozen!");

  if ((cctk_nghostzones[0]<3 || cctk_nghostzones[1]<3 || cctk_nghostzones[2]<3))
    CCTK_VError(VERR_DEF_PARAMS,"ERROR: The version of PPM in this thorn requires 3 ghostzones. You only have (%d,%d,%d) ghostzones!",cctk_nghostzones[0],cctk_nghostzones[1],cctk_nghostzones[2]);

  if(cctk_iteration==0) {
    CCTK_VInfo(CCTK_THORNSTRING,"Setting Symmetry = %s... at iteration = %d",Symmetry,cctk_iteration);

    int sym[3];

    if(CCTK_EQUALS(Symmetry,"none")) {
      /* FIRST SET NO SYMMETRY OPTION */
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      SetCartSymGN(cctkGH,sym,"IllinoisGRMHD::grmhd_conservatives");
      SetCartSymGN(cctkGH,sym,"IllinoisGRMHD::A_x_tilde");
      SetCartSymGN(cctkGH,sym,"IllinoisGRMHD::A_y_tilde");
      SetCartSymGN(cctkGH,sym,"IllinoisGRMHD::A_z_tilde");
      SetCartSymGN(cctkGH,sym,"IllinoisGRMHD::Phi_tilde");
      SetCartSymGN(cctkGH,sym,"IllinoisGRMHD::grmhd_velocities");
    } else if(CCTK_EQUALS(Symmetry,"equatorial")) {
      /* THEN SET EQUATORIAL SYMMETRY OPTION */
      // Set default to no symmetry, which is correct for scalars and most vectors:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      SetCartSymGN(cctkGH,sym,"IllinoisGRMHD::grmhd_conservatives");
      // Don't worry about the wrong sym values since A_{\mu} is staggered
      // and we're going to impose the symmetry separately
      SetCartSymGN(cctkGH,sym,"IllinoisGRMHD::A_x_tilde");
      SetCartSymGN(cctkGH,sym,"IllinoisGRMHD::A_y_tilde");
      SetCartSymGN(cctkGH,sym,"IllinoisGRMHD::A_z_tilde");
      SetCartSymGN(cctkGH,sym,"IllinoisGRMHD::Phi_tilde");
      SetCartSymGN(cctkGH,sym,"IllinoisGRMHD::grmhd_velocities");

      // Then set unstaggered B field variables
      sym[2] = -Sym_Bz;
      SetCartSymVN(cctkGH, sym,"IllinoisGRMHD::Bx_center");
      SetCartSymVN(cctkGH, sym,"IllinoisGRMHD::By_center");
      sym[2] = Sym_Bz;
      SetCartSymVN(cctkGH, sym,"IllinoisGRMHD::Bz_center");

      sym[2] = -1;
      SetCartSymVN(cctkGH, sym,"IllinoisGRMHD::S_z_tilde");
      SetCartSymVN(cctkGH, sym,"IllinoisGRMHD::vz");
    } else {
      CCTK_VError(VERR_DEF_PARAMS,"IllinoisGRMHD_initsymbound: Should not be here; picked an impossible symmetry.");
    }
  }
}


