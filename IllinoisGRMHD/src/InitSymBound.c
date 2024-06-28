/*
  Set the symmetries for the IllinoisGRMHD variables
*/

#include "IllinoisGRMHD.h"
#include "Symmetry.h"

void IllinoisGRMHD_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_IllinoisGRMHD_InitSymBound;
  DECLARE_CCTK_PARAMETERS;

  if( ( CCTK_EQUALS(Matter_BC,"frozen") && !CCTK_EQUALS(EM_BC,"frozen") ) ||
      ( !CCTK_EQUALS(Matter_BC,"frozen") && CCTK_EQUALS(EM_BC,"frozen") ) )
    CCTK_VERROR("If Matter_BC or EM_BC is set to FROZEN, BOTH must be set to frozen!");

  if ((cctk_nghostzones[0]<3 || cctk_nghostzones[1]<3 || cctk_nghostzones[2]<3))
    CCTK_VERROR("ERROR: The version of PPM in this thorn requires 3 ghostzones. You only have (%d,%d,%d) ghostzones!",cctk_nghostzones[0],cctk_nghostzones[1],cctk_nghostzones[2]);

  if(cctk_iteration==0) {
    CCTK_VINFO("Setting Symmetry = %s... at iteration = %d",Symmetry,cctk_iteration);

    int sym[3];

    sym[0] = 1; sym[1] = 1; sym[2] = 1;
    SetCartSymGN(cctkGH, sym, "IllinoisGRMHD::grmhd_conservatives");
    SetCartSymGN(cctkGH, sym, "IllinoisGRMHD::grmhd_primitives_allbutBi");
    SetCartSymGN(cctkGH, sym, "IllinoisGRMHD::grmhd_B_center");
    SetCartSymGN(cctkGH, sym, "IllinoisGRMHD::Ax");
    SetCartSymGN(cctkGH, sym, "IllinoisGRMHD::Ay");
    SetCartSymGN(cctkGH, sym, "IllinoisGRMHD::Az");
    SetCartSymGN(cctkGH, sym, "IllinoisGRMHD::phitilde");
    if(ghl_params->evolve_entropy)
      SetCartSymGN(cctkGH, sym, "IllinoisGRMHD::ent_star");
    if(ghl_eos->eos_type == ghl_eos_tabulated)
      SetCartSymGN(cctkGH, sym, "IllinoisGRMHD::Ye_star");

    if(CCTK_EQUALS(Symmetry, "equatorial")) {
      sym[2] = -Sym_Bz;
      SetCartSymVN(cctkGH, sym, "IllinoisGRMHD::Bx_center");
      SetCartSymVN(cctkGH, sym, "IllinoisGRMHD::By_center");

      sym[2] = Sym_Bz;
      SetCartSymVN(cctkGH, sym, "IllinoisGRMHD::Bz_center");

      sym[2] = -1;
      SetCartSymVN(cctkGH, sym, "IllinoisGRMHD::Stildez");
      SetCartSymVN(cctkGH, sym, "IllinoisGRMHD::vz");
    } else if (!CCTK_EQUALS(Symmetry, "none")) {
      CCTK_ERROR("IllinoisGRMHD_initsymbound: Should not be here; picked an impossible symmetry.");
    }
  }
}
