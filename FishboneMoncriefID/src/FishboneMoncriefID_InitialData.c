#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "FishboneMoncriefID.h"

void FishboneMoncriefID_InitialData(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Leo says: get EOS key based on new parameter eos_type
  int eos_key;
  if( CCTK_Equals(eos_type,"polytrope") ) {
    eos_key = POLYTROPE_EOS;
  }
  else if( CCTK_Equals(eos_type,"tabulated") ) {
    eos_key = TABULATED_EOS;
  }
  else {
    CCTK_VERROR("Unsupported value for parameter eos_type: %s",eos_type);
  }

  if( verbosity_level > 0 ) {
    CCTK_VINFO("Fishbone-Moncrief Disk Initial data.");
    CCTK_VINFO("Input parameters:");
    CCTK_VINFO("  - a                = %e",a);
    CCTK_VINFO("  - M                = %e",M);
    CCTK_VINFO("  - r_in             = %e",r_in);
    CCTK_VINFO("  - r_at_max_density = %e",r_at_max_density);
  }

  if( eos_key == POLYTROPE_EOS ) {
    if( verbosity_level > 0 ) CCTK_VINFO("  - EOS type: polytrope - kappa = %e | gamma = %e",kappa,gamma);
    FishboneMoncriefID_initial_polytrope_eos(CCTK_PASS_CTOC);
  }
  else if( eos_key == TABULATED_EOS ) {
    if( verbosity_level > 0 ) CCTK_VINFO("  - EOS type: tabulated");
    if( !CCTK_IsThornActive("EOS_Omni") ) CCTK_ERROR("Cannot set tabulated EOS initial data without EOS_Omni. Please activate EOS_Omni and set the nuc_eos parameters in your parfile.");
    FishboneMoncriefID_initial_tabulated_eos(CCTK_PASS_CTOC);
  }

}
