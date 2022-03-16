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
    CCTK_VERROR("Unknown parameter type for eos_type: %s",eos_type);
  }

  CCTK_VINFO("Fishbone-Moncrief Disk Initial data.");
  CCTK_VINFO("Using input parameters of\n a = %e,\n M = %e,\n r_in = %e,\n r_at_max_density = %e",a,M,r_in,r_at_max_density);

  if( eos_key == POLYTROPE_EOS ) {
    CCTK_VINFO("  - kappa = %e\n  - gamma = %e",kappa,gamma);
    FishboneMoncriefID_initial_polytrope_eos(CCTK_PASS_CTOC);
  }
  else if( eos_key == TABULATED_EOS ) {
    CCTK_VINFO("  - Using tabulated equation of state");
    if( !CCTK_IsThornActive("EOS_Omni") ) CCTK_ERROR("Cannot set tabulated EOS initial data without EOS_Omni. Please activate EOS_Omni and set the nuc_eos parameters in your parfile.");
    FishboneMoncriefID_initial_tabulated_eos(CCTK_PASS_CTOC);
  }

}
