// Thorn      : IllinoisGRMHD
// File       : EOS_get_key.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: In this file we set up a function which
//              retrieves the appropriate EOS key/handle
//              from the EOS_Omni thorn, setting the
//              igm_eos_key parameter appropriately.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void IllinoisGRMHD_EOS_get_key( CCTK_ARGUMENTS ) {

  // This gives us the pointer igm_eos_key
  DECLARE_CCTK_ARGUMENTS_IllinoisGRMHD_EOS_get_key;
  // This gives us the parameter igm_eos_type
  DECLARE_CCTK_PARAMETERS;

  // Now check what we have
  if( CCTK_EQUALS(igm_eos_type,"Hybrid") ) {
    *igm_eos_key = EOS_Omni_GetHandle("Hybrid");
    CCTK_VInfo(CCTK_THORNSTRING,"Hybrid EOS selected. igm_eos_key set to %d",*igm_eos_key);
  }
  else if( CCTK_EQUALS(igm_eos_type,"Tabulated") || CCTK_EQUALS(igm_eos_type,"nuc_eos") ) {
    *igm_eos_key = EOS_Omni_GetHandle("nuc_eos");
    CCTK_VInfo(CCTK_THORNSTRING,"Tabulated (nuc_eos) EOS selected. igm_eos_key set to %d",*igm_eos_key);
  }
  else {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Unknown EOS type: %s. ABORTING!",igm_eos_type);
  }

}
