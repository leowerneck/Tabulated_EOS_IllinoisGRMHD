// Thorn      : IllinoisGRMHD
// File       : EOS_initialize_parameters.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: In this file we provide a function to initialize
//              the EOS parameters struct based on the parameter
//              file configuration.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "IllinoisGRMHD_headers.h"
#include "EOS_headers.hh"
#include "con2prim_headers.h"

void initialize_igm_eos_parameters_from_input( const CCTK_INT* igm_eos_key,const CCTK_REAL cctk_time,igm_eos_parameters &eos ) {

  // Set the EOS key
  // IllinoisGRMHD_EOS_get_key(&eos.key);
  eos.key         = *igm_eos_key;
  eos.c2p_routine = Noble2D;

  if( eos.key == EOS_Omni_GetHandle("Hybrid") ) {
    initialize_Hybrid_EOS_parameters_from_input(eos);
  }
  else if( eos.key == EOS_Omni_GetHandle("nuc_eos") ) {
    initialize_Tabulated_EOS_parameters_from_input(cctk_time,eos);
  }
  else {
    CCTK_VError(VERR_DEF_PARAMS,"Unknown EOS key: %d. ABORTING!",eos.key);
  }
  
}
