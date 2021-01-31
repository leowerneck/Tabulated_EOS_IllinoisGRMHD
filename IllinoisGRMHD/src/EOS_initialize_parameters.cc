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

  DECLARE_CCTK_PARAMETERS;
  
  // Set the EOS key
  eos.key           = *igm_eos_key;

  // Get con2prim key
  eos.c2p_routine   = con2prim_get_key(igm_con2prim_routine);

  // Get con2prim backup keys
  eos.c2p_backup[0] = con2prim_get_key(igm_con2prim_backup_routine[0]);
  eos.c2p_backup[1] = con2prim_get_key(igm_con2prim_backup_routine[1]);
  eos.c2p_backup[2] = con2prim_get_key(igm_con2prim_backup_routine[2]);

  // Whether or not to evolve the entropy
  eos.evolve_entropy = igm_evolve_entropy;

  // Maximum Lorentz factor
  eos.W_max = GAMMA_SPEED_LIMIT;
  // Inverse of Lorentz factor squared (used by Palenzuela con2prim routines)
  eos.inv_W_max_squared = 1.0/(SQR(eos.W_max));

  // EOS specific initialization
  eos.is_Hybrid      = false;
  eos.is_Tabulated   = false;
  if( eos.key == EOS_Omni_GetHandle("Hybrid") ) {
    eos.is_Hybrid = true;
    initialize_Hybrid_EOS_parameters_from_input(eos);
  }
  else if( eos.key == EOS_Omni_GetHandle("nuc_eos") ) {
    eos.is_Tabulated = true;
    initialize_Tabulated_EOS_parameters_from_input(cctk_time,eos);
  }
  else {
    CCTK_VError(VERR_DEF_PARAMS,"Unknown EOS key: %d. ABORTING!",eos.key);
  }
  
}
