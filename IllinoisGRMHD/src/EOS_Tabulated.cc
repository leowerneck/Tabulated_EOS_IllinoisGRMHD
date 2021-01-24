// Thorn      : IllinoisGRMHD
// File       : EOS_Tabulated.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: In this file we provide tabulated EOS functions.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "IllinoisGRMHD_headers.h"

//----------- EOS_Omni stuff -----------
// The following extern variables are declared in
// the EOS_Omni thorn file readtable.cc. They are
// initialized by the function nuc_eos_C_ReadTable(),
// in the same file.
namespace nuc_eos {
  extern double eos_rhomin,eos_yemin,eos_tempmin;
  extern double eos_rhomax,eos_yemax,eos_tempmax;
  extern double energy_shift;
}
namespace nuc_eos_private {
  extern double nrho,nye,ntemp;
  extern double *restrict alltables;
}

// This is defined in the nuc_eos.hh file
#ifndef NTABLES
#define NTABLES (19)
#endif

// These are useful for us, though not defined
// explicitly by EOS_Omni
const CCTK_INT table_key_pressure = 0;
const CCTK_INT table_key_epsilon  = 1;
const CCTK_INT table_key_entropy  = 2;

// These are also very useful
const CCTK_INT have_eps  = 0;
const CCTK_INT have_temp = 1;
const CCTK_INT have_ent  = 2;
//--------------------------------------

// EOS_Omni does not provide functions to obtain the
// maximum table value of e.g. the specific internal
// energy. This function does that.
CCTK_REAL get_EOS_table_max( const int which_var ) {

  // It's simply too annoying to keep appending the
  // name to all the variables. This namespace is
  // also very small, so we should be okay.
  using namespace nuc_eos_private;

  // Loop over the table, searching for the maximum value
  CCTK_INT  totalsize     = nrho * nye * ntemp;
  CCTK_REAL var_max_value = alltables[which_var];

  for(int i=0;i<totalsize;i++) {
    CCTK_REAL var_aux = alltables[which_var + NTABLES*i];
    if( var_aux > var_max_value ) var_max_value = var_aux;
  }
  return var_max_value;
  
}

// EOS_Omni does not provide functions to obtain the
// maximum table value of e.g. the specific internal
// energy. This function does that.
CCTK_REAL get_EOS_table_min( const int which_var ) {

  // It's simply too annoying to keep appending the
  // name to all the variables. This namespace is
  // also very small, so we should be okay.
  using namespace nuc_eos_private;

  // Loop over the table, searching for the maximum value
  CCTK_INT  totalsize     = nrho * nye * ntemp;
  CCTK_REAL var_min_value = alltables[which_var];

  for(int i=0;i<totalsize;i++) {
    CCTK_REAL var_aux = alltables[which_var + NTABLES*i];
    if( var_aux < var_min_value ) var_min_value = var_aux;
  }
  return var_min_value;
  
}

void initialize_Tabulated_EOS_parameters_from_input( const CCTK_REAL cctk_time,igm_eos_parameters& eos ) {

  DECLARE_CCTK_PARAMETERS;

  // It's simply too annoying to keep appending the
  // name to all the variables. This namespace is
  // also very small, so we should be okay.
  using namespace nuc_eos;

  // Which variable do we want to reconstruct during PPM?
  // The tabulated EOS case actually allows for all the
  // possibilities.
  if( CCTK_EQUALS(igm_PPM_reconstructed_variable,"pressure") ) {
    eos.PPM_reconstructed_var = PRESSURE;
  }
  else if( CCTK_EQUALS(igm_PPM_reconstructed_variable,"epsilon") ) {
    eos.PPM_reconstructed_var = EPSILON;
  }
  else if( CCTK_EQUALS(igm_PPM_reconstructed_variable,"entropy") ) {
    eos.PPM_reconstructed_var = ENTROPY;
  }
  else {
    CCTK_VError(VERR_DEF_PARAMS,
                "PPM reconstruction of variable \"%s\" not supported with Tabulated EOS. "
                "Can only reconstruct: pressure, epsilon, entropy. ABORTING!",
                igm_PPM_reconstructed_variable);
  }

  // Initialize tabulated EOS parameters

  // Whether or not to evolve the entropy
  eos.evolve_entropy = igm_evolve_entropy;

  // Check if it's time to begin temperature evolution
  if( igm_evolve_temperature && (cctk_time >= igm_freeze_T_evolution_until_cctk_time) ) {
    eos.evolve_T = true;
  }
  else {
    eos.evolve_T = false;
  }

  // Root-finding precision (for table inversions)
  eos.root_finding_precision = igm_eos_root_finding_precision;

  // --------- Atmospheric values ---------
  // Atmospheric rho
  eos.rho_atm = rho_b_atm;
  // Atmospheric electron fraction
  eos.Ye_atm  = igm_Ye_atm;
  // Atmospheric temperature fraction
  eos.T_atm   = igm_T_atm;
  // Compute P, eps, and S in the atmosphere
  if( eos.evolve_entropy ) {
    get_P_eps_and_S_from_rho_Ye_and_T( eos,
                                       eos.rho_atm,eos.Ye_atm,eos.T_atm,
                                       &eos.P_atm,&eos.eps_atm,&eos.S_atm );
  }
  else {
    get_P_and_eps_from_rho_Ye_and_T( eos,
                                     eos.rho_atm,eos.Ye_atm,eos.T_atm,
                                     &eos.P_atm,&eos.eps_atm );
  }
  // Atmospheric tau
  eos.tau_atm = tau_atm;
  // --------------------------------------

  // -------------- Ceilings --------------
  // Get the maximum pressure. Remember that the alltables
  // array actually constains ln(press), so we must adjust
  // appropriately.
  const CCTK_REAL eos_prsmax = exp(get_EOS_table_max( table_key_pressure ));
  // Then get the maximum value of eps. Remember that the 
  // alltables array actually contains ln(eps + energy_shift),
  // so we must adjust appropriately.
  const CCTK_REAL eos_epsmax = exp(get_EOS_table_max( table_key_epsilon )) - energy_shift;
  // Finally, get the maximum entropy
  const CCTK_REAL eos_entmax = get_EOS_table_max( table_key_entropy );
  // Now set the EOS struct variables
  eos.rho_max = eos_rhomax  * igm_eos_table_ceiling_safety_factor;
  eos.Ye_max  = eos_yemax   * igm_eos_table_ceiling_safety_factor;
  eos.T_max   = eos_tempmax * igm_eos_table_ceiling_safety_factor;
  eos.P_max   = eos_prsmax  * igm_eos_table_ceiling_safety_factor;
  eos.eps_max = eos_epsmax  * igm_eos_table_ceiling_safety_factor;
  eos.S_max   = eos_entmax  * igm_eos_table_ceiling_safety_factor;
  eos.W_max   = GAMMA_SPEED_LIMIT;
  // --------------------------------------

  // --------------- Floors ---------------
  // Get the miniimum pressure. Remember that the alltables
  // array actually constains ln(press), so we must adjust
  // appropriately.
  const CCTK_REAL eos_prsmin = exp(get_EOS_table_min( table_key_pressure ));
  // Then get the minimum value of eps. Remember that the 
  // alltables array actually contains ln(eps + energy_shift),
  // so we must adjust appropriately.
  const CCTK_REAL eos_epsmin = exp(get_EOS_table_min( table_key_epsilon )) - energy_shift;
  // Finally, get the minimum entropy
  const CCTK_REAL eos_entmin = get_EOS_table_min( table_key_entropy );
  // Now set the EOS struct variables
  eos.rho_min = eos_rhomin  * igm_eos_table_floor_safety_factor;
  eos.Ye_min  = eos_yemin   * igm_eos_table_floor_safety_factor;
  eos.T_min   = eos_tempmin * igm_eos_table_floor_safety_factor;
  eos.P_min   = eos_prsmin  * igm_eos_table_floor_safety_factor;
  eos.eps_min = eos_epsmin  * igm_eos_table_floor_safety_factor;
  eos.S_min   = eos_entmin  * igm_eos_table_floor_safety_factor;
  // --------------------------------------

  // All done!

}

void get_P_and_eps_from_rho_Ye_and_T( const igm_eos_parameters eos,
                                      const CCTK_REAL rho,
                                      const CCTK_REAL Y_e,
                                      const CCTK_REAL T,
                                      CCTK_REAL *restrict P,
                                      CCTK_REAL *restrict eps ) {
  // Set up for table call
  CCTK_INT  npoints  = 1;
  CCTK_INT  keyerr   = 0;
  CCTK_INT  anyerr   = 0;
  // The EOS_Omni function does not like some of its arguments
  // being constant (even when they are not supposed to change).
  // We declare auxiliary variables here to avoid errors.
  CCTK_REAL rho_in   = rho;
  CCTK_REAL Y_e_in   = Y_e;
  CCTK_REAL T_in     = T;
  CCTK_REAL P_out    = 0.0;
  CCTK_REAL eps_out  = 0.0;

  // Perform the table interpolations
  EOS_Omni_press( eos.key,have_temp,eos.root_finding_precision,npoints,
                  &rho_in,&eps_out,&T_in,&Y_e_in,&P_out,
                  &keyerr,&anyerr );

  // Now update P and eps
  *P   = P_out;
  *eps = eps_out;

  // FIXME: add error handling!  
}


void get_P_eps_and_S_from_rho_Ye_and_T( const igm_eos_parameters eos,
                                        const CCTK_REAL rho,
                                        const CCTK_REAL Y_e,
                                        const CCTK_REAL T,
                                        CCTK_REAL *restrict P,
                                        CCTK_REAL *restrict eps,
                                        CCTK_REAL *restrict S ) {
  // Set up for table call
  CCTK_INT  npoints  = 1;
  CCTK_INT  keyerr   = 0;
  CCTK_INT  anyerr   = 0;
  // The EOS_Omni function does not like some of its arguments
  // being constant (even when they are not supposed to change).
  // We declare auxiliary variables here to avoid errors.
  CCTK_REAL rho_in   = rho;
  CCTK_REAL Y_e_in   = Y_e;
  CCTK_REAL T_in     = T;
  CCTK_REAL P_out    = 0.0;
  CCTK_REAL eps_out  = 0.0;
  CCTK_REAL S_out    = 0.0;
  CCTK_REAL dummy    = 0.0;

  // Perform the table interpolations
  EOS_Omni_short( eos.key,have_temp,eos.root_finding_precision,npoints,
                  &rho_in,&eps_out,&T_in,&Y_e_in,&P_out,&S_out,
                  &dummy,&dummy,&dummy,&dummy,&dummy,
                  &keyerr,&anyerr );

  // Now update P, eps, and S
  *P   = P_out;
  *eps = eps_out;
  *S   = S_out;
  
  // FIXME: add error handling!
}

void get_P_eps_and_T_from_rho_Ye_and_S( const igm_eos_parameters eos,
                                        const CCTK_REAL rho,
                                        const CCTK_REAL Y_e,
                                        const CCTK_REAL S,
                                        CCTK_REAL *restrict P,
                                        CCTK_REAL *restrict eps,
                                        CCTK_REAL *restrict T ) {
  // Set up for table call
  CCTK_INT  npoints  = 1;
  CCTK_INT  keyerr   = 0;
  CCTK_INT  anyerr   = 0;
  // The EOS_Omni function does not like some of its arguments
  // being constant (even when they are not supposed to change).
  // We declare auxiliary variables here to avoid errors.
  CCTK_REAL rho_in   = rho;
  CCTK_REAL Y_e_in   = Y_e;
  CCTK_REAL S_in     = S;
  CCTK_REAL P_out    = 0.0;
  CCTK_REAL eps_out  = 0.0;
  CCTK_REAL T_out    = *T;
  CCTK_REAL dummy    = 0.0;

  // Perform the table interpolations
  EOS_Omni_short( eos.key,have_ent,eos.root_finding_precision,npoints,
                  &rho_in,&eps_out,&T_out,&Y_e_in,&P_out,&S_in,
                  &dummy,&dummy,&dummy,&dummy,&dummy,
                  &keyerr,&anyerr );

  // Now update P, eps, and S
  *P   = P_out;
  *eps = eps_out;
  *T   = T_out;
  
  // FIXME: add error handling!
}
