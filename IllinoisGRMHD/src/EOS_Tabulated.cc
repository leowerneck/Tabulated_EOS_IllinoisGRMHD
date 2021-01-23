// Thorn      : IllinoisGRMHD
// File       : EOS_Tabulated.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: In this file we provide tabulated EOS functions.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "IllinoisGRMHD_headers.h"

void initialize_Tabulated_EOS_parameters_from_input( igm_eos_parameters& eos ) {

  DECLARE_CCTK_PARAMETERS;

  // Initialize tabulated EOS parameters

  // Which variable do we want to reconstruct during PPM?
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
                "PPM reconstruction of variable \"%s\" not supported. "
                "Can only reconstruct: pressure, epsilon, entropy. ABORTING!",
                igm_PPM_reconstructed_variable);
  }

  // Whether or not to evolve the entropy
  eos.evolve_entropy = igm_evolve_entropy;

  // Root-finding precision (for table inversions)
  eos.root_finding_precision = igm_eos_root_finding_precision;

  // Atmospheric rho
  eos.rho_b_atm = rho_b_atm;
  // Atmospheric electron fraction
  eos.Ye_atm    = igm_Ye_atm;
  // Atmospheric temperature fraction
  eos.T_atm     = igm_T_atm;

  // Compute P, eps, and S in the atmosphere
  if( eos.evolve_entropy ) {
    get_P_eps_and_S_from_rho_Ye_and_T( eos,
                                       eos.rho_b_atm,eos.Ye_atm,eos.T_atm,
                                       &eos.P_atm,&eos.eps_atm,&eos.S_atm );
  }
  else {
    get_P_and_eps_from_rho_Ye_and_T( eos,
                                     eos.rho_b_atm,eos.Ye_atm,eos.T_atm,
                                     &eos.P_atm,&eos.eps_atm );
  }

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
  CCTK_INT  havetemp = 1;
  CCTK_INT  keyerr   = 0;
  CCTK_INT  anyerr   = 0;
  // The temperature cannot be constant in the following call
  CCTK_REAL T_in     = T;

  // Perform the table interpolations
  EOS_Omni_press( eos.key,havetemp,eos.root_finding_precision,npoints,
                  &rho,eps,&T_in,&Y_e,P,
                  &keyerr,&anyerr );
  
  // TODO: add error handling
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
  CCTK_INT  havetemp = 1;
  CCTK_INT  keyerr   = 0;
  CCTK_INT  anyerr   = 0;
  CCTK_REAL dummy    = 0.0;
  // The temperature cannot be constant in the following call
  CCTK_REAL T_in     = T;

  // Perform the table interpolations
  EOS_Omni_short( eos.key,havetemp,eos.root_finding_precision,npoints,
                  &rho,eps,&T_in,&Y_e,P,S,
                  &dummy,&dummy,&dummy,&dummy,&dummy,
                  &keyerr,&anyerr );

  // TODO: add error handling
}
