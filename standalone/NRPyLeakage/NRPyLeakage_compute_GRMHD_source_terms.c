#include "Basic_defines.h"
#include "NRPyEOS.h"
#include "NRPyLeakage.h"

/*
 * (c) Leo Werneck
 * Compute GRMHD source terms following Ruffert et al. (1996)
 * https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
void NRPyLeakage_compute_GRMHD_source_terms(const int which_constants_to_use,
                                            const NRPyEOS_params *restrict eos_params,
                                            const REAL rho_b,
                                            const REAL Y_e,
                                            const REAL T,
                                            const REAL tau_nue,
                                            const REAL tau_anue,
                                            const REAL tau_nux,
                                            REAL *restrict R_source,
                                            REAL *restrict Q_source) {


  switch (which_constants_to_use) {
    case USE_NRPY_CONSTANTS:
      NRPyLeakage_compute_GRMHD_source_terms_nrpy_constants(eos_params,rho_b,Y_e,T,
                                                            tau_nue,tau_anue,tau_nux,
                                                            R_source,Q_source);
      break;
    case USE_HARM_CONSTANTS:
      NRPyLeakage_compute_GRMHD_source_terms_harm_constants(eos_params,rho_b,Y_e,T,
                                                            tau_nue,tau_anue,tau_nux,
                                                            R_source,Q_source);
      break;
    default:
      fprintf(stderr,"(NRPyLeakage) ERROR: Unknown constant type (%d) in NRPyLeakage_compute_GRMHD_source_terms().\n",which_constants_to_use);
      fprintf(stderr,"(NRPyLeakage) Options are: USE_NRPY_CONSTANTS (%d) and USE_HARM_CONSTANTS (%d)\n",USE_NRPY_CONSTANTS,USE_HARM_CONSTANTS);
      fprintf(stderr,"(NRPyLeakage) Aborting!\n");
      exit(1);
  }
}
