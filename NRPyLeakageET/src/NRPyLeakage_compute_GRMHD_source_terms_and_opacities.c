#include "NRPyLeakage.h"

/*
 * (c) Leo Werneck
 * Compute GRMHD source terms following Ruffert et al. (1996)
 * https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
void NRPyLeakage_compute_GRMHD_source_terms_and_opacities(const int which_constants_to_use,
                                                          const REAL rho_b,
                                                          const REAL Y_e,
                                                          const REAL T,
                                                          const REAL *restrict tau_nue,
                                                          const REAL *restrict tau_anue,
                                                          const REAL *restrict tau_nux,
                                                          REAL *restrict R_source,
                                                          REAL *restrict Q_source,
                                                          REAL *restrict kappa_nue,
                                                          REAL *restrict kappa_anue,
                                                          REAL *restrict kappa_nux) {


  switch (which_constants_to_use) {
    case USE_NRPY_CONSTANTS:
      NRPyLeakage_compute_GRMHD_source_terms_and_opacities_nrpy_constants(rho_b,Y_e,T,
                                                                          tau_nue,tau_anue,tau_nux,
                                                                          R_source,Q_source,
                                                                          kappa_nue,kappa_anue,kappa_nux);
      break;
  case USE_HARM_CONSTANTS:
    NRPyLeakage_compute_GRMHD_source_terms_and_opacities_harm_constants(rho_b,Y_e,T,
                                                                        tau_nue,tau_anue,tau_nux,
                                                                        R_source,Q_source,
                                                                        kappa_nue,kappa_anue,kappa_nux);
      break;
    default:
      fprintf(stderr,"(NRPyLeakage) ERROR: Unknown constant type (%d) in NRPyLeakage_compute_GRMHD_source_terms().\n",which_constants_to_use);
      fprintf(stderr,"(NRPyLeakage) Options are: USE_NRPY_CONSTANTS (%d) and USE_HARM_CONSTANTS (%d)\n",USE_NRPY_CONSTANTS,USE_HARM_CONSTANTS);
      fprintf(stderr,"(NRPyLeakage) Aborting!\n");
      exit(1);
  }
}
