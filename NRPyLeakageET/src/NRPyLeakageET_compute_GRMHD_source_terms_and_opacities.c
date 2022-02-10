#include "NRPyLeakageET.h"

/*
 * (c) Leo Werneck
 * Compute GRMHD source terms following Ruffert et al. (1996)
 * https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
void NRPyLeakageET_compute_GRMHD_source_terms_and_opacities(const int which_constants_to_use,
                                                            const CCTK_REAL rho_b,
                                                            const CCTK_REAL Y_e,
                                                            const CCTK_REAL T,
                                                            const CCTK_REAL *restrict tau_nue,
                                                            const CCTK_REAL *restrict tau_anue,
                                                            const CCTK_REAL *restrict tau_nux,
                                                            CCTK_REAL *restrict R_source,
                                                            CCTK_REAL *restrict Q_source,
                                                            CCTK_REAL *restrict kappa_nue,
                                                            CCTK_REAL *restrict kappa_anue,
                                                            CCTK_REAL *restrict kappa_nux) {


  switch (which_constants_to_use) {
    case USE_NRPY_CONSTANTS:
      NRPyLeakageET_compute_GRMHD_source_terms_and_opacities_nrpy_constants(rho_b,Y_e,T,
                                                                            tau_nue,tau_anue,tau_nux,
                                                                            R_source,Q_source,
                                                                            kappa_nue,kappa_anue,kappa_nux);
      break;
  case USE_HARM_CONSTANTS:
    NRPyLeakageET_compute_GRMHD_source_terms_and_opacities_harm_constants(rho_b,Y_e,T,
                                                                          tau_nue,tau_anue,tau_nux,
                                                                          R_source,Q_source,
                                                                          kappa_nue,kappa_anue,kappa_nux);
      break;
    default:
      fprintf(stderr,"(NRPyLeakageET) ERROR: Unknown constant type (%d) in NRPyLeakage_compute_GRMHD_source_terms().\n",which_constants_to_use);
      fprintf(stderr,"(NRPyLeakageET) Options are: USE_NRPY_CONSTANTS (%d) and USE_HARM_CONSTANTS (%d)\n",USE_NRPY_CONSTANTS,USE_HARM_CONSTANTS);
      fprintf(stderr,"(NRPyLeakageET) Aborting!\n");
      exit(1);
  }
}
