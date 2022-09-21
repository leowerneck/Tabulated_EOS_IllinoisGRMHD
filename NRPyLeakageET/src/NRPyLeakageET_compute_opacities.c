#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "NRPyLeakageET.h"

/*
 * (c) Leo Werneck
 * Compute GRMHD source terms following Ruffert et al. (1996)
 * https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
void NRPyLeakageET_compute_opacities(CCTK_ARGUMENTS) {
  #ifdef DECLARE_CCTK_ARGUMENTS_NRPyLeakageET_Initialize
  DECLARE_CCTK_ARGUMENTS_CHECKED(NRPyLeakageET_Initialize);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;

  if(!NRPyLeakageET_ProcessOwnsData()) return;

  if(verbosity_level>1) CCTK_VINFO("Computing opacities at ref. lvl. %d...",GetRefinementLevel(cctkGH));

  switch (constants_key) {
  case USE_NRPY_CONSTANTS:
#pragma omp parallel
    for(int k=0;k<cctk_lsh[2];k++) {
      for(int j=0;j<cctk_lsh[1];j++) {
        for(int i=0;i<cctk_lsh[0];i++) {
          // Step 1: Set local gridfunction index
          const CCTK_INT index = CCTK_GFINDEX3D(cctkGH,i,j,k);

          // Step 2: Declare variables
          const CCTK_REAL rhoL = rho[index];
          CCTK_REAL kappa_nue[2], kappa_anue[2], kappa_nux[2];
          CCTK_REAL tau_nue[2], tau_anue[2], tau_nux[2];

          // Step 3: Check density threshold, compute opacities
          if( rhoL < rho_min_threshold || rhoL > rho_max_threshold ) {
            // Step 3.a: Below density threshold; set opacities and optical
            // depths to zero
            kappa_nue[0] = kappa_anue[0] = kappa_nux[0] = 0.0;
            kappa_nue[1] = kappa_anue[1] = kappa_nux[1] = 0.0;
            tau_nue  [0] = tau_anue  [0] = tau_nux  [0] = 0.0;
            tau_nue  [1] = tau_anue  [1] = tau_nux  [1] = 0.0;
          }
          else {
            const CCTK_REAL gxxL  = gxx[index];
            const CCTK_REAL gxyL  = gxy[index];
            const CCTK_REAL gxzL  = gxz[index];
            const CCTK_REAL gyyL  = gyy[index];
            const CCTK_REAL gyzL  = gyz[index];
            const CCTK_REAL gzzL  = gzz[index];
            const CCTK_REAL gdet  = fabs(gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + gxzL * gxyL * gyzL
                                       - gxzL * gyyL * gxzL - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL);
            const CCTK_REAL phiL  = (1.0/12.0) * log(gdet);
            const CCTK_REAL psiL  = exp(phiL);
            const CCTK_REAL psi2L = psiL *psiL;
            const CCTK_REAL psi4L = psi2L*psi2L;
            const CCTK_REAL psi6L = psi4L*psi2L;
            if( psi6L > psi6_threshold ) {
              kappa_nue[0] = kappa_anue[0] = kappa_nux[0] = 0.0;
              kappa_nue[1] = kappa_anue[1] = kappa_nux[1] = 0.0;
              tau_nue  [0] = tau_anue  [0] = tau_nux  [0] = 0.0;
              tau_nue  [1] = tau_anue  [1] = tau_nux  [1] = 0.0;
            }
            else {
              // Step 3.b: Above density threshold, compute opacities
              // Step 3.b.i: Read from main memory
              const CCTK_REAL Y_eL         = Y_e[index];
              const CCTK_REAL temperatureL = temperature[index];
              tau_nue [0]                  = tau_0_nue [index];
              tau_nue [1]                  = tau_1_nue [index];
              tau_anue[0]                  = tau_0_anue[index];
              tau_anue[1]                  = tau_1_anue[index];
              tau_nux [0]                  = tau_0_nux [index];
              tau_nux [1]                  = tau_1_nux [index];

              // Step 3.b.ii: Compute opacities
              NRPyLeakageET_compute_opacities_nrpy_constants(rhoL,Y_eL,temperatureL,
                                                             tau_nue,tau_anue,tau_nux,
                                                             kappa_nue,kappa_anue,kappa_nux);
            }
          }

          // Step 4: Write to main memory
          kappa_0_nue [index] = kappa_nue [0];
          kappa_1_nue [index] = kappa_nue [1];
          kappa_0_anue[index] = kappa_anue[0];
          kappa_1_anue[index] = kappa_anue[1];
          kappa_0_nux [index] = kappa_nux [0];
          kappa_1_nux [index] = kappa_nux [1];

          // Step 5: Copy the optical depths to the previous time level as a placeholder
          tau_0_nue_p [index] = tau_nue [0];
          tau_1_nue_p [index] = tau_nue [1];
          tau_0_anue_p[index] = tau_anue[0];
          tau_1_anue_p[index] = tau_anue[1];
          tau_0_nux_p [index] = tau_nux [0];
          tau_1_nux_p [index] = tau_nux [1];
        }
      }
    }
    break;
  case USE_HARM_CONSTANTS:
#pragma omp parallel
    for(int k=0;k<cctk_lsh[2];k++) {
      for(int j=0;j<cctk_lsh[1];j++) {
        for(int i=0;i<cctk_lsh[0];i++) {
          // Step 1: Set local gridfunction index
          const CCTK_INT index = CCTK_GFINDEX3D(cctkGH,i,j,k);

          // Step 2: Declare variables
          const CCTK_REAL rhoL = rho[index];
          CCTK_REAL kappa_nue[2], kappa_anue[2], kappa_nux[2];
          CCTK_REAL tau_nue[2], tau_anue[2], tau_nux[2];

          // Step 3: Check density threshold, compute opacities
          if( rhoL < rho_min_threshold || rhoL > rho_max_threshold ) {
            // Step 3.a: Below density threshold; set opacities and optical
            // depths to zero
            kappa_nue[0] = kappa_anue[0] = kappa_nux[0] = 0.0;
            kappa_nue[1] = kappa_anue[1] = kappa_nux[1] = 0.0;
            tau_nue  [0] = tau_anue  [0] = tau_nux  [0] = 0.0;
            tau_nue  [1] = tau_anue  [1] = tau_nux  [1] = 0.0;
          }
          else {
            const CCTK_REAL gxxL  = gxx[index];
            const CCTK_REAL gxyL  = gxy[index];
            const CCTK_REAL gxzL  = gxz[index];
            const CCTK_REAL gyyL  = gyy[index];
            const CCTK_REAL gyzL  = gyz[index];
            const CCTK_REAL gzzL  = gzz[index];
            const CCTK_REAL gdet  = fabs(gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + gxzL * gxyL * gyzL
                                       - gxzL * gyyL * gxzL - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL);
            const CCTK_REAL phiL  = (1.0/12.0) * log(gdet);
            const CCTK_REAL psiL  = exp(phiL);
            const CCTK_REAL psi2L = psiL *psiL;
            const CCTK_REAL psi4L = psi2L*psi2L;
            const CCTK_REAL psi6L = psi4L*psi2L;
            if( psi6L > psi6_threshold ) {
              kappa_nue[0] = kappa_anue[0] = kappa_nux[0] = 0.0;
              kappa_nue[1] = kappa_anue[1] = kappa_nux[1] = 0.0;
              tau_nue  [0] = tau_anue  [0] = tau_nux  [0] = 0.0;
              tau_nue  [1] = tau_anue  [1] = tau_nux  [1] = 0.0;
            }
            else {
              // Step 3.b: Above density threshold, compute opacities
              // Step 3.b.i: Read from main memory
              const CCTK_REAL Y_eL         = Y_e[index];
              const CCTK_REAL temperatureL = temperature[index];
              tau_nue [0]                  = tau_0_nue [index];
              tau_nue [1]                  = tau_1_nue [index];
              tau_anue[0]                  = tau_0_anue[index];
              tau_anue[1]                  = tau_1_anue[index];
              tau_nux [0]                  = tau_0_nux [index];
              tau_nux [1]                  = tau_1_nux [index];

              // Step 3.b.ii: Compute opacities
              NRPyLeakageET_compute_opacities_harm_constants(rhoL,Y_eL,temperatureL,
                                                             tau_nue,tau_anue,tau_nux,
                                                             kappa_nue,kappa_anue,kappa_nux);
            }
          }

          // Step 4: Write to main memory
          kappa_0_nue [index] = kappa_nue [0];
          kappa_1_nue [index] = kappa_nue [1];
          kappa_0_anue[index] = kappa_anue[0];
          kappa_1_anue[index] = kappa_anue[1];
          kappa_0_nux [index] = kappa_nux [0];
          kappa_1_nux [index] = kappa_nux [1];

          // Step 5: Copy the optical depths to the previous time level as a placeholder
          tau_0_nue_p [index] = tau_nue [0];
          tau_1_nue_p [index] = tau_nue [1];
          tau_0_anue_p[index] = tau_anue[0];
          tau_1_anue_p[index] = tau_anue[1];
          tau_0_nux_p [index] = tau_nux [0];
          tau_1_nux_p [index] = tau_nux [1];
        }
      }
    }
    break;
  default:
    CCTK_VINFO("Unknown constant type (%d) in NRPyLeakageET_compute_opacities()",constants_key);
    CCTK_VINFO("Options are: USE_NRPY_CONSTANTS (%d) and USE_HARM_CONSTANTS (%d)",USE_NRPY_CONSTANTS,USE_HARM_CONSTANTS);
    CCTK_ERROR("Aborting!");
  }

  if(verbosity_level>1) CCTK_VINFO("Finished computing opacities at ref. lvl. %d",GetRefinementLevel(cctkGH));
}
