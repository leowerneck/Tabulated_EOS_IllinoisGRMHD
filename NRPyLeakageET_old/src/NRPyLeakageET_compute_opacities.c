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
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(!NRPyLeakageET_ProcessOwnsData()) return;

  if(verbosity_level>1) CCTK_VInfo(CCTK_THORNSTRING,"Computing opacities at ref. lvl. %d...",GetRefinementLevel(cctkGH));

  switch (constants_key) {
  case USE_NRPY_CONSTANTS:
#pragma omp parallel
    for(int k=0;k<cctk_lsh[2];k++) {
      for(int j=0;j<cctk_lsh[1];j++) {
        for(int i=0;i<cctk_lsh[0];i++) {
          // Step 1: Set local gridfunction index
          const CCTK_INT index = CCTK_GFINDEX3D(cctkGH,i,j,k);

          // Step 2: Declare variables, read in density
          const CCTK_REAL rhoL = rho[index];
          CCTK_REAL kappa_nue[2], kappa_anue[2], kappa_nux[2];
          CCTK_REAL tau_nue[2], tau_anue[2], tau_nux[2];

          // Step 3: Check density threshold, compute opacities
          if( rhoL < rho_threshold ) {
            // Step 3.a: Below density threshold; set opacities and optical
            // depths to zero
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

          // Step 2: Declare variables, read in density
          const CCTK_REAL rhoL = rho[index];
          CCTK_REAL kappa_nue[2], kappa_anue[2], kappa_nux[2];
          CCTK_REAL tau_nue[2], tau_anue[2], tau_nux[2];

          // Step 3: Check density threshold, compute opacities
          if( rhoL < rho_threshold ) {
            // Step 3.a: Below density threshold; set opacities and optical
            // depths to zero
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
    fprintf(stderr,"(NRPyLeakageET) ERROR: Unknown constant type (%d) in NRPyLeakageET_compute_GRMHD_source_terms().\n",constants_key);
    fprintf(stderr,"(NRPyLeakageET) Options are: USE_NRPY_CONSTANTS (%d) and USE_HARM_CONSTANTS (%d)\n",USE_NRPY_CONSTANTS,USE_HARM_CONSTANTS);
    fprintf(stderr,"(NRPyLeakageET) Aborting!\n");
    exit(1);
  }

  if(verbosity_level>1) CCTK_VInfo(CCTK_THORNSTRING,"Finished computing opacities at ref. lvl. %d",GetRefinementLevel(cctkGH));
}
