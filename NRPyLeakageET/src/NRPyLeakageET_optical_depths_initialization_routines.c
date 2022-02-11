#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "NRPyLeakageET.h"

void NRPyLeakageET_optical_depths_initialize_to_zero(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Step 1: Initialize all time levels of the optical
  //         depth gridfunctions to zero.
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) {
    for(int j=0;j<cctk_lsh[1];j++) {
      for(int i=0;i<cctk_lsh[0];i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        tau_0_nue [index] = 0.0;
        tau_1_nue [index] = 0.0;
        tau_0_anue[index] = 0.0;
        tau_1_anue[index] = 0.0;
        tau_0_nux [index] = 0.0;
        tau_1_nux [index] = 0.0;

        tau_0_nue_p [index] = 0.0;
        tau_1_nue_p [index] = 0.0;
        tau_0_anue_p[index] = 0.0;
        tau_1_anue_p[index] = 0.0;
        tau_0_nux_p [index] = 0.0;
        tau_1_nux_p [index] = 0.0;

        tau_0_nue_p_p [index] = 0.0;
        tau_1_nue_p_p [index] = 0.0;
        tau_0_anue_p_p[index] = 0.0;
        tau_1_anue_p_p[index] = 0.0;
        tau_0_nux_p_p [index] = 0.0;
        tau_1_nux_p_p [index] = 0.0;
      }
    }
  }

  CCTK_INFO("Initialized all optical depths gridfunctions to zero");
}

void NRPyLeakageET_initialization_driver(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Step 1: Compute opacities and store results of previous
  //         iteration on previous time level
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) {
    for(int j=0;j<cctk_lsh[1];j++) {
      for(int i=0;i<cctk_lsh[0];i++) {
        // Step 1.a: Set gridpoint index
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // Step 1.b: Read from main memory
        const CCTK_REAL rhoL         = rho[index];
        const CCTK_REAL Y_eL         = Y_e[index];
        const CCTK_REAL temperatureL = temperature[index];
        const CCTK_REAL tau_nueL [2] = {tau_0_nue [index],tau_1_nue [index]};
        const CCTK_REAL tau_anueL[2] = {tau_0_anue[index],tau_1_anue[index]};
        const CCTK_REAL tau_nuxL [2] = {tau_0_nux [index],tau_1_nux [index]};

        // Step 1.c: Compute opacities
        CCTK_REAL kappa_nueL[2], kappa_anueL[2], kappa_nuxL[2];
        NRPyLeakageET_compute_opacities(NRPyLeakageET_constants_key,
                                        rhoL,Y_eL,temperatureL,
                                        tau_nueL,tau_anueL,tau_nuxL,
                                        kappa_nueL,kappa_anueL,kappa_nuxL);

        // Step 1.d: Write opacities to main memory
        kappa_0_nue [index] = kappa_nueL [0];
        kappa_1_nue [index] = kappa_nueL [1];
        kappa_0_anue[index] = kappa_anueL[0];
        kappa_1_anue[index] = kappa_anueL[1];
        kappa_0_nux [index] = kappa_nuxL [0];
        kappa_1_nux [index] = kappa_nuxL [1];

        // Step 1.e: Copy optical depth to previous time level
        tau_0_nue_p [index] = tau_nueL [0];
        tau_1_nue_p [index] = tau_nueL [1];
        tau_0_anue_p[index] = tau_anueL[0];
        tau_1_anue_p[index] = tau_anueL[1];
        tau_0_nux_p [index] = tau_nuxL [0];
        tau_1_nux_p [index] = tau_nuxL [1];
      }
    }
  }

  // Step 2: Update optical depths
  NRPyLeakageET_optical_depths_PathOfLeastResistance(CCTK_PASS_CTOC);

  // All done!
}

void NRPyLeakageET_copy_opacities_and_optical_depths_to_previous_time_levels(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Step 1: Copy data from the current time level to previous time levels
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) {
    for(int j=0;j<cctk_lsh[1];j++) {
      for(int i=0;i<cctk_lsh[0];i++) {
        // Step 1.a: Set gridpoint index
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // Step 1.b: Read from main memory
        const CCTK_REAL kappa_0_nueL  = kappa_0_nue [index];
        const CCTK_REAL kappa_1_nueL  = kappa_1_nue [index];
        const CCTK_REAL kappa_0_anueL = kappa_0_anue[index];
        const CCTK_REAL kappa_1_anueL = kappa_1_anue[index];
        const CCTK_REAL kappa_0_nuxL  = kappa_0_nux [index];
        const CCTK_REAL kappa_1_nuxL  = kappa_1_nux [index];
        const CCTK_REAL tau_0_nueL    = tau_0_nue   [index];
        const CCTK_REAL tau_1_nueL    = tau_1_nue   [index];
        const CCTK_REAL tau_0_anueL   = tau_0_anue  [index];
        const CCTK_REAL tau_1_anueL   = tau_1_anue  [index];
        const CCTK_REAL tau_0_nuxL    = tau_0_nux   [index];
        const CCTK_REAL tau_1_nuxL    = tau_1_nux   [index];

        // Step 1.c: Write opacities to main memory
        kappa_0_nue_p [index] = kappa_0_nueL;
        kappa_1_nue_p [index] = kappa_1_nueL;
        kappa_0_anue_p[index] = kappa_0_anueL;
        kappa_1_anue_p[index] = kappa_1_anueL;
        kappa_0_nux_p [index] = kappa_0_nuxL;
        kappa_1_nux_p [index] = kappa_1_nuxL;

        kappa_0_nue_p_p [index] = kappa_0_nueL;
        kappa_1_nue_p_p [index] = kappa_1_nueL;
        kappa_0_anue_p_p[index] = kappa_0_anueL;
        kappa_1_anue_p_p[index] = kappa_1_anueL;
        kappa_0_nux_p_p [index] = kappa_0_nuxL;
        kappa_1_nux_p_p [index] = kappa_1_nuxL;

        // Step 1.d: Write optical depths to main memory
        tau_0_nue_p [index] = tau_0_nueL;
        tau_1_nue_p [index] = tau_1_nueL;
        tau_0_anue_p[index] = tau_0_anueL;
        tau_1_anue_p[index] = tau_1_anueL;
        tau_0_nux_p [index] = tau_0_nuxL;
        tau_1_nux_p [index] = tau_1_nuxL;

        tau_0_nue_p_p [index] = tau_0_nueL;
        tau_1_nue_p_p [index] = tau_1_nueL;
        tau_0_anue_p_p[index] = tau_0_anueL;
        tau_1_anue_p_p[index] = tau_1_anueL;
        tau_0_nux_p_p [index] = tau_0_nuxL;
        tau_1_nux_p_p [index] = tau_1_nuxL;
      }
    }
  }
}

void NRPyLeakageET_copy_optical_depths_from_previous_time_level(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
    // Step 1: Copy optical depths from the previous time level to the current time level
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) {
    for(int j=0;j<cctk_lsh[1];j++) {
      for(int i=0;i<cctk_lsh[0];i++) {
        // Step 1.a: Set gridpoint index
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // Step 1.d: Copy optical depths
        tau_0_nue [index] = tau_0_nue_p [index];
        tau_1_nue [index] = tau_1_nue_p [index];
        tau_0_anue[index] = tau_0_anue_p[index];
        tau_1_anue[index] = tau_1_anue_p[index];
        tau_0_nux [index] = tau_0_nux_p [index];
        tau_1_nux [index] = tau_1_nux_p [index];
      }
    }
  }
}
