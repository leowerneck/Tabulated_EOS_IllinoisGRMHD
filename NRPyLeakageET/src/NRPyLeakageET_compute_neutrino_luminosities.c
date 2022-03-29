#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "NRPyLeakageET.h"

/*
 * (c) Leo Werneck
 * Compute GRMHD source terms following Ruffert et al. (1996)
 * https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
void NRPyLeakageET_compute_neutrino_luminosities(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(!NRPyLeakageET_ProcessOwnsData()) return;

  if(verbosity_level>1) CCTK_VINFO("Computing luminosities at ref. lvl. %d...",GetRefinementLevel(cctkGH));

  switch (constants_key) {
  case USE_NRPY_CONSTANTS:
#pragma omp parallel for
    for(int k=cctk_nghostzones[2];k<cctk_lsh[2]-cctk_nghostzones[2];k++) {
      for(int j=cctk_nghostzones[1];j<cctk_lsh[1]-cctk_nghostzones[1];j++) {
        for(int i=cctk_nghostzones[0];i<cctk_lsh[0]-cctk_nghostzones[0];i++) {
          // Step 1: Set gridpoint index
          const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

          const CCTK_REAL rhoL         = rho[index];
          if( rhoL < rho_threshold ) {
            lum_nue[index] = lum_anue[index] = lum_nux[index] = 0.0;
          }
          else {
            // Step 2: Read from main memory
            const CCTK_REAL alpL         = alp[index];
            const CCTK_REAL gxxL         = gxx[index];
            const CCTK_REAL gxyL         = gxy[index];
            const CCTK_REAL gxzL         = gxz[index];
            const CCTK_REAL gyyL         = gyy[index];
            const CCTK_REAL gyzL         = gyz[index];
            const CCTK_REAL gzzL         = gzz[index];
            const CCTK_REAL Y_eL         = Y_e[index];
            const CCTK_REAL temperatureL = temperature[index];
            const CCTK_REAL wL           = w_lorentz[index];
            const CCTK_REAL tau_nueL [2] = {tau_0_nue [index],tau_1_nue [index]};
            const CCTK_REAL tau_anueL[2] = {tau_0_anue[index],tau_1_anue[index]};
            const CCTK_REAL tau_nuxL [2] = {tau_0_nux [index],tau_1_nux [index]};

            // Step 3: Compute neutrino luminosities
            CCTK_REAL lum_nueL,lum_anueL,lum_nuxL;
            NRPyLeakageET_compute_neutrino_luminosities_nrpy_constants(alpL,gxxL,gxyL,gxzL,gyyL,gyzL,gzzL,
                                                                       rhoL,Y_eL,temperatureL,wL,
                                                                       tau_nueL,tau_anueL,tau_nuxL,
                                                                       &lum_nueL,&lum_anueL,&lum_nuxL);

            // Step 4: Write to main memory
            lum_nue [index] = lum_nueL;
            lum_anue[index] = lum_anueL;
            lum_nux [index] = lum_nuxL;
          }
        }
      }
    }
    break;
  case USE_HARM_CONSTANTS:
#pragma omp parallel for
    for(int k=cctk_nghostzones[2];k<cctk_lsh[2]-cctk_nghostzones[2];k++) {
      for(int j=cctk_nghostzones[1];j<cctk_lsh[1]-cctk_nghostzones[1];j++) {
        for(int i=cctk_nghostzones[0];i<cctk_lsh[0]-cctk_nghostzones[0];i++) {
          // Step 1: Set gridpoint index
          const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

          const CCTK_REAL rhoL         = rho[index];
          if( rhoL < rho_threshold ) {
            lum_nue[index] = lum_anue[index] = lum_nux[index] = 0.0;
          }
          else {
            // Step 2: Read from main memory
            const CCTK_REAL alpL         = alp[index];
            const CCTK_REAL gxxL         = gxx[index];
            const CCTK_REAL gxyL         = gxy[index];
            const CCTK_REAL gxzL         = gxz[index];
            const CCTK_REAL gyyL         = gyy[index];
            const CCTK_REAL gyzL         = gyz[index];
            const CCTK_REAL gzzL         = gzz[index];
            const CCTK_REAL Y_eL         = Y_e[index];
            const CCTK_REAL temperatureL = temperature[index];
            const CCTK_REAL wL           = w_lorentz[index];
            const CCTK_REAL tau_nueL [2] = {tau_0_nue [index],tau_1_nue [index]};
            const CCTK_REAL tau_anueL[2] = {tau_0_anue[index],tau_1_anue[index]};
            const CCTK_REAL tau_nuxL [2] = {tau_0_nux [index],tau_1_nux [index]};

            // Step 3: Compute neutrino luminosities
            CCTK_REAL lum_nueL,lum_anueL,lum_nuxL;
            NRPyLeakageET_compute_neutrino_luminosities_harm_constants(alpL,gxxL,gxyL,gxzL,gyyL,gyzL,gzzL,
                                                                       rhoL,Y_eL,temperatureL,wL,
                                                                       tau_nueL,tau_anueL,tau_nuxL,
                                                                       &lum_nueL,&lum_anueL,&lum_nuxL);

            // Step 4: Write to main memory
            lum_nue [index] = lum_nueL;
            lum_anue[index] = lum_anueL;
            lum_nux [index] = lum_nuxL;
          }
        }
      }
    }
    break;
    default:
      CCTK_VWARN(CCTK_WARN_ALERT,"Unknown constant type (%d) in NRPyLeakageET_compute_neutrino_luminosities()",constants_key);
      CCTK_VWARN(CCTK_WARN_ALERT,"Options are: USE_NRPY_CONSTANTS (%d) and USE_HARM_CONSTANTS (%d)",USE_NRPY_CONSTANTS,USE_HARM_CONSTANTS);
      CCTK_VWARN(CCTK_WARN_ABORT,"Aborting!");
  }
}
