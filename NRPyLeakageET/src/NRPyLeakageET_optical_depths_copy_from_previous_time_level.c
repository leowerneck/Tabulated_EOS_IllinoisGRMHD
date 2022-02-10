#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void NRPyLeakageET_optical_depths_copy_from_previous_time_level(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Copy data from previous time level
#pragma omp parallel for
  for(int k=cctk_nghostzones[2];k<cctk_lsh[2]-cctk_nghostzones[2];k++) {
    for(int j=cctk_nghostzones[1];j<cctk_lsh[1]-cctk_nghostzones[1];j++) {
      for(int i=cctk_nghostzones[0];i<cctk_lsh[0]-cctk_nghostzones[0];i++) {
        const int index   = CCTK_GFINDEX3D(cctkGH,i,j,k);
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
