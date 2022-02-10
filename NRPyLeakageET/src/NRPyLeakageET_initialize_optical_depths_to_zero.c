#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void NRPyLeakageET_initialize_optical_depths_to_zero(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

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
