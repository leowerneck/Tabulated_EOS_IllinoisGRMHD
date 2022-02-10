#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


void NRPyLeakage_test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VInfo(CCTK_THORNSTRING,"test   : %e %e %e %e %e %e",tau_0_nue[10],tau_1_nue[10],tau_0_anue[10],tau_1_anue[10],tau_0_nux[10],tau_1_nux[10]);
  CCTK_VInfo(CCTK_THORNSTRING,"test p : %e %e %e %e %e %e",tau_0_nue_p[10],tau_1_nue_p[10],tau_0_anue_p[10],tau_1_anue_p[10],tau_0_nux_p[10],tau_1_nux_p[10]);
  CCTK_VInfo(CCTK_THORNSTRING,"test pp: %e %e %e %e %e %e",tau_0_nue_p_p[10],tau_1_nue_p_p[10],tau_0_anue_p_p[10],tau_1_anue_p_p[10],tau_0_nux_p_p[10],tau_1_nux_p_p[10]);
}
