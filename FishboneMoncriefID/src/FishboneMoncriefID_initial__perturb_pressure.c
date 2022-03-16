#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "FishboneMoncriefID.h"

void FishboneMoncriefID_initial__perturb_pressure(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  for(CCTK_INT k=0;k<cctk_lsh[2];k++) {
    for(CCTK_INT j=0;j<cctk_lsh[1];j++) {
      for(CCTK_INT i=0;i<cctk_lsh[0];i++) {
        CCTK_INT idx = CCTK_GFINDEX3D(cctkGH,i,j,k);
        // Generate random number in range [0,1),
        // snippet courtesy http://daviddeley.com/random/crandom.htm
        CCTK_REAL random_number_between_0_and_1 = ( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );

        CCTK_REAL random_number_between_min_and_max = random_min + (random_max - random_min)*random_number_between_0_and_1;
        press[idx] = press[idx]*(1.0 + random_number_between_min_and_max);
        // Add 1e-300 to rho to avoid division by zero when density is zero.
        eps[idx] = press[idx] / ((rho[idx] + 1e-300) * (gamma - 1.0));
      }
    }
  }
}
