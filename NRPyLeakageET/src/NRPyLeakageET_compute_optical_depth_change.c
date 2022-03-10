#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "NRPyLeakageET.h"

void NRPyLeakageET_CopyOpticalDepthsToAux(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(!NRPyLeakageET_ProcessOwnsData()) return;

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) {
    for(int j=0;j<cctk_lsh[1];j++) {
      for(int i=0;i<cctk_lsh[0];i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        tau_0_nue_aux [index] = tau_0_nue [index];
        tau_1_nue_aux [index] = tau_1_nue [index];
        tau_0_anue_aux[index] = tau_0_anue[index];
        tau_1_anue_aux[index] = tau_1_anue[index];
        tau_0_nux_aux [index] = tau_0_nux [index];
        tau_1_nux_aux [index] = tau_1_nux [index];
      }
    }
  }
}

void NRPyLeakageET_compute_optical_depth_change(CCTK_ARGUMENTS, const int it) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(!NRPyLeakageET_ProcessOwnsData()) return;

  if(verbosity_level>1) CCTK_VInfo(CCTK_THORNSTRING,"Beginning to compute changes in optical depths at Ref. Lev. %d",GetRefinementLevel(cctkGH));

  char filename[256];
  sprintf(filename,"output/optical_depths_rl%d_it%03d",GetRefinementLevel(cctkGH),it);
  FILE *fp = fopen(filename,"w");

#pragma omp parallel for
  for(int k=cctk_nghostzones[2];k<cctk_lsh[2]-cctk_nghostzones[2];k++) {
    for(int j=cctk_nghostzones[1];j<cctk_lsh[1]-cctk_nghostzones[1];j++) {
      for(int i=cctk_nghostzones[0];i<cctk_lsh[0]-cctk_nghostzones[0];i++) {
        // Step 1: Set local index
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        if( rho[index] < rho_threshold ) {
          tau_0_nue_aux [index] = 0.0;
          tau_1_nue_aux [index] = 0.0;
          tau_0_anue_aux[index] = 0.0;
          tau_1_anue_aux[index] = 0.0;
          tau_0_nux_aux [index] = 0.0;
          tau_1_nux_aux [index] = 0.0;
        }
        else {
          // Step 2: Read from main memory
          const CCTK_REAL tau_0_nueL      = tau_0_nue     [index];
          const CCTK_REAL tau_1_nueL      = tau_1_nue     [index];
          const CCTK_REAL tau_0_anueL     = tau_0_anue    [index];
          const CCTK_REAL tau_1_anueL     = tau_1_anue    [index];
          const CCTK_REAL tau_0_nuxL      = tau_0_nux     [index];
          const CCTK_REAL tau_1_nuxL      = tau_1_nux     [index];
          const CCTK_REAL tau_0_nue_auxL  = tau_0_nue_aux [index];
          const CCTK_REAL tau_1_nue_auxL  = tau_1_nue_aux [index];
          const CCTK_REAL tau_0_anue_auxL = tau_0_anue_aux[index];
          const CCTK_REAL tau_1_anue_auxL = tau_1_anue_aux[index];
          const CCTK_REAL tau_0_nux_auxL  = tau_0_nux_aux [index];
          const CCTK_REAL tau_1_nux_auxL  = tau_1_nux_aux [index];

          if( j==cctk_lsh[1]/2 && k==cctk_lsh[2]/2 ) {
            fprintf(fp,"%d %e %e %e\n",i,x[index],tau_0_nueL,tau_0_nue_auxL);
          }

          // Step 3: Compute differences in tau
          const CCTK_REAL diff_tau_0_nue  = tau_0_nueL -tau_0_nue_auxL;
          const CCTK_REAL diff_tau_1_nue  = tau_1_nueL -tau_1_nue_auxL;
          const CCTK_REAL diff_tau_0_anue = tau_0_anueL-tau_0_anue_auxL;
          const CCTK_REAL diff_tau_1_anue = tau_1_anueL-tau_1_anue_auxL;
          const CCTK_REAL diff_tau_0_nux  = tau_0_nuxL -tau_0_nux_auxL;
          const CCTK_REAL diff_tau_1_nux  = tau_1_nuxL -tau_1_nux_auxL;

          // Step 4: Write to main memory
          tau_0_nue_aux [index] = diff_tau_0_nue *diff_tau_0_nue;
          tau_1_nue_aux [index] = diff_tau_1_nue *diff_tau_1_nue;
          tau_0_anue_aux[index] = diff_tau_0_anue*diff_tau_0_anue;
          tau_1_anue_aux[index] = diff_tau_1_anue*diff_tau_1_anue;
          tau_0_nux_aux [index] = diff_tau_0_nux *diff_tau_0_nux;
          tau_1_nux_aux [index] = diff_tau_1_nux *diff_tau_1_nux;
        }
      }
    }
  }
  fclose(fp);

  if(verbosity_level>1) CCTK_VInfo(CCTK_THORNSTRING,"Finished computing changes in optical depths at Ref. Lev. %d",GetRefinementLevel(cctkGH));
}
