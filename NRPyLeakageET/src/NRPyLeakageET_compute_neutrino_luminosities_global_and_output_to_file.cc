#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "carpet.hh"
#include "NRPyLeakageET.h"

void NRPyLeakageET_set_neutrino_luminosities_to_zero(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(!NRPyLeakageET_ProcessOwnsData()) return;

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) {
    for(int j=0;j<cctk_lsh[1];j++) {
      for(int i=0;i<cctk_lsh[0];i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        lum_nue [index] = 0.0;
        lum_anue[index] = 0.0;
        lum_nux [index] = 0.0;
      }
    }
  }
}

void NRPyLeakageET_compute_neutrino_luminosities_global_and_output_to_file(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if( cctk_iteration%compute_luminosities_every == 0 ) {
    // Step 1: Compute the neutrino luminosities
    BEGIN_REFLEVEL_LOOP(cctkGH) {
      BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
        BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
          NRPyLeakageET_set_neutrino_luminosities_to_zero(CCTK_PASS_CTOC);
          NRPyLeakageET_compute_neutrino_luminosities(CCTK_PASS_CTOC);
        } END_COMPONENT_LOOP;
      } END_MAP_LOOP;
    } END_REFLEVEL_LOOP;

    // Step 2: Perform the global integration
    // Step 2.a: Get sum operation handle
    const int op_sum = CCTK_ReductionHandle("sum");

    // Step 2.b: Set variables we want to reduce
    const char varnames[3][24] = {"NRPyLeakageET::lum_nue ",
                                  "NRPyLeakageET::lum_anue",
                                  "NRPyLeakageET::lum_nux "};


    // Step 2.c: Get indices for the variables we want to reduce
    const int varindices[3] = {CCTK_VarIndex("NRPyLeakageET::lum_nue"),
                               CCTK_VarIndex("NRPyLeakageET::lum_anue"),
                               CCTK_VarIndex("NRPyLeakageET::lum_nux")};

    // Step 2.d: Compute the global sum of the luminosities
    const CCTK_REAL d3x = CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2);
    CCTK_REAL global_luminosities[3];
    for(int i=0;i<3;i++) {
      int varindex = varindices[i];
      assert(varindex>=0);
      CCTK_REAL global_luminosity = 0.0;
      const int ierr = CCTK_Reduce(cctkGH,
                                   -1, // target processors; -1 -> all
                                   op_sum,
                                   1,  // Number of outputs
                                   CCTK_VARIABLE_REAL,
                                   &global_luminosity,
                                   1,  // Number of inputs
                                   varindex);
      if( ierr ) CCTK_VERROR("Error in reduction of %s",varnames[i]);
      global_luminosities[i] = d3x * global_luminosity;
    }

    // Step 3: Now output the luminosities to file
    char filename[512];
    sprintf(filename,"%s/%s",out_dir,luminosities_outfile);
    FILE *fp;
    if( cctk_iteration == 0 ) {
      fp = fopen(filename,"w");
      if( !fp ) CCTK_VERROR("Could not open file %s",luminosities_outfile);
      fprintf(fp,"# NRPyLeakageET output: Neutrino luminosities integrated over entire grid\n");
      fprintf(fp,"# Column 1: cctk_iteration\n");
      fprintf(fp,"# Column 2: cctk_time\n");
      fprintf(fp,"# Column 3: Electron neutrino luminosity\n");
      fprintf(fp,"# Column 4: Electron antineutrino luminosity\n");
      fprintf(fp,"# Column 5: Heavy lepton neutrinos luminosity\n");
    }
    else {
      fp = fopen(filename,"a");
      if( !fp ) CCTK_VERROR("Could not open file %s",luminosities_outfile);
    }

    fprintf(fp,"%d %.15e %.15e %.15e %.15e\n",cctk_iteration,cctk_time,global_luminosities[0],global_luminosities[1],global_luminosities[2]);
    fclose(fp);
  }

}
