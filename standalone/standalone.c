#include "Basic_defines.h"

#include "NRPyEOS.h"
#include "NRPyLeakage.h"
#include "Leakage_tests.h"
#include "harm_neutrinos.h"
#include "harm_units.h"

void calc_taus_(const NRPyEOS_params *eos_params,
                const REAL *eos_tempmin,
                const REAL *atm_temp,
                const REAL *eos_yemax,
                const REAL *eos_yemin,
                const REAL *rho_gf,
                const REAL *kappa_gf,
                REAL *rho,
                REAL *temp,
                REAL *ye,
                REAL oldtauruff[3][1],
                REAL tauruff[3][1],
                REAL chiross[3][1],
                REAL heatflux[3][1],
                REAL heaterms[3],
                REAL heateave[3],
                const int *nzones,
                REAL *rad,
                REAL *ds,
                REAL compos[4][1],
                REAL *xentropy);

void opacities_dependence_on_density_and_electron_fraction(const NRPyEOS_params *eos_params, const REAL T) {

  // Step 1: Set limits
  const REAL lrmin = log(eos_params->eos_rhomin);
  const REAL yemin = eos_params->eos_yemin;

  // Step 2: Set optical depths
  const REAL tau_nue [2] = {0.0,0.0};
  const REAL tau_anue[2] = {0.0,0.0};
  const REAL tau_nux [2] = {0.0,0.0};

  // Step 3: Loop over the grid rho-Y_e, computing opacities
  FILE *fp = fopen("opacities_rho_ye.asc","w");
  for(int j=0;j<eos_params->nye;j++) {
    const REAL Y_e = yemin + eos_params->dye*j;
    for(int i=0;i<eos_params->nrho;i++) {
      const REAL rho = exp(lrmin + eos_params->drho*i);
      REAL kappa_nue[2], kappa_anue[2], kappa_nux[2];
      NRPyLeakage_compute_opacities(USE_HARM_CONSTANTS,eos_params,rho,Y_e,T,tau_nue,tau_anue,tau_nux,
                                    kappa_nue,kappa_anue,kappa_nux);
      fprintf(fp,"%e %e %e %e %e %e %e %e\n",rho,Y_e,kappa_nue[0],kappa_nue[1],kappa_anue[0],kappa_anue[1],kappa_nux[0],kappa_nux[1]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}

int main(int argc, char **argv) {

  if( argc < 3 ) {
    fprintf(stderr,"(NRPyLeakage) Correct usage is ./Leakage_standalone eos_file_path test_type [parameters]\n");
    fprintf(stderr,"(NRPyLeakage) Test types are:\n");
    fprintf(stderr,"(NRPyLeakage)     0: Optically thin gas (tests number/cooling rates)\n");
    fprintf(stderr,"(NRPyLeakage)     1: Constant density sphere (tests optical depth)\n");
    fprintf(stderr,"(NRPyLeakage)     2: 2D data for opacities as a function of (rho,Y_e) for given temperature in MeV\n");
    fprintf(stderr,"(NRPyLeakage)     3: Equation of state table interpolation benchmark\n");
    fprintf(stderr,"(NRPyLeakage)     4: Table slicer - constant Y_e and S\n");
    exit(1);
  }

  const int test_type = atoi(argv[2]);

  if( test_type == 0 && argc != 6 ) {
    fprintf(stderr,"(NRPyLeakage) For the optically thin gas test one must also provide 3 additional arguments:\n");
    fprintf(stderr,"(NRPyLeakage)     The gas density in code units\n");
    fprintf(stderr,"(NRPyLeakage)     The gas electron fraction\n");
    fprintf(stderr,"(NRPyLeakage)     The gas temperature in MeV\n");
    exit(1);
  }

  if( test_type == 2 && argc != 4 ) {
    fprintf(stderr,"(NRPyLeakage) For the 2D data for opacities test one must also provide 1 additional argument:\n");
    fprintf(stderr,"(NRPyLeakage)     The temperature in MeV\n");
    exit(1);
  }

  if( test_type == 4 && argc != 5 ) {
    fprintf(stderr,"(NRPyLeakage) For the table slicer one must also provide 2 additional arguments:\n");
    fprintf(stderr,"(NRPyLeakage)     The electron fraction\n");
    fprintf(stderr,"(NRPyLeakage)     The specific entropy in kb/baryon\n");
    exit(1);
  }

  // Read EOS table
  NRPyEOS_params eos_params;
  NRPyEOS_readtable_set_EOS_params(argv[1],&eos_params);

  if( test_type == 0 ) {
    // Optically thin gas test
    const REAL rho_b = strtod(argv[3],NULL);
    const REAL Y_e   = strtod(argv[4],NULL);
    const REAL T     = strtod(argv[5],NULL);
    OpticallyThinGas_NRPyLeakage(USE_HARM_CONSTANTS,&eos_params,rho_b,Y_e,T);
    OpticallyThinGas_harm_leakage(USE_HARM_CONSTANTS,&eos_params,rho_b,Y_e,T);
  }
  else if( test_type == 1 ) {
    ConstantDensitySphere_NRPyLeakage(&eos_params);
  }
  else if( test_type == 2 ) {
    const REAL T = strtod(argv[3],NULL);
    opacities_dependence_on_density_and_electron_fraction(&eos_params,T);
  }
  else if( test_type == 3 ) {
    NRPyEOS_benchmark(&eos_params);
  }
  else if( test_type == 4 ) {
    const REAL Y_e = strtod(argv[3],NULL);
    const REAL S   = strtod(argv[4],NULL);
    NRPyEOS_constant_Ye_and_S_slice(&eos_params,Y_e,S);
  }
  else{
    fprintf(stderr,"(NRPyLeakage) Unknown test type %d\n",test_type);
    exit(1);
  }

  // Free memory allocated for EOS table
  NRPyEOS_free_memory(&eos_params);

  return 0;
}
