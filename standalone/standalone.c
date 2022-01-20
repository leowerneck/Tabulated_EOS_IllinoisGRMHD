#include "Basic_defines.h"

#include "NRPyEOS.h"
#include "NRPyLeakage.h"
#include "Leakage_tests.h"

int main(int argc, char **argv) {

  if( argc != 2 ) {
    fprintf(stderr,"(Leakage) Correct usage is ./Leakage_standalone eos_file_path\n");
    exit(1);
  }

  // Read EOS table
  NRPyEOS_params eos_params;
  NRPyEOS_readtable_set_EOS_params(argv[1],&eos_params);

  // Test Leakage
  const REAL rho_b    = 1e-12;
  const REAL Y_e      = 0.05;
  const REAL T        = 1.0;
  const REAL tau_nue  = 0.0;
  const REAL tau_anue = 0.0;
  const REAL tau_nux  = 0.0;

  REAL P,epsilon;
  NRPyEOS_P_and_eps_from_rho_Ye_T(&eos_params,rho_b,Y_e,T,&P,&epsilon);
  
  // Print info
  printf("(Leakage) Density           = %.15e\n",rho_b);
  printf("(Leakage) Electron fraction = %.15e\n",Y_e);
  printf("(Leakage) Temperature       = %.15e\n",T);
  printf("(Leakage) Pressure          = %.15e\n",P);
  printf("(Leakage) Energy            = %.15e\n",epsilon);

  REAL R_source,Q_source;
  NRPyLeakage_compute_GRMHD_source_terms(&eos_params,rho_b,Y_e,T,tau_nue,tau_anue,tau_nux,&R_source,&Q_source);
  printf("(Leakage) R_source          = %.15e\n",R_source);
  printf("(Leakage) Q_source          = %.15e\n",Q_source);
  printf("(Leakage) rhs_Y_e           = %.15e\n",R_source/rho_b);
  printf("(Leakage) rhs_eps           = %.15e\n",Q_source/rho_b);

  // Now compute optically thin
  // Leakage_optically_thin_regime_semi_analytic(rho_b,Y_e,T);

  // Free memory allocated for EOS table
  NRPyEOS_free_memory(&eos_params);

  return 0;
}
