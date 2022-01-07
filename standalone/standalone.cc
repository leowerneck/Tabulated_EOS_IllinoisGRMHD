#include "Basic_defines.hh"
#include "Leakage.hh"
#include "WVU_EOS_Tabulated_headers.hh"

#include "Leakage_tests.hh"

void Leakage_compute_free_rates(const REAL rho_b,
                const REAL Y_e,
                const REAL T,
                const REAL tau_nue,
                const REAL tau_anue,
                REAL *restrict R_free_total_nue ,
                REAL *restrict R_free_total_anue,
                REAL *restrict R_free_total_nux ,
                REAL *restrict Q_free_total_nue ,
                REAL *restrict Q_free_total_anue,
                REAL *restrict Q_free_total_nux);

int main(int argc, char **argv) {

  if( argc != 2 ) {
    fprintf(stderr,"(Leakage) Correct usage is ./Leakage_standalone eos_file_path\n");
    exit(1);
  }

  // Read EOS table
  WVU_EOS_ReadTable(argv[1]);

  // Test Leakage
  const REAL rho_b    = 1e-12;
  const REAL Y_e      = 0.05;
  const REAL T        = 1.0;
  const REAL tau_nue  = 0.0;
  const REAL tau_anue = 0.0;
  const REAL tau_nux  = 0.0;

  REAL P,epsilon;
  WVU_EOS_P_and_eps_from_rho_Ye_T(rho_b,Y_e,T,&P,&epsilon);
  
  // Print info
  printf("(Leakage) Density           = %.15e\n",rho_b);
  printf("(Leakage) Electron fraction = %.15e\n",Y_e);
  printf("(Leakage) Temperature       = %.15e\n",T);
  printf("(Leakage) Pressure          = %.15e\n",P);
  printf("(Leakage) Energy            = %.15e\n",epsilon);

  REAL R_source,Q_source;
  Leakage_compute_GRMHD_source_terms(rho_b,Y_e,T,tau_nue,tau_anue,tau_nux,&R_source,&Q_source);
  printf("(Leakage) R_source          = %.15e\n",R_source);
  printf("(Leakage) Q_source          = %.15e\n",Q_source);
  printf("(Leakage) rhs_Y_e           = %.15e\n",R_source/rho_b);
  printf("(Leakage) rhs_eps           = %.15e\n",Q_source/rho_b);

  // Now compute optically thin
  Leakage_optically_thin_regime_semi_analytic(rho_b,Y_e,T);

  // Free memory allocated for EOS table
  WVU_EOS_free_memory();

  return 0;
}
