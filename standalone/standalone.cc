#include "Basic_defines.hh"
#include "Leakage.hh"
#include "WVU_EOS_Tabulated_headers.hh"

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
  const REAL Y_e      = 0.5;
  const REAL T        = 1.0;
  const REAL tau_nue  = 1e100;
  const REAL tau_anue = 1e100;
  REAL R_free_total_nue,R_free_total_anue,R_free_total_nux;
  REAL Q_free_total_nue,Q_free_total_anue,Q_free_total_nux;
  Leakage_compute_free_rates(rho_b,Y_e,T,tau_nue,tau_anue,
                             &R_free_total_nue,&R_free_total_anue,&R_free_total_nux,
                             &Q_free_total_nue,&Q_free_total_anue,&Q_free_total_nux);

  printf("(Leakage) R_free_total_nue  = %.15e\n",R_free_total_nue);
  printf("(Leakage) R_free_total_anue = %.15e\n",R_free_total_anue);
  printf("(Leakage) R_free_total_nux  = %.15e\n",R_free_total_nux);
  printf("(Leakage) Q_free_total_nue  = %.15e\n",Q_free_total_nue);
  printf("(Leakage) Q_free_total_anue = %.15e\n",Q_free_total_anue);
  printf("(Leakage) Q_free_total_nux  = %.15e\n",Q_free_total_nux);

  printf("(Leakage) R_source          = %.15e\n",(R_free_total_anue-R_free_total_nue)*Leakage::units_cgs_to_geom_density*Leakage::amu/rho_b);

  // printf("(Leakage) R_beta_nue         = %.15e\n",R_beta_nue);
  // printf("(Leakage) R_beta_anue        = %.15e\n",R_beta_anue);
  // printf("(Leakage) R_pair_nue_anue    = %.15e\n",R_pair_nue_anue);
  // printf("(Leakage) R_pair_nux_anux    = %.15e\n",R_pair_nux_anux);
  // printf("(Leakage) R_plasmon_nue_anue = %.15e\n",R_plasmon_nue_anue);
  // printf("(Leakage) R_plasmon_nux_anux = %.15e\n",R_plasmon_nux_anux);
  // printf("(Leakage) R_Brems_nui_anui   = %.15e\n",R_Brems_nui_anui);

  // printf("(Leakage) Q_beta_nue         = %.15e\n",Q_beta_nue);
  // printf("(Leakage) Q_beta_anue        = %.15e\n",Q_beta_anue);
  // printf("(Leakage) Q_pair_nue_anue    = %.15e\n",Q_pair_nue_anue);
  // printf("(Leakage) Q_pair_nux_anux    = %.15e\n",Q_pair_nux_anux);
  // printf("(Leakage) Q_plasmon_nue_anue = %.15e\n",Q_plasmon_nue_anue);
  // printf("(Leakage) Q_plasmon_nux_anux = %.15e\n",Q_plasmon_nux_anux);
  // printf("(Leakage) Q_Brems_nui_anui   = %.15e\n",Q_Brems_nui_anui);

  // Free memory allocated for EOS table
  WVU_EOS_free_memory();

  return 0;
}
