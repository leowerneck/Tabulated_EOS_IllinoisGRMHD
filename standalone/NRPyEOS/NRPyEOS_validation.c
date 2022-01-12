#include "NRPyEOS.h"

void NRPyEOS_validation(const NRPyEOS_params *restrict eos_params, const int N, const double Ye, const char *restrict filename) {
  // Step 1: Set test parameters
  const int keys[NRPyEOS_ntablekeys] = {
    NRPyEOS_press_key,NRPyEOS_eps_key,NRPyEOS_entropy_key,
    NRPyEOS_munu_key,NRPyEOS_cs2_key,NRPyEOS_depsdT_key,
    NRPyEOS_dPdrho_key,NRPyEOS_dPdeps_key,NRPyEOS_muhat_key,
    NRPyEOS_mu_e_key,NRPyEOS_mu_p_key,NRPyEOS_mu_n_key,
    NRPyEOS_Xa_key,NRPyEOS_Xh_key,NRPyEOS_Xn_key,NRPyEOS_Xp_key,
    NRPyEOS_Abar_key,NRPyEOS_Zbar_key,NRPyEOS_Gamma_key
  };
  NRPyEOS_error_report report;
  double outvars[NRPyEOS_ntablekeys];

  // Step 2.a: Print information about the test
  printf("(NRPyEOS) Performing test for Ye = %.2lf\n",Ye);
  printf("(NRPyEOS) The following quantities will be interpolated:\n");
  for(int i=0;i<NRPyEOS_ntablekeys;i++)
    printf("(NRPyEOS) key %2d: %s\n",i,table_var_names[i]);

  // Step 3: Perform EOS test
  FILE *fp = fopen(filename,"w");
  for(int j=0;j<N;j++) {
    // Step 3.a: Set local temperature
    const double T = eos_params->eos_tempmin + j*eos_params->dtemp;
    for(int i=0;i<N;i++) {
      // Step 3.b: Set local density
      const double rho = eos_params->eos_rhomin + i*eos_params->drho;

      // Step 3.c: Interpolate all EOS quantities
      NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, NRPyEOS_ntablekeys,rho,Ye,T, keys,outvars, &report );

      // Step 3.d: Output results to file
      fprintf(fp,"%d %d %.15e %.15e",i,j,rho,T);
      for(int nn=0;nn<NRPyEOS_ntablekeys;nn++)
        fprintf(fp," %.15e",outvars[nn]);
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}

int main(int argc, char **argv) {

  // Step 0: Check correct usage
  if( argc != 2 ) {
    fprintf(stderr,"(NRPyEOS_validation) Correct usage is ./NRPyEOS_validation eos_file_path\n");
    exit(1);
  }

  // Step 1: Read the EOS table; initialize the EOS parameters struct
  NRPyEOS_params eos_params;
  NRPyEOS_readtable_set_EOS_params(argv[1],&eos_params);

  // Step 2: Perform validation tests
  // Step 2.a: Y_e = 0.5
  NRPyEOS_validation(&eos_params,100,0.5,"NRPyEOS_test_Y_e_0.5.txt");
  // Step 2.b: Y_e = 0.05
  NRPyEOS_validation(&eos_params,100,0.05,"NRPyEOS_test_Y_e_0.05.txt");

  // Step 3: Free memory
  NRPyEOS_free_memory(&eos_params);

  return 0;
}
