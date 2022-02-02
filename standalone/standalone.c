#include "Basic_defines.h"

#include "NRPyEOS.h"
#include "NRPyLeakage.h"
#include "Leakage_tests.h"
#include "harm_neutrinos.h"
#include "harm_units.h"

int main(int argc, char **argv) {

  if( argc != 5 ) {
    fprintf(stderr,"(Leakage) Correct usage is ./Leakage_standalone eos_file_path rho Ye T\n");
    exit(1);
  }

  // Read EOS table
  NRPyEOS_params eos_params;
  NRPyEOS_readtable_set_EOS_params(argv[1],&eos_params);

  NRPyEOS_benchmark(&eos_params);

  // // Test Leakage
  // const REAL rho_b       = strtod(argv[2],NULL);
  // const REAL Y_e         = strtod(argv[3],NULL);
  // const REAL T           = strtod(argv[4],NULL);
  // const REAL tau_nue[2]  = {0.0,0.0};
  // const REAL tau_anue[2] = {0.0,0.0};
  // const REAL tau_nux[2]  = {0.0,0.0};

  // REAL P,eps;
  // NRPyEOS_P_and_eps_from_rho_Ye_T(&eos_params,rho_b,Y_e,T,&P,&eps);

  // const int which_constants = USE_HARM_CONSTANTS;

  // // Print info
  // printf("(NRPyLeakage) ******************************************\n");
  // printf("(NRPyLeakage) Basic information:\n");
  // printf("(NRPyLeakage) Density           = %.15e\n",rho_b);
  // printf("(NRPyLeakage) Electron fraction = %.15e\n",Y_e);
  // printf("(NRPyLeakage) Temperature       = %.15e\n",T);
  // printf("(NRPyLeakage) Pressure          = %.15e\n",P);
  // printf("(NRPyLeakage) Energy            = %.15e\n",eps);

  // REAL R_source,Q_source,kappa_nue[2],kappa_anue[2],kappa_nux[2];
  // NRPyLeakage_compute_GRMHD_source_terms_and_opacities(which_constants,&eos_params,rho_b,Y_e,T,tau_nue,tau_anue,tau_nux,
  //                                                      &R_source,&Q_source,kappa_nue,kappa_anue,kappa_nux);
  // printf("(NRPyLeakage) ******************************************\n");
  // printf("(NRPyLeakage) Values obtained using NRPyLeakage:\n");
  // printf("(NRPyLeakage) R_source     = %.15e\n",R_source);
  // printf("(NRPyLeakage) Q_source     = %.15e\n",Q_source);
  // printf("(NRPyLeakage) kappa_0_nue  = %.15e\n",kappa_nue[0]);
  // printf("(NRPyLeakage) kappa_1_nue  = %.15e\n",kappa_nue[1]);
  // printf("(NRPyLeakage) kappa_0_anue = %.15e\n",kappa_anue[0]);
  // printf("(NRPyLeakage) kappa_1_anue = %.15e\n",kappa_anue[1]);
  // printf("(NRPyLeakage) kappa_0_nux  = %.15e\n",kappa_nux[0]);
  // printf("(NRPyLeakage) kappa_1_nux  = %.15e\n",kappa_nux[1]);
  // printf("(NRPyLeakage) rhs_Y_e      = %.15e\n",R_source/rho_b);
  // printf("(NRPyLeakage) rhs_eps      = %.15e\n",Q_source/rho_b);

  // double prims[NUMPRIMS];
  // prims[RHO ] = rho_b;
  // prims[YE  ] = Y_e;
  // prims[TEMP] = T;
  // prims[UU  ] = eps*prims[RHO];
  // for(int i=0;i<3;i++) {
  //   prims[U1+i] = 0.0;
  //   prims[B1+i] = 0.0;
  // }

  // double optical_depth[N_OPTICAL_DEPTHS];
  // for(int i=0;i<N_OPTICAL_DEPTHS;i++) optical_depth[i] = 0.0;

  // double R_function, Q_function;
  // neutrino_absorption_heating_rate(&eos_params, prims, optical_depth, &R_function, &Q_function);
  // printf("(NRPyLeakage) ******************************************\n");
  // printf("(NRPyLeakage) Values obtained using HARM:\n");
  // printf("(NRPyLeakage) R_source = %.15e\n",R_function);
  // printf("(NRPyLeakage) Q_source = %.15e\n",Q_function);
  // printf("(NRPyLeakage) rhs_Y_e  = %.15e\n",R_function/rho_b);
  // printf("(NRPyLeakage) rhs_eps  = %.15e\n",Q_function/rho_b);
  // printf("(NRPyLeakage) ******************************************\n");
  // printf("(NRPyLeakage) Relative errors:\n");
  // printf("(NRPyLeakage) R_source: %e\n",fabs(1.0-R_source/R_function));
  // printf("(NRPyLeakage) Q_source: %e\n",fabs(1.0-Q_source/Q_function));
  // printf("(NRPyLeakage) ******************************************\n");
  

  // Now compute optically thin
  // OpticallyThinGas_NRPyLeakage(which_constants,&eos_params,rho_b,Y_e,T);
  // OpticallyThinGas_harm_leakage(which_constants,&eos_params,rho_b,Y_e,T);

  // Free memory allocated for EOS table
  NRPyEOS_free_memory(&eos_params);

  return 0;
}
