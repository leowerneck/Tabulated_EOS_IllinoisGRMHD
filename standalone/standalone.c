#include "Basic_defines.h"

#include "NRPyEOS.h"
#include "NRPyLeakage.h"
#include "Leakage_tests.h"
#include "harm_neutrinos.h"
#include "harm_units.h"

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
  const REAL Y_e      = 0.5;
  const REAL T        = 1.0;
  const REAL tau_nue  = 0.0;
  const REAL tau_anue = 0.0;
  const REAL tau_nux  = 0.0;

  REAL P,eps;
  NRPyEOS_P_and_eps_from_rho_Ye_T(&eos_params,rho_b,Y_e,T,&P,&eps);

  const int which_constants = USE_HARM_CONSTANTS;

  REAL cgs_to_geom_density, cgs_to_geom_R, cgs_to_geom_Q;
  if( which_constants == USE_NRPY_CONSTANTS ) {
    cgs_to_geom_density = NRPyLeakage_units_cgs_to_geom_density;
    cgs_to_geom_R       = NRPyLeakage_units_cgs_to_geom_R;
    cgs_to_geom_Q       = NRPyLeakage_units_cgs_to_geom_Q;
  }
  else if( which_constants == USE_HARM_CONSTANTS ) {
    cgs_to_geom_density = NRPyLeakage_harm_units_cgs_to_geom_density;
    cgs_to_geom_R       = NRPyLeakage_harm_units_cgs_to_geom_R;
    cgs_to_geom_Q       = NRPyLeakage_harm_units_cgs_to_geom_Q;
  }
  else {
    fprintf(stderr,"(NRPyLeakage) Unsupported constants %d\n",which_constants);
    exit(1);
  }

  printf("(NRPyLeakage) **********************************************\n");
  printf("(NRPyLeakage) Units conversion factors:\n");
  printf("(NRPyLeakage) Density [cgs to  code]: %.15e\n",cgs_to_geom_density);
  printf("(NRPyLeakage) R       [cgs to  code]: %.15e\n",cgs_to_geom_R);
  printf("(NRPyLeakage) Q       [cgs to  code]: %.15e\n",cgs_to_geom_Q);
  printf("(NRPyLeakage) **********************************************\n");
  printf("(HARMLeakage) **********************************************\n");
  printf("(HARMLeakage) Units conversion factors:\n");
  printf("(HARMLeakage) Density [cgs to  code]: %.15e\n",RHOGF);
  printf("(HARMLeakage) R       [cgs to  code]: %.15e\n",R_CONV);
  printf("(HARMLeakage) Q       [cgs to  code]: %.15e\n",Q_CONV);
  printf("(HARMLeakage) **********************************************\n");
  
  // Print info
  printf("(NRPyLeakage) Density           = %.15e\n",rho_b);
  printf("(NRPyLeakage) Electron fraction = %.15e\n",Y_e);
  printf("(NRPyLeakage) Temperature       = %.15e\n",T);
  printf("(NRPyLeakage) Pressure          = %.15e\n",P);
  printf("(NRPyLeakage) Energy            = %.15e\n",eps);

  REAL R_source,Q_source;
  NRPyLeakage_compute_GRMHD_source_terms(which_constants,&eos_params,rho_b,Y_e,T,tau_nue,tau_anue,tau_nux,&R_source,&Q_source);
  printf("(NRPyLeakage) R_source          = %.15e\n",R_source);
  printf("(NRPyLeakage) Q_source          = %.15e\n",Q_source);
  printf("(NRPyLeakage) rhs_Y_e           = %.15e\n",R_source/rho_b);
  printf("(NRPyLeakage) rhs_eps           = %.15e\n",Q_source/rho_b);


  double prims[NUMPRIMS];
  prims[RHO ] = rho_b;
  prims[YE  ] = Y_e;
  prims[TEMP] = T;
  prims[UU  ] = eps*prims[RHO];
  for(int i=0;i<3;i++) {
    prims[U1+i] = 0.0;
    prims[B1+i] = 0.0;
  }

  double optical_depth[N_OPTICAL_DEPTHS];
  for(int i=0;i<N_OPTICAL_DEPTHS;i++) optical_depth[i] = 0.0;

  double R_function, Q_function;
  neutrino_absorption_heating_rate(&eos_params, prims, optical_depth, &R_function, &Q_function);
  printf("(HARMLeakage) **********************************************\n");
  printf("(HARMLeakage) R_source          = %.15e\n",R_function);
  printf("(HARMLeakage) Q_source          = %.15e\n",Q_function);
  printf("(HARMLeakage) rhs_Y_e           = %.15e\n",R_function/rho_b);
  printf("(HARMLeakage) rhs_eps           = %.15e\n",Q_function/rho_b);
  

  // Now compute optically thin
  // Leakage_optically_thin_regime_semi_analytic(rho_b,Y_e,T);

  // Free memory allocated for EOS table
  NRPyEOS_free_memory(&eos_params);

  return 0;
}
