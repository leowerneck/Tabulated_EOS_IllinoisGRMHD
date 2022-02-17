#include "Basic_defines.h"

#include "NRPyEOS.h"
#include "NRPyLeakage.h"
#include "Leakage_tests.h"
#include "harm_neutrinos.h"
#include "harm_units.h"

void nrpyeos_fortran_interface_(const NRPyEOS_params *eos_params,
                                const REAL *rho,
                                const REAL *Y_e,
                                const REAL *T,
                                REAL *P,
                                REAL *eps);

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

int main(int argc, char **argv) {

  if( argc != 5 ) {
    fprintf(stderr,"(Leakage) Correct usage is ./Leakage_standalone eos_file_path rho Ye T\n");
    exit(1);
  }

  // Read EOS table
  NRPyEOS_params eos_params;
  NRPyEOS_readtable_set_EOS_params(argv[1],&eos_params);

  // NRPyEOS_benchmark(&eos_params);

  // ConstantDensitySphere_NRPyLeakage(&eos_params);

  // Test Leakage
  const REAL rho_b = strtod(argv[2],NULL); // 0.00141722994401086;//6.0e7*NRPyLeakage_units_cgs_to_geom_D;//
  const REAL Y_e   = strtod(argv[3],NULL); // 0.057773852636472;//0.5;//
  const REAL T     = strtod(argv[4],NULL); // 0.01;//
  REAL tau_nue[2]  = {0.0,0.0};
  REAL tau_anue[2] = {0.0,0.0};
  REAL tau_nux[2]  = {0.0,0.0};

  REAL P,eps;
  NRPyEOS_P_and_eps_from_rho_Ye_T(&eos_params,rho_b,Y_e,T,&P,&eps);

  const int which_constants = USE_HARM_CONSTANTS;

  // Print info
  printf("(NRPyLeakage) ******************************************\n");
  printf("(NRPyLeakage) Basic information:\n");
  printf("(NRPyLeakage) Density           = %.15e\n",rho_b);
  printf("(NRPyLeakage) Electron fraction = %.15e\n",Y_e);
  printf("(NRPyLeakage) Temperature       = %.15e\n",T);
  printf("(NRPyLeakage) Pressure          = %.15e\n",P);
  printf("(NRPyLeakage) Energy            = %.15e\n",eps);

  P = 0.0/0.0;
  eps = 0.0/0.0;
  nrpyeos_fortran_interface_(&eos_params,&rho_b,&Y_e,&T,&P,&eps);
  // simplest_fortran_interface_(&rho_b,&Y_e,&T,&P,&eps);
  printf("(NRPyLeakage) ******************************************\n");
  printf("(NRPyLeakage) Basic information (FORTRAN interface):\n");
  printf("(NRPyLeakage) Density           = %.15e\n",rho_b);
  printf("(NRPyLeakage) Electron fraction = %.15e\n",Y_e);
  printf("(NRPyLeakage) Temperature       = %.15e\n",T);
  printf("(NRPyLeakage) Pressure          = %.15e\n",P);
  printf("(NRPyLeakage) Energy            = %.15e\n",eps);

  printf("(NRPyLeakage) ******************************************\n");
  printf("(NRPyLeakage) Values obtained using ZelmaniLeak:\n");
  const int nzones       = 1;
  const REAL eos_tempmin = eos_params.eos_tempmin;
  const REAL atm_temp    = eos_params.eos_tempmin;
  const REAL eos_yemin   = eos_params.eos_yemin;
  const REAL eos_yemax   = eos_params.eos_yemax;
  const REAL rho_gf      = NRPyLeakage_units_geom_to_cgs_D;
  const REAL kappa_gf    = NRPyLeakage_units_geom_to_cgs_L;
  REAL rho               = rho_b * NRPyLeakage_units_cgs_to_geom_D;
  REAL temp              = T;
  REAL ye                = Y_e;
  REAL r                 = 1;
  REAL ds                = 1;
  REAL oldtauruff[3][nzones], tauruff[3][nzones], chiross[3][nzones], heatflux[3][nzones], compos[4][nzones];;
  for(int i=0;i<nzones;i++) {
    for(int n=0;n<3;n++) {
      oldtauruff[n][i] = 0;
      tauruff   [n][i] = 0;
      chiross   [n][i] = 0;
      heatflux  [n][i] = 0;
      compos    [n][i] = 0;
    }
    compos[3][i] = 0;
  }
  REAL heaterms[3]={0,0,0}, heateave[3]={0,0,0};
  REAL S=0;
  calc_taus_(&eos_params,&eos_tempmin,&atm_temp,&eos_yemax,&eos_yemin,&rho_gf,&kappa_gf,
             &rho,&temp,&ye,oldtauruff,tauruff,chiross,heatflux,heaterms,heateave,
             &nzones,
             &r,&ds,compos,&S);

  // printf("compos: %e %e %e %e\n",compos[0][0],compos[1][0],compos[2][0],compos[3][0]);

  REAL R_source,Q_source,kappa_nue[2],kappa_anue[2],kappa_nux[2];
  NRPyLeakage_compute_GRMHD_source_terms_and_opacities(which_constants,&eos_params,rho_b,Y_e,T,tau_nue,tau_anue,tau_nux,
                                                       &R_source,&Q_source,kappa_nue,kappa_anue,kappa_nux);
  printf("(NRPyLeakage) ******************************************\n");
  printf("(NRPyLeakage) Values obtained using NRPyLeakage:\n");
  printf("(NRPyLeakage) R_source     = %.15e\n",R_source);
  printf("(NRPyLeakage) Q_source     = %.15e\n",Q_source);
  printf("(NRPyLeakage) kappa_0_nue  = %.15e\n",kappa_nue[0]);
  printf("(NRPyLeakage) kappa_1_nue  = %.15e\n",kappa_nue[1]);
  printf("(NRPyLeakage) kappa_0_anue = %.15e\n",kappa_anue[0]);
  printf("(NRPyLeakage) kappa_1_anue = %.15e\n",kappa_anue[1]);
  printf("(NRPyLeakage) kappa_0_nux  = %.15e\n",kappa_nux[0]);
  printf("(NRPyLeakage) kappa_1_nux  = %.15e\n",kappa_nux[1]);
  printf("(NRPyLeakage) rhs_Y_e      = %.15e\n",R_source/rho_b);
  printf("(NRPyLeakage) rhs_eps      = %.15e\n",Q_source/rho_b);

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
