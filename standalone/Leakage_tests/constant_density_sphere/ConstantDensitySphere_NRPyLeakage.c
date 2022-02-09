#include "Basic_defines.h"
#include "NRPyEOS.h"
#include "NRPyLeakage.h"

void ConstantDensitySphere_NRPyLeakage(const NRPyEOS_params *restrict eos_params) {

  // Step 1: Set basic parameters
  const int N0     = 32;
  const int N1     = N0;
  const int N2     = N0;
  const int Ng0    =  2;
  const int Ng1    = Ng0;
  const int Ng2    = Ng0;
  const int Nt0    = N0 + 2*Ng0;
  const int Nt1    = N1 + 2*Ng1;
  const int Nt2    = N2 + 2*Ng2;
  const int Ntotal = Nt0 * Nt1 * Nt2;
  const REAL xmax  = +2;
  const REAL xmin  = -2;
  const REAL ymax  = xmax;
  const REAL ymin  = xmin;
  const REAL zmax  = xmax;
  const REAL zmin  = xmin;
  const REAL dx    = (xmax-xmin)/N0;
  const REAL dy    = (xmax-xmin)/N1;
  const REAL dz    = (xmax-xmin)/N2;
  const REAL rSph  = 1;

  // Step 2: Allocate memory for the metric, opacities, and optical depths
  REAL *gammaDD00 = (REAL *)malloc(sizeof(REAL)*Ntotal);
  REAL *gammaDD11 = (REAL *)malloc(sizeof(REAL)*Ntotal);
  REAL *gammaDD22 = (REAL *)malloc(sizeof(REAL)*Ntotal);
  REAL *kappa_nue[2], *kappa_anue[2], *kappa_nux[2];
  REAL *tau_nue[2], *tau_anue[2], *tau_nux[2];
  for(int i=0;i<2;i++) {
   kappa_nue [i] = (REAL *)malloc(sizeof(REAL)*Ntotal);
   kappa_anue[i] = (REAL *)malloc(sizeof(REAL)*Ntotal);
   kappa_nux [i] = (REAL *)malloc(sizeof(REAL)*Ntotal);
   tau_nue   [i] = (REAL *)malloc(sizeof(REAL)*Ntotal);
   tau_anue  [i] = (REAL *)malloc(sizeof(REAL)*Ntotal);
   tau_nux   [i] = (REAL *)malloc(sizeof(REAL)*Ntotal);
  }

  // Step 3: Define hydro quantities at the sphere interior and exterior
  //         The parameters below are extracted from https://arxiv.org/abs/2106.05356
  // Step 3.a: Interior of the sphere (sec. 4.2.1 of above referece)
  const REAL rho_interior = 9.8e13 * NRPyLeakage_units_cgs_to_geom_D;
  const REAL Y_e_interior = 0.1;
  const REAL T_interior   = 8.0;
  // Step 3.b: Exterior of the sphere (sec. 4.2.1 of above referece)
  const REAL rho_exterior = 6.0e7 * NRPyLeakage_units_cgs_to_geom_D;
  const REAL Y_e_exterior = 0.5;
  const REAL T_exterior   = 0.01;

  // Step 4: Compute opacities in the interior and exterior
  const REAL tau_nue_in [2] = {0.0,0.0};
  const REAL tau_anue_in[2] = {0.0,0.0};
  const REAL tau_nux_in [2] = {0.0,0.0};
  // Step 4.a: Interior
  REAL __attribute__((unused)) R_source,Q_source;
  REAL kappa_nue_interior[2], kappa_anue_interior[2], kappa_nux_interior[2];
  NRPyLeakage_compute_GRMHD_source_terms_and_opacities(USE_NRPY_CONSTANTS,
                                                       eos_params,rho_interior,Y_e_interior,T_interior,
                                                       tau_nue_in,tau_anue_in,tau_nux_in,
                                                       &R_source,&Q_source,
                                                       kappa_nue_interior,kappa_anue_interior,kappa_nux_interior);
  // Step 4.b: Exterior
  REAL kappa_nue_exterior[2], kappa_anue_exterior[2], kappa_nux_exterior[2];
  NRPyLeakage_compute_GRMHD_source_terms_and_opacities(USE_NRPY_CONSTANTS,
                                                       eos_params,rho_exterior,Y_e_exterior,T_exterior,
                                                       tau_nue_in,tau_anue_in,tau_nux_in,
                                                       &R_source,&Q_source,
                                                       kappa_nue_exterior,kappa_anue_exterior,kappa_nux_exterior);

  // Step 5: Print basic information
  fprintf(stderr,"(ConstantDensitySphere) Test information:\n");
  fprintf(stderr,"(ConstantDensitySphere)     Domain properties:\n");
  fprintf(stderr,"(ConstantDensitySphere)         - Sphere radius = %22.15e\n",rSph);
  fprintf(stderr,"(ConstantDensitySphere)         - xmin          = %22.15e\n",xmin);
  fprintf(stderr,"(ConstantDensitySphere)         - xmax          = %22.15e\n",xmax);
  fprintf(stderr,"(ConstantDensitySphere)         - ymin          = %22.15e\n",ymin);
  fprintf(stderr,"(ConstantDensitySphere)         - ymax          = %22.15e\n",ymax);
  fprintf(stderr,"(ConstantDensitySphere)         - zmin          = %22.15e\n",zmin);
  fprintf(stderr,"(ConstantDensitySphere)         - zmax          = %22.15e\n",zmax);
  fprintf(stderr,"(ConstantDensitySphere)         - Nx            = (%d) + 2x(%d)\n",N0,Ng0);
  fprintf(stderr,"(ConstantDensitySphere)         - Ny            = (%d) + 2x(%d)\n",N1,Ng1);
  fprintf(stderr,"(ConstantDensitySphere)         - Nz            = (%d) + 2x(%d)\n",N2,Ng2);
  fprintf(stderr,"(ConstantDensitySphere)     Hydro quantities at sphere interior:\n");
  fprintf(stderr,"(ConstantDensitySphere)         - rho           = %22.15e\n",rho_interior);
  fprintf(stderr,"(ConstantDensitySphere)         - Y_e           = %22.15e\n",Y_e_interior);
  fprintf(stderr,"(ConstantDensitySphere)         -  T            = %22.15e\n",T_interior);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_0_nue   = %22.15e\n",kappa_nue_interior[0]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_1_nue   = %22.15e\n",kappa_nue_interior[1]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_0_anue  = %22.15e\n",kappa_anue_interior[0]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_1_anue  = %22.15e\n",kappa_anue_interior[1]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_0_nux   = %22.15e\n",kappa_nux_interior[0]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_1_nux   = %22.15e\n",kappa_nux_interior[1]);
  fprintf(stderr,"(ConstantDensitySphere)     Hydro quantities at sphere exterior:\n");
  fprintf(stderr,"(ConstantDensitySphere)         - rho           = %22.15e\n",rho_exterior);
  fprintf(stderr,"(ConstantDensitySphere)         - Y_e           = %22.15e\n",Y_e_exterior);
  fprintf(stderr,"(ConstantDensitySphere)         -  T            = %22.15e\n",T_exterior);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_0_nue   = %22.15e\n",kappa_nue_exterior[0]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_1_nue   = %22.15e\n",kappa_nue_exterior[1]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_0_anue  = %22.15e\n",kappa_anue_exterior[0]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_1_anue  = %22.15e\n",kappa_anue_exterior[1]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_0_nux   = %22.15e\n",kappa_nux_exterior[0]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_1_nux   = %22.15e\n",kappa_nux_exterior[1]);

  // Step 3: Loop over the grid and compute the opacities
// #pragma omp parallel for
//   for(int k=0;k<Nt2;k++) {
//     for(int j=0;j<Nt1;j++) {
//       for(int i=0;i<Nt0;i++) {
        
//       }
//     }
//   }
  

  free(gammaDD00);
  free(gammaDD11);
  free(gammaDD22);
  for(int i=0;i<2;i++) {
    free(kappa_nue [i]);
    free(kappa_anue[i]);
    free(kappa_nux [i]);
    free(tau_nue   [i]);
    free(tau_anue  [i]);
    free(tau_nux   [i]);
  }
}
