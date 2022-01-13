#include <stdio.h>
#include <stdlib.h>

#include "NRPyEOS.h"
#include "harm3d_neutrinos.h"

int main(int argc, char **argv) {

  // Step 0: Check correct usage
  if( argc != 2 ) {
    fprintf(stderr,"(harm3d neutrinos standalone) Correct usage is ./harm3d_neutrinos_standalone eos_file_path\n");
    exit(1);
  }

  // Step 1: Read in the EOS
  NRPyEOS_params eos_params;
  NRPyEOS_readtable_set_EOS_params(argv[1],&eos_params);

  // Step 2: Variables needed by harm3d_neutrinos
  // Step 2.a: Grid point indices
  int ii=0,jj=0,kk=0;

  // Step 2.b: Set metric quantities to flat space
  // Step 2.b.i: Four-metric
  struct of_geom geom;
  geom.gcon[0][0] = -1.0;
  geom.gcov[0][0] = -1.0;
  for(int i=1;i<4;i++) {
    geom.gcon[0][i] = 0.0;
    geom.gcon[i][0] = 0.0;
    geom.gcon[i][i] = 1.0;
    geom.gcov[0][i] = 0.0;
    geom.gcov[i][0] = 0.0;
    geom.gcov[i][i] = 1.0;
    for(int j=i+1;j<4;j++) {
      geom.gcon[i][j] = 0.0;
      geom.gcon[j][i] = 0.0;
      geom.gcov[i][j] = 0.0;
      geom.gcov[j][i] = 0.0;
    }
  }
  // Step 2.b.ii: Four-metric determinant
  geom.g     = 1.0;
  geom.g_inv = 1.0;
  // Step 2.b.iii: Lapse function
  geom.alpha = 1.0;
  // Step 2.b.iv: Shift vector
  for(int i=0;i<3;i++) geom.beta[i] = 0.0;
  // Step 2.b.v: Unit normal vector
  geom.ncon[0] = 1.0;
  for(int i=1;i<4;i++) geom.ncon[i] = 0.0;

  // Step 2.c: Set primitives
  double prims[NUMPRIMS];
  // Step 2.c.i: First set the base primitives
  prims[RHO ] = 1e-12;
  prims[YE  ] = 0.5;
  prims[TEMP] = 1.0;
  for(int i=0;i<3;i++) {
    prims[U1+i] = 0.0;
    prims[B1+i] = 0.0;
  }
  // Step 2.c.ii: Now set u = eps*rho
  double eps;
  NRPyEOS_eps_from_rho_Ye_T(&eos_params,prims[RHO],prims[YE],prims[TEMP],&eps);
  prims[UU] = eps*prims[RHO];
  // Step 2.d: Print primitive information
  printf("(harm3d) Density    : %.15e g/cm^3\n",prims[RHO ] * INVRHOGF);
  printf("(harm3d) Temperature: %.15e\n"       ,prims[TEMP]);
  printf("(harm3d) e- fraction: %.15e MeV\n"   ,prims[YE  ]);
  printf("(harm3d) Energy     : %.15e\n"       ,prims[UU  ]);
  printf("(harm3d) Velocity x : %.15e\n"       ,prims[U1  ]);
  printf("(harm3d) Velocity y : %.15e\n"       ,prims[U2  ]);
  printf("(harm3d) Velocity z : %.15e\n"       ,prims[U3  ]);
  printf("(harm3d) B-field  x : %.15e\n"       ,prims[B1  ]);
  printf("(harm3d) B-field  y : %.15e\n"       ,prims[B2  ]);
  printf("(harm3d) B-field  z : %.15e\n"       ,prims[B3  ]);
  
  // double *optical_depth;
  // double *dU;
  // struct of_state q;
  

  // neutrino_source_func(&eos_params,ph,optical_depth,&q,&geom,ii,jj,kk,dU);

  // Free memory
  NRPyEOS_free_memory(&eos_params);

  return 0;
}
