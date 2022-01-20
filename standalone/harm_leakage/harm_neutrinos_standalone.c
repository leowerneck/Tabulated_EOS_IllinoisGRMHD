#include <stdio.h>
#include <stdlib.h>

#include "NRPyEOS.h"
#include "harm_neutrinos.h"

int main(int argc, char **argv) {

  // Step 0: Check correct usage
  if( argc != 2 ) {
    fprintf(stderr,"(harm3d neutrinos standalone) Correct usage is ./harm3d_neutrinos_standalone eos_file_path\n");
    exit(1);
  }

  // Step 1: Read in the EOS
  NRPyEOS_params eos_params;
  NRPyEOS_readtable_set_EOS_params(argv[1],&eos_params);

  // Step 2: Set variables needed by HARM
  // Step 2.a: Primitives
  double prims[NUMPRIMS];
  // Step 2.a.i: First set the base primitives
  prims[RHO ] = 1e-12;
  prims[YE  ] = 0.5;
  prims[TEMP] = 1.0;
  for(int i=0;i<3;i++) {
    prims[U1+i] = 0.0;
    prims[B1+i] = 0.0;
  }
  // Step 2.a.ii: Now set u = eps*rho
  double eps;
  NRPyEOS_eps_from_rho_Ye_T(&eos_params,prims[RHO],prims[YE],prims[TEMP],&eps);
  prims[UU] = eps*prims[RHO];

  // Step 2.b: Set optical depth to zero
  double optical_depth[N_OPTICAL_DEPTHS];
  for(int i=0;i<N_OPTICAL_DEPTHS;i++) optical_depth[i] = 0.0;

  // Step 3: Compute leakage source terms
  double R_function, Q_function;
  neutrino_absorption_heating_rate(&eos_params, prims, optical_depth, &R_function, &Q_function);

  // Step 4: Print information about the test
  printf("(harm) Density    : %.15e g/cm^3\n",prims[RHO ] * INVRHOGF);
  printf("(harm) Temperature: %.15e\n"       ,prims[TEMP]);
  printf("(harm) e- fraction: %.15e MeV\n"   ,prims[YE  ]);
  printf("(harm) Energy     : %.15e\n"       ,prims[UU  ]);
  printf("(harm) Velocity x : %.15e\n"       ,prims[U1  ]);
  printf("(harm) Velocity y : %.15e\n"       ,prims[U2  ]);
  printf("(harm) Velocity z : %.15e\n"       ,prims[U3  ]);
  printf("(harm) B-field  x : %.15e\n"       ,prims[B1  ]);
  printf("(harm) B-field  y : %.15e\n"       ,prims[B2  ]);
  printf("(harm) B-field  z : %.15e\n"       ,prims[B3  ]);
  printf("(harm) R_function : %.15e\n"       ,R_function/prims[RHO]);
  printf("(harm) Q_function : %.15e\n"       ,Q_function/prims[RHO]);

  // Step 5: Free memory
  NRPyEOS_free_memory(&eos_params);

  return 0;
}
