#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NRPyEOS.h"

// This function slices a 3D EOS table assuming
// a constant electron fraction and entropy. The
// output is a file containing a 1D table with
// three columns:
//
//   1)  h       - the specific enthalpy
//   2) log(rho) - natural log of density
//   3) log(T)   - natural log of temperature
//
// Using this file one can determine (rho,T) as
// a function of the specific enthalpy through
// interpolations.

void NRPyEOS_constant_Ye_and_S_slice( const NRPyEOS_params *restrict eos_params,
                                      const double Y_e,
                                      const double S ) {

  printf("(NRPyEOS) Slicing table at:\n");
  printf("(NRPyEOS)    S  = %g\n",S);
  printf("(NRPyEOS)   Y_e = %g\n",Y_e);

  // Step 1: Open the output file
  char filename[256];
  sprintf(filename,"constant_Ye_and_S_slice__%g_%g.txt",Y_e,S);
  FILE *fp = fopen(filename,"w");
  fprintf(fp,"# NRPyEOS table slicer\n");
  fprintf(fp,"# Slicing table at:\n");
  fprintf(fp,"#    S  = %g\n",S);
  fprintf(fp,"#   Y_e = %g\n",Y_e);
  fprintf(fp,"# All quantities below in geometrized (G=c=Msun=1) units:\n");
  fprintf(fp,"#   Column 1: log10 specific enthalpy\n");
  fprintf(fp,"#   Column 2: log10 baryonic density\n");
  fprintf(fp,"#   Column 3: log10 temperature (in MeV)\n");

  // Step 3: Loop over the density and output local quantities to file
  for(int ir=0;ir<eos_params->nrho;ir++) {

    // Step 3.a: Set the local density
    const double rho = exp(eos_params->logrho[ir]);

    // Step 3.b: Set the temperature guess to 1.0
    double T = 0.01;

    // Step 3.c: Given (rho,Y_e,S), determine T, P, and eps
    double P, eps;
    NRPyEOS_P_eps_and_T_from_rho_Ye_S(eos_params, rho, Y_e, S, &P, &eps, &T);

    // Step 3.d: Compute the local enthalpy
    const double h = 1 + eps + P/rho;

    // Step 3.e: Output to file
    fprintf(fp,"%.15e %.15e %.15e\n",log10(h),log10(rho),log10(T));
  }

  fclose(fp);

  printf("(NRPyEOS) Finished slicing table. Output file: %s\n",filename);
}
