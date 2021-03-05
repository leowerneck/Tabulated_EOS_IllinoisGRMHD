#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "util_Table.h"
#include <assert.h>
#include <math.h>
#include <string.h>

#include "integrands.h"

void VolumeIntegrals_vacuum_ComputeIntegrand(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int which_integral = NumIntegrals - *IntegralCounter + 1;

  // ---------------- Get BSSN constraints gridfunctions ----------------
  const CCTK_INT timelevel = 0;
  // Get pointer to Hamiltonian constraint gridfunction
  CCTK_REAL* H_gf   = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel, HamiltonianVarString));
  if(  !H_gf  ) CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING, "Couldn't get data pointer of input array variable '%s'", HamiltonianVarString);

  // Get pointer to Momentum constraint gridfunction in the x-direction
  CCTK_REAL* MU0_gf = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel, Momentum0VarString));;
  if( !MU0_gf ) CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING, "Couldn't get data pointer of input array variable '%s'", Momentum0VarString);

  // Get pointer to Momentum constraint gridfunction in the y-direction
  CCTK_REAL* MU1_gf = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel, Momentum1VarString));;
  if( !MU1_gf ) CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING, "Couldn't get data pointer of input array variable '%s'", Momentum1VarString);

  // Get pointer to Momentum constraint gridfunction in the z-direction
  CCTK_REAL* MU2_gf = (CCTK_REAL*)(CCTK_VarDataPtr(cctkGH,timelevel, Momentum2VarString));;
  if( !MU2_gf ) CCTK_VError (__LINE__, __FILE__, CCTK_THORNSTRING, "Couldn't get data pointer of input array variable '%s'", Momentum2VarString);
  // --------------------------------------------------------------------

  /* Note: Must extend this if/else statement if adding a new integrand! */
  if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"H_M_CnstraintsL2")) { 
#pragma omp parallel for collapse(3)
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  VolumeIntegrals_vacuum_L2_integrand(VolIntegrand1, index,H_gf,x,y,z);
	  VolumeIntegrals_vacuum_L2_integrand(VolIntegrand2, index,MU0_gf,x,y,z);
	  VolumeIntegrals_vacuum_L2_integrand(VolIntegrand3, index,MU1_gf,x,y,z);
	  VolumeIntegrals_vacuum_L2_integrand(VolIntegrand4, index,MU2_gf,x,y,z); }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"centeroflapse")) { 
#pragma omp parallel for collapse(3)
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  VolumeIntegrals_vacuum_CoL_integrand(VolIntegrand1,VolIntegrand2,VolIntegrand3,VolIntegrand4, index,alp,x,y,z); }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"usepreviousintegrands")) {
    /* Do Nothing; the action for this is below. */
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"one")) {
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { int index = CCTK_GFINDEX3D(cctkGH,i,j,k); VolIntegrand1[index]=1.0; }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"ADM_Mass")) {
    double idx = 1.0/CCTK_DELTA_SPACE(0);
    double idy = 1.0/CCTK_DELTA_SPACE(1);
    double idz = 1.0/CCTK_DELTA_SPACE(2);
#pragma omp parallel for collapse(3)
    for(int k=1;k<cctk_lsh[2]-1;k++) for(int j=1;j<cctk_lsh[1]-1;j++) for(int i=1;i<cctk_lsh[0]-1;i++) {
          VolumeIntegrals_vacuum_ADM_Mass_integrand_eval_derivs(VolIntegrand2,VolIntegrand3,VolIntegrand4, cctkGH, i,j,k, idx,idy,idz, alp,gxx,gxy,gxz,gyy,gyz,gzz);
        }
#pragma omp parallel for collapse(3)
    for(int k=2;k<cctk_lsh[2]-2;k++) for(int j=2;j<cctk_lsh[1]-2;j++) for(int i=2;i<cctk_lsh[0]-2;i++) {
          VolumeIntegrals_vacuum_ADM_Mass_integrand(VolIntegrand1, cctkGH, i,j,k, idx,idy,idz, VolIntegrand2, VolIntegrand3, VolIntegrand4);
        }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"ADM_Momentum")) {
    double idx = 1.0/CCTK_DELTA_SPACE(0);
    double idy = 1.0/CCTK_DELTA_SPACE(1);
    double idz = 1.0/CCTK_DELTA_SPACE(2);
#pragma omp parallel for collapse(3)
    for(int k=1;k<cctk_lsh[2]-1;k++) for(int j=1;j<cctk_lsh[1]-1;j++) for(int i=1;i<cctk_lsh[0]-1;i++) {
          VolumeIntegrals_vacuum_ADM_Momentum_integrand_eval_derivs(VolIntegrand2,VolIntegrand3,VolIntegrand4, cctkGH, i,j,k, idx,idy,idz, alp,gxx,gxy,gxz,gyy,gyz,gzz,kxx,kxy,kxz,kyy,kyz,kzz);
        }
#pragma omp parallel for collapse(3)
    for(int k=2;k<cctk_lsh[2]-2;k++) for(int j=2;j<cctk_lsh[1]-2;j++) for(int i=2;i<cctk_lsh[0]-2;i++) {
          VolumeIntegrals_vacuum_ADM_Momentum_integrand(VolIntegrand1, cctkGH, i,j,k, idx,idy,idz, VolIntegrand2,VolIntegrand3,VolIntegrand4);
        }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"ADM_Angular_Momentum")) {
    double idx = 1.0/CCTK_DELTA_SPACE(0);
    double idy = 1.0/CCTK_DELTA_SPACE(1);
    double idz = 1.0/CCTK_DELTA_SPACE(2);
#pragma omp parallel for collapse(3)
    for(int k=1;k<cctk_lsh[2]-1;k++) for(int j=1;j<cctk_lsh[1]-1;j++) for(int i=1;i<cctk_lsh[0]-1;i++) {
          VolumeIntegrals_vacuum_ADM_Angular_Momentum_integrand_eval_derivs(VolIntegrand2,VolIntegrand3,VolIntegrand4, cctkGH, i,j,k, idx,idy,idz, x,y,z, alp,gxx,gxy,gxz,gyy,gyz,gzz,kxx,kxy,kxz,kyy,kyz,kzz);
        }
#pragma omp parallel for collapse(3)
    for(int k=2;k<cctk_lsh[2]-2;k++) for(int j=2;j<cctk_lsh[1]-2;j++) for(int i=2;i<cctk_lsh[0]-2;i++) {
          VolumeIntegrals_vacuum_ADM_Angular_Momentum_integrand(VolIntegrand1, cctkGH, i,j,k, idx,idy,idz, VolIntegrand2,VolIntegrand3,VolIntegrand4);
        }
  } else {
    /* Print a warning if no integrand is computed because Integration_quantity_keyword unrecognized. */
    printf("VolumeIntegrals: WARNING: Integrand not computed. Did not understand Integration_quantity_keyword[%d] = %s\n",which_integral,Integration_quantity_keyword[which_integral]);
  }

  if(cctk_iteration==0) {
    volintegral_inside_sphere__center_x[which_integral]=volintegral_sphere__center_x_initial[which_integral];
    volintegral_inside_sphere__center_y[which_integral]=volintegral_sphere__center_y_initial[which_integral];
    volintegral_inside_sphere__center_z[which_integral]=volintegral_sphere__center_z_initial[which_integral];

    volintegral_outside_sphere__center_x[which_integral]=volintegral_sphere__center_x_initial[which_integral];
    volintegral_outside_sphere__center_y[which_integral]=volintegral_sphere__center_y_initial[which_integral];
    volintegral_outside_sphere__center_z[which_integral]=volintegral_sphere__center_z_initial[which_integral];
  }

  if(volintegral_sphere__tracks__amr_centre[which_integral]!=-1) {
    int which_centre = volintegral_sphere__tracks__amr_centre[which_integral];

    volintegral_inside_sphere__center_x[which_integral] = position_x[which_centre];
    volintegral_inside_sphere__center_y[which_integral] = position_y[which_centre];
    volintegral_inside_sphere__center_z[which_integral] = position_z[which_centre];

    volintegral_outside_sphere__center_x[which_integral] = position_x[which_centre];
    volintegral_outside_sphere__center_y[which_integral] = position_y[which_centre];
    volintegral_outside_sphere__center_z[which_integral] = position_z[which_centre];
  }

  /* ZERO OUT INTEGRATION REGIONS */

  /* The below code supports zeroing out of arbitrary spherical shells. 
     In the case of integration INSIDE a full sphere, this code also supports
     moving spheres if AMR centre tracking is enabled. This can be used to 
     track compact objects. 

     Here's one way:
     Generally the lapse at the center of these objects is minimized. 
     If the object has a length scale of R and is known to be centered at x,y,z, compute 
     X^i = Integral [ (1-lapse) x^i dV] / Integral [(1-lapse) dV]
     over the spherical volume centered at x,y,z with radius R. 
     This should yield X,Y,Z = x,y,z to good approximation. 
     If you enable AMR centre tracking, it should work just as well
     as any other method for tracking compact objects, if not better. */

  /* Set integrands to zero outside a sphere centered at x,y,z. 
     I.e., this results in the integral being restricted INSIDE sphere */
  if(volintegral_inside_sphere__radius[which_integral]>0.0) {
    double radius=volintegral_inside_sphere__radius[which_integral];
    double xprime=volintegral_sphere__center_x_initial[which_integral];
    double yprime=volintegral_sphere__center_y_initial[which_integral];
    double zprime=volintegral_sphere__center_z_initial[which_integral];
    if(cctk_iteration>0 && (amr_centre__tracks__volintegral_inside_sphere[which_integral]!=-1 || volintegral_sphere__tracks__amr_centre[which_integral]!=-1) ) {
      xprime=volintegral_inside_sphere__center_x[which_integral];
      yprime=volintegral_inside_sphere__center_y[which_integral];
      zprime=volintegral_inside_sphere__center_z[which_integral];
    }
#pragma omp parallel for collapse(3)
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  double x_minus_xprime = x[index] - xprime;
	  double y_minus_xprime = y[index] - yprime;
	  double z_minus_xprime = z[index] - zprime;
	  if(sqrt(x_minus_xprime*x_minus_xprime + y_minus_xprime*y_minus_xprime + z_minus_xprime*z_minus_xprime)>radius) 
	    VolIntegrand1[index]=VolIntegrand2[index]=VolIntegrand3[index]=VolIntegrand4[index]=0.0;
	}
  }

  /* Set integrands to zero inside a sphere centered at x,y,z. 
     I.e., this results in the integral being restricted OUTSIDE sphere.
     Combine this with above to get spherical shell. 
  
     Note that volume integrals outside a sphere are fixed at the 
     original sphere center position for all time, unlike volume 
     integrals inside a sphere. We do this because the latter are 
     generally used for tracking moving compact objects or other things. */
  if(volintegral_outside_sphere__radius[which_integral]>0.0) {
    double radius=volintegral_outside_sphere__radius[which_integral];
    double xprime=volintegral_outside_sphere__center_x[which_integral];
    double yprime=volintegral_outside_sphere__center_y[which_integral];
    double zprime=volintegral_outside_sphere__center_z[which_integral];
#pragma omp parallel for collapse(3)
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  double x_minus_xprime = x[index] - xprime;
	  double y_minus_xprime = y[index] - yprime;
	  double z_minus_xprime = z[index] - zprime;
	  if(sqrt(x_minus_xprime*x_minus_xprime + y_minus_xprime*y_minus_xprime + z_minus_xprime*z_minus_xprime)<=radius) 
	    VolIntegrand1[index]=VolIntegrand2[index]=VolIntegrand3[index]=VolIntegrand4[index]=0.0;
	}
  }
}
