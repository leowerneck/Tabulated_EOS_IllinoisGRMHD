#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "IllinoisGRMHD_headers.h"
#include "con2prim_headers.h"

/*****************************************************************************/
/*********** IMPORTED FROM A STRIPPED VERSION OF THE HEADER FILES ************/
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>

void NR_2D_WT( const igm_eos_parameters eos,
               const int safe_guess,
               const double tol_x,
               const double S_squared,
               const double BdotS,
               const double B_squared,
               const double *restrict S_con,
               const double *restrict con,
               double *restrict prim,
               bool *restrict c2p_failed );

/*****************************************************************************/
/************** ILLINOISGRMHD INTERFACE FOR CERDA-DURAN CON2PRIM *************/
/*****************************************************************************/
//
// This routine expects 8 conservs as input:
//
// -> D     = rho_star      / sqrt(gamma) = W * rho (W is the Lorentz factor)
// -> DYe   = Ye_Star       / sqrt(gamma)
// -> B^{i} = \tilde{B}^{i} / sqrt(gamma) / sqrt(4pi)
// -> S_{i} = \tilde{S}_{i} / sqrt(gamma)
//
// From the input quantities, we compute B_{i} and S^{i}
int con2prim_CerdaDuran2D( const igm_eos_parameters eos,
                           const CCTK_REAL *restrict adm_quantities,
                           const CCTK_REAL *restrict con,
                           CCTK_REAL *restrict prim,
                           output_stats& stats ) {

  // Set gamma_{ij}
  CCTK_REAL gammaDD[3][3];
  gammaDD[0][0] = adm_quantities[GXX];
  gammaDD[0][1] = gammaDD[1][0] = adm_quantities[GXY];
  gammaDD[0][2] = gammaDD[2][0] = adm_quantities[GXZ];
  gammaDD[1][1] = adm_quantities[GYY];
  gammaDD[1][2] = gammaDD[2][1] = adm_quantities[GYZ];
  gammaDD[2][2] = adm_quantities[GZZ];

  // Set gamma^{ij}
  CCTK_REAL gammaUU[3][3];
  gammaUU[0][0] = adm_quantities[GUPXX];
  gammaUU[0][1] = gammaUU[1][0] = adm_quantities[GUPXY];
  gammaUU[0][2] = gammaUU[2][0] = adm_quantities[GUPXZ];
  gammaUU[1][1] = adm_quantities[GUPYY];
  gammaUU[1][2] = gammaUU[2][1] = adm_quantities[GUPYZ];
  gammaUU[2][2] = adm_quantities[GUPZZ];

  // Read in B^{i}
  CCTK_REAL BU[3];
  for(int i=0;i<3;i++) BU[i] = con[B1_con+i];

  // Compute B_{i} = gamma_{ij}B^{j}
  CCTK_REAL BD[3];
  for(int i=0;i<3;i++) {
    BD[i] = 0;
    for(int j=0;j<3;j++) {
      BD[i] += gammaDD[i][j] * BU[j];
    }
  }

  // Read in S_{i}.
  CCTK_REAL SD[3];
  for(int i=0;i<3;i++) SD[i] = con[S1_cov+i];

  // Compute S^{i} = gamma^{ij}S_{j}
  CCTK_REAL SU[3];
  for(int i=0;i<3;i++) {
    SU[i] = 0;
    for(int j=0;j<3;j++) {
      SU[i] += gammaUU[i][j] * SD[j];
    }
  }

  // Need to calculate for (21) and (22) in Cerda-Duran 2008
  // B * S = B^i * S_i
  CCTK_REAL BdotS = 0.0;
  for(int i=0;i<3;i++) BdotS += BU[i] * SD[i];
  
  // B^2 = B^i * B_i
  CCTK_REAL B_squared = 0.0;
  for(int i=0;i<3;i++) B_squared += BU[i] * BD[i];
  
  // S^2 = S^i * S_i
  CCTK_REAL S_squared = 0.0;
  for(int i=0;i<3;i++) S_squared += SU[i] * SD[i];

  bool c2p_failed = false;
  int safe_guess  = 0;
  double tol_x    = 5e-9;
  NR_2D_WT( eos, safe_guess,tol_x, S_squared,BdotS,B_squared, SU,con,prim, &c2p_failed );

  if( c2p_failed ) {
    // If failed to recover the prims, try again with safe guesses
    int safe_guess=1;
    NR_2D_WT( eos, safe_guess,tol_x, S_squared,BdotS,B_squared, SU,con,prim, &c2p_failed );
  }

  return c2p_failed;

}

/*****************************************************************************/
/*********************** CERDA-DURAN CON2PRIM FUNCTIONS **********************/
/*****************************************************************************/
void calc_WT_max( const igm_eos_parameters eos,
                  const double B_squared,
                  const double *restrict con,
                  double *restrict xmax ) {

  // Calculate maximum values for x = (rho, T) ("safe guess" initial values)
  // cf. Cerda-Duran et al. 2008, Eq. (39)-(42)

  double rhomax = con[DD];
  double epsmax = (con[TAU] - B_squared/2.0) / con[DD];

  // ensure that rhomax and epsmax are in validity range of EOS
  // Note that in setting the IllinoisGRMHD EOS parameters, we
  // already impose a safety factor on the table floors and
  // ceilings, so we don't need to use the prefactors of 95%
  if(rhomax > eos.rho_max) rhomax = eos.rho_max;
  if(epsmax > eos.eps_max) epsmax = eos.eps_max;

  // Now compute P max and T max
  double xye   = con[YE]/con[DD];
  double xtemp = eos.T_max; // initial guess, choose large enough
  double xprs  = 0.0;
  WVU_EOS_P_and_T_from_rho_Ye_eps( rhomax,xye,epsmax, &xprs,&xtemp );

  // Now set W_max and T_max
  xmax[0] = 1.0e4;
  xmax[1] = xtemp;

}


void calc_prim_from_x_2D_WT( const igm_eos_parameters eos,
                             const double BdotS, 
                             const double B_squared,
                             const double *restrict S_con,
                             const double *restrict con,
                             double *restrict prim,
                             double *restrict x ) {
  
  // Recover the primitive variables from the scalars (W,Z) 
  // and conserved variables, Eq. (23)-(25) in Cerdá-Durán et al. 2008
  double W = x[0];
  double T = x[1];
  
  // Calculate press, eps etc. from (rho, temp, Ye) using EOS,
  // required for consistency
  double xrho  = con[DD]/W;
  double xye   = con[YE]/con[DD];
  double xtemp = T;
  double xeps  = 0.0;
  double xprs  = 0.0;  
  WVU_EOS_P_and_eps_from_rho_Ye_T( xrho,xye,xtemp, &xprs,&xeps );
  
  double Z = con[DD] * (1.0 + xeps + xprs/xrho) * W;

  // Eq. (24) in Siegel et al. 2018, with S^{i} := gamma^{ij} S_{j}
  // The extra factor of W converts v^{i} to tilde(u)^{i}.
  prim[UTCON1  ] = W*(S_con[0] + (BdotS)*con[B1_con]/Z)/(Z+B_squared);
  prim[UTCON2  ] = W*(S_con[1] + (BdotS)*con[B2_con]/Z)/(Z+B_squared);
  prim[UTCON3  ] = W*(S_con[2] + (BdotS)*con[B3_con]/Z)/(Z+B_squared);
  prim[RHO     ] = xrho;
  prim[YE      ] = xye;
  prim[TEMP    ] = T;
  prim[PRESS   ] = xprs;
  prim[EPS     ] = xeps;
  prim[B1_con  ] = con[B1_con];
  prim[B2_con  ] = con[B2_con];
  prim[B3_con  ] = con[B3_con];
  prim[WLORENTZ] = W;
  
}

void NR_step_2D_WT( const igm_eos_parameters eos,
                    const double S_squared,
                    const double BdotS,
                    const double B_squared,
                    const double *restrict con, 
                    double *restrict x,
                    double *restrict dx,
                    double *restrict f ) {
  // Finding the roots of f(x):
  //
  // x_{n+1} = x_{n} - f(x)/J = x_{n} + dx_{n}
  // 
  // where J is the Jacobian matrix J_{ij} = df_i/dx_j
  //
  // Here, compute dx = [dW, dT]
  double W = x[0];
  double T = x[1];

  // Need partial derivatives of specific internal energy and pressure wrt density and 
  // temperature. Those need to be based on primitives computed from Newton-Raphson state
  // vector x and conservatives
  double rho      = con[DD]/W;
  double ye       = con[YE]/con[DD];
  double P        = 0.0;
  double eps      = 0.0;
  double dPdrho   = 0.0;
  double dPdT     = 0.0;
  double depsdrho = 0.0;
  double depsdT   = 0.0;
  WVU_EOS_P_eps_dPdrho_dPdT_depsdrho_and_depsdT_from_rho_Ye_T(rho,ye,T,&P,&eps,&dPdrho,&dPdT,&depsdrho,&depsdT);

  // h = 1 + eps + P/rho
  double h     = 1.0 + eps + P / rho;

  // z = rho*h*W^{2} = ( rho + eps*rho + P )*W^{2}
  double W_sqr = W*W;
  double z     = rho * h * W_sqr;

  // Some useful auxiliary variables
  double BdotSsqr       = BdotS*BdotS;
  double z_plus_Bsq     = z + B_squared;
  double z_plus_Bsq_sqr = z_plus_Bsq*z_plus_Bsq;
  double z_sqr          = z*z;

  // dz/dW = 2 * rho * h * W = 2 * D * h
  double dz_dW = con[DD]*h;

  // dz/dT = ( rho*deps/dT + dP/dT )*W^{2}
  double dz_dT = (rho * depsdT + dPdT)*W_sqr;

  // f1 = ( (z+B^{2})^{2} - S^{2} - (2z+B^{2})(B.S)^{2}/z^{2} )W^{2} - (z+B^{2})^{2}
  f[0] = (z_plus_Bsq_sqr - S_squared - (z+z_plus_Bsq)*BdotSsqr/z_sqr)*W_sqr - z_plus_Bsq_sqr;

  // df1/dz = ( 2(z+B^{2}) + 2(B.S)^{2}/z^{2} + 2(B^{2})(B.S)^{2}/z^{3} )W^{2} - 2(z+B^{2})
  double df1_dz = (2.0*z_plus_Bsq + (2.0/(z_sqr) + 2.0*B_squared/(z_sqr*z))*BdotSsqr)*W_sqr - 2.0*z_plus_Bsq;

  // Ignoring the fact that z = z(W) for now, we have:
  // df1/dW = 2W( (z+B^{2})^{2} - S^{2} - (2z+B^{2})(B.S)^{2}/z^{2} )
  double df1_dW = 2.0*(z_plus_Bsq_sqr - S_squared - (z+z_plus_Bsq)*(BdotSsqr)/z_sqr)*W;

  // f2 = (tau + D - z - B^{2} + (B.S)^{2}/(2z^{2}) + P)W^{2} + B^{2}/2
  f[1] = (con[TAU] + con[DD] - z - B_squared + (BdotSsqr)/(2.0*z_sqr) + P)*W_sqr + 0.5*B_squared;

  // df2/dz = ( -1 - (B.S)^{2}/z^{3} ) W^{2}
  double df2_dz = (-1.0-BdotSsqr/(z_sqr*z))*W_sqr;

  // df2/dW = 2(tau + D - z - B^{2} + (B.S)^{2}/(2z^{2}) + P)W
  double df2_dW = 2.0*(con[TAU] + con[DD] - z - B_squared + BdotSsqr/(2.0*z_sqr) + P)*W;

  // df2/dP = W^{2}
  double df2_dP = W_sqr;

  // Now compute the Jacobian matrix
  double J[2][2];
  // df1_dW
  J[0][0] = df1_dz * dz_dW + df1_dW;
  double a = J[0][0];

  // df1_dT
  J[0][1] = df1_dz * dz_dT;
  double b = J[0][1];

  // df2_dW
  J[1][0] = df2_dz * dz_dW + df2_dW;
  double c = J[1][0];

  // df2_dT
  J[1][1] = df2_dz * dz_dT + df2_dP * dPdT;
  double d = J[1][1];

  // Then the inverse Jacobian Matrix
  double Ji[2][2];
  double detJ = a*d - b*c; 
  Ji[0][0]    = +d/detJ;
  Ji[0][1]    = -b/detJ;
  Ji[1][0]    = -c/detJ;
  Ji[1][1]    = +a/detJ;

  // Compute the step size
  dx[0] = -Ji[0][0]* f[0] - Ji[0][1]* f[1];
  dx[1] = -Ji[1][0]* f[0] - Ji[1][1]* f[1];

}


void NR_2D_WT( const igm_eos_parameters eos,
               const int safe_guess,
               const double tol_x,
               const double S_squared,
               const double BdotS,
               const double B_squared,
               const double *restrict S_con,
               const double *restrict con,
               double *restrict prim,
               bool *restrict c2p_failed ) {

  // 2D Newton-Raphson scheme, using state vector x = (W, T) and 2D function
  // f(x) = (f1(x), f2(x)) given by Eqs. (27), (28) of Siegel et al. 2018

  // initialize Newton-Raphson state vector x  = (W, T)
  double x[2];
  for(int i=0;i<2;i++) {
    x[i] = 0.0;
  }

  int count = 0;      // count number of iterations
  double f[2];        // root finding function f
  double x_lowlim[2]; // lower limits on W, T
  double dx[2];       // displacement vector
  double x_old[2];    // old state vector
  double error[2];    // error vector abs(dx/x)

  x_lowlim[0] = 1.0;
  x_lowlim[1] = eos.T_min;

  // set initial guess for Newton-Raphson state vector x = (W, T)
  if( safe_guess==0 ) {
    x[0] = prim[WLORENTZ];
    x[1] = prim[TEMP    ];
  }
  else if( safe_guess==1 ) {
    calc_WT_max(eos,B_squared,con,x);
  }

  // initialize variables
  for(int i=0;i<2;i++) {
    x_old[i] = 0.0;
    dx[i]    = 0.0;
    f[i]     = 0.0;
    error[i] = 0.0;
    //check for NaNs
    if( x[i] != x[i] ) {
      *c2p_failed = true;
      return;
    }
  }

  int i_extra     = 0;
  int doing_extra = 0;

  if (EXTRA_NEWT_ITER==0) {
    i_extra=-1;
  }

  bool keep_iterating=true;
  double maxerror;

  while(keep_iterating) {

    // do Newton-Raphson step
    NR_step_2D_WT(eos,S_squared, BdotS, B_squared, con, x, dx, f);

    // Update x vector and compute error
    for(int i=0;i<2;i++) {

      x_old[i] = x[i];
      // update x and exclude unphysical regime
      x[i]     = fmax(x[i]+ dx[i], x_lowlim[i]);
      error[i] = fabs((x[i]-x_old[i])/x[i]);

      // Check for NaNs
      if( x[i] != x[i] ) {
        *c2p_failed = true;
        return;
      }

    }
    maxerror = MAX(error[0], error[1]);
    count++;

    // termination criterion
    if( (fabs(maxerror) <= tol_x) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++;

    if( ((fabs(maxerror) <= tol_x)&&(doing_extra == 0))
        || (i_extra >= EXTRA_NEWT_ITER) || (count >= (MAX_NEWT_ITER)) ) {
      keep_iterating = false;
    }

  } // END of while(keep_iterating)

  //  Check for bad untrapped divergences
  if( (!CCTK_isfinite(f[0])) || (!CCTK_isfinite(f[1])) ) {
    *c2p_failed = true;
  }

  if( fabs(maxerror) <= tol_x ){
    *c2p_failed = false;
  }
  else if( (fabs(maxerror) <= tol_x) && (fabs(maxerror) > tol_x) ){
    *c2p_failed = false;
  }
  else {
    *c2p_failed = true;
  }

  // Recover the primitive variables from the final scalars x = (W, T)
  calc_prim_from_x_2D_WT( eos,BdotS,B_squared,S_con, con,prim,x );
}
