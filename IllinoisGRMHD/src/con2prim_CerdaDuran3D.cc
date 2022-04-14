#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "IllinoisGRMHD_headers.h"
#include "con2prim_headers.h"
#include "con2prim_helpers.h"
#include "inlined_functions.h"

/*****************************************************************************/
/*********** IMPORTED FROM A STRIPPED VERSION OF THE HEADER FILES ************/
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>

void NR_3D_WZT( const igm_eos_parameters eos,
                const CCTK_REAL tol_x,
                const CCTK_REAL S_squared,
                const CCTK_REAL BdotS,
                const CCTK_REAL B_squared,
                const CCTK_REAL *restrict S_con,
                const CCTK_REAL *restrict con,
                CCTK_REAL *restrict prim,
                bool *restrict c2p_failed );

static inline CCTK_REAL compute_W_from_cons( const igm_eos_parameters eos,
					     const CCTK_REAL *restrict con,
					     const CCTK_REAL S_squared,
					     const CCTK_REAL B_squared,
					     const CCTK_REAL BdotS ) {

  // Use the same as the Palenzuela routine
  const CCTK_REAL q = con[TAU]/con[DD];
  const CCTK_REAL r = S_squared/(con[DD]*con[DD]);
  const CCTK_REAL s = B_squared/con[DD];
  const CCTK_REAL t = BdotS/(pow(con[DD],1.5));
  const CCTK_REAL x = 2.0+2.0*q-s;

  CCTK_REAL Wminus2 = 1.0 - ( x*x*r + (2*x+s)*t*t ) / ( x*x*(x+s)*(x+s) );
  Wminus2           = fmin(fmax(Wminus2,eos.inv_W_max_squared ), 1.0);
  const CCTK_REAL W = pow(Wminus2, -0.5);
  return W;

}

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
int con2prim_CerdaDuran3D( const igm_eos_parameters eos,
                           const CCTK_REAL *restrict adm_quantities,
                           const CCTK_REAL *restrict con,
                           CCTK_REAL *restrict prim,
                           output_stats& stats ) {

  // Set gamma_{ij} and gamma^{ij}
  CCTK_REAL gammaDD[3][3],gammaUU[3][3];
  set_gammaDD_and_gammaUU_from_ADM_quantities(adm_quantities,gammaDD,gammaUU);

  // Read in B^{i}
  CCTK_REAL BU[3];
  for(int i=0;i<3;i++) BU[i] = con[B1_con+i];

  // Read in S_{i}
  CCTK_REAL SD[3];
  for(int i=0;i<3;i++) SD[i] = con[S1_cov+i];

  // Compute B_{i} = gamma_{ij}B^{j}
  CCTK_REAL BD[3];
  raise_or_lower_indices_3d( BU,gammaDD, BD );

  // Compute S^{i} = gamma^{ij}S_{j}
  CCTK_REAL SU[3];
  raise_or_lower_indices_3d( SD,gammaUU, SU );

  // S^2 = S^i * S_i
  CCTK_REAL S_squared = 0.0;
  for(int i=0;i<3;i++) S_squared += SU[i] * SD[i];

  // Enforce ceiling on S^{2} (A5 of Palenzuela et al. https://arxiv.org/pdf/1505.01607.pdf)
  CCTK_REAL S_squared_max = SQR( con[DD] + con[TAU] );
  if( S_squared > 0.9999 * S_squared_max ) {
    // Compute rescaling factor
    CCTK_REAL rescale_factor_must_be_less_than_one = sqrt(0.9999*S_squared_max/S_squared);
    // Rescale S_{i}
    for(int i=0;i<3;i++) SD[i] *= rescale_factor_must_be_less_than_one;
    // S_{i} has been rescaled. Recompute S^{i}.
    raise_or_lower_indices_3d( SD,gammaUU, SU );
    // Now recompute S^{2} := gamma^{ij}S^{i}S_{j}.
    S_squared = 0.0;
    for(int i=0;i<3;i++) S_squared += SU[i] * SD[i];
    // Check if the fix was successful
    if( simple_rel_err(S_squared,0.9999*S_squared_max) > 1e-12 ) CCTK_VError(VERR_DEF_PARAMS,"Incompatible values of S_squared after rescaling: %.15e %.15e\n",S_squared,0.9999*S_squared_max);
  }

  // Need to calculate for (21) and (22) in Cerda-Duran 2008
  // B * S = B^i * S_i
  CCTK_REAL BdotS = 0.0;
  for(int i=0;i<3;i++) BdotS += BU[i] * SD[i];

  // B^2 = B^i * B_i
  CCTK_REAL B_squared = 0.0;
  for(int i=0;i<3;i++) B_squared += BU[i] * BD[i];

  bool c2p_failed  = false;
  CCTK_REAL tol_x    = 5e-9;
  NR_3D_WZT( eos, tol_x, S_squared,BdotS,B_squared, SU,con,prim, &c2p_failed );

  return c2p_failed;

}

/*****************************************************************************/
/*********************** CERDA-DURAN CON2PRIM FUNCTIONS **********************/
/*****************************************************************************/
void calc_WZT_guesses( const igm_eos_parameters eos,
		       const CCTK_REAL S_squared,
		       const CCTK_REAL B_squared,
		       const CCTK_REAL BdotS,
		       const CCTK_REAL *restrict con,
		       CCTK_REAL *restrict xmax ) {

  // First compute W following Palenzuela
  const CCTK_REAL W = compute_W_from_cons(eos, con, S_squared, B_squared, BdotS);

  // Compute P and eps
  const CCTK_REAL xrho  = con[DD]/W;
  const CCTK_REAL xye   = con[YE]/con[DD];
  const CCTK_REAL xtemp = eos.T_max; // initial guess, choose large enough
  CCTK_REAL xprs        = 0.0;
  CCTK_REAL xeps        = 0.0;
  WVU_EOS_P_and_eps_from_rho_Ye_T( xrho, xye, xtemp, &xprs, &xeps );

  // Compute z = rho * h * W^2
  const CCTK_REAL h = 1.0 + xeps + xprs/xrho;
  const CCTK_REAL z = xrho * h * W * W;

  // Now set W, z, and T
  xmax[0] = W;
  xmax[1] = z;
  xmax[2] = xtemp;

}


void calc_prim_from_x_3D_WZT( const igm_eos_parameters eos,
                              const CCTK_REAL BdotS,
                              const CCTK_REAL B_squared,
                              const CCTK_REAL *restrict S_con,
                              const CCTK_REAL *restrict con,
                              CCTK_REAL *restrict prim,
                              CCTK_REAL *restrict x ) {

  // Recover the primitive variables from the scalars (W,Z)
  // and conserved variables, Eq. (23)-(25) in Cerdá-Durán et al. 2008
  CCTK_REAL W = x[0];
  CCTK_REAL z = x[1];
  CCTK_REAL T = x[2];

  // Calculate press, eps etc. from (rho, temp, Ye) using EOS,
  // required for consistency
  CCTK_REAL xrho  = con[DD]/W;
  CCTK_REAL xye   = con[YE]/con[DD];
  CCTK_REAL xtemp = T;
  CCTK_REAL xeps  = 0.0;
  CCTK_REAL xprs  = 0.0;
  WVU_EOS_P_and_eps_from_rho_Ye_T( xrho,xye,xtemp, &xprs,&xeps );

  // Eq. (24) in Siegel et al. 2018, with S^{i} := gamma^{ij} S_{j}
  // The extra factor of W converts v^{i} to tilde(u)^{i}.
  prim[UTCON1  ] = W*(S_con[0] + (BdotS)*con[B1_con]/z)/(z+B_squared);
  prim[UTCON2  ] = W*(S_con[1] + (BdotS)*con[B2_con]/z)/(z+B_squared);
  prim[UTCON3  ] = W*(S_con[2] + (BdotS)*con[B3_con]/z)/(z+B_squared);
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

void NR_step_3D_eps( const igm_eos_parameters eos,
                     const CCTK_REAL S_squared,
                     const CCTK_REAL BdotS,
                     const CCTK_REAL B_squared,
                     const CCTK_REAL *restrict con,
                     CCTK_REAL *restrict x,
                     CCTK_REAL *restrict dx,
                     CCTK_REAL *restrict f ) {
  // Finding the roots of f(x):
  //
  // x_{n+1} = x_{n} - f(x)/J = x_{n} + dx_{n}
  //
  // where J is the Jacobian matrix J_{ij} = df_i/dx_j
  //
  // Here, compute dx = [dW, dT]
  CCTK_REAL W = x[0];
  CCTK_REAL z = x[1];
  CCTK_REAL T = x[2];

  // Need partial derivatives of specific internal energy and pressure wrt density and
  // temperature. Those need to be based on primitives computed from Newton-Raphson state
  // vector x and conservatives
  CCTK_REAL xrho      = con[DD]/W;
  CCTK_REAL xye       = con[YE]/con[DD];
  CCTK_REAL xprs      = 0.0;
  CCTK_REAL xeps      = 0.0;
  CCTK_REAL xdPdrho   = 0.0;
  CCTK_REAL xdPdT     = 0.0;
  CCTK_REAL xdepsdrho = 0.0;
  CCTK_REAL xdepsdT   = 0.0;
  WVU_EOS_P_eps_dPdrho_dPdT_depsdrho_and_depsdT_from_rho_Ye_T(xrho,xye,T,&xprs,&xeps,&xdPdrho,&xdPdT,&xdepsdrho,&xdepsdT);

  // Some useful auxiliary variables
  CCTK_REAL BdotSsqr       = BdotS*BdotS;
  CCTK_REAL z_plus_Bsq     = z + B_squared;
  CCTK_REAL z_plus_Bsq_sqr = z_plus_Bsq*z_plus_Bsq;
  CCTK_REAL z_sqr          = z*z;
  CCTK_REAL W_sqr          = W*W;
  CCTK_REAL inv_W          = 1.0/W;
  CCTK_REAL inv_W_sqr      = 1.0/W_sqr;
  CCTK_REAL inv_D          = 1.0/con[DD];

  // Now compute eps(W,z)
  CCTK_REAL eps = - 1.0 + (z*inv_W - xprs*W)*inv_D;
  eps = fmax(eps,eos.eps_min);

  //---------------------------------------------------------------
  //------------ Equation (29) in Siegel et al. (2017) ------------
  //---------------------------------------------------------------
  // f0 = ( (z+B^{2})^{2} - S^{2} - (2z+B^{2})(B.S)^{2}/z^{2} )W^{2} - (z+B^{2})^{2}
  CCTK_REAL f0 = (z_plus_Bsq_sqr - S_squared - (z+z_plus_Bsq)*BdotSsqr/z_sqr)*W_sqr - z_plus_Bsq_sqr;

  // df0/dW = 2W( (z+B^{2})^{2} - S^{2} - (2z+B^{2})(B.S)^{2}/z^{2} )
  CCTK_REAL a  = 2.0*(z_plus_Bsq_sqr - S_squared - (z+z_plus_Bsq)*(BdotSsqr)/z_sqr)*W;

  // df0/dz = ( 2(z+B^{2}) + 2(B.S)^{2}/z^{2} + 2(B^{2})(B.S)^{2}/z^{3} )W^{2} - 2(z+B^{2})
  CCTK_REAL b  = (2.0*z_plus_Bsq + (2.0/(z_sqr) + 2.0*B_squared/(z_sqr*z))*BdotSsqr)*W_sqr - 2.0*z_plus_Bsq;

  // df0/dT = 0
  CCTK_REAL c  = 0.0;
  //---------------------------------------------------------------

  //---------------------------------------------------------------
  //------------ Equation (30) in Siegel et al. (2017) ------------
  //---------------------------------------------------------------
  // f1 = (tau + D - z - B^{2} + (B.S)^{2}/(2z^{2}) + P)W^{2} + B^{2}/2
  CCTK_REAL f1 = (con[TAU] + con[DD] - z - B_squared + (BdotSsqr)/(2.0*z_sqr) + xprs)*W_sqr + 0.5*B_squared;

  // Recall that D = W rho => rho = D / W. Thus
  //
  // dP/dW = (dP/drho)(drho/dW) = -W^{2}D(dP/drho)
  //
  // df1/dW = 2(tau + D - z - B^{2} + (B.S)^{2}/(2z^{2}) + P)W + dPdW * W^{2}
  CCTK_REAL d  = 2.0*(con[TAU] + con[DD] - z - B_squared + BdotSsqr/(2.0*z_sqr) + xprs)*W - con[DD]*xdPdrho;

  // df1/dz = ( -1 - (B.S)^{2}/z^{3} ) W^{2}
  CCTK_REAL e  = (-1.0-BdotSsqr/(z_sqr*z))*W_sqr;

  // df1/dT = W^{2}dPdT
  CCTK_REAL fc = W_sqr * xdPdT;

  //---------------------------------------------------------------

  //---------------------------------------------------------------
  //------------ Equation (31) in Siegel et al. (2017) ------------
  //---------------------------------------------------------------
  // Note that we use the *pressure* instead of eps
  // f2 = P(W,z) - P(rho,Ye,T)
  CCTK_REAL f2 = eps - xeps;

  // df2/dW = deps(W,z)/dW - (deps(rho,Ye,T)/drho)*(drho/dW)
  //        = -z/(D W^{2}) - P/W + (1/W)(dP/drho) + (D/W^{2})(deps/drho)
  CCTK_REAL g  = (con[DD]*xdepsdrho - z*inv_D)*inv_W_sqr + (xdPdrho - xprs)*inv_W;

  // df2/dz = deps(W,z)/dz = 1/(D W)
  CCTK_REAL h  = inv_W*inv_D;

  // df2/dT = -(W/D)(dP/dT) - deps(rho,Ye,T)/dT
  CCTK_REAL k  = -W*inv_D*xdPdT - xdepsdT;

  // Set the vector of equations
  f[0] = f0;
  f[1] = f1;
  f[2] = f2;

  // Compute the determinant
  CCTK_REAL A    = e*k-fc*h;
  CCTK_REAL B    = fc*g-d*k;
  CCTK_REAL C    = d*h-e*g;
  CCTK_REAL detJ = a*(A) + b*(B) + c*(C);

  //Compute the matrix inverse
  CCTK_REAL Ji[3][3];
  Ji[0][0] = A/detJ;
  Ji[1][0] = B/detJ;
  Ji[2][0] = C/detJ;
  Ji[0][1] = (c*h-b*k)/detJ;
  Ji[1][1] = (a*k-c*g)/detJ;
  Ji[2][1] = (g*b-a*h)/detJ;
  Ji[0][2] = (b*fc-c*e)/detJ;
  Ji[1][2] = (c*d-a*fc)/detJ;
  Ji[2][2] = (a*e-b*d)/detJ;

  // Compute the step size
  dx[0] = -Ji[0][0]* f[0] - Ji[0][1]* f[1] - Ji[0][2]* f[2];
  dx[1] = -Ji[1][0]* f[0] - Ji[1][1]* f[1] - Ji[1][2]* f[2];
  dx[2] = -Ji[2][0]* f[0] - Ji[2][1]* f[1] - Ji[2][2]* f[2];

}

void NR_3D_WZT( const igm_eos_parameters eos,
                const CCTK_REAL tol_x,
                const CCTK_REAL S_squared,
                const CCTK_REAL BdotS,
                const CCTK_REAL B_squared,
                const CCTK_REAL *restrict S_con,
                const CCTK_REAL *restrict con,
                CCTK_REAL *restrict prim,
                bool *restrict c2p_failed ) {

  // 2D Newton-Raphson scheme, using state vector x = (W, T) and 2D function
  // f(x) = (f1(x), f2(x)) given by Eqs. (27), (28) of Siegel et al. 2018

  // initialize Newton-Raphson state vector x  = (W, T)
  CCTK_REAL x[3];
  for(int i=0;i<3;i++) {
    x[i] = 0.0;
  }

  int count = 0;      // count number of iterations
  CCTK_REAL f[3];        // root finding function f
  CCTK_REAL x_lowlim[3]; // lower limits on W, T
  CCTK_REAL dx[3];       // displacement vector
  CCTK_REAL x_old[3];    // old state vector
  CCTK_REAL error[3];    // error vector abs(dx/x)

  x_lowlim[0] = 1.0;
  x_lowlim[1] = eos.rho_min;
  x_lowlim[2] = eos.T_atm;//exp( ( log(eos.T_max)+log(eos.T_min) )/2.0 );

  // set initial guess for Newton-Raphson state vector x = (W, T)
  calc_WZT_guesses(eos,S_squared,B_squared,BdotS,con,x);

  // printf("Guesses: %e %e %e\n",x[0],x[1],x[2]);

  // initialize variables
  for(int i=0;i<3;i++) {
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
  CCTK_REAL maxerror;

  while(keep_iterating) {

    // do Newton-Raphson step
    NR_step_3D_eps(eos,S_squared, BdotS, B_squared, con, x, dx, f);

    // Update x vector and compute error
    for(int i=0;i<3;i++) {

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

    // printf("Iter %d: %e %e %e | error: %e %e %e\n",count,x[0],x[1],x[2],error[0],error[1],error[2]);
    // getchar();

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
  if( (!robust_isfinite(f[0])) || (!robust_isfinite(f[1])) ) {
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
  calc_prim_from_x_3D_WZT( eos,BdotS,B_squared,S_con, con,prim,x );
}
