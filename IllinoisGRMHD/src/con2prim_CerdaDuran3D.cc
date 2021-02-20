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

void NR_3D_WZT( const igm_eos_parameters eos,
                const int safe_guess,
                const int use_eps_or_P,
                const CCTK_REAL tol_x,
                const CCTK_REAL S_squared,
                const CCTK_REAL BdotS,
                const CCTK_REAL B_squared,
                const CCTK_REAL *restrict S_con,
                const CCTK_REAL *restrict con,
                CCTK_REAL *restrict prim,
                bool *restrict c2p_failed );

static const int use_eps = 0;
static const int use_P   = 1;

static inline void compute_W_from_cons( const CCTK_REAL *restrict con,
                                        CCTK_REAL *restrict W ) {

  // Use the same as the Palenzuela routine
  
  
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

  bool c2p_failed  = false;
  int safe_guess   = 0;
  int use_eps_or_P = use_eps;
  CCTK_REAL tol_x    = 5e-9;
  NR_3D_WZT( eos, safe_guess,use_eps_or_P,tol_x, S_squared,BdotS,B_squared, SU,con,prim, &c2p_failed );

  if( c2p_failed ) {
    // If failed to recover the prims, try again with safe guesses
    int safe_guess=1;
    NR_3D_WZT( eos, safe_guess,use_eps_or_P,tol_x, S_squared,BdotS,B_squared, SU,con,prim, &c2p_failed );
  }

  return c2p_failed;

}

/*****************************************************************************/
/*********************** CERDA-DURAN CON2PRIM FUNCTIONS **********************/
/*****************************************************************************/
void calc_WZT_from_prim( const igm_eos_parameters eos,
                         const CCTK_REAL *restrict prim,
                         const CCTK_REAL *restrict con,
                         CCTK_REAL *restrict x ) {

  // Start by computing P and eps using the EOS
  CCTK_REAL W   = prim[WLORENTZ];
  CCTK_REAL rho = con[DD]/W;
  CCTK_REAL Ye  = con[YE]/con[DD];
  CCTK_REAL T   = prim[TEMP];
  CCTK_REAL P   = 0.0;
  CCTK_REAL eps = 0.0;
  WVU_EOS_P_and_eps_from_rho_Ye_T( rho,Ye,T, &P,&eps );

  // Then compute the enthalpy
  CCTK_REAL h = 1.0 + eps + P / rho;

  // Then compute z
  CCTK_REAL z = rho * h * W * W;

  // Then set the x vector
  x[0] = W;
  x[1] = z;
  x[2] = T;

}

void calc_WZT_max( const igm_eos_parameters eos,
                   const CCTK_REAL B_squared,
                   const CCTK_REAL *restrict con,
                   CCTK_REAL *restrict xmax ) {

  // Calculate maximum values for x = (rho, T) ("safe guess" initial values)
  // cf. Cerda-Duran et al. 2008, Eq. (39)-(42)

  CCTK_REAL rhomax = con[DD];
  CCTK_REAL epsmax = (con[TAU] - 0.5*B_squared) / con[DD];

  // ensure that rhomax and epsmax are in validity range of EOS
  // Note that in setting the IllinoisGRMHD EOS parameters, we
  // already impose a safety factor on the table floors and
  // ceilings, so we don't need to use the prefactors of 95%
  if(rhomax > eos.rho_max) rhomax = eos.rho_max;
  if(epsmax > eos.eps_max) epsmax = eos.eps_max;

  // Now compute P max and T max
  CCTK_REAL xye   = con[YE]/con[DD];
  CCTK_REAL xtemp = eos.T_max; // initial guess, choose large enough
  CCTK_REAL xprs  = 0.0;
  WVU_EOS_P_and_T_from_rho_Ye_eps( rhomax,xye,epsmax, &xprs,&xtemp );

  // Now set W_max and T_max
  xmax[0] = 1.0e4;
  xmax[1] = con[TAU] + xprs + con[DD] - 0.5*B_squared;
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

void NR_step_3D_P( const igm_eos_parameters eos,
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

  // Now compute P(W,z)
  CCTK_REAL P = z*inv_W_sqr - con[DD] * ( 1.0 + xeps) * inv_W ;

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
  CCTK_REAL f2 = P - xprs;

  // df2/dW = dP(W,z)/dW - (dP(rho,Ye,T)/drho)*(drho/dW)
  //        = D( 1.0 + eps + (D/W)(deps/drho) )/W^{2} -2z/W^{3} + (D/W^{2})(dP/drho)
  CCTK_REAL g  = (con[DD]*(1.0 + xeps + con[DD]*xdepsdrho/W) - 2.0*z*inv_W + con[DD]*xdPdrho)*inv_W_sqr;

  // df2/dz = dP(W,z)/dz = 1/W^{2}
  CCTK_REAL h  = inv_W_sqr;

  // df2/dT = -(D/W)(deps/dT) - dP(rho,Ye,T)/dT
  CCTK_REAL k  = -con[DD]*inv_W*xdepsdT - xdPdT;

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
                const int safe_guess,
                const int use_eps_or_P,
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
  if( safe_guess==0 ) {
    calc_WZT_from_prim(eos,prim,con,x);
  }
  else if( safe_guess==1 ) {
    calc_WZT_max(eos,B_squared,con,x);
  }

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
    if( use_eps_or_P == use_eps ) {
      NR_step_3D_eps(eos,S_squared, BdotS, B_squared, con, x, dx, f);
    }
    else if( use_eps_or_P == use_P ) {
      NR_step_3D_P(eos,S_squared, BdotS, B_squared, con, x, dx, f);
    }
    else {
      CCTK_ERROR("(CerdaDuran3D) use_eps_or_P can only be \"use_eps\" or \"use_P\". ABORTING.");
    }

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
  calc_prim_from_x_3D_WZT( eos,BdotS,B_squared,S_con, con,prim,x );
}
