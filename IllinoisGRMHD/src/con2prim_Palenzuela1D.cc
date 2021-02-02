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

void palenzuela( const igm_eos_parameters eos,
                 bool *restrict c2p_failed, const double S_squared, const double BdotS,
                 const double B_squared, const double *restrict con, double *restrict prim,
                 const double *restrict SU, const double tol_x, bool use_epsmin );

double zbrent( double (*func)(const igm_eos_parameters, double, double *restrict, bool, double *restrict),
               const igm_eos_parameters eos, double *restrict param,
               double *restrict temp_guess, double x1, double x2, double tol_x, bool *restrict c2p_failed, bool use_epsmin );

void calc_prim( const igm_eos_parameters eos,
                const double x, const double *restrict con, const double * param, const double temp_guess, double *restrict prim,
                const double S_squared, const double BdotS,
                const double B_squared, const double *restrict SU );

double func_root( const igm_eos_parameters eos, double x, double *restrict param, bool use_epsmin, double *restrict temp_guess );


/*****************************************************************************/
/************** ILLINOISGRMHD INTERFACE FOR PALENZUELA CON2PRIM **************/
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
int con2prim_Palenzuela1D( const igm_eos_parameters eos,
                           const CCTK_REAL *restrict adm_quantities,
                           const CCTK_REAL *restrict con,
                           CCTK_REAL *restrict prim ) {

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
  palenzuela( eos, &c2p_failed, S_squared,BdotS,B_squared, con, prim, SU, 1e-12, true );

  return c2p_failed;

}

/*****************************************************************************/
/*********************** PALENZUELA CON2PRIM FUNCTIONS ***********************/
/*****************************************************************************/
void palenzuela(const igm_eos_parameters eos,
                bool *restrict c2p_failed, const double S_squared, const double BdotS,
                const double B_squared, const double * con, double *restrict prim,
                const CCTK_REAL *restrict SU, const double tol_x, bool use_epsmin)
{

  // main root finding routine
  // using scheme of Palenzuela et al. 2015, PRD, 92, 044045

  // some quantities computed from the conservatives,
  // defined as order unity quantities
  double q = con[TAU]/con[DD];
  double r = S_squared/(con[DD]*con[DD]);
  double s = B_squared/con[DD];
  double t = BdotS/(pow(con[DD],1.5));

  double param[6];
  param[par_q] = q;
  param[par_r] = r;
  param[par_s] = s;
  param[par_t] = t;
  param[conDD] = con[DD];
  param[conYE] = con[YE];

  // bracket for x
  double xlow = 1.0+q-s;
  double xup  = 2.0+2.0*q-s;

  // initial guess for temperature
  double temp_guess=prim[TEMP];

  // find x, this is the recovery process
  double x = zbrent(*func_root, eos, param, &temp_guess, xlow, xup, tol_x, c2p_failed, use_epsmin);

  // calculate final set of primitives
  calc_prim(eos, x, con, param, temp_guess, prim, S_squared, BdotS, B_squared, SU);

}

void calc_prim(const igm_eos_parameters eos,
               const double x, const double *restrict con, const double *restrict param, const double temp_guess, double *restrict prim,
               const double S_squared, const double BdotS,
               const double B_squared, const double *restrict SU ) {
  
  // Recover the primitive variables prim from x, q, r, s, t, con

  double q = param[par_q];
  double r = param[par_r];
  double s = param[par_s];
  double t = param[par_t];

  double Wminus2 = 1.0 - ( x*x*r + (2*x+s)*t*t ) / ( x*x*(x+s)*(x+s) );
  Wminus2 = fmin(fmax(Wminus2,eos.inv_W_max_squared ), 1.0);
  double W= pow(Wminus2, -0.5);

  double rho   = con[DD]/W;
  double ye    = con[YE]/con[DD];
  double temp  = temp_guess;
  double press = 0.0;
  double eps   = 0.0;
  double ent   = 0.0;
  if( eos.evolve_T ) {
    eps = W - 1.0 + (1.0-W*W)*x/W + W*(q - s + t*t/(2*x*x) + s/(2*W*W)  );
    eps=fmax(eps, eos.eps_min);
    get_P_S_and_T_from_rho_Ye_and_eps( eos, rho,ye,eps, &press,&ent,&temp );
  }
  else {
    get_P_eps_and_S_from_rho_Ye_and_T( eos, rho,ye,temp, &press,&eps,&ent );
  }

  const CCTK_REAL Z = x*rho*W;

  // Now we compute the velocity using eq. (24) in https://arxiv.org/pdf/1712.07538.pdf.
  // Note, however, that IllinoisGRMHD expects the velocity tilde(u)^{i} := W v^{i},
  // and wo we already perform the conversion here (see the discussion above eq. (26) in
  // https://arxiv.org/pdf/astro-ph/0512420.pdf for the definition of tilde(u)^{i}).
  prim[UTCON1  ] = W*(SU[0] + (BdotS)*con[B1_con]/Z)/(Z+B_squared);
  prim[UTCON2  ] = W*(SU[1] + (BdotS)*con[B2_con]/Z)/(Z+B_squared);
  prim[UTCON3  ] = W*(SU[2] + (BdotS)*con[B3_con]/Z)/(Z+B_squared);
  prim[RHO     ] = rho;
  prim[EPS     ] = eps;
  prim[B1_con  ] = con[B1_con];
  prim[B2_con  ] = con[B2_con];
  prim[B3_con  ] = con[B3_con];
  prim[TEMP    ] = temp;
  prim[YE      ] = ye;
  prim[PRESS   ] = press;
  prim[WLORENTZ] = W;
  prim[ENT     ] = ent;
}

double func_root(const igm_eos_parameters eos, double x, double *restrict param, bool use_epsmin, double *restrict temp_guess) {

  // computes f(x) from x and q,r,s,t

  const double q  = param[par_q];
  const double r  = param[par_r];
  const double s  = param[par_s];
  const double t  = param[par_t];

  // (i)
  double Wminus2 = 1.0 - ( x*x*r + (2*x+s)*t*t  )/ (x*x*(x+s)*(x+s));
  Wminus2 = fmin(fmax(Wminus2,eos.inv_W_max_squared ), 1.0);
  const double W= pow(Wminus2, -0.5);

  // (ii)
  double rho  = param[conDD]/W;
  double ye   = param[conYE]/param[conDD];
  double temp = *temp_guess;
  double P    = 0.0;
  double eps  = 0.0;
  if( eos.evolve_T ) {
    eps = W - 1.0 + (1.0-W*W)*x/W + W*(q - s + t*t/(2*x*x) + s/(2*W*W));
    eps = fmax(eps, eos.eps_min);
    get_P_and_T_from_rho_Ye_and_eps( eos,rho,ye,eps, &P, &temp );
  }
  else {
    get_P_and_eps_from_rho_Ye_and_T( eos,rho,ye,eps, &P, &temp );
  }

  // (iv)
  double ans = x- (1.0 + eps + P/rho)*W;

  return ans;
}

/****************************************************************************/
/************************* BRENT's METHOD FUNCTIONS *************************/
/****************************************************************************/

#define ITMAX 300   //Maximum allowed number of iterations.
#define PREC 3.0e-15
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


double zbrent(double (*func)(const igm_eos_parameters, double, double *, bool, double *restrict),
              const igm_eos_parameters eos, double *restrict param,
              double * temp_guess, double x1, double x2, double tol_x, bool *restrict c2p_failed, bool use_epsmin)
{

  // Implementation of Brent’s method to find the root of a function func to accuracy tol_x. The root
  // must be bracketed by x1 and x2.
  // Based on W. H. Press et al. 2007, Numerical Recipes in C: The Art of Scientific Computing, 3rd edition, 
  // Cambridge University Press, ISBN 978-0-521-88068-8

  int iter = 0;

  //termination critera definitions
  // int MAX_NEWT_ITER   = c2p.max_iterations;
  // int EXTRA_NEWT_ITER = c2p.extra_iterations;
  double ans = 0.0;

  double a=x1,b=x2,c=x2,d=0.0,e=0.0,min1=0.0,min2=0.0;
  double temp_guess_a = *temp_guess;
  double temp_guess_b = *temp_guess;
  double temp_guess_c = *temp_guess;
  double fa=(*func)(eos, a, param, use_epsmin, &temp_guess_a);
  double fb=(*func)(eos, b, param, use_epsmin, &temp_guess_b);
  double fc,p,q,r,s,tol1;

  // root must be bracketed
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    *c2p_failed = true;
    return b;
  }
  fc=fb;
  temp_guess_c = temp_guess_b;

  int i_extra = 0;
  int doing_extra = 0;

  if (EXTRA_NEWT_ITER==0) {
    i_extra=-1;
  }

  bool keep_iterating = 1;
  double maxerror = 0.0;

  while (keep_iterating) {
 
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      // reset c, adjust bracket interval d
      c=a;
      fc=fa;
      e=d=b-a;
      temp_guess_c = temp_guess_a;
    }
    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
      temp_guess_a = temp_guess_b;
      temp_guess_b = temp_guess_c;
      temp_guess_c = temp_guess_a;
    }

    tol1=2.0*PREC*fabs(b)+0.5*tol_x;
    maxerror=0.5*(c-b);
    // check convergence
    if (fabs(maxerror) <= tol1 || fb == 0.0){
      // we are done
      *c2p_failed = false;
      keep_iterating = 0;
      ans = b;
    }
    else if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      // inverse quadratic interpolation
      s=fb/fa;   
    
      if (a == c) {
        // fake quad interpolation
        p=2.0*maxerror*s;
        q=1.0-s;
      } else {
        // actual quad interpolation
        q=fa/fc;
        r=fb/fc;
        p=s*(2.0*maxerror*q*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
      }

      if (p > 0.0) q = -q;  // check whether in bounds
      p=fabs(p);
      min1=3.0*maxerror*q-fabs(tol1*q);
      min2=fabs(e*q);
      
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        // interpolation successful
        e=d;
        d=p/q;
      } else {
        // interpolation failed, use bisection
        d=maxerror;
        e=d;
      }

      // update a to be last best guess
      a=b;
      fa=fb;
      temp_guess_a = temp_guess_b;

      // update trial root
      if (fabs(d) > tol1)  b += d;
      else                 b += SIGN(tol1,maxerror);
      fb=(*func)(eos, b, param, use_epsmin, &temp_guess_b);


    } else {  

      // Bounds decreasing too slowly, use bisection

      d=maxerror;
      e=d;

      // update a to be last best guess
      a=b;
      fa=fb;
      temp_guess_a = temp_guess_b;

      // update trial root
      if (fabs(d) > tol1)  b += d;
      else                 b += SIGN(tol1, maxerror);
      fb=(*func)(eos, b, param, use_epsmin, &temp_guess_b);
    }

    ++iter;

    // termination criterion
    if( (fabs(maxerror) <= tol1 ) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(maxerror) <= tol1)&&(doing_extra == 0)) 
        || (i_extra > EXTRA_NEWT_ITER) || (iter >= (ITMAX)) ) {
      keep_iterating = 0;
    }

  } //end for loop

  if( (!CCTK_isfinite(fb)) ) {
    *c2p_failed = true;
    return b;
  } if( (fabs(maxerror) <= tol1 || fb == 0.0)){
    *c2p_failed = false;
    return ans;
  } else if( (fabs(maxerror) <= tol1) && (fabs(maxerror) > tol1) ){
    *c2p_failed = false;
    return ans;
  } else {
    *c2p_failed = true;
    return b;
  }

}
