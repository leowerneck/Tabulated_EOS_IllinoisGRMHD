#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "IllinoisGRMHD_headers.h"
#include "harm_primitives_headers.h"

/*****************************************************************************/
/*********** IMPORTED FROM A STRIPPED VERSION OF THE HEADER FILES ************/
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>

#define par_q (0)
#define par_r (1)
#define par_s (2)
#define par_t (3)
#define conDD (4)
#define conYE (5)
#define conDS (6) // Entropy

struct c2p_report {
  bool failed;
  bool adjust_cons;
  char err_msg[200];
  int count;
  bool retry;
  int c2p_keyerr;
  int nEOScalls;
};

void palenzuela_entropy(const igm_eos_parameters eos,
                        struct c2p_report * c2p_rep, const double S_squared, const double BdotS, 
                        const double B_squared, const double * con, double * prim, const double g_con[4][4], 
                        const double g_cov[4][4], const double tol_x, bool use_epsmin);

double zbrent_entropy(double (*func)(const igm_eos_parameters, double, double *, struct c2p_report *, bool, double *),
                      const igm_eos_parameters eos, double * param, 
                      double * temp_guess, double x1, double x2, double tol_x, struct c2p_report * c2p_rep, bool use_epsmin);

void calc_prim_entropy(const igm_eos_parameters eos,
                       const double x, const double * con, const double * param, const double temp_guess, double * prim, 
                       const double S_squared, const double BdotS, 
                       const double B_squared, const double g4up[4][4], struct c2p_report * c2p_rep );

double func_root_entropy(const igm_eos_parameters eos, double x, double * param, struct c2p_report * c2p_rep, bool use_epsmin, double * temp_guess);


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
int con2prim_palenzuela_entropy( const igm_eos_parameters eos,
                                 const CCTK_REAL g4dn[4][4],
                                 const CCTK_REAL g4up[4][4],
                                 CCTK_REAL *con,
                                 CCTK_REAL *prim ) {

  // Compute scalars needed by the con2prim routine      
  // Lower indices - covariant
  CCTK_REAL B1_cov = g4dn[1][1]*con[B1_con] + g4dn[1][2]*con[B2_con] + g4dn[1][3]*con[B3_con];
  CCTK_REAL B2_cov = g4dn[2][1]*con[B1_con] + g4dn[2][2]*con[B2_con] + g4dn[2][3]*con[B3_con];
  CCTK_REAL B3_cov = g4dn[3][1]*con[B1_con] + g4dn[3][2]*con[B2_con] + g4dn[3][3]*con[B3_con];
  
  // Raise indices - contravariant
  CCTK_REAL S1_con = g4up[1][1]*con[S1_cov] + g4up[1][2]*con[S2_cov] + g4up[1][3]*con[S3_cov];
  CCTK_REAL S2_con = g4up[2][1]*con[S1_cov] + g4up[2][2]*con[S2_cov] + g4up[2][3]*con[S3_cov];
  CCTK_REAL S3_con = g4up[3][1]*con[S1_cov] + g4up[3][2]*con[S2_cov] + g4up[3][3]*con[S3_cov];
 
  // Need to calculate for (21) and (22) in Cerda-Duran 2008
  // B * S = B^i * S_i
  CCTK_REAL BdotS = con[B1_con]*con[S1_cov] + con[B2_con]*con[S2_cov] + con[B3_con]*con[S3_cov];
  
  // B^2 = B^i * B_i
  CCTK_REAL B_squared = con[B1_con]*B1_cov + con[B2_con]*B2_cov + con[B3_con]*B3_cov;
  
  // S^2 = S^i * S_i
  CCTK_REAL S_squared = S1_con*con[S1_cov] + S2_con*con[S2_cov] + S3_con*con[S3_cov];

  struct c2p_report report;
  palenzuela_entropy( eos, &report, S_squared,BdotS,B_squared, con, prim, g4up, g4dn, 1e-12, true );

  return report.failed;

}

/*****************************************************************************/
/*********************** PALENZUELA CON2PRIM FUNCTIONS ***********************/
/*****************************************************************************/
void palenzuela_entropy(const igm_eos_parameters eos,
                        struct c2p_report * c2p_rep, const double S_squared, const double BdotS, 
                        const double B_squared, const double * con, double * prim, const double g_con[4][4],
                        const double g_cov[4][4], const double tol_x, bool use_epsmin)
{

  // main root finding routine
  // using scheme of Palenzuela et al. 2015, PRD, 92, 044045

  // some quantities computed from the conservatives,
  // defined as order unity quantities
  double q = con[TAU]/con[DD];
  double r = S_squared/(con[DD]*con[DD]);
  double s = B_squared/con[DD];
  double t = BdotS/(pow(con[DD],1.5));

  double param[7];
  param[par_q] = q;
  param[par_r] = r;
  param[par_s] = s;
  param[par_t] = t;
  param[conDD] = con[DD];
  param[conYE] = con[YE];
  param[conDS] = con[DS];

  // bracket for x
  double xlow = 1.0+q-s;
  double xup  = 2.0+2.0*q-s;

  // initial guess for temperature
  double temp_guess=prim[TEMP];

  // reset #EOS calls counter
  c2p_rep->nEOScalls = 0;

  // find x, this is the recovery process
  double x = zbrent_entropy(*func_root_entropy, eos, param, &temp_guess, xlow, xup, tol_x, c2p_rep, use_epsmin);

  // calculate final set of primitives
  calc_prim_entropy(eos, x, con, param, temp_guess, prim, S_squared, BdotS, B_squared, g_con, c2p_rep);

}

void calc_prim_entropy(const igm_eos_parameters eos,
                       const double x, const double * con, const double * param, const double temp_guess, double * prim, 
                       const double S_squared, const double BdotS, 
                       const double B_squared, const double g4up[4][4], struct c2p_report * c2p_rep ) {
  
  // Recover the primitive variables prim from x, r, s, t, con
  double r = param[par_r];
  double s = param[par_s];
  double t = param[par_t];

  double Wminus2 = 1.0 - ( x*x*r + (2*x+s)*t*t ) / ( x*x*(x+s)*(x+s) );
  Wminus2 = fmin(fmax(Wminus2 ,1e-2 ), 1.0);
  double W= pow(Wminus2, -0.5);

  double rho   = con[DD]/W;
  double ent   = con[DS]/W;
  double ye    = con[YE]/con[DD];
  double temp  = temp_guess;
  double press = 0.0;
  double eps   = 0.0;
  if( eos.evolve_T ) {
    get_P_eps_and_T_from_rho_Ye_and_S( eos,rho,ye,ent, &press,&eps,&temp );
  }
  else {
    get_P_eps_and_S_from_rho_Ye_and_T( eos,rho,ye,temp, &press,&eps,&ent );
  }

  const CCTK_REAL Z = x*rho*W;
  // Lower indices - covariant
  // CCTK_REAL B1_cov = g_cov[1][1]*con[B1_con]+g_cov[1][2]*con[B2_con]+g_cov[1][3]*con[B3_con];
  // CCTK_REAL B2_cov = g_cov[1][2]*con[B1_con]+g_cov[2][2]*con[B2_con]+g_cov[2][3]*con[B3_con];
  // CCTK_REAL B3_cov = g_cov[1][3]*con[B1_con]+g_cov[2][3]*con[B2_con]+g_cov[3][3]*con[B3_con];

  // Leo's note: The original routine output the *covariant* 3-velocity, v_{i}. This
  // is not what we need in IGM. So instead of lowering the indices of B^{i} and computing
  //
  // v_{i} = S_{i}/(z+B^{2}) + (B.S)B_{i}/(z(z+B^{2})) ,
  //
  // which was what was originally done (see commented lines of code above), we will instead
  // *raise* the indices of S_{i}, so that we may evaluate v^{i} using equation (24) of
  // Siegel et al., i.e.
  // .---------------------------------------------------.
  // | v^{i} = S^{i}/(z+B^{2}) + (B.S)B^{i}/(z(z+B^{2})) |
  // .---------------------------------------------------.
  // Raise indices - contravariant
  CCTK_REAL S1_con = g4up[1][1]*con[S1_cov] + g4up[1][2]*con[S2_cov] + g4up[1][3]*con[S3_cov];
  CCTK_REAL S2_con = g4up[1][2]*con[S1_cov] + g4up[2][2]*con[S2_cov] + g4up[2][3]*con[S3_cov];
  CCTK_REAL S3_con = g4up[1][3]*con[S1_cov] + g4up[2][3]*con[S2_cov] + g4up[3][3]*con[S3_cov];

  // These lines of code were used to evaluate v_{i}
  // prim[v1_cov  ] = (con[S1_cov] + (BdotS)*B1_cov/Z)/(Z+B_squared);
  // prim[v2_cov  ] = (con[S2_cov] + (BdotS)*B2_cov/Z)/(Z+B_squared);
  // prim[v3_cov  ] = (con[S3_cov] + (BdotS)*B3_cov/Z)/(Z+B_squared);
  // We replaced them with the ones below, used to compute v^{i} instead
  prim[UTCON1  ] = W*(S1_con + (BdotS)*con[B1_con]/Z)/(Z+B_squared);
  prim[UTCON2  ] = W*(S2_con + (BdotS)*con[B2_con]/Z)/(Z+B_squared);
  prim[UTCON3  ] = W*(S3_con + (BdotS)*con[B3_con]/Z)/(Z+B_squared);

  // These lines of code are kept the same
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

double func_root_entropy(const igm_eos_parameters eos, double x, double * param, struct c2p_report * c2p_rep, bool use_epsmin, double * temp_guess) {

  // computes f(x) from x and r,s,t
  const double r  = param[par_r];
  const double s  = param[par_s];
  const double t  = param[par_t];

  double Wminus2 = 1.0 - ( x*x*r + (2*x+s)*t*t  )/ (x*x*(x+s)*(x+s));  
  Wminus2 = fmin(fmax(Wminus2 ,1e-2 ), 1.0);
  const double W = pow(Wminus2, -0.5);

  double rho  = param[conDD]/W;
  double ye   = param[conYE]/param[conDD];
  double ent  = param[conDS]/W;
  double P    = 0.0;
  double temp = *temp_guess;
  double eps  = 0.0;

  if( eos.evolve_T ) {
    get_P_eps_and_T_from_rho_Ye_and_S( eos,rho,ye,ent, &P,&eps,&temp );
  }
  else {
    get_P_eps_and_S_from_rho_Ye_and_T( eos,rho,ye,temp, &P,&eps,&ent );
  }
  
  double ans = x - (1.0 + eps + P/rho)*W;

  return ans;

}

/****************************************************************************/
/************************* BRENT's METHOD FUNCTIONS *************************/
/****************************************************************************/

#define ITMAX 300   //Maximum allowed number of iterations.
#define PREC 3.0e-15
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


double zbrent_entropy(double (*func)(const igm_eos_parameters, double, double *, struct c2p_report *, bool, double *),
                      const igm_eos_parameters eos, double * param, 
                      double * temp_guess, double x1, double x2, double tol_x, struct c2p_report * c2p_rep, bool use_epsmin)
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
  double fa=(*func)(eos, a, param, c2p_rep, use_epsmin, &temp_guess_a);
  double fb=(*func)(eos, b, param, c2p_rep, use_epsmin, &temp_guess_b);
  double fc,p,q,r,s,tol1;

  // root must be bracketed
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    c2p_rep->failed = true;
    c2p_rep->c2p_keyerr = 1;
    c2p_rep->count = iter;

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
      c2p_rep->failed = false;
      c2p_rep->c2p_keyerr = 0;
      c2p_rep->count = iter;
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
      fb=(*func)(eos, b, param, c2p_rep, use_epsmin, &temp_guess_b);


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
      fb=(*func)(eos, b, param, c2p_rep, use_epsmin, &temp_guess_b);
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

  c2p_rep->count = iter;
 
  if( (!CCTK_isfinite(fb)) ) {
    c2p_rep->failed = true;
    c2p_rep->c2p_keyerr = 2;
    return b;
  } if( (fabs(maxerror) <= tol1 || fb == 0.0)){
    c2p_rep->failed = false;
    return ans;
  } else if( (fabs(maxerror) <= tol1) && (fabs(maxerror) > tol1) ){
    c2p_rep->failed = false;
    return ans;
  } else {
    c2p_rep->failed = true;
    c2p_rep->c2p_keyerr = 1;
    return b;
  }

}
