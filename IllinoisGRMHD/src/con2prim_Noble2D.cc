#include "cctk.h"
#include "cctk_Parameters.h"

#include "IllinoisGRMHD_headers.h"
#include "EOS_headers.hh"
#include "con2prim_headers.h"
#include "harm_u2p_util.h"

/***********************************************************************************
    Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble,
                   Gabor Toth, and Luca Del Zanna

                        HARM  version 1.0   (released May 1, 2006)

    This file is part of HARM.  HARM is a program that solves hyperbolic
    partial differential equations in conservative form using high-resolution
    shock-capturing techniques.  This version of HARM has been configured to
    solve the relativistic magnetohydrodynamic equations of motion on a
    stationary black hole spacetime in Kerr-Schild coordinates to evolve
    an accretion disk model.

    You are morally obligated to cite the following two papers in his/her
    scientific literature that results from use of any part of HARM:

    [1] Gammie, C. F., McKinney, J. C., \& Toth, G.\ 2003,
        Astrophysical Journal, 589, 444.

    [2] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006,
        Astrophysical Journal, 641, 626.


    Further, we strongly encourage you to obtain the latest version of
    HARM directly from our distribution website:
    http://rainman.astro.uiuc.edu/codelib/


    HARM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    HARM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HARM; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************

utoprim_2d.c:
---------------

    Uses the 2D method:
       -- solves for two independent variables (W,v^2) via a 2D
          Newton-Raphson method
       -- can be used (in principle) with a general equation of state.

  -- Currently returns with an error state (>0) if a negative rest-mass
      density or internal energy density is calculated.  You may want
      to change this aspect of the code so that it still calculates the
      velocity and so that you can floor the densities.  If you want to
      change this aspect of the code please comment out the "return(retval)"
      statement after "retval = 5;" statement in Utoprim_new_body();

******************************************************************************/

#ifdef  NEWT_DIM
#undef  NEWT_DIM
#endif
#define NEWT_DIM (2)

// Declarations:
CCTK_REAL vsq_calc(harm_aux_vars_struct& harm_aux,CCTK_REAL W);

int general_newton_raphson( igm_eos_parameters eos, harm_aux_vars_struct& harm_aux,
                            CCTK_REAL x[], int n,
                            void (*funcd)(igm_eos_parameters, harm_aux_vars_struct&,CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL [][NEWT_DIM], CCTK_REAL *, CCTK_REAL *, int) );

void func_vsq( igm_eos_parameters eos, harm_aux_vars_struct& harm_aux,CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], CCTK_REAL [][NEWT_DIM], CCTK_REAL *f, CCTK_REAL *df, int n);

CCTK_REAL x1_of_x0(harm_aux_vars_struct& harm_aux,CCTK_REAL x0 ) ;

CCTK_REAL pressure_W_vsq(igm_eos_parameters eos, CCTK_REAL W, CCTK_REAL vsq, CCTK_REAL &D) ;
CCTK_REAL dpdW_calc_vsq(igm_eos_parameters eos, CCTK_REAL W, CCTK_REAL vsq);
CCTK_REAL dpdvsq_calc(igm_eos_parameters eos, CCTK_REAL W, CCTK_REAL vsq, CCTK_REAL &D);
int Utoprim_new_body(const igm_eos_parameters eos,
                     const CCTK_REAL *restrict U,
                     const CCTK_REAL gcov[NDIM][NDIM],
                     const CCTK_REAL gcon[NDIM][NDIM],
                     CCTK_REAL *restrict prim);

/**********************************************************************/
/******************************************************************

  Utoprim_2d():

  -- Driver for new prim. var. solver.  The driver just translates
     between the two sets of definitions for U and P.  The user may
     wish to alter the translation as they see fit.  Note that Greek
     indices run 0,1,2,3 and Latin indices run 1,2,3 (spatial only).


             /     rho u^t     \
         U = | T^t_t + rho u^t |  sqrt(-det(g_{\mu\nu}))
             |      T^t_i      |
             \       B^i       /

             /     rho     \
         P = |     uu      |
             | \tilde{u}^i |
             \     B^i     /


   Arguments:
       U[NPR]           = conserved variables (current values on input/output);
       gcov[NDIM][NDIM] = covariant form of the metric ;
       gcon[NDIM][NDIM] = contravariant form of the metric ;
       gdet             = sqrt( - determinant of the metric) ;
       prim[NPR]        = primitive variables (guess on input, calculated values on
                                                output if there are no problems);

   -- NOTE: for those using this routine for special relativistic MHD and are
            unfamiliar with metrics, merely set
              gcov = gcon = diag(-1,1,1,1)  and gdet = 1.  ;

******************************************************************/


int con2prim_Noble2D( const igm_eos_parameters eos,
                      const CCTK_REAL g4dn[4][4],
                      const CCTK_REAL g4up[4][4],
                      const CCTK_REAL *restrict cons,
                      CCTK_REAL *restrict prim ) {

  return( Utoprim_new_body(eos, cons, g4dn, g4up, prim) );

}



/**********************************************************************/
/**********************************************************************************

  Utoprim_new_body():

     -- Attempt an inversion from U to prim using the initial guess prim.

     -- This is the main routine that calculates auxiliary quantities for the
        Newton-Raphson routine.

  -- assumes that
             /   rho gamma   \
         U = | alpha T^t_\mu |
             \   alpha B^i   /



             /     rho     \
      prim = |     uu      |
             | \tilde{u}^i |
             \  alpha B^i  /


return:  (i*100 + j)  where
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used)
                   or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the
                   given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence
                   (occurrence of "nan" or "+/-inf" ;

         j = 0 -> success
             1 -> failure: some sort of failure in Newton-Raphson;
             2 -> failure: utsq<0 w/ initial p[] guess;
             3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: v^2 > 1
             5 -> failure: rho,uu <= 0 ;

**********************************************************************************/

int Utoprim_new_body( const igm_eos_parameters eos,
                      const CCTK_REAL *restrict U,
                      const CCTK_REAL gcov[NDIM][NDIM],
                      const CCTK_REAL gcon[NDIM][NDIM],
                      CCTK_REAL *restrict prim )
{

  CCTK_REAL x_2d[NEWT_DIM];
  CCTK_REAL Bcon[NDIM],Bcov[NDIM],Qcov[NDIM],Qcon[NDIM],ncov[NDIM],ncon[NDIM],Qtcon[NDIM];
  CCTK_REAL rho0,u,p,w,gammasq,gtmp,W_last,W,utsq,vsq;
  int i,j, n, retval, i_increase;

  // Contains Bsq,QdotBsq,Qsq,Qtsq,Qdotn,QdotB,D,gamma,gamma_times_S,ye
  harm_aux_vars_struct harm_aux;

  n = NEWT_DIM ;

  // Assume ok initially:
  retval = 0;

  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
  Bcon[0] = 0. ;
  for(i=1;i<4;i++) Bcon[i] = U[BCON1+i-1] ;

  lower_g(Bcon,gcov,Bcov) ;

  for(i=0;i<4;i++) Qcov[i] = U[QCOV0+i] ;
  raise_g(Qcov,gcon,Qcon) ;


  harm_aux.Bsq = 0. ;
  for(i=1;i<4;i++) harm_aux.Bsq += Bcon[i]*Bcov[i] ;

  harm_aux.QdotB = 0. ;
  for(i=0;i<4;i++) harm_aux.QdotB += Qcov[i]*Bcon[i] ;
  harm_aux.QdotBsq = harm_aux.QdotB*harm_aux.QdotB ;

  ncov_calc(gcon,ncov) ;
  // FIXME: The exact form of n^{\mu} can be found
  //        in eq. (2.116) and implementing it
  //        directly is a lot more efficient than
  //        performing n^{\mu} = g^{\mu\nu}n_{nu}
  raise_g(ncov,gcon,ncon);

  harm_aux.Qdotn = Qcon[0]*ncov[0] ;

  harm_aux.Qsq = 0. ;
  for(i=0;i<4;i++) harm_aux.Qsq += Qcov[i]*Qcon[i] ;

  harm_aux.Qtsq = harm_aux.Qsq + harm_aux.Qdotn*harm_aux.Qdotn ;

  harm_aux.D    = U[RHO];

  /* calculate W from last timestep and use for guess */
  utsq = 0. ;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++) utsq += gcov[i][j]*prim[UTCON1+i-1]*prim[UTCON1+j-1] ;


  if( (utsq < 0.) && (fabs(utsq) < 1.0e-13) ) {
    utsq = fabs(utsq);
  }
  if(utsq < 0. || utsq > UTSQ_TOO_BIG) {
    retval = 2;
    return(retval) ;
  }

  gammasq = 1. + utsq ;
  harm_aux.gamma  = sqrt(gammasq);

  // Always calculate rho from D and gamma so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq) .
  rho0 = harm_aux.D / harm_aux.gamma ;
  u = prim[UU];

  p = 0;
  if( eos.is_Hybrid ) {
    p = pressure_rho0_u(eos, rho0,u);
  }
  else if( eos.is_Tabulated ) {
    harm_aux.ye            = U[YE]/U[RHO];
    harm_aux.gamma_times_S = U[WS];
    harm_aux.use_entropy   = false;
    CCTK_REAL xrho         = rho0;
    CCTK_REAL xye          = harm_aux.ye;
    CCTK_REAL xtemp        = prim[TEMP];
    CCTK_REAL xprs         = 0.0;
    CCTK_REAL xeps         = 0.0;
    CCTK_REAL xdepsdT      = 0.0;

    // Now compute P and eps from (rho,Ye,T). Note that
    // at this point we do not know W, so we do not
    // use the entropy in this function call.
    get_P_eps_and_depsdT_from_rho_Ye_and_T( eos,xrho,xye,xtemp, &xprs,&xeps,&xdepsdT );
    p = xprs;
    u = xeps*xrho;
    if( xdepsdT < eos.depsdT_threshold ) harm_aux.use_entropy = true;
  }

  w = rho0 + u + p ;
  W_last = w*gammasq ;

  // Make sure that W is large enough so that v^2 < 1 :
  i_increase = 0;
  while( (( W_last*W_last*W_last * ( W_last + 2.*harm_aux.Bsq )
            - harm_aux.QdotBsq*(2.*W_last + harm_aux.Bsq) ) <= W_last*W_last*(harm_aux.Qtsq-harm_aux.Bsq*harm_aux.Bsq))
         && (i_increase < 10) ) {
    W_last *= 10.;
    i_increase++;
  }

  // Calculate W and vsq:
  x_2d[0] =  fabs( W_last );
  x_2d[1] = x1_of_x0( harm_aux,W_last ) ;

  retval = general_newton_raphson( eos, harm_aux,x_2d, n, func_vsq) ;

  W = x_2d[0];
  vsq = x_2d[1];

  /* Problem with solver, so return denoting error before doing anything further */
  if( (retval != 0) || (W == FAIL_VAL) ) {
    retval = retval*100+1;
    return(retval);
  }
  else{
    if(W <= 0. || W > W_TOO_BIG) {
      retval = 3;
      return(retval) ;
    }
  }

  // Calculate v^2:
  if( vsq >= 1. ) {
    vsq = 1.-2.e-16;
    //retval = 4;
    //return(retval) ;
  }


  // Recover the primitive variables from the scalars and conserved variables:
  gtmp = sqrt(1. - vsq);
  harm_aux.gamma = 1./gtmp ;
  rho0 = harm_aux.D * gtmp;

  w = W * (1. - vsq) ;

  if( eos.is_Hybrid ) {
    p = pressure_rho0_w(eos, rho0,w) ;
    u = w - (rho0 + p) ; // u = rho0 eps, w = rho0 h
    prim[RHO] = rho0 ;
    prim[UU ] = u ;
  }
  else {
    CCTK_REAL xrho  = rho0;
    CCTK_REAL xye   = harm_aux.ye;
    CCTK_REAL xtemp = prim[TEMP];
    CCTK_REAL xent  = harm_aux.gamma_times_S / W;
    CCTK_REAL xprs=0,xuu=0,xeps=0;
    if( harm_aux.use_entropy ) {
      get_P_eps_and_T_from_rho_Ye_and_S( eos,xrho,xye,xent, &xprs,&xeps,&xtemp );
    }
    else {
      xprs  = -0.5*harm_aux.Bsq/(harm_aux.gamma*harm_aux.gamma)+harm_aux.Qdotn+W+harm_aux.Bsq-0.5*harm_aux.QdotBsq/(W*W);;
      xuu   = (W-harm_aux.D*harm_aux.gamma-xprs*harm_aux.gamma*harm_aux.gamma)/(harm_aux.D*harm_aux.gamma) * rho0;
      xeps  = xuu/xrho;
      get_P_S_and_T_from_rho_Ye_and_eps( eos,xrho,xye,xeps, &xprs,&xent,&xtemp );      
    }
    
    // Update P and T in the prim array
    prim[RHO  ] = xrho;
    prim[YE   ] = harm_aux.ye;
    prim[TEMP ] = MIN(MAX(xtemp,eos.T_atm),eos.T_max);
    prim[PRESS] = xprs;
    prim[EPS  ] = xeps;
    prim[ENT  ] = xent;
  }

  if( (rho0 <= 0.) || (u <= 0.) ) {
    // User may want to handle this case differently, e.g. do NOT return upon
    // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:

    retval = 5;
    //return(retval) ;
  }

  for(i=1;i<4;i++) Qtcon[i] = Qcon[i] + ncon[i] * harm_aux.Qdotn;
  for(i=1;i<4;i++) prim[UTCON1+i-1] = harm_aux.gamma/(W+harm_aux.Bsq) * ( Qtcon[i] + harm_aux.QdotB*Bcon[i]/W ) ;

  /* set field components */
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;

  /* done! */
  return(retval) ;

}



/**********************************************************************/
/****************************************************************************
   vsq_calc():

      -- evaluate v^2 (spatial, normalized velocity) from
            W = \gamma^2 w

****************************************************************************/
CCTK_REAL vsq_calc(harm_aux_vars_struct& harm_aux,CCTK_REAL W)
{
  CCTK_REAL Wsq = W*W ;
  CCTK_REAL Xsq = (harm_aux.Bsq + W) * (harm_aux.Bsq + W);

  return(  ( Wsq * harm_aux.Qtsq  + harm_aux.QdotBsq * (harm_aux.Bsq + 2.*W)) / (Wsq*Xsq) );
}



/********************************************************************

  x1_of_x0():

    -- calculates v^2 from W  with some physical bounds checking;
    -- asumes x0 is already physical
    -- makes v^2 physical  if not;

*********************************************************************/

CCTK_REAL x1_of_x0(harm_aux_vars_struct& harm_aux,CCTK_REAL x0 )
{
  CCTK_REAL dv  = 1.e-15;
  CCTK_REAL vsq = fabs(vsq_calc(harm_aux,x0)) ; // guaranteed to be positive

  return( ( vsq > 1. ) ? (1.0 - dv) : vsq   );
}


/********************************************************************

  validate_x():

    -- makes sure that x[0,1] have physical values, based upon
       their definitions:

*********************************************************************/

void validate_x(CCTK_REAL x[2], CCTK_REAL x0[2] )
{

  CCTK_REAL dv = 1.e-15;

  /* Always take the absolute value of x[0] and check to see if it's too big:  */
  x[0] = fabs(x[0]);
  x[0] = (x[0] > W_TOO_BIG) ?  x0[0] : x[0];


  x[1] = (x[1] < 0.) ?   0.       : x[1];  /* if it's too small */
  x[1] = (x[1] > 1.) ?  (1. - dv) : x[1];  /* if it's too big   */

  return;

}


/************************************************************

  general_newton_raphson():

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

*****************************************************************/
int general_newton_raphson( igm_eos_parameters eos, harm_aux_vars_struct& harm_aux,
                            CCTK_REAL x[], int n,
                            void (*funcd)(igm_eos_parameters, harm_aux_vars_struct&,CCTK_REAL [], CCTK_REAL [], CCTK_REAL [],
                                          CCTK_REAL [][NEWT_DIM], CCTK_REAL *, CCTK_REAL *, int) )
{
  CCTK_REAL f, df, dx[NEWT_DIM], x_old[NEWT_DIM];
  CCTK_REAL resid[NEWT_DIM], jac[NEWT_DIM][NEWT_DIM];
  CCTK_REAL errx, x_orig[NEWT_DIM];
  int id, i_extra, doing_extra;
  int n_iter;
  int keep_iterating;


  // Initialize various parameters and variables:
  n_iter = 0;
  errx = 1. ;
  df = f = 1.;
  i_extra = doing_extra = 0;
  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;

  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) {

    (*funcd) (eos, harm_aux, x, dx, resid, jac, &f, &df, n);  /* returns with new dx, f, df */


    /* Save old values before calculating the new: */
    errx = 0.;
    for( id = 0; id < n ; id++) {
      x_old[id] = x[id] ;
    }

    /* Make the newton step: */
    for( id = 0; id < n ; id++) {
      x[id] += dx[id]  ;
    }

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    validate_x( x, x_old ) ;


    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/

    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0))
        || (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)

  /*  Check for bad untrapped divergences : */
  if( (finite(f)==0) ||  (finite(df)==0) ) {
    return(2);
  }


  if( fabs(errx) > MIN_NEWT_TOL){
    //CCTK_VInfo(CCTK_THORNSTRING,"%d %e %e %e %e",n_iter,f,df,errx,MIN_NEWT_TOL);
    return(1);
  }
  if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    return(0);
  }
  if( fabs(errx) <= NEWT_TOL ){
    return(0);
  }

  return(0);

}




/**********************************************************************/
/*********************************************************************************
   func_vsq():

        -- calculates the residuals, and Newton step for general_newton_raphson();
        -- for this method, x=W,vsq here;

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
*********************************************************************************/

void func_vsq(igm_eos_parameters eos, harm_aux_vars_struct& harm_aux, CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[],
              CCTK_REAL jac[][NEWT_DIM], CCTK_REAL *f, CCTK_REAL *df, int n)
{

  CCTK_REAL W   = x[0];
  CCTK_REAL vsq = x[1];
  CCTK_REAL Wsq = W*W;

  CCTK_REAL p_tmp, dPdvsq, dPdW;

  if( eos.is_Hybrid ) {
    p_tmp  = pressure_W_vsq( eos, W, vsq , harm_aux.D);
    dPdW   = dpdW_calc_vsq( eos, W, vsq );
    dPdvsq = dpdvsq_calc( eos, W, vsq, harm_aux.D );
  }
  else {
    // Here we must compute dPdW and dPdvsq for tabulated EOS.
    // We will be using the expressions 
    CCTK_REAL gamma_sq = 1.0/(1.0-vsq);
    harm_aux.gamma     = sqrt(gamma_sq);
    if( harm_aux.gamma > eos.W_max ) {
      harm_aux.gamma = eos.W_max;
      gamma_sq       = harm_aux.gamma * harm_aux.gamma;
    }
    
    const CCTK_REAL rho = MAX(harm_aux.D / harm_aux.gamma,eos.rho_min);
    const CCTK_REAL xye = harm_aux.ye;
    const CCTK_REAL h   = fabs(W /(rho*gamma_sq)); // W := rho*h*gamma^{2}
    const CCTK_REAL ent = harm_aux.gamma_times_S / harm_aux.gamma;
    CCTK_REAL T         = eos.T_atm;
    CCTK_REAL prs       = 0.0;
    CCTK_REAL eps       = 0.0;
    CCTK_REAL dPdrho    = 0.0;
    CCTK_REAL dPdeps    = 0.0;

    // Now compute the pressure and its derivatives with respect to rho and eps
    if( harm_aux.use_entropy ) {
      get_P_eps_T_dPdrho_and_dPdeps_from_rho_Ye_and_S(eos,rho,xye,ent, &prs,&T,&eps,&dPdrho,&dPdeps);
    }
    else {
      get_P_eps_T_dPdrho_and_dPdeps_from_rho_Ye_and_h(eos,rho,xye,h, &prs,&eps,&T,&dPdrho,&dPdeps);
    }

    // Set P
    p_tmp = prs;

    // Now compute dP/dW
    const CCTK_REAL dPdeps_o_rho = dPdeps/rho;
    dPdW = ( dPdeps_o_rho/(1.+dPdeps_o_rho) )/gamma_sq;

    // And finally dP/d(v^{2})
    const CCTK_REAL dPdvsq_1 = -0.5*harm_aux.D*harm_aux.gamma*dPdrho;
    const CCTK_REAL dPdvsq_2 = -0.5*(W + prs*gamma_sq)/rho;
    dPdvsq = (dPdvsq_1 + dPdeps*dPdvsq_2)/(1+dPdeps_o_rho);
  }

  // These expressions were calculated using Mathematica, but made into efficient
  // code using Maple.  Since we know the analytic form of the equations, we can
  // explicitly calculate the Newton-Raphson step:

  CCTK_REAL t2  = -0.5*harm_aux.Bsq+dPdvsq;
  CCTK_REAL t3  = harm_aux.Bsq+W;
  CCTK_REAL t4  = t3*t3;
  CCTK_REAL t9  = 1/Wsq;
  CCTK_REAL t11 = harm_aux.Qtsq-vsq*t4+harm_aux.QdotBsq*(harm_aux.Bsq+2.0*W)*t9;
  CCTK_REAL t16 = harm_aux.QdotBsq*t9;
  CCTK_REAL t18 = -harm_aux.Qdotn-0.5*harm_aux.Bsq*(1.0+vsq)+0.5*t16-W+p_tmp;
  CCTK_REAL t21 = 1/t3;
  CCTK_REAL t23 = 1/W;
  CCTK_REAL t24 = t16*t23;
  CCTK_REAL t25 = -1.0+dPdW-t24;
  CCTK_REAL t35 = t25*t3+(harm_aux.Bsq-2.0*dPdvsq)*(harm_aux.QdotBsq+vsq*Wsq*W)*t9*t23;
  CCTK_REAL t36 = 1/t35;
  CCTK_REAL t40 = (vsq+t24)*t3;
  dx[0] = -(t2*t11+t4*t18)*t21*t36;
  dx[1] = -(-t25*t11-2.0*t40*t18)*t21*t36;
  //detJ = t3*t35; // <- set but not used...
  jac[0][0] = -2.0*t40;
  jac[0][1] = -t4;
  jac[1][0] = t25;
  jac[1][1] = t2;
  resid[0] = t11;
  resid[1] = t18;



  *df = -resid[0]*resid[0] - resid[1]*resid[1];

  *f = -0.5 * ( *df );

}



/**********************************************************************
 **********************************************************************

 The following routines specify the equation of state.  All routines
  above here should be indpendent of EOS.  If the user wishes
  to use another equation of state, the below functions must be replaced
  by equivalent routines based upon the new EOS.

  **********************************************************************
  **********************************************************************/


/**********************************************************************/
/**********************************************************************
  pressure_W_vsq():

        -- Hybrid single and piecewise polytropic equation of state;
        -- pressure as a function of P_cold, eps_cold, W, vsq, and D:
**********************************************************************/
CCTK_REAL pressure_W_vsq(igm_eos_parameters eos, CCTK_REAL W, CCTK_REAL vsq, CCTK_REAL &D)
{

  // Compute gamma^{-2} = 1 - v^{2} and gamma^{-1}
  CCTK_REAL inv_gammasq = 1.0 - vsq;
  CCTK_REAL inv_gamma   = sqrt(inv_gammasq);

  // Compute rho_b = D / gamma
  CCTK_REAL rho_b = D*inv_gamma;

  // Compute P_cold and eps_cold
  CCTK_REAL P_cold, eps_cold;
  compute_P_cold__eps_cold(eos,rho_b, P_cold,eps_cold);

  // Compute p = P_{cold} + P_{th}
  return( ( P_cold + (eos.Gamma_th - 1.0)*( W*inv_gammasq - D*inv_gamma*( 1.0 + eps_cold ) ) )/eos.Gamma_th );

}



/**********************************************************************/
/**********************************************************************
  dpdW_calc_vsq():

      -- partial derivative of pressure with respect to W;
**********************************************************************/
CCTK_REAL dpdW_calc_vsq(igm_eos_parameters eos,CCTK_REAL W, CCTK_REAL vsq)
{

  return( (eos.Gamma_th - 1.0) * (1.0 - vsq) /  eos.Gamma_th  ) ;

}


/**********************************************************************/
/**********************************************************************
  dpdvsq_calc():

      -- partial derivative of pressure with respect to vsq
**********************************************************************/
CCTK_REAL dpdvsq_calc(igm_eos_parameters eos, CCTK_REAL W, CCTK_REAL vsq, CCTK_REAL &D)
{

  // Set gamma and rho
  CCTK_REAL gamma = 1.0/sqrt(1.0 - vsq);
  CCTK_REAL rho_b = D/gamma;

  // Compute P_cold and eps_cold
  CCTK_REAL P_cold, eps_cold;
  compute_P_cold__eps_cold(eos,rho_b, P_cold,eps_cold);

  // Set basic polytropic quantities
  int polytropic_index = find_polytropic_K_and_Gamma_index(eos,rho_b);
  CCTK_REAL Gamma_ppoly_tab = eos.Gamma_ppoly_tab[polytropic_index];


  /* Now we implement the derivative of P_cold with respect
   * to v^{2}, given by
   *  ----------------------------------------------------
   * | dP_cold/dvsq = gamma^{2 + Gamma_{poly}/2} P_{cold} |
   *  ----------------------------------------------------
   */
  CCTK_REAL dPcold_dvsq = P_cold * pow(gamma,2.0 + 0.5*Gamma_ppoly_tab);


  /* Now we implement the derivative of eps_cold with respect
   * to v^{2}, given by
   *  -----------------------------------------------------------------------------------
   * | deps_cold/dvsq = gamma/(D*(Gamma_ppoly_tab-1)) * (dP_cold/dvsq + gamma^{2} P_cold / 2) |
   *  -----------------------------------------------------------------------------------
   */
  CCTK_REAL depscold_dvsq = ( gamma/(D*(Gamma_ppoly_tab-1.0)) ) * ( dPcold_dvsq + 0.5*gamma*gamma*P_cold );

  /* Now we implement the derivative of p_hybrid with respect
   * to v^{2}, given by
   *  -----------------------------------------------------------------------------
   * | dp/dvsq = Gamma_th^{-1}( dP_cold/dvsq                                       |
   * |                          + (Gamma_{th}-1)*(-W                               |
   * |                                            + D gamma (1 + eps_cold)/2       |
   * |                                            - (D/gamma) * deps_cold/dvsq) )  |
   *  -----------------------------------------------------------------------------
   */
  return( ( dPcold_dvsq + (eos.Gamma_th-1.0)*( -W + D*gamma*(1+eps_cold)/2.0 - D*depscold_dvsq/gamma ) )/eos.Gamma_th );
}


/******************************************************************************
             END   OF   UTOPRIM_2D.C
******************************************************************************/
