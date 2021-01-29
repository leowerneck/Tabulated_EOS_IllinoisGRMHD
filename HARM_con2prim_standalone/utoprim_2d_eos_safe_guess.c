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

utoprim_2d_fast.c: 
---------------

-- optimized version;  changes made near "//-fast" 

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

#include "u2p_defs.h"
#include "options.h"

#define NEWT_DIM 2


static void raise(double *ucov, struct of_geom *geom, double *ucon)
{

        ucon[0] = geom->gcon[0][0]*ucov[0]
                + geom->gcon[0][1]*ucov[1]
                + geom->gcon[0][2]*ucov[2]
                + geom->gcon[0][3]*ucov[3] ;
        ucon[1] = geom->gcon[1][0]*ucov[0]
                + geom->gcon[1][1]*ucov[1]
                + geom->gcon[1][2]*ucov[2]
                + geom->gcon[1][3]*ucov[3] ;
        ucon[2] = geom->gcon[2][0]*ucov[0]
                + geom->gcon[2][1]*ucov[1]
                + geom->gcon[2][2]*ucov[2]
                + geom->gcon[2][3]*ucov[3] ;
        ucon[3] = geom->gcon[3][0]*ucov[0]
                + geom->gcon[3][1]*ucov[1]
                + geom->gcon[3][2]*ucov[2]
                + geom->gcon[3][3]*ucov[3] ;

        return ;
}


/* these variables need to be shared between the functions
   Utoprim_1D, residual, and utsq */

FTYPE Bsq,QdotBsq,Qtsq,Qdotn,D,alpha,half_Bsq ;

extern void EOS_press_with_T(double rho, double * eps,
        double * temp, double ye, double * prs);

extern void EOS_P_from_hrho_dPdrho_dPdeps(const double rho, const double enth,
          double * temp_guess, double ye, double * eps, double * press, double * dPdrho,
          double * dPdeps, double * entr);

extern void EOS_press(double rho, double * eps,
        double * temp, double ye, double * prs);
// Declarations: 
static FTYPE vsq_calc(FTYPE W);
static int Utoprim_new_body(FTYPE *U, struct of_geom *geom, FTYPE *prim, FTYPE *gamma, FTYPE *bsq);
static int general_newton_raphson( FTYPE x[], int n, 
            FTYPE D, FTYPE ye0, FTYPE * temp_guess, 
				   void (*funcd) (FTYPE [], FTYPE [], 
						  FTYPE [], FTYPE [][NEWT_DIM], 
						  FTYPE *, FTYPE *, FTYPE , FTYPE, FTYPE *) );

static  void func_vsq( FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], 
		       FTYPE *f, FTYPE *df, FTYPE D, FTYPE ye0, FTYPE * temp_guess);
static FTYPE x1_of_x0(FTYPE x0 ) ;
static FTYPE eos_info(FTYPE W, FTYPE vsq, FTYPE *dpdw, FTYPE *dpdvsq);
static  FTYPE eos_quantities_general(FTYPE W, FTYPE vsq, FTYPE *dpdw, FTYPE *dpdvsq, 
         FTYPE D, FTYPE ye0, FTYPE * temp_guess,  FTYPE *p_tmp);

/**********************************************************************/
/******************************************************************

  Utoprim_2d_fast():
  
  -- Driver for new prim. var. solver.  The driver just translates
     between the two sets of definitions for U and P.  The user may 
     wish to alter the translation as they see fit.  

     It assumes that on input/output:

              /  rho u^t           \
         U =  |  T^t_t   + rho u^t |  sqrt(-det(g_{\mu\nu}))  
              |  T^t_\mu           |  
              \   B^i              /

             /    rho        \
	 P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /


   Arguments:
       U[NPR]    = conserved variables (current values on input/output);
       gcov[NDIM][NDIM] = covariant form of the metric ;
       gcon[NDIM][NDIM] = contravariant form of the metric ;
       gdet             = sqrt( - determinant of the metric) ;
       prim[NPR] = primitive variables (guess on input, calculated values on 
                                        output if there are no problems);
  
   -- NOTE: for those using this routine for special relativistic MHD and are
            unfamiliar with metrics, merely set 
              gcov = gcon = diag(-1,1,1,1)  and gdet = 1.  ;

******************************************************************/
int Utoprim_2d_eos_safe_guess(FTYPE *U, struct of_geom *geom, FTYPE *prim, FTYPE *gamma, FTYPE *bsq)
{

  FTYPE U_tmp[NPR], prim_tmp[NPR];
  int i, j, ret; 
  FTYPE geomfactor,inv_gdet;

  
  /* Set the geometry variables: */
  inv_gdet = geom->g_inv; 
  alpha = geom->alpha;
  geomfactor = alpha * inv_gdet;

  /* First update the primitive B-fields */
  prim[BCON1] = U[BCON1] * inv_gdet ;
  prim[BCON2] = U[BCON2] * inv_gdet ;
  prim[BCON3] = U[BCON3] * inv_gdet ;

  if( U[0] <= 0. ) {
    fprintf(stderr, "Negative U[0] found!!  We encourage you to figure out how this weird thing happened!! \n"); 
    return(-100);
  }
  
  /* Transform the CONSERVED variables into the new system */
  U_tmp[RHO]    = geomfactor * U[RHO] ;
  U_tmp[UU]     = geomfactor * (U[UU] - U[RHO]) ;
  U_tmp[UTCON1] = geomfactor * U[UTCON1] ;
  U_tmp[UTCON2] = geomfactor * U[UTCON2] ;
  U_tmp[UTCON3] = geomfactor * U[UTCON3] ;
  U_tmp[BCON1 ] = geomfactor * U[BCON1] ;
  U_tmp[BCON2 ] = geomfactor * U[BCON2] ;
  U_tmp[BCON3 ] = geomfactor * U[BCON3] ;
  U_tmp[YE]     = geomfactor * U[YE] ;

  
  /* Transform the PRIMITIVE variables into the new system */

  prim_tmp[RHO   ] = prim[RHO   ];
  prim_tmp[UU    ] = prim[UU    ];
  prim_tmp[UTCON1] = prim[UTCON1];
  prim_tmp[UTCON2] = prim[UTCON2];
  prim_tmp[UTCON3] = prim[UTCON3];

  prim_tmp[BCON1] =  U_tmp[BCON1 ] ;
  prim_tmp[BCON2] =  U_tmp[BCON2 ] ;
  prim_tmp[BCON3] =  U_tmp[BCON3 ] ;

  prim_tmp[YE    ] = prim[YE    ];
  prim_tmp[TEMP    ] = prim[TEMP    ];

  ret = Utoprim_new_body(U_tmp, geom, prim_tmp, gamma, bsq);

  /* Use new primitive variables if there was no problem : */ 
  if( ret == 0 ) {
    prim[RHO   ] = prim_tmp[RHO   ];
    prim[UU    ] = prim_tmp[UU    ];
    prim[UTCON1] = prim_tmp[UTCON1];
    prim[UTCON2] = prim_tmp[UTCON2];
    prim[UTCON3] = prim_tmp[UTCON3];
    prim[YE    ] = prim_tmp[YE    ];
    prim[TEMP    ] = prim_tmp[TEMP    ];
  }

  return( ret ) ;

}


/**********************************************************************/
/**********************************************************************************

  Utoprim_new_body():

     -- Attempt an inversion from U to prim using the initial guess prim.

     -- This is the main routine that calculates auxiliary quantities for the 
        Newton-Raphson routine. 

  -- assumes that 
             /  rho gamma        \
         U = |  alpha T^t_\mu    |
             \  alpha B^i        /



               /    rho        \
	prim = |    uu         |
               | \tilde{u}^i   |
               \  alpha B^i   /


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
static int Utoprim_new_body(FTYPE *U, struct of_geom *geom, FTYPE *prim, FTYPE *gamma_out, FTYPE *bsq)
{

  FTYPE x_2d[NEWT_DIM];
  FTYPE QdotB,Bcon[NDIM],*Qcov,Qcon[NDIM],Qsq;
  FTYPE rho0,ye0,u,p,w,gammasq,gamma,gtmp,W_last,W,utsq,vsq;
  int i,j, retval, i_increase ;


  // Assume ok initially:
  retval = 0;

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
  Bcon[0] = 0. ;
  Bcon[1] = U[BCON1];
  Bcon[2] = U[BCON2];
  Bcon[3] = U[BCON3];
  //-faster   lower(Bcon,geom,Bcov);
  //-faster   Bsq = Bcon[1]*Bcov[1] + Bcon[2]*Bcov[2] + Bcon[3]*Bcov[3] ;
  Bsq =   geom->gcov[1][1]*Bcon[1]*Bcon[1]
        + geom->gcov[2][2]*Bcon[2]*Bcon[2]
        + geom->gcov[3][3]*Bcon[3]*Bcon[3]
     + 2*(geom->gcov[1][2]*Bcon[1]*Bcon[2]
	+ geom->gcov[1][3]*Bcon[1]*Bcon[3]
	+ geom->gcov[2][3]*Bcon[2]*Bcon[3]) ;


  Qcov    = &(U[QCOV0]) ;  
  raise(Qcov,geom,Qcon) ;

  QdotB = Qcov[1]*Bcon[1] + Qcov[2]*Bcon[2] + Qcov[3]*Bcon[3] ;
  QdotBsq = QdotB*QdotB ;

  //-fast  ncov_calc(gcon,ncov) ;
//-faster  ncov[0] = -geom->alpha;
//-faster  ncov[1] = ncov[2] = ncov[3] = 0.;

  //-fast   raise_g(ncov,gcon,ncon);
//-faster   ncon[0] = -alpha * geom->gcon[0][0]; 
//-faster   ncon[1] = -alpha * geom->gcon[0][1]; 
//-faster   ncon[2] = -alpha * geom->gcon[0][2]; 
//-faster   ncon[3] = -alpha * geom->gcon[0][3]; 
  
  //-faster  Qdotn = Qcon[0]*ncov[0] ;
  Qdotn = -alpha*Qcon[0];

  //-fast Qsq = 0. ;
  //-fast  for(i=0;i<4;i++) Qsq += Qcov[i]*Qcon[i] ;
  Qtsq = DOT(Qcov,Qcon);
  Qtsq += Qdotn*Qdotn ;

  D = U[RHO] ;
  ye0=U[YE]/U[RHO];
  
  half_Bsq = 0.5*Bsq;

  /* calculate W from last timestep and use for guess */
  //-fast   utsq = 0. ;
  //-fast   for(i=1;i<4;i++)
  //-fast     for(j=1;j<4;j++) utsq += gcov[i][j]*prim[UTCON1+i-1]*prim[UTCON1+j-1] ;


  gammasq = GAMMAMAX*GAMMAMAX ;
  gamma  = sqrt(gammasq);
	
  // Always calculate rho from D and gamma so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq) . 
  //getting safe guesses, Cerd'a-Dur'an et al. (2008)
  rho0 = D ;

  double tau;
  tau=-Qdotn-D;
  u=tau-Bsq/2;


  //fprintf(stderr, "prim[UU] %e\n", prim[UU]);
  //fprintf(stderr, "TEMP %d\n", TEMP);

   double temp_guess;
   temp_guess=prim[TEMP];
  
#if( TABULATED_EOS_TABLES )
  if(rho0 < eos_rhomin) {
    rho0=eos_rhomin;
  }
   
  double xeps = u/rho0;
  double xtemp = temp_guess;
  double xP = 0.;
  //EOS_press_with_T(rho0,&xeps,&xtemp,ye0, &xP);
  EOS_press(rho0,&xeps,&xtemp,ye0, &xP);
  p = xP;

#else 
 
  p = pressure_rho0_u(rho0,u) ;

#endif



  W_last = tau+p+D-Bsq/2. ;

//  fprintf(stdout,"Conserved = %26.16e %26.16e %26.16e %26.16e %26.16e \n", U[0], U[1], U[2], U[3], U[4]); 
//  fprintf(stdout,"Prims     = %26.16e %26.16e %26.16e %26.16e %26.16e \n", prim[0], prim[1], prim[2], prim[3], prim[4]); 
//  fprintf(stdout,"newPrims  = %26.16e %26.16e %26.16e %26.16e %26.16e \n", rho0,u,p,w,W_last);
//  fprintf(stdout,"OTHERS    = %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e \n", Qdotn,Qsq,D,half_Bsq,utsq,gamma); 
//  fflush(stdout);

  // Make sure that W is large enough so that v^2 < 1 : 
  i_increase = 0;
  while( (( W_last*W_last*W_last * ( W_last + 2.*Bsq ) 
	    - QdotBsq*(2.*W_last + Bsq) ) <= W_last*W_last*(Qtsq-Bsq*Bsq))
	 && (i_increase < 10) ) {
    W_last *= 10.;
    i_increase++;
  }
  
  // Calculate W and vsq: 
  x_2d[0] =  fabs( W_last );
  x_2d[1] = 1.-1./(gammasq);//x1_of_x0( W_last ) ;
  retval = general_newton_raphson( x_2d, NEWT_DIM, D, ye0, &temp_guess, func_vsq ) ;  

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
    retval = 4;
    return(retval) ;
  }

     double gamma_sq_nr, gamma_nr;
     gamma_sq_nr = 1.0/(1. -vsq);
     gamma_nr=sqrt(gamma_sq_nr);

     if (gamma_nr>GAMMAMAX){
           retval = 5;
           return(retval) ;
       }


  // Recover the primitive variables from the scalars and conserved variables:
  gtmp = sqrt(1. - vsq);
  gamma = 1./gtmp ;
  rho0 = D * gtmp;
  *gamma_out = gamma;

  w = W * (1. - vsq) ;


  
#if( TABULATED_EOS_TABLES )
//   FTYPE p_tmp,dPdvsq, dPdW;;
//    double xtemp2;
//    xtemp2=temp_guess;
//    eos_quantities_general(W,vsq, &dPdW, &dPdvsq, D, ye0, &xtemp2, &p_tmp);
//    u=(W-D*gamma-p_tmp*gamma*gamma)/(D*gamma) * rho0;
    double xtemp2;
    xtemp2=temp_guess;
    double p_tmp;  
    p_tmp=-0.5*Bsq/(gamma*gamma) +Qdotn+W+Bsq-0.5*QdotBsq/(W*W);;
    u=(W-D*gamma-p_tmp*gamma*gamma)/(D*gamma) * rho0;
    double xeps2=u/rho0;
    double xP2 = p_tmp;

    EOS_press(rho0,&xeps2,&xtemp2,ye0, &xP2);

#else
  p = pressure_rho0_w(rho0,w) ;
  u = w - (rho0 + p) ;
#endif

  // User may want to handle this case differently, e.g. do NOT return upon 
  // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:
  if( treat_floor_as_failure && ((rho0 <= 0.) || (u <= 0.)) ) { 
    retval = 5;
    return(retval) ;
  }
  prim[RHO] = rho0 ;
  prim[UU] = u ;

#if( TABULATED_EOS_TABLES )
  if (xtemp2<eos_tempmin){
   xtemp2=eos_tempmin;
  }


  prim[TEMP] = xtemp2;

#endif

  //-fast  for(i=1;i<4;i++)  Qtcon[i] = Qcon[i] + ncon[i] * Qdotn;

  //-fast  for(i=1;i<4;i++) prim[UTCON1+i-1] = gamma/(W+Bsq) * ( Qtcon[i] + QdotB*Bcon[i]/W ) ;
  double g_o_WBsq = gamma/(W+Bsq);
  double QdB_o_W  = QdotB / W; 
  *bsq = Bsq * (1.-vsq) + QdB_o_W*QdB_o_W;
  prim[UTCON1] = g_o_WBsq * ( Qcon[1] + geom->ncon[1] * Qdotn + QdB_o_W*Bcon[1] ) ;
  prim[UTCON2] = g_o_WBsq * ( Qcon[2] + geom->ncon[2] * Qdotn + QdB_o_W*Bcon[2] ) ;
  prim[UTCON3] = g_o_WBsq * ( Qcon[3] + geom->ncon[3] * Qdotn + QdB_o_W*Bcon[3] ) ;
	
  /* set field components */
  //-fast  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;


  /* done! */
  return(retval) ;

}


/**********************************************************************/ 
/****************************************************************************
   vsq_calc(): 
    
      -- evaluate v^2 (spatial, normalized velocity) from 
            W = \gamma^2 w 

****************************************************************************/
static FTYPE vsq_calc(FTYPE W)
{
	FTYPE Wsq,Xsq,Bsq_W;
	
	Wsq = W*W ;
	Bsq_W = (Bsq + W);
	Xsq = Bsq_W * Bsq_W;

	return(  ( Wsq * Qtsq  + QdotBsq * (Bsq_W + W)) / (Wsq*Xsq) );
}


/********************************************************************

  x1_of_x0(): 
           
    -- calculates v^2 from W  with some physical bounds checking;
    -- asumes x0 is already physical
    -- makes v^2 physical  if not;

*********************************************************************/

static FTYPE x1_of_x0(FTYPE x0 ) 
{
  FTYPE x1,vsq;
  FTYPE dv = 1.e-15;
  

  vsq = fabs(vsq_calc(x0)) ; // guaranteed to be positive 


  return( ( vsq > 1. ) ? (1.0 - dv) : vsq   ); 

}

/********************************************************************

  validate_x(): 
           
    -- makes sure that x[0,1] have physical values, based upon 
       their definitions:
    
*********************************************************************/

static void validate_x(FTYPE x[2], FTYPE x0[2] ) 
{
  
  FTYPE const dv = 1.e-15;

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
static int general_newton_raphson( FTYPE x[], int n, 
            FTYPE D, FTYPE ye0, FTYPE * temp_guess, 
			    void (*funcd) (FTYPE [], FTYPE [], FTYPE [], 
					   FTYPE [][NEWT_DIM], FTYPE *, 
					   FTYPE *,FTYPE , FTYPE, FTYPE *) )
{
  FTYPE f, df, dx[NEWT_DIM], x_old[NEWT_DIM];
  FTYPE resid[NEWT_DIM], jac[NEWT_DIM][NEWT_DIM];
  FTYPE errx, x_orig[NEWT_DIM];
  int    n_iter, id, jd, i_extra, doing_extra;
  FTYPE dW,dvsq,vsq_old,vsq,W,W_old;
  FTYPE const dv = (1.-1.e-15);

  int   keep_iterating;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;
  //-fast  for( id = 0; id < NEWT_DIM ; id++)  x_old[id] = x_orig[id] = x[id] ;
  x_old[0] = x_orig[0] = x[0] ;
  x_old[1] = x_orig[1] = x[1] ;

  vsq_old = vsq = W = W_old = 0.;
  n_iter = 0;

  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, &df, D,ye0,temp_guess);  /* returns with new dx, f, df */
      

    /* Save old values before calculating the new: */
    errx = 0.;
    //-fast    for( id = 0; id < NEWT_DIM ; id++) {   x_old[id] = x[id] ;   }
    x_old[0] = x[0] ; 
    x_old[1] = x[1] ; 

    /* Make the newton step: */
    //-fast    for( id = 0; id < NEWT_DIM ; id++) {    x[id] += dx[id]  ;    }
    x[0] += dx[0]  ; 
    x[1] += dx[1]  ; 

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    //--faster  validate_x( x, x_old ) ;
    if( x[0] < 0. ) {  x[0] = fabs(x[0]);  } 
    else { 
     if(x[0] > W_TOO_BIG)  { x[0] = x_old[0] ; }
    }

    if( x[1] < 0. ) {  x[1] = 0.; } 
    else { 
      if( x[1] > 1. ) { x[1] = dv; }
    }

    double gamma_nr, gamma_sq_nr;
    
     gamma_sq_nr = 1.0/(1. - x[1]);
     gamma_nr=sqrt(gamma_sq_nr);

     if (gamma_nr>GAMMAMAX){
      gamma_nr=GAMMAMAX;
      gamma_sq_nr=gamma_nr*gamma_nr;
      x[1] = 1.-1./gamma_sq_nr;
       }




      
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
  if( (!isfinite(f)) ||  (!isfinite(df)) ) {
    fprintf(stderr, "infinite \n");
    fprintf(stderr, "x[0] %e\n", x[0] );
    fprintf(stderr, "x[1] %e\n", x[1] );
    fprintf(stderr, "dx0 %e\n", dx[0]);
    fprintf(stderr, "dx1 %e\n", dx[1]);
    fprintf(stderr, "residuals %e %e\n", resid[0], resid[1]);
    return(2);
  }


  if( fabs(errx) <= NEWT_TOL ){
    return(0);
  }
  else if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    return(0);
  }
  else {
    fprintf(stderr,"error bad in 2deos safe_guess %e \n",errx);
    fprintf(stderr, "x[0] %e\n", x[0] );
    fprintf(stderr, "x[1] %e\n", x[1] );
    //fprintf(stderr, "error bad %e\n", errx);
    fprintf(stderr, "dx0 %e\n", dx[0]);
    fprintf(stderr, "dx1 %e\n", dx[1]);
    fprintf(stderr, "residuals %e %e\n", resid[0], resid[1]);
    fprintf(stderr,"original 0, 1 %26.16e %26.16e\n",x_orig[0],x_orig[1]);
    fprintf(stderr, "D,ye0,temp_guess %26.16e %26.16e %26.16e\n",D,ye0,*temp_guess);
//  //  fprintf(stdout,"newPrims  = %26.16e %26.16e %26.16e %26.16e %26.16e \n", rho0,u,p,w,W_last);
//  //  fprintf(stdout,"OTHERS    = %26.16e %26.16e %26.16e %26.16e %26.16e %26.16e \n", Qdotn,Qsq,D,half_Bsq,utsq,gamma); 
    return(1);
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

static void func_vsq(FTYPE x[], FTYPE dx[], FTYPE resid[], 
		     FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df,FTYPE D, FTYPE ye0, FTYPE * temp_guess)
{

  FTYPE gamma_sq, gamma;
  FTYPE  W, vsq, Wsq, p_tmp, dPdvsq, dPdW;
  FTYPE t11, t16,t18,t2,t21,   t23,   t24,   t25,   t3,   t35,   t36,   t4,   t40;


  W = x[0];
  vsq = x[1];
  
  Wsq = W*W;
  gamma_sq = 1.0/(1. - vsq);
  gamma=sqrt(gamma_sq);
  
  if (gamma>GAMMAMAX){
   gamma=GAMMAMAX;
   gamma_sq=gamma*gamma;
   vsq = 1.-1./gamma_sq;
   }


#if( TABULATED_EOS_TABLES )
  eos_quantities_general(W,vsq, &dPdW, &dPdvsq, D, ye0, temp_guess, &p_tmp);
#else
  p_tmp = eos_info(W, vsq, &dPdW, &dPdvsq);
#endif




  // These expressions were calculated using Mathematica, but made into efficient 
  // code using Maple.  Since we know the analytic form of the equations, we can 
  // explicitly calculate the Newton-Raphson step: 

  t2 = -half_Bsq+dPdvsq;
  t3 = Bsq+W;
  t4 = t3*t3;
  t23 = 1/W;
  t16 = QdotBsq*t23*t23;
  t11 = Qtsq-vsq*t4+t16*(Bsq+W+W);
  t18 = -Qdotn-half_Bsq*(1.0+vsq)+0.5*t16-W+p_tmp;
  t24 = t16*t23;
  t25 = -1.0+dPdW-t24;
  t35 = t25*t3+(Bsq-2.0*dPdvsq)*(t16+vsq*W)*t23;
  //  t21 = 1/t3;
  //  t36 = 1/t35;
  t21 = 1/(t3*t35);
  dx[0] = -(t2*t11+t4*t18)*t21;
  t40 = -2*(vsq+t24)*t3;
  dx[1] = -(-t25*t11+t40*t18)*t21;
  //  detJ = t3*t35;
  jac[0][0] = t40;
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
  eos_info(): 
 
      -- partial derivative of pressure with respect to vsq
 **********************************************************************/
static FTYPE eos_info(FTYPE W, FTYPE vsq, FTYPE *dpdw, FTYPE *dpdvsq)
{
  register double ftmp,gtmp;

  ftmp = 1. - vsq;
  gtmp = sqrt(ftmp);

  *dpdw =  GAMMA_M1_O_GAMMA * ftmp ;
  *dpdvsq =  GAMMA_M1_O_GAMMA * ( 0.5 * D/gtmp  -  W ) ;

  return( GAMMA_M1_O_GAMMA * ( W * ftmp  -  D * gtmp )  );  // p 

}


/**********************************************************************/
/********************************************************************** 
  eos_quantities_general():
  - provides p, dp/dW, dp/d(v^2) for given W, v^2, conservatives
  - performs an inversion from specific enthalpy to temperature and
    computes the other quantities
 **********************************************************************/
static double eos_quantities_general(FTYPE W, FTYPE vsq, FTYPE *dpdw, FTYPE *dpdvsq, 
         FTYPE D, FTYPE ye0, FTYPE * temp_guess, FTYPE *p_tmp)
{ 
  double gamma_sq = 1.0/(1. - vsq);
  double gamma = sqrt(gamma_sq);

   if (gamma>GAMMAMAX){
      gamma=GAMMAMAX;
      gamma_sq=GAMMAMAX*GAMMAMAX;
    }

  double dpdeps = 0.0;
  double dpdrho = 0.0;
  double p = 0.0;
  double eps = 0.0;
  double entr = 0.0;
  double rho = D /gamma;
  double enth = W/gamma_sq/rho; // specific enthalpy h, W = rho*h*gamma^2
    
    if(rho < eos_rhomin) {
      rho=eos_rhomin;
    }

    if (enth<0){
    enth=fabs(enth);
     }

  EOS_P_from_hrho_dPdrho_dPdeps(rho, enth, temp_guess, ye0, &eps, &p, &dpdrho, 
            &dpdeps, &entr); 

  *p_tmp=p;

  double dpdeps_o_rho= dpdeps/rho;
  *dpdw = ( dpdeps_o_rho/(1.+dpdeps_o_rho) )/gamma_sq;

  double dpdvsq_1 = -0.5*D*gamma*dpdrho;
  double dpdvsq_2 = -0.5*(W + p*gamma_sq)/rho;

  *dpdvsq = (dpdvsq_1 + dpdeps*dpdvsq_2)/(1+dpdeps_o_rho);

  return 0;
}


/****************************************************************************** 
             END   OF   UTOPRIM_2D_fast.C
 ******************************************************************************/



