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

#include "u2p_defs.h"
#include "options.h"

#define NEWT_DIM 3

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
FTYPE Bsq,QdotBsq,Qtsq,Qdotn,QdotB,D,alpha,half_Bsq;

extern void EOS_press_with_T(double rho, double * eps,
        double * temp, double ye, double * prs);

extern void EOS_P_from_hrho_dPdrho_dPdeps(const double rho, const double enth,
          double * temp_guess, double ye, double * eps, double * press, double * dPdrho,
          double * dPdeps, double * entr);

extern void EP_dEdgamma_dEdW_dEdT(double * Eprim, double * Pprim, double * dEdgamma, double * dEdW, double * dEdT, 
       double * dpdrho, double * dpdT, double * x, const double D, double ye);

extern void EOS_press(double rho, double * eps,
        double * temp, double ye, double * prs);

// Declarations: 
static int Utoprim_new_body(FTYPE *U, struct of_geom *geom, FTYPE *prim, FTYPE *gamma, FTYPE *bsq);
static int general_newton_raphson3d( FTYPE x[], int n, 
            FTYPE D, FTYPE ye,
				   void (*funcd) (FTYPE [], FTYPE [], 
						  FTYPE [], FTYPE [][NEWT_DIM], 
						  FTYPE *, FTYPE *, FTYPE , FTYPE));

static  void func_vsq( FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], 
		       FTYPE *f, FTYPE *df, FTYPE D, FTYPE ye);
static FTYPE x1_of_x0(FTYPE x0 ) ;
static FTYPE vsq_calc(FTYPE W);
static FTYPE pressure_W_vsq(FTYPE W, FTYPE vsq) ;
static FTYPE dpdW_calc_vsq(FTYPE W, FTYPE vsq);
static FTYPE dpdvsq_calc(FTYPE W, FTYPE vsq);
static  FTYPE eos_quantities_general(FTYPE W, FTYPE vsq, FTYPE *dpdw, FTYPE *dpdvsq,
         FTYPE D, FTYPE ye, FTYPE * temp_guess,FTYPE *p_tmp);


/**********************************************************************/
/******************************************************************

  Utoprim_2d():
  
  -- Driver for new prim. var. solver.  The driver just translates
     between the two sets of definitions for U and P.  The user may 
     wish to alter the translation as they see fit.  


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
int Utoprim_3d_eos(FTYPE *U, struct of_geom *geom, FTYPE *prim, FTYPE *gamma, FTYPE *bsq)
{


  //fprintf(stderr, "first U[RHO]/prim[RHO]%e\n", U[RHO]/prim[RHO]);
  FTYPE U_tmp[NPR], prim_tmp[NPR];
  int i, j, ret; 
  FTYPE geomfactor,inv_gdet;


  /* Set the geometry variables: */
  inv_gdet = geom->g_inv; 
  alpha = geom->alpha;
  geomfactor = alpha * inv_gdet;

  //fprintf(stderr, "temp_guess%15.6E\n",temp_guess);
  //fprintf(stderr, "g_inv%15.6E\n",inv_gdet);

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
  prim_tmp[TEMP  ] = prim[TEMP  ];

 // fprintf(stderr, "prim_tmp[TEMP  ] %e\n", prim_tmp[TEMP  ]);


  ret = Utoprim_new_body(U_tmp, geom, prim_tmp, gamma, bsq);

  /* Use new primitive variables if there was no problem : */ 
  if( ret == 0 ) {
    prim[RHO   ] = prim_tmp[RHO   ];
    prim[UU    ] = prim_tmp[UU    ];
    prim[UTCON1] = prim_tmp[UTCON1];
    prim[UTCON2] = prim_tmp[UTCON2];
    prim[UTCON3] = prim_tmp[UTCON3];
    prim[YE    ] = prim_tmp[YE    ];
    prim[TEMP  ] = prim_tmp[TEMP  ];
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


  FTYPE x_3d[NEWT_DIM];
  FTYPE Bcon[NDIM],*Qcov,Qcon[NDIM],Qsq;
  FTYPE rho0,u,p,w,gammasq,gamma,gtmp,W_last,W,utsq,vsq, ye0, Temp;
  FTYPE rhomax, umax, epsmax, tempmax;
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

  if(ye0 >1.) {
    fprintf(stderr," in 3d eos ye0 %e \n", ye0);
    fprintf(stderr," U[rho] %e \n", U[RHO]);
    fprintf(stderr," U[YE] %e \n", U[YE]);

    ye0=1;
  }


  half_Bsq = 0.5*Bsq;
  //fprintf(stderr, "half_Bsq %e\n",half_Bsq );

  /* calculate W from last timestep and use for guess */
  //-fast   utsq = 0. ;
  //-fast   for(i=1;i<4;i++)
  //-fast     for(j=1;j<4;j++) utsq += gcov[i][j]*prim[UTCON1+i-1]*prim[UTCON1+j-1] ;

  utsq =  geom->gcov[1][1]*prim[U1]*prim[U1]
        + geom->gcov[2][2]*prim[U2]*prim[U2]
        + geom->gcov[3][3]*prim[U3]*prim[U3]
     + 2*(geom->gcov[1][2]*prim[U1]*prim[U2]
	+ geom->gcov[1][3]*prim[U1]*prim[U3]
	+ geom->gcov[2][3]*prim[U2]*prim[U3]) ;

  if( (utsq < 0.) && (fabs(utsq) < 1.0e-13) ) { 
    utsq = fabs(utsq);
  }
  if(utsq < 0. || utsq > UTSQ_TOO_BIG) {
    retval = 2;
    return(retval) ;
  }

  gammasq = 1. + utsq ;
  gamma  = sqrt(gammasq);
	

  rho0 = D/gamma;
  u=prim[UU];//epsmax*rhomax;
  //fprintf(stderr, "eos_epsmin %e\n", eos_epsmin);
  //epsmax =umax/rhomax; //(-Qdotn-D-Bsq/2.)/D  ;//
  //umax=epsmax*rhomax;

  double temp_guess;

  temp_guess=prim[TEMP];
  tempmax=prim[TEMP];

#if( TABULATED_EOS_TABLES )
  if(rho0 < eos_rhomin) {
    rho0=eos_rhomin;
  }   
  double xeps = u/rho0;
  double xtemp = prim[TEMP];
  double xP = 0.;
  EOS_press_with_T(rho0,&xeps,&xtemp,ye0, &xP);
  //EOS_press(rho0,&xeps,&xtemp,ye0, &xP);
  p = xP;

#else 
 
  p = pressure_rho0_u(rho0,u) ;

#endif

  w = rho0 + u + p ;

  W_last = gamma*gamma*(rho0+u+p);//-Qdotn+p-Bsq/2. ;

  // Make sure that W is large enough so that v^2 < 1 : 
  i_increase = 0;
  while( (( W_last*W_last*W_last * ( W_last + 2.*Bsq ) 
	    - QdotBsq*(2.*W_last + Bsq) ) <= W_last*W_last*(Qtsq-Bsq*Bsq))
	 && (i_increase < 10) ) {
    W_last *= 10.;
    i_increase++;
  }


  
  // Calculate gamma,W and T: 
  x_3d[0] = gamma; //Lorentz factor
  x_3d[1] = fabs( W_last );
  x_3d[2] =  tempmax;

  //fprintf(stderr, "x_3d %26.16e %26.16e   %26.16e \n", x_3d[0],x_3d[1],x_3d[2] ); 

  //temp_guess=1.;
  // printf("3D NR in : %.15e %.15e %.15e\n",x_3d[0],x_3d[1],x_3d[2]);
  retval = general_newton_raphson3d( x_3d, NEWT_DIM,  D, ye0,func_vsq) ;
  // printf("3D NR out: %.15e %.15e %.15e (retval = %d)\n",x_3d[0],x_3d[1],x_3d[2],retval);
 
  //fprintf(stderr, "x_3d sfter nr %26.16e %26.16e   %26.16e \n", x_3d[0],x_3d[1],x_3d[2] );

  gamma = x_3d[0];
  W = x_3d[1];
  Temp=x_3d[2];

  vsq=x1_of_x0(W);//1.-1./(gamma*gamma);



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

  // Recover the primitive variables from the scalars and conserved variables:

  rho0 = D / gamma;

  gammasq=gamma*gamma;


  w = W * (1. - vsq) ;


#if( TABULATED_EOS_TABLES )
   FTYPE p_tmp,dPdvsq, dPdW;;
    double xtemp2;
    xtemp2=Temp;
    eos_quantities_general(W,vsq, &dPdW, &dPdvsq, D, ye0, &xtemp2, &p_tmp);
    u=(W-D*gamma-p_tmp*gamma*gamma)/(D*gamma) * rho0;
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


  prim[TEMP] = xtemp2;


  //
  //

 // prim[TEMP] = Temp;
  
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

  //fprintf(stderr, "prim_inside[UU]%e\n", prim[UU    ]);
  /* done! */
  return(retval) ;

}



/************************************************************

  general_newton_raphson(): 

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

*****************************************************************/
static int general_newton_raphson3d( FTYPE x[], int n, 
          FTYPE D, FTYPE ye,
			    void (*funcd) (FTYPE [], FTYPE [], FTYPE [], 
					   FTYPE [][NEWT_DIM], FTYPE *, 
					   FTYPE *, FTYPE, FTYPE) )
{


  FTYPE f, df, dx[NEWT_DIM], x_old[NEWT_DIM];
  FTYPE resid[NEWT_DIM], jac[NEWT_DIM][NEWT_DIM];
  FTYPE errx, erry, errz, errmax, x_orig[NEWT_DIM];
  int    n_iter, id, jd, i_extra, doing_extra;
  FTYPE dW,dgamma,gamma_old,gamma,W,W_old, dT, T, T_old;
  FTYPE const dv = (1.-1.e-15);

  int   keep_iterating;
  
  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;

  //-fast  for( id = 0; id < NEWT_DIM ; id++)  x_old[id] = x_orig[id] = x[id] ;
  x_old[0] = x_orig[0] = x[0] ;
  x_old[1] = x_orig[1] = x[1] ;
  x_old[2] = x_orig[2] = x[2] ;

  //fprintf(stderr, "x_old[0] %e\n", x_old[0]);
  //fprintf(stderr, "x_old[1] %e\n", x_old[1]);
  //fprintf(stderr, "x_old[2] %e\n", x_old[2]);

  //gamma_old = gamma = W = W_old =T=T_old=0.;
  n_iter = 0;

  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) {

    // printf("NR iter: %02d: in   x_{i} = %23.15e %23.15e %23.15e\n",n_iter+1,x[0],x[1],x[2]);

    (*funcd) (x, dx, resid, jac, &f, &df, D,ye);  /* returns with new dx, f, df */

    // printf("NR iter: %02d: out dx_{i} = %23.15e %23.15e %23.15e\n",n_iter+1,dx[0],dx[1],dx[2]);
      

    /* Save old values before calculating the new: */
    errx = 0.;
    erry = 0.;
    errz = 0.;
    errmax = 0.;
    //-fast    for( id = 0; id < NEWT_DIM ; id++) {   x_old[id] = x[id] ;   }
    x_old[0] = x[0] ; 
    x_old[1] = x[1] ; 
    x_old[2] = x[2] ;


    /* Make the newton step: */
    //-fast    for( id = 0; id < NEWT_DIM ; id++) {    x[id] += dx[id]  ;    }
    x[0] += dx[0]  ; 
    x[1] += dx[1]  ; 
    x[2] += dx[2]  ; 

    //fprintf(stderr, "x[0] %e\n", x[0] );
    //fprintf(stderr, "x[1] %e\n", x[1] );
    //fprintf(stderr, "x[2] %e\n", x[2] );



    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);
    erry = (x[1]==0.) ?  fabs(dx[1]) : fabs(dx[1]/x[1]);
    errz = (x[2]==0.) ?  fabs(dx[2]) : fabs(dx[2]/x[2]);
    //fprintf(stderr, "errx %e\n", errx);
    //fprintf(stderr, "fabs(dx[0]) %e\n", fabs(dx[0]));
    //fprintf(stderr, "x[0] %e\n", x[0]);

    errmax=MAX(errx,MAX(erry,errz));
    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    //--faster  validate_x( x, x_old ) ;

    if( x[0] < 1. ) {  x[0] = 1.; } 

    if ( x[0] > GAMMAMAX ) {  x[0] = x_old[0]; }

    if( x[1] < 0. ) {  x[1] = fabs(x[1]);  } 
    else { 
      if(x[1] > W_TOO_BIG)  { x[1] = x_old[1] ; }
    }

    if( x[2] < eos_tempmin ) {  x[2] = x_old[2]; }
    if( x[2] > eos_tempmax ) {  x[2] = x_old[2]; } 

    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    if( (fabs(errmax) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }
    if( doing_extra == 1 ) i_extra++ ;
    if( ((fabs(errmax) <= NEWT_TOL)&&(doing_extra == 0)) 
	|| (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }
    n_iter++;

  }   // END of while(keep_iterating)

  /*  Check for bad untrapped divergences : */
  if( (!isfinite(f)) ||  (!isfinite(df)) ) {
    return(2);
  }

  if( fabs(errmax) <= NEWT_TOL ){
    return(0);
  }
  else if( (fabs(errmax) <= MIN_NEWT_TOL) && (fabs(errmax) > NEWT_TOL) ){
    return(0);
  }
  else {
    fprintf(stderr,"error bad in 3deos%e \n",errmax);
    fprintf(stderr, "x[0] %e\n", x[0] );
    fprintf(stderr, "x[1] %e\n", x[1] );
    fprintf(stderr, "x[2] %e\n", x[2] );
    //fprintf(stderr, "error bad %e\n", errx);
    fprintf(stderr, "dx0 %e\n", dx[0]);
    fprintf(stderr, "dx1 %e\n", dx[1]);
    fprintf(stderr, "dx2 %e\n", dx[2]);
    fprintf(stderr, "residuals %e %e %e \n", resid[0], resid[1], resid[2]);
    fprintf(stderr,"original 0, 1 2 %26.16e %26.16e   %26.16e \n",x_orig[0],x_orig[1], x_orig[2]);
    fprintf(stderr, "D,ye0,temp_guess %26.16e %26.16e %26.16e\n",D,ye,temp_guess);

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
		     FTYPE J[][NEWT_DIM], FTYPE *f, FTYPE *df, 
         FTYPE D, FTYPE ye)
{

  
  FTYPE  W, gamma, T, vsq, Wsq, p_tmp,eps_tmp, temp, tmp2,tmp3;

  FTYPE E,E_EOS,P_EOS,dEdgamma,dEdW,dEdT,dpEOSdrho,dpEOSdT;
  gamma = x[0];
  W = x[1];
  T= x[2];
  vsq=x1_of_x0(W);//1.-1./(gamma*gamma);

  
  Wsq = W*W;


  //fprintf(stderr, "before EP_dEdgamma_dEdW_dEdT\n");

  // printf("EOS call: %.15e %.15e %.15e %.15e %.15e\n",D/gamma,ye,T,W,gamma);

  EP_dEdgamma_dEdW_dEdT(&E_EOS,&P_EOS,&dEdgamma,&dEdW,&dEdT,&dpEOSdrho,&dpEOSdT,x,D,ye);

  // printf("EOS out : %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",P_EOS,E_EOS,dpEOSdrho,dpEOSdT,dEdgamma,dEdW,dEdT);

  //fprintf(stderr, "after EP_dEdgamma_dEdW_dEdT\n");

  // Eq. (25) in Cerda-Duran et al. 2008
  E = -1.0 + W/(D*gamma) - P_EOS*gamma/D;
  E=fmax(E,eos_epsmin);

  // printf("E = %.15e\n",E);
  
  // dequation2/dgamma 
  J[0][0] = 2.0*((W + Bsq)*(W + Bsq)-Qtsq-(2.0*W + Bsq)*((QdotB)*(QdotB))/(W*W))*gamma;
  double a = J[0][0];
  
  // dequation2/dW 
  J[0][1] = (2.0*(W + Bsq) + (2.0/(W*W) + 2.0*Bsq/(W*W*W))*(QdotB*QdotB))*gamma*gamma - 2.0*(W + Bsq);
  double b = J[0][1];

  // dequation2/dT 
  J[0][2] = 0;
  double c = J[0][2];

  // dequation1/dgamma 
  J[1][0] = 2.0*(-Qdotn - W - Bsq + (QdotB*QdotB)/(2.0*W*W) + P_EOS)*gamma - D*dpEOSdrho;
  double d = J[1][0];
  
  // dequation1/dW
  J[1][1] = (-1.0-(QdotB*QdotB)/(W*W*W))*gamma*gamma;
  double e = J[1][1];

// dequation1/dT 
  J[1][2] = gamma*gamma*dpEOSdT;
  double fc = J[1][2];
  
  // d/dgamma (E-E(rho,T,Y_e))
  J[2][0] = dEdgamma;
  double g = J[2][0];
  
  // d/dW (E-E(rho,T,Y_e))
  J[2][1] = dEdW;
  double h = J[2][1];
  
  // d/dT (E-E(rho,T,Y_e))
  J[2][2] = dEdT;
  double k = J[2][2];

  // printf("a,b,c,d,e,fc,g,h,k = %.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n",a,b,c,d,e,fc,g,h,k);


 // compute f(x) from (21) (equation1), (22)(equation2) and (28)(equation3) in Cerda-Duran et al. 2008

  // printf("%.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",Qdotn,W,Bsq,QdotB,P_EOS,gamma,Bsq);
  
  resid[0] = ((W+Bsq)*(W+Bsq) - Qtsq - (2.0*W +Bsq)*(QdotB*QdotB)/(W*W))*gamma*gamma - ((W + Bsq)*(W + Bsq));
  resid[1] = (-Qdotn - W - Bsq + (QdotB*QdotB)/(2.0*W*W) + P_EOS)*gamma*gamma + 0.5*Bsq;
  resid[2] = E - E_EOS;

  // printf("resid[0] = %.15e | resid[1] = %.15e | resid[2] = %.15e\n",resid[0],resid[1],resid[2]);
  //fprintf(stderr, "E %e\n", E);
  //fprintf(stderr, "E_EOS %e\n", E_EOS);

  // Compute the determinant
  double A = e*k-fc*h;
  double B = fc*g-d*k;
  double C = d*h-e*g;
  double detJ = a*(A) + b*(B) + c*(C);


  //Compute the matrix inverse
  double Ji[NEWT_DIM][NEWT_DIM];
  Ji[0][0] = A/detJ;
  Ji[1][0] = B/detJ;
  Ji[2][0] = C/detJ;
  Ji[0][1] = (c*h-b*k)/detJ;
  Ji[1][1] = (a*k-c*g)/detJ;
  Ji[2][1] = (g*b-a*h)/detJ;
  Ji[0][2] = (b*fc-c*e)/detJ;
  Ji[1][2] = (c*d-a*fc)/detJ;
  Ji[2][2] = (a*e-b*d)/detJ;

  // printf("Ji = \n%23.15e %23.15e %23.15e\n%23.15e %23.15e %23.15e\n%23.15e %23.15e %23.15e\n",
         // J[0][0],J[0][1],J[0][2],
         // J[1][0],J[1][1],J[1][2],
         // J[2][0],J[2][1],J[2][2]);

  // Compute the step size
  dx[0] = (-Ji[0][0]* resid[0] - Ji[0][1]* resid[1] - Ji[0][2]* resid[2]);
  dx[1] = (-Ji[1][0]* resid[0] - Ji[1][1]* resid[1] - Ji[1][2]* resid[2]);
  dx[2] = (-Ji[2][0]* resid[0] - Ji[2][1]* resid[1] - Ji[2][2]* resid[2]);

  // printf("dx_{i} = %23.15e %23.15e %23.15e\n",dx[0],dx[1],dx[2]);


  *df = -resid[0]*resid[0] - resid[1]*resid[1]-resid[2]*resid[2];

  *f = -0.5 * ( *df );

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


/**********************************************************************/
/********************************************************************** 
 *   eos_quantities_general():
 *     - provides p, dp/dW, dp/d(v^2) for given W, v^2, conservatives
 *       - works for 3-parameter EOS
 *         - performs an inversion from specific enthalpy to temperature and
 *             computes the other quantities
 *              **********************************************************************/
static double eos_quantities_general(FTYPE W, FTYPE vsq, FTYPE *dpdw, FTYPE *dpdvsq,
         FTYPE D, FTYPE ye, FTYPE * temp_guess, FTYPE *p_tmp)
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
  double rhow = D /gamma;
  double enth = W/gamma_sq/rhow; // specific enthalpy h, W = rho*h*gamma^2
if(rhow < eos_rhomin) {
      rhow=eos_rhomin;
    }

    if (enth<0){
    enth=fabs(enth);
     }

  EOS_P_from_hrho_dPdrho_dPdeps(rhow, enth, temp_guess, ye, &eps, &p, &dpdrho,
            &dpdeps, &entr);
  *p_tmp=p;




  double dpdeps_o_rho= dpdeps/rhow;
  *dpdw = ( dpdeps_o_rho/(1.+dpdeps_o_rho) )/gamma_sq;

  double dpdvsq_1 = -0.5*D*gamma*dpdrho;
  double dpdvsq_2 = -0.5*(W + p*gamma_sq)/rhow;

  *dpdvsq = (dpdvsq_1 + dpdeps*dpdvsq_2)/(1+dpdeps_o_rho);


  return 0;


}
/****************************************************************************** 
             END   OF   UTOPRIM_3D.C
 ******************************************************************************/
