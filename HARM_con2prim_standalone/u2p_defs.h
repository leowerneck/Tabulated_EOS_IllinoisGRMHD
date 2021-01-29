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

#include "decs.h"

 /* your choice of floating-point data type */
#define FTYPE double    

#define NPR  (NP)

/* Adiabatic index used for the state equation */
#define GAMMA             (gam)
#define GAMMA_M1_O_GAMMA  (gam_m1_o_gam)
#define G_ISOTHERMAL      (1.)
double eos_rhomax, eos_rhomin;
double eos_tempmax, eos_tempmin;
double eos_epsmax, eos_epsmin;

/* use K(s)=K(r)=const. (G_ATM = gam) of time or  T = T(r) = const. of time (G_ATM = 1.) */
#define USE_ISENTROPIC 1

#if( USE_ISENTROPIC ) 
#define G_ATM gam
#else
#define G_ATM G_ISOTHERMAL
#endif


#define MAX_NEWT_ITER 40     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
#define NEWT_TOL   3.5e-9    /* Min. of tolerance allowed for Newton-Raphson iterations */
#define MIN_NEWT_TOL  5.0e-10    /* Max. of tolerance allowed for Newton-Raphson iterations */
#define EXTRA_NEWT_ITER 2

#define NEWT_TOL2     1.0e-15      /* TOL of new 1D^*_{v^2} gnr2 method */
#define MIN_NEWT_TOL2 1.0e-10  /* TOL of new 1D^*_{v^2} gnr2 method */

#define W_TOO_BIG	1.e20	/* \gamma^2 (\rho_0 + u + p) is assumed
                                  to always be smaller than this.  This
				  is used to detect solver failures */
#define UTSQ_TOO_BIG	1.e20    /* \tilde{u}^2 is assumed to be smaller
                                  than this.  Used to detect solver
				  failures */

#define FAIL_VAL  1.e30    /* Generic value to which we set variables when a problem arises */

#define NUMEPSILON (2.2204460492503131e-16)


/* some mnemonics */
/* for primitive variables */
#define UTCON1 	(U1)
#define UTCON2 	(U2)
#define UTCON3 	(U3)
#define BCON1	(B1)
#define BCON2	(B2)
#define BCON3	(B3)

/* for conserved variables */
#define QCOV0	1
#define QCOV1	2
#define QCOV2	3
#define QCOV3	4

/**************************************************/
/**************************************************
  The following functions assume a Gamma-law EOS:
***************************************************/
/* 
pressure as a function of rho0 and u 
this is used by primtoU and Utoprim_?D
*/
#define pressure_rho0_u(__rho0, __u)  ( ((GAMMA - 1.)*__u)  )

/* 
pressure as a function of rho0 and w = rho0 + u + p 
this is used by primtoU and Utoprim_1D
*/
#define pressure_rho0_w(__rho0, __w)  ( (GAMMA_M1_O_GAMMA*(w - rho0)) )


