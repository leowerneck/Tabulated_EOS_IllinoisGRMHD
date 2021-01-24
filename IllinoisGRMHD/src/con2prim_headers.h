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
#ifndef HARM_PRIMITIVES_HEADERS_H_
#define HARM_PRIMITIVES_HEADERS_H_

// These are used to select the con2prim routine
static const int Noble2D              = 0;
static const int Palenzuela1D         = 1;
static const int Palenzuela1D_entropy = 2;

static const int NPR =8;
static const int NDIM=4;

/* Adiabatic index used for the state equation */
//#define GAMMA    (2.0)

static const CCTK_REAL G_ISOTHERMAL = 1.0;

/* use K(s)=K(r)=const. (G_ATM = GAMMA) of time or  T = T(r) = const. of time (G_ATM = 1.) */
/*
  #define USE_ISENTROPIC 1

  #if( USE_ISENTROPIC )
  #define G_ATM GAMMA
  #else
  #define G_ATM G_ISOTHERMAL
  #endif
*/

static const int MAX_NEWT_ITER=30;     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
//#define MAX_NEWT_ITER 300     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
static const CCTK_REAL NEWT_TOL    =1.0e-10;    /* Min. of tolerance allowed for Newton-Raphson iterations */
static const CCTK_REAL MIN_NEWT_TOL=1.0e-10;    /* Max. of tolerance allowed for Newton-Raphson iterations */
static const int EXTRA_NEWT_ITER=0; /* ZACH SAYS: Original value = 2. But I don't think this parameter > 0 is warranted. Just slows the code for no reason, since our tolerances are fine. */

static const CCTK_REAL NEWT_TOL2    =1.0e-15;      /* TOL of new 1D^*_{v^2} gnr2 method */
static const CCTK_REAL MIN_NEWT_TOL2=1.0e-10;  /* TOL of new 1D^*_{v^2} gnr2 method */

static const CCTK_REAL W_TOO_BIG    =1.e20;    /* \gamma^2 (\rho_0 + u + p) is assumed
                                                  to always be smaller than this.  This
                                                  is used to detect solver failures */
static const CCTK_REAL UTSQ_TOO_BIG =1.e20;    /* \tilde{u}^2 is assumed to be smaller
                                                  than this.  Used to detect solver
                                                  failures */

static const CCTK_REAL FAIL_VAL     =1.e30;    /* Generic value to which we set variables when a problem arises */

static const CCTK_REAL NUMEPSILON=(2.2204460492503131e-16);

/* some mnemonics */
/* for primitive variables */
static const int UU       =1;
static const int UTCON1   =2;
static const int UTCON2   =3;
static const int UTCON3   =4;
static const int BCON1    =5;
static const int BCON2    =6;
static const int BCON3    =7;

/* for conserved variables */
static const int QCOV0    =1;
static const int QCOV1    =2;
static const int QCOV2    =3;
static const int QCOV3    =4;

/* Also for primitive variables */
static const int RHO      =0;
static const int EPS      =1;
static const int v1_con   =2;
static const int v2_con   =3;
static const int v3_con   =4;
static const int B1_con   =5;
static const int B2_con   =6;
static const int B3_con   =7;
static const int YE       =8;
static const int TEMP     =9;
static const int PRESS    =10;
static const int WLORENTZ =11;
static const int ENT      =12;
static const int numprims =13; // rho, eps, v^{x,y,z}, B^{x,y,z}, Ye, T, P, W, S

/* Also for conservative variables */
static const int DD       =0;
static const int TAU      =1;
static const int S1_cov   =2;
static const int S2_cov   =3;
static const int S3_cov   =4;
static const int DS       =9;
static const int numcons  =10; // D, tau, S_{x,y,z}, B^{x,y,z}, DYe, DS

/********************************************************************************************/
// Function prototype declarations:

void con2prim_select( int* c2p_key,
                      int (*con2prim)( const igm_eos_parameters,
                                       const CCTK_REAL[4][4],const CCTK_REAL[4][4],
                                       const CCTK_REAL *restrict,CCTK_REAL *restrict ) );

int con2prim_Noble2D( const igm_eos_parameters eos,
                      const CCTK_REAL g4dn[4][4],
                      const CCTK_REAL g4up[4][4],
                      const CCTK_REAL *restrict cons,
                      CCTK_REAL *restrict prim );

int con2prim_Palenzuela1D( const igm_eos_parameters eos,
                           const CCTK_REAL g4dn[4][4],
                           const CCTK_REAL g4up[4][4],
                           const CCTK_REAL *restrict cons,
                           CCTK_REAL *restrict prim );

int con2prim_Palenzuela1D_entropy( const igm_eos_parameters eos,
                                   const CCTK_REAL g4dn[4][4],
                                   const CCTK_REAL g4up[4][4],
                                   const CCTK_REAL *restrict cons,
                                   CCTK_REAL *restrict prim );

inline int IllinoisGRMHD_conservative_to_primitive(const int index,const int i,const int j,const int k,CCTK_REAL *X,CCTK_REAL *Y,CCTK_REAL *Z,
                                                   CCTK_REAL *METRIC,CCTK_REAL *METRIC_PHYS,CCTK_REAL *METRIC_LAP_PSI4,
                                                   CCTK_REAL *CONSERVS,CCTK_REAL *PRIMS,
                                                   CCTK_REAL g4dn[NDIM][NDIM],CCTK_REAL g4up[NDIM][NDIM],
                                                   output_stats& stats,igm_eos_parameters& eos);

inline int font_fix__hybrid_EOS(CCTK_REAL &u_x, CCTK_REAL &u_y, CCTK_REAL &u_z,CCTK_REAL *CONSERVS,CCTK_REAL *PRIMS,CCTK_REAL *METRIC_PHYS,CCTK_REAL *METRIC_LAP_PSI4, igm_eos_parameters eos);
void eigenvalues_3by3_real_sym_matrix(CCTK_REAL & lam1, CCTK_REAL & lam2, CCTK_REAL & lam3,
                                      CCTK_REAL M11, CCTK_REAL M12, CCTK_REAL M13, CCTK_REAL M22, CCTK_REAL M23, CCTK_REAL M33);

void raise_g(CCTK_REAL vcov[NDIM], CCTK_REAL gcon[NDIM][NDIM], CCTK_REAL vcon[NDIM]);
void lower_g(CCTK_REAL vcon[NDIM], CCTK_REAL gcov[NDIM][NDIM], CCTK_REAL vcov[NDIM]);
void ncov_calc(CCTK_REAL gcon[NDIM][NDIM],CCTK_REAL ncov[NDIM]);
CCTK_REAL pressure_rho0_u(igm_eos_parameters eos, CCTK_REAL rho0, CCTK_REAL u);
CCTK_REAL pressure_rho0_w(igm_eos_parameters eos, CCTK_REAL rho0, CCTK_REAL w);

void set_cons_from_PRIMS_and_CONSERVS( const igm_eos_parameters eos,
                                       const CCTK_REAL *restrict METRIC,
                                       const CCTK_REAL *restrict METRIC_LAP_PSI4,
                                       const CCTK_REAL *restrict PRIMS,
                                       const CCTK_REAL *restrict CONSERVS,
                                       CCTK_REAL *restrict cons);

void set_prim_from_PRIMS_and_CONSERVS( const igm_eos_parameters eos,
                                       const int which_guess,
                                       const CCTK_REAL *restrict METRIC,
                                       const CCTK_REAL *restrict METRIC_LAP_PSI4,
                                       const CCTK_REAL *restrict PRIMS,
                                       const CCTK_REAL *restrict CONSERVS,
                                       const CCTK_REAL *restrict cons,
                                       CCTK_REAL *restrict prim );

/********************************************************************************************/

#endif /* HARM_PRIMITIVES_HEADERS_H_ */

