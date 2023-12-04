/*********************************************
 * Evaluate RHS of GRMHD & induction equations
 * (vector potential prescription), using the
 * generalized Lorenz gauge condition for the
 * EM gauge.
 *
 * Based originally on the Illinois GRMHD code,
 * written by Matt Duez, Yuk Tung Liu, and Branson
 * Stephens (original version), and then developed
 * primarily by Zachariah Etienne, Yuk Tung Liu,
 * and Vasileios Paschalidis.
 *
 * Rewritten for public release in 2013
 *      by Zachariah B. Etienne
 *
 * References:
 * Original unigrid GRMHD evolution prescription:
 *    http://arxiv.org/abs/astro-ph/0503420
 * Vector potential formulation in full GR:
 *    http://arxiv.org/abs/1007.2848
 * Improved EM gauge conditions for AMR grids:
 *    http://arxiv.org/abs/1110.4633
 * Generalized Lorenz gauge prescription:
 *    http://arxiv.org/abs/1207.3354
 *
 * Note that the Generalized Lorenz gauge strength
 *  parameter has units of 1/M, just like the \eta
 *  parameter in the gamma-driving shift condition,
 *  so setting it too large will result in violation
 *  of the CFL condition.
 *
 * This version of PPM implements the standard
 * Colella & Woodward PPM, though modified as in GRHydro
 * to have 3 ghostzones instead of 4.
 *********************************************/


#include "cctk.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "IllinoisGRMHD_headers.h" /* Generic #define's and function prototypes */
#include "driver_evaluate_MHD_rhs.h" /* Function prototypes for this file only */
#include "inlined_functions.h"
#include "GRHayLib.h"

#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

// The inner two points of the interpolation function use
// the value of A_in, and the outer two points use A_out.
#define A_out -0.0625
#define A_in  0.5625
//Interpolates to the -1/2 face of point Var
#define COMPUTE_FCVAL(Varm2,Varm1,Var,Varp1) (A_out*(Varm2) + A_in*(Varm1) + A_in*(Var) + A_out*(Varp1))

// Computes 4th-order derivative
#define B_out -1.0/12.0
#define B_in  2.0/3.0
#define COMPUTE_DERIV(Varm2,Varm1,Varp1,Varp2) (B_in*(Varp1 - Varm1) + B_out*(Varp2 - Varm2))

static void GRHayLMHD_tabulated_entropy_evaluate_sources_rhs(CCTK_ARGUMENTS);
static void GRHayLMHD_tabulated_entropy_calculate_flux_dir_rhs(
      const cGH *restrict cctkGH,
      const int flux_dir,
      const CCTK_REAL **B_center,
      const CCTK_REAL *restrict B_stagger,
      CCTK_REAL **vel_r,
      CCTK_REAL **vel_l,
      CCTK_REAL *restrict cmin,
      CCTK_REAL *restrict cmax);
static void GRHayLMHD_evaluate_psi6phi_and_A_gauge_rhs(CCTK_ARGUMENTS);
static void GRHayLMHD_reconstruction_loop(
      const cGH *restrict cctkGH,
      const int flux_dir, const int num_vars,
      const int *restrict var_indices,
      const double *rho_b,
      const double *pressure,
      const double *v_flux,
      const double **in_prims,
      double **out_prims_r,
      double **out_prims_l);
static void GRHayLMHD_A_flux_rhs(
      const cGH *restrict cctkGH,
      const int flux_dir,
      /*const*/ double **in_prims_r,
      /*const*/ double **in_prims_l,
      /*const*/ double **cmin,
      /*const*/ double **cmax,
      double *restrict A_rhs);
static void GRHayLMHD_compute_metric_derivs(
      const cGH *cctkGH,
      const int i, const int j, const int k,
      const int flux_dir,
      const CCTK_REAL dxi,
      const CCTK_REAL *restrict lapse,
      const CCTK_REAL *restrict betax,
      const CCTK_REAL *restrict betay,
      const CCTK_REAL *restrict betaz,
      const CCTK_REAL *restrict gxx,
      const CCTK_REAL *restrict gxy,
      const CCTK_REAL *restrict gxz,
      const CCTK_REAL *restrict gyy,
      const CCTK_REAL *restrict gyz,
      const CCTK_REAL *restrict gzz,
      ghl_metric_quantities *restrict metric_derivs);

static void GRHayLMHD_interpolate_metric_to_face(
      const cGH *cctkGH,
      const int i, const int j, const int k,
      const int flux_dir,
      const CCTK_REAL *restrict lapse,
      const CCTK_REAL *restrict betax,
      const CCTK_REAL *restrict betay,
      const CCTK_REAL *restrict betaz,
      const CCTK_REAL *restrict gxx,
      const CCTK_REAL *restrict gxy,
      const CCTK_REAL *restrict gxz,
      const CCTK_REAL *restrict gyy,
      const CCTK_REAL *restrict gyz,
      const CCTK_REAL *restrict gzz,
      ghl_metric_quantities *restrict metric);

static double get_Gamma_eff_tabulated(
      const double rho_in,
      const double press_in) {
  return 1.0;
}

extern "C" void IllinoisGRMHD_driver_evaluate_MHD_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  int levelnumber = GetRefinementLevel(cctkGH);

  if(CCTK_Equals(verbose, "essential+iteration output")) {
    CCTK_VInfo(CCTK_THORNSTRING,"***** Iter. # %d, Lev: %d, Integrating to time: %e *****",cctk_iteration,levelnumber,cctk_delta_time/cctk_levfac[0]+cctk_time);
  }

  if( sizeof(CCTK_REAL) < 8 ) CCTK_VError(VERR_DEF_PARAMS,"Error: IllinoisGRMHD assumes that CCTK_REAL is a double precision number. Setting otherwise will likely cause havoc with the conserv_to_prims solver.");

  if(cctk_nghostzones[0]<3 || cctk_nghostzones[1]<3 || cctk_nghostzones[2]<3) { CCTK_VError(VERR_DEF_PARAMS,"ERROR. Need at least 3 ghostzones for IllinoisGRMHD evolutions."); }

  CCTK_REAL dX[3] = { CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1), CCTK_DELTA_SPACE(2) };


  /**********************************
   * Piecewise Polytropic EOS Patch *
   *   Setting up the EOS struct    *
   **********************************/
  /*
   * The short piece of code below takes care
   * of initializing the EOS parameters.
   * Please refer to the "inlined_functions.h"
   * source file for the documentation on the
   * function.
   */
  igm_eos_parameters eos;
  initialize_igm_eos_parameters_from_input(igm_eos_key,cctk_time,eos);


  // in_prims,out_prims_r, and out_prims_l are arrays of pointers to the actual gridfunctions.
  gf_and_gz_struct in_prims[MAXNUMVARS],out_prims_r[MAXNUMVARS],out_prims_l[MAXNUMVARS];
  int which_prims_to_reconstruct[MAXNUMVARS],num_prims_to_reconstruct;

  /* SET POINTERS TO GRMHD GRIDFUNCTIONS */
  // The order here MATTERS, and must be consistent with the global variable declarations in
  //   evaluate_MHD_rhs_headers.h (look for RHOB=0, etc.)
  //   For example, in_prims[0] _must_ be rho_b.
  in_prims[RHOB       ].gf=rho_b;           out_prims_r[RHOB       ].gf=rho_br;      out_prims_l[RHOB       ].gf=rho_bl;
  in_prims[PRESSURE   ].gf=P;               out_prims_r[PRESSURE   ].gf=Pr;          out_prims_l[PRESSURE   ].gf=Pl;
  in_prims[VX         ].gf=vx;              out_prims_r[VX         ].gf=vxr;         out_prims_l[VX         ].gf=vxl;
  in_prims[VY         ].gf=vy;              out_prims_r[VY         ].gf=vyr;         out_prims_l[VY         ].gf=vyl;
  in_prims[VZ         ].gf=vz;              out_prims_r[VZ         ].gf=vzr;         out_prims_l[VZ         ].gf=vzl;
  in_prims[BX_CENTER  ].gf=Bx;              out_prims_r[BX_CENTER  ].gf=Bxr;         out_prims_l[BX_CENTER  ].gf=Bxl;
  in_prims[BY_CENTER  ].gf=By;              out_prims_r[BY_CENTER  ].gf=Byr;         out_prims_l[BY_CENTER  ].gf=Byl;
  in_prims[BZ_CENTER  ].gf=Bz;              out_prims_r[BZ_CENTER  ].gf=Bzr;         out_prims_l[BZ_CENTER  ].gf=Bzl;
  in_prims[BX_STAGGER ].gf=Bx_stagger;      out_prims_r[BX_STAGGER ].gf=Bx_staggerr; out_prims_l[BX_STAGGER ].gf=Bx_staggerl;
  in_prims[BY_STAGGER ].gf=By_stagger;      out_prims_r[BY_STAGGER ].gf=By_staggerr; out_prims_l[BY_STAGGER ].gf=By_staggerl;
  in_prims[BZ_STAGGER ].gf=Bz_stagger;      out_prims_r[BZ_STAGGER ].gf=Bz_staggerr; out_prims_l[BZ_STAGGER ].gf=Bz_staggerl;
  in_prims[VXR        ].gf=vxr;             out_prims_r[VXR        ].gf=vxrr;        out_prims_l[VXR        ].gf=vxrl;
  in_prims[VYR        ].gf=vyr;             out_prims_r[VYR        ].gf=vyrr;        out_prims_l[VYR        ].gf=vyrl;
  in_prims[VZR        ].gf=vzr;             out_prims_r[VZR        ].gf=vzrr;        out_prims_l[VZR        ].gf=vzrl;
  in_prims[VXL        ].gf=vxl;             out_prims_r[VXL        ].gf=vxlr;        out_prims_l[VXL        ].gf=vxll;
  in_prims[VYL        ].gf=vyl;             out_prims_r[VYL        ].gf=vylr;        out_prims_l[VYL        ].gf=vyll;
  in_prims[VZL        ].gf=vzl;             out_prims_r[VZL        ].gf=vzlr;        out_prims_l[VZL        ].gf=vzll;
  in_prims[YEPRIM     ].gf=igm_Ye;          out_prims_r[YEPRIM     ].gf=Yer;         out_prims_l[YEPRIM     ].gf=Yel;
  in_prims[TEMPERATURE].gf=igm_temperature; out_prims_r[TEMPERATURE].gf=Tr;          out_prims_l[TEMPERATURE].gf=Tl;
  in_prims[EPSILON    ].gf=igm_eps;         out_prims_r[EPSILON    ].gf=epsr;        out_prims_l[EPSILON    ].gf=epsl;
  in_prims[ENTROPY    ].gf=igm_entropy;     out_prims_r[ENTROPY    ].gf=Sr;          out_prims_l[ENTROPY    ].gf=Sl;

  // Prims are defined AT ALL GRIDPOINTS, so we set the # of ghostzones to zero:
  for(int i=0;i<MAXNUMVARS;i++) for(int j=1;j<=3;j++) { in_prims[i].gz_lo[j]=0; in_prims[i].gz_hi[j]=0; }
  // Left/right variables are not yet defined, yet we set the # of gz's to zero by default:
  for(int i=0;i<MAXNUMVARS;i++) for(int j=1;j<=3;j++) { out_prims_r[i].gz_lo[j]=0; out_prims_r[i].gz_hi[j]=0; }
  for(int i=0;i<MAXNUMVARS;i++) for(int j=1;j<=3;j++) { out_prims_l[i].gz_lo[j]=0; out_prims_l[i].gz_hi[j]=0; }


  // Convert ADM variables (from ADMBase) to the BSSN-based variables expected by this routine.
  IllinoisGRMHD_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij(cctkGH,cctk_lsh,  gxx,gxy,gxz,gyy,gyz,gzz,alp,
                                                                              gtxx,gtxy,gtxz,gtyy,gtyz,gtzz,
                                                                              gtupxx,gtupxy,gtupxz,gtupyy,gtupyz,gtupzz,
                                                                              phi_bssn,psi_bssn,lapm1);

  /* SET POINTERS TO METRIC GRIDFUNCTIONS */
  CCTK_REAL *metric[NUMVARS_FOR_METRIC_FACEVALS]; // "metric" here is array of pointers to the actual gridfunctions.
  int ww=0;
  metric[ww]=phi_bssn;ww++;
  metric[ww]=psi_bssn;ww++;
  metric[ww]=gtxx;    ww++;
  metric[ww]=gtxy;    ww++;
  metric[ww]=gtxz;    ww++;
  metric[ww]=gtyy;    ww++;
  metric[ww]=gtyz;    ww++;
  metric[ww]=gtzz;    ww++;
  metric[ww]=lapm1;   ww++;
  metric[ww]=betax;   ww++;
  metric[ww]=betay;   ww++;
  metric[ww]=betaz;   ww++;
  metric[ww]=gtupxx;  ww++;
  metric[ww]=gtupyy;  ww++;
  metric[ww]=gtupzz;  ww++;

  /* SET POINTERS TO STRESS-ENERGY TENSOR GRIDFUNCTIONS */
  CCTK_REAL *TUPmunu[10];// "TUPmunu" here is array of pointers to the actual gridfunctions.
  ww=0;
  TUPmunu[ww]=TUPtt; ww++;
  TUPmunu[ww]=TUPtx; ww++;
  TUPmunu[ww]=TUPty; ww++;
  TUPmunu[ww]=TUPtz; ww++;
  TUPmunu[ww]=TUPxx; ww++;
  TUPmunu[ww]=TUPxy; ww++;
  TUPmunu[ww]=TUPxz; ww++;
  TUPmunu[ww]=TUPyy; ww++;
  TUPmunu[ww]=TUPyz; ww++;
  TUPmunu[ww]=TUPzz; ww++;


  // 1) First initialize {rho_star_rhs,tau_rhs,st_x_rhs,st_y_rhs,st_z_rhs} to zero
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
        int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
        Ax_rhs      [index] = 0.0;
        Ay_rhs      [index] = 0.0;
        Az_rhs      [index] = 0.0;
        psi6phi_rhs [index] = 0.0;

        rho_star_rhs[index] = 0.0;
        st_x_rhs    [index] = 0.0;
        st_y_rhs    [index] = 0.0;
        st_z_rhs    [index] = 0.0;
        tau_rhs     [index] = 0.0;

        Ye_star_rhs [index] = 0.0;
        S_star_rhs  [index] = 0.0;

	s_tau       [index] = 0.0;
	s_sx        [index] = 0.0;
	s_sy        [index] = 0.0;
	s_sz        [index] = 0.0;
      }

  // Here, we:
  // 1) Compute tau_rhs extrinsic curvature terms, and
  // 2) Compute TUPmunu.
  // This function is housed in the file: "compute_tau_rhs_extrinsic_curvature_terms_and_TUPmunu.C"
  //compute_tau_rhs_extrinsic_curvature_terms_and_TUPmunu(eos,cctkGH,cctk_lsh,cctk_nghostzones,dX,metric,in_prims,TUPmunu,
  //                                                      gtupxy,gtupxz,gtupyz,
  //                                                      kxx,kxy,kxz,kyy,kyz,kzz,
  //                                                      tau_rhs, s_tau);
  GRHayLMHD_tabulated_entropy_evaluate_sources_rhs(cctkGH);

  int flux_dirn;
  flux_dirn=1;
  // First compute ftilde, which is used for flattening left and right face values
  // This function is housed in the file: "reconstruct_set_of_prims_PPM.C"
  ftilde_gf_compute(cctkGH,cctk_lsh,flux_dirn, in_prims, ftilde_gf);



  /* There are two stories going on here:
   * 1) Computation of \partial_x on RHS of \partial_t {rho_star,tau,mhd_st_{x,y,z}},
   *    via PPM reconstruction onto (i-1/2,j,k), so that
   *    \partial_x F = [ F(i+1/2,j,k) - F(i-1/2,j,k) ] / dx
   * 2) Computation of \partial_t A_i, where A_i are *staggered* gridfunctions,
   *    where A_x is defined at (i,j+1/2,k+1/2), A_y at (i+1/2,j,k+1/2), etc.
   *    Ai_rhs = \partial_t A_i = \epsilon_{ijk} \psi^{6} v^j B^k,
   *    where \epsilon_{ijk} is the flat-space antisymmetric operator.
   * 2A) Az_rhs is defined at (i+1/2,j+1/2,k), and it depends on {Bx,By,vx,vy},
   *     so the trick is to reconstruct {Bx,By,vx,vy} cleverly to get to these
   *     staggered points. For example:
   * 2Aa) vx and vy are at (i,j,k), and we reconstruct them to (i-1/2,j,k) below. After
   *      this, we'll reconstruct again in the y-dir'n to get {vx,vy} at (i-1/2,j-1/2,k)
   * 2Ab) By_stagger is at (i,j+1/2,k), and we reconstruct below to (i-1/2,j+1/2,k). */
  ww=0;
  which_prims_to_reconstruct  [ww]=RHOB;                      ww++;
  which_prims_to_reconstruct  [ww]=eos.PPM_reconstructed_var; ww++;
  if( eos.is_Tabulated ) {
    which_prims_to_reconstruct[ww]=YEPRIM;                    ww++;
    // which_prims_to_reconstruct[ww]=TEMPERATURE;               ww++;
  }
  which_prims_to_reconstruct  [ww]=VX;                        ww++;
  which_prims_to_reconstruct  [ww]=VY;                        ww++;
  which_prims_to_reconstruct  [ww]=VZ;                        ww++;
  which_prims_to_reconstruct  [ww]=BY_CENTER;                 ww++;
  which_prims_to_reconstruct  [ww]=BZ_CENTER;                 ww++;
  which_prims_to_reconstruct  [ww]=BY_STAGGER;                ww++;
  num_prims_to_reconstruct=ww;
  // This function is housed in the file: "reconstruct_set_of_prims_PPM.C"
  reconstruct_set_of_prims_PPM(eos,cctkGH,cctk_lsh,flux_dirn,num_prims_to_reconstruct,which_prims_to_reconstruct,
                               in_prims,out_prims_r,out_prims_l,ftilde_gf,temporary);

  if( eos.is_Tabulated ) {
    // check_temperature_reconstruction(eos,cctkGH,cctk_lsh,in_prims,out_prims_r,out_prims_l);
    compute_remaining_prims_on_right_and_left_face(eos,cctkGH,cctk_lsh,in_prims,out_prims_r,out_prims_l);
  }

  //Right and left face values of BI_CENTER are used in mhdflux computation (first to compute b^a).
  //   Instead of reconstructing, we simply set B^x face values to be consistent with BX_STAGGER.
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
        int index=CCTK_GFINDEX3D(cctkGH,i,j,k), indexim1=CCTK_GFINDEX3D(cctkGH,i-1+(i==0),j,k); /* indexim1=0 when i=0 */
        out_prims_r[BX_CENTER].gf[index]=out_prims_l[BX_CENTER].gf[index]=in_prims[BX_STAGGER].gf[indexim1]; }
  // Then add fluxes to RHS for hydro variables {rho_b,P,vx,vy,vz}:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_hydro_rhss.C"
  add_fluxes_and_source_terms_to_hydro_rhss(eos,flux_dirn,cctkGH,cctk_lsh,cctk_nghostzones,dX,
                                            metric,TUPmunu,
                                            num_prims_to_reconstruct,
                                            in_prims,out_prims_r,out_prims_l,
                                            cmax_x,cmin_x,
                                            rho_star_flux,tau_flux,st_x_flux,st_y_flux,st_z_flux,Ye_star_flux,S_star_flux,
                                            rho_star_rhs,tau_rhs,st_x_rhs,st_y_rhs,st_z_rhs,Ye_star_rhs,S_star_rhs,
                                            s_tau, s_sx, s_sy, s_sz);

  // Note that we have already reconstructed vx and vy along the x-direction,
  //   at (i-1/2,j,k). That result is stored in v{x,y}{r,l}.  Bx_stagger data
  //   are defined at (i+1/2,j,k).
  // Next goal: reconstruct Bx, vx and vy at (i+1/2,j+1/2,k).
  flux_dirn=2;
  // First compute ftilde, which is used for flattening left and right face values
  // This function is housed in the file: "reconstruct_set_of_prims_PPM.C"
  ftilde_gf_compute(cctkGH,cctk_lsh,flux_dirn, in_prims, ftilde_gf);

  // in_prims[{VXR,VXL,VYR,VYL}].gz_{lo,hi} ghostzones are set to all zeros, which
  //    is incorrect. We fix this below.
  // [Note that this is a cheap operation, copying only 8 integers and a pointer.]
  in_prims[VXR]=out_prims_r[VX];
  in_prims[VXL]=out_prims_l[VX];
  in_prims[VYR]=out_prims_r[VY];
  in_prims[VYL]=out_prims_l[VY];

  /* There are two stories going on here:
   * 1) Computation of \partial_y on RHS of \partial_t {rho_star,tau,mhd_st_{x,y,z}},
   *    via PPM reconstruction onto (i,j-1/2,k), so that
   *    \partial_y F = [ F(i,j+1/2,k) - F(i,j-1/2,k) ] / dy
   * 2) Computation of \partial_t A_i, where A_i are *staggered* gridfunctions,
   *    where A_x is defined at (i,j+1/2,k+1/2), A_y at (i+1/2,j,k+1/2), etc.
   *    Ai_rhs = \partial_t A_i = \epsilon_{ijk} \psi^{6} v^j B^k,
   *    where \epsilon_{ijk} is the flat-space antisymmetric operator.
   * 2A) Az_rhs is defined at (i+1/2,j+1/2,k), and it depends on {Bx,By,vx,vy},
   *     so the trick is to reconstruct {Bx,By,vx,vy} cleverly to get to these
   *     staggered points. For example:
   * 2Aa) VXR = [right-face of vx reconstructed along x-direction above] is at (i-1/2,j,k),
   *      and we reconstruct it to (i-1/2,j-1/2,k) below. Similarly for {VXL,VYR,VYL}
   * 2Ab) Bx_stagger is at (i+1/2,j,k), and we reconstruct to (i+1/2,j-1/2,k) below
   * 2Ac) By_stagger is at (i-1/2,j+1/2,k) already for Az_rhs, from the previous step.
   * 2B) Ax_rhs is defined at (i,j+1/2,k+1/2), and it depends on {By,Bz,vy,vz}.
   *     Again the trick is to reconstruct these onto these staggered points.
   * 2Ba) Bz_stagger is at (i,j,k+1/2), and we reconstruct to (i,j-1/2,k+1/2) below */
  ww=0;
  // NOTE! The order of variable reconstruction is important here,
  //   as we don't want to overwrite {vxr,vxl,vyr,vyl}!
  which_prims_to_reconstruct[ww]=VXR;                       ww++;
  which_prims_to_reconstruct[ww]=VYR;                       ww++;
  which_prims_to_reconstruct[ww]=VXL;                       ww++;
  which_prims_to_reconstruct[ww]=VYL;                       ww++;
  num_prims_to_reconstruct=ww;
  // This function is housed in the file: "reconstruct_set_of_prims_PPM.C"
  reconstruct_set_of_prims_PPM(eos,cctkGH,cctk_lsh,flux_dirn,num_prims_to_reconstruct,which_prims_to_reconstruct,
                               in_prims,out_prims_r,out_prims_l,ftilde_gf,temporary);
  ww=0;
  // Reconstruct other primitives last!
  which_prims_to_reconstruct  [ww]=RHOB;                      ww++;
  which_prims_to_reconstruct  [ww]=eos.PPM_reconstructed_var; ww++;
  if( eos.is_Tabulated ) {
    which_prims_to_reconstruct[ww]=YEPRIM;                    ww++;
    // which_prims_to_reconstruct[ww]=TEMPERATURE;               ww++;
  }
  which_prims_to_reconstruct  [ww]=VX;                        ww++;
  which_prims_to_reconstruct  [ww]=VY;                        ww++;
  which_prims_to_reconstruct  [ww]=VZ;                        ww++;
  which_prims_to_reconstruct  [ww]=BX_CENTER;                 ww++;
  which_prims_to_reconstruct  [ww]=BZ_CENTER;                 ww++;
  which_prims_to_reconstruct  [ww]=BX_STAGGER;                ww++;
  which_prims_to_reconstruct  [ww]=BZ_STAGGER;                ww++;
  num_prims_to_reconstruct=ww;
  // This function is housed in the file: "reconstruct_set_of_prims_PPM.C"
  reconstruct_set_of_prims_PPM(eos,cctkGH,cctk_lsh,flux_dirn,num_prims_to_reconstruct,which_prims_to_reconstruct,
                               in_prims,out_prims_r,out_prims_l,ftilde_gf,temporary);

  if( eos.is_Tabulated ) {
    // check_temperature_reconstruction(eos,cctkGH,cctk_lsh,in_prims,out_prims_r,out_prims_l);
    compute_remaining_prims_on_right_and_left_face(eos,cctkGH,cctk_lsh,in_prims,out_prims_r,out_prims_l);
  }

  //Right and left face values of BI_CENTER are used in mhdflux computation (first to compute b^a).
  //   Instead of reconstructing, we simply set B^y face values to be consistent with BY_STAGGER.
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
        int index=CCTK_GFINDEX3D(cctkGH,i,j,k), indexjm1=CCTK_GFINDEX3D(cctkGH,i,j-1+(j==0),k); /* indexjm1=0 when j=0 */
        out_prims_r[BY_CENTER].gf[index]=out_prims_l[BY_CENTER].gf[index]=in_prims[BY_STAGGER].gf[indexjm1]; }
  // Then add fluxes to RHS for hydro variables {rho_b,P,vx,vy,vz}:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_hydro_rhss.C"
  add_fluxes_and_source_terms_to_hydro_rhss(eos,flux_dirn,cctkGH,cctk_lsh,cctk_nghostzones,dX,
                                            metric,TUPmunu,
                                            num_prims_to_reconstruct,
                                            in_prims,out_prims_r,out_prims_l,
                                            cmax_y,cmin_y,
                                            rho_star_flux,tau_flux,st_x_flux,st_y_flux,st_z_flux,Ye_star_flux,S_star_flux,
                                            rho_star_rhs,tau_rhs,st_x_rhs,st_y_rhs,st_z_rhs,Ye_star_rhs,S_star_rhs,
                                            s_tau, s_sx, s_sy, s_sz);

  /*****************************************
   * COMPUTING RHS OF A_z, BOOKKEEPING NOTE:
   * We want to compute
   * \partial_t A_z - [gauge terms] = \psi^{6} (v^x B^y - v^y B^x).
   * A_z is defined at (i+1/2,j+1/2,k).
   * ==========================
   * Where defined  | Variables
   * (i-1/2,j-1/2,k)| {vxrr,vxrl,vxlr,vxll,vyrr,vyrl,vylr,vyll}
   * (i+1/2,j-1/2,k)| {Bx_stagger_r,Bx_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i-1/2,j+1/2,k)| {By_stagger_r,By_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j,k)        | {phi}
   * ==========================
   ******************************************/
  // Interpolates to i+1/2
#define IPH(METRICm1,METRICp0,METRICp1,METRICp2) (-0.0625*((METRICm1) + (METRICp2)) + 0.5625*((METRICp0) + (METRICp1)))
  // Next compute phi at (i+1/2,j+1/2,k):
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=1;j<cctk_lsh[1]-2;j++) for(int i=1;i<cctk_lsh[0]-2;i++) {
        temporary[CCTK_GFINDEX3D(cctkGH,i,j,k)]=
          IPH(IPH(phi_bssn[CCTK_GFINDEX3D(cctkGH,i-1,j-1,k)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j-1,k)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i+1,j-1,k)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i+2,j-1,k)]),
              IPH(phi_bssn[CCTK_GFINDEX3D(cctkGH,i-1,j  ,k)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j  ,k)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i+1,j  ,k)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i+2,j  ,k)]),
              IPH(phi_bssn[CCTK_GFINDEX3D(cctkGH,i-1,j+1,k)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j+1,k)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i+1,j+1,k)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i+2,j+1,k)]),
              IPH(phi_bssn[CCTK_GFINDEX3D(cctkGH,i-1,j+2,k)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j+2,k)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i+1,j+2,k)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i+2,j+2,k)]));
      }


  int A_directionz=3;
  A_i_rhs_no_gauge_terms(A_directionz,cctkGH,cctk_lsh,cctk_nghostzones,out_prims_r,out_prims_l,temporary,cmax_x,cmin_x,cmax_y,cmin_y, Az_rhs);


  // in_prims[{VYR,VYL,VZR,VZL}].gz_{lo,hi} ghostzones are not correct, so we fix
  //    this below.
  // [Note that this is a cheap operation, copying only 8 integers and a pointer.]
  in_prims[VYR]=out_prims_r[VY];
  in_prims[VYL]=out_prims_l[VY];
  in_prims[VZR]=out_prims_r[VZ];
  in_prims[VZL]=out_prims_l[VZ];

  flux_dirn=3;
  // First compute ftilde, which is used for flattening left and right face values
  // This function is housed in the file: "reconstruct_set_of_prims_PPM.C"
  ftilde_gf_compute(cctkGH,cctk_lsh,flux_dirn, in_prims, ftilde_gf);

  /* There are two stories going on here:
   * 1) Single reconstruction to (i,j,k-1/2) for {rho,P,vx,vy,vz,Bx,By,Bz} to compute
   *    z-dir'n advection terms in \partial_t {rho_star,tau,mhd_st_{x,y,z}} at (i,j,k)
   * 2) Multiple reconstructions for *staggered* gridfunctions A_i:
   *    Ai_rhs = \partial_t A_i = \epsilon_{ijk} \psi^{6} v^j B^k,
   *    where \epsilon_{ijk} is the flat-space antisymmetric operator.
   * 2A) Ax_rhs is defined at (i,j+1/2,k+1/2), depends on v{y,z} and B{y,z}
   * 2Aa) v{y,z}{r,l} are at (i,j-1/2,k), so we reconstruct here to (i,j-1/2,k-1/2)
   * 2Ab) Bz_stagger{r,l} are at (i,j-1/2,k+1/2) already.
   * 2Ac) By_stagger is at (i,j+1/2,k), and below we reconstruct its value at (i,j+1/2,k-1/2)
   * 2B) Ay_rhs is defined at (i+1/2,j,k+1/2), depends on v{z,x} and B{z,x}.
   * 2Ba) v{x,z} are reconstructed to (i,j,k-1/2). Later we'll reconstruct again to (i-1/2,j,k-1/2).
   * 2Bb) Bz_stagger is at (i,j,k+1/2). Later we will reconstruct to (i-1/2,j,k+1/2).
   * 2Bc) Bx_stagger is at (i+1/2,j,k), and below we reconstruct its value at (i+1/2,j,k-1/2)
   */
  ww=0;
  // NOTE! The order of variable reconstruction is important here,
  //   as we don't want to overwrite {vxr,vxl,vyr,vyl}!
  which_prims_to_reconstruct  [ww]=VYR;                       ww++;
  which_prims_to_reconstruct  [ww]=VZR;                       ww++;
  which_prims_to_reconstruct  [ww]=VYL;                       ww++;
  which_prims_to_reconstruct  [ww]=VZL;                       ww++;
  num_prims_to_reconstruct=ww;
  // This function is housed in the file: "reconstruct_set_of_prims_PPM.C"
  reconstruct_set_of_prims_PPM(eos,cctkGH,cctk_lsh,flux_dirn,num_prims_to_reconstruct,which_prims_to_reconstruct,
                               in_prims,out_prims_r,out_prims_l,ftilde_gf,temporary);
  // Reconstruct other primitives last!
  ww=0;
  which_prims_to_reconstruct  [ww]=RHOB;                      ww++;
  which_prims_to_reconstruct  [ww]=eos.PPM_reconstructed_var; ww++;
  if( eos.is_Tabulated ) {
    which_prims_to_reconstruct[ww]=YEPRIM;                    ww++;
    // which_prims_to_reconstruct[ww]=TEMPERATURE;               ww++;
  }
  which_prims_to_reconstruct  [ww]=VX;                        ww++;
  which_prims_to_reconstruct  [ww]=VY;                        ww++;
  which_prims_to_reconstruct  [ww]=VZ;                        ww++;
  which_prims_to_reconstruct  [ww]=BX_CENTER;                 ww++;
  which_prims_to_reconstruct  [ww]=BY_CENTER;                 ww++;
  which_prims_to_reconstruct  [ww]=BX_STAGGER;                ww++;
  which_prims_to_reconstruct  [ww]=BY_STAGGER;                ww++;
  num_prims_to_reconstruct=ww;
  // This function is housed in the file: "reconstruct_set_of_prims_PPM.C"
  reconstruct_set_of_prims_PPM(eos,cctkGH,cctk_lsh,flux_dirn,num_prims_to_reconstruct,which_prims_to_reconstruct,
                               in_prims,out_prims_r,out_prims_l,ftilde_gf,temporary);

  if( eos.is_Tabulated ) {
    // check_temperature_reconstruction(eos,cctkGH,cctk_lsh,in_prims,out_prims_r,out_prims_l);
    compute_remaining_prims_on_right_and_left_face(eos,cctkGH,cctk_lsh,in_prims,out_prims_r,out_prims_l);
  }

  //Right and left face values of BI_CENTER are used in mhdflux computation (first to compute b^a).
  //   Instead of reconstructing, we simply set B^z face values to be consistent with BZ_STAGGER.
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
        int index=CCTK_GFINDEX3D(cctkGH,i,j,k), indexkm1=CCTK_GFINDEX3D(cctkGH,i,j,k-1+(k==0)); /* indexkm1=0 when k=0 */
        out_prims_r[BZ_CENTER].gf[index]=out_prims_l[BZ_CENTER].gf[index]=in_prims[BZ_STAGGER].gf[indexkm1]; }

  // Then add fluxes to RHS for hydro variables {rho_b,P,vx,vy,vz}:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_hydro_rhss.C"
  add_fluxes_and_source_terms_to_hydro_rhss(eos,flux_dirn,cctkGH,cctk_lsh,cctk_nghostzones,dX,
                                            metric,TUPmunu,
                                            num_prims_to_reconstruct,
                                            in_prims,out_prims_r,out_prims_l,
                                            cmax_z,cmin_z,
                                            rho_star_flux,tau_flux,st_x_flux,st_y_flux,st_z_flux,Ye_star_flux,S_star_flux,
                                            rho_star_rhs,tau_rhs,st_x_rhs,st_y_rhs,st_z_rhs,Ye_star_rhs,S_star_rhs,
                                            s_tau, s_sx, s_sy, s_sz);

  // in_prims[{VYR,VYL,VZR,VZL}].gz_{lo,hi} ghostzones are not set correcty.
  //    We fix this below.
  // [Note that this is a cheap operation, copying only 8 integers and a pointer.]
  in_prims[VXR]=out_prims_r[VX];
  in_prims[VZR]=out_prims_r[VZ];
  in_prims[VXL]=out_prims_l[VX];
  in_prims[VZL]=out_prims_l[VZ];
  // FIXME: lines above seem to be inconsistent with lines below.... Possible bug, not major enough to affect evolutions though.
  in_prims[VZR].gz_lo[1]=in_prims[VZR].gz_hi[1]=0;
  in_prims[VXR].gz_lo[1]=in_prims[VXR].gz_hi[1]=0;
  in_prims[VZL].gz_lo[1]=in_prims[VZL].gz_hi[1]=0;
  in_prims[VXL].gz_lo[1]=in_prims[VXL].gz_hi[1]=0;

  /*****************************************
   * COMPUTING RHS OF A_x, BOOKKEEPING NOTE:
   * We want to compute
   * \partial_t A_x - [gauge terms] = \psi^{6} (v^y B^z - v^z B^y).
   * A_x is defined at (i,j+1/2,k+1/2).
   * ==========================
   * Where defined  | Variables
   * (i,j-1/2,k-1/2)| {vyrr,vyrl,vylr,vyll,vzrr,vzrl,vzlr,vzll}
   * (i,j+1/2,k-1/2)| {By_stagger_r,By_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j-1/2,k+1/2)| {Bz_stagger_r,Bz_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j,k)        | {phi}
   * ==========================
   ******************************************/
  // Next compute phi at (i,j+1/2,k+1/2):
#pragma omp parallel for
  for(int k=1;k<cctk_lsh[2]-2;k++) for(int j=1;j<cctk_lsh[1]-2;j++) for(int i=0;i<cctk_lsh[0];i++) {
        temporary[CCTK_GFINDEX3D(cctkGH,i,j,k)]=
          IPH(IPH(phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j-1,k-1)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j,k-1)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j+1,k-1)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j+2,k-1)]),
              IPH(phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j-1,k  )],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j,k  )],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j+1,k  )],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j+2,k  )]),
              IPH(phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j-1,k+1)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j,k+1)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j+1,k+1)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j+2,k+1)]),
              IPH(phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j-1,k+2)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j,k+2)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j+1,k+2)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j+2,k+2)]));
      }


  int A_directionx=1;
  A_i_rhs_no_gauge_terms(A_directionx,cctkGH,cctk_lsh,cctk_nghostzones,out_prims_r,out_prims_l,temporary,cmax_y,cmin_y,cmax_z,cmin_z, Ax_rhs);


  // We reprise flux_dirn=1 to finish up computations of Ai_rhs's!
  flux_dirn=1;
  // First compute ftilde, which is used for flattening left and right face values
  // This function is housed in the file: "reconstruct_set_of_prims_PPM.C"
  ftilde_gf_compute(cctkGH,cctk_lsh,flux_dirn, in_prims, ftilde_gf);

  ww=0;
  // NOTE! The order of variable reconstruction is important here,
  //   as we don't want to overwrite {vxr,vxl,vyr,vyl}!
  which_prims_to_reconstruct[ww]=VXR;       ww++;
  which_prims_to_reconstruct[ww]=VZR;       ww++;
  which_prims_to_reconstruct[ww]=VXL;       ww++;
  which_prims_to_reconstruct[ww]=VZL;       ww++;
  which_prims_to_reconstruct[ww]=BZ_STAGGER;ww++;
  num_prims_to_reconstruct=ww;
  // This function is housed in the file: "reconstruct_set_of_prims_PPM.C"
  reconstruct_set_of_prims_PPM(eos,cctkGH,cctk_lsh,flux_dirn,num_prims_to_reconstruct,which_prims_to_reconstruct,
                               in_prims,out_prims_r,out_prims_l,ftilde_gf,temporary);

  /*****************************************
   * COMPUTING RHS OF A_y, BOOKKEEPING NOTE:
   * We want to compute
   * \partial_t A_y - [gauge terms] = \psi^{6} (v^z B^x - v^x B^z).
   * A_y is defined at (i+1/2,j,k+1/2).
   * ==========================
   * Where defined  | Variables
   * (i-1/2,j,k-1/2)| {vxrr,vxrl,vxlr,vxll,vzrr,vzrl,vzlr,vzll}
   * (i+1/2,j,k-1/2)| {Bx_stagger_r,Bx_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i-1/2,j,k+1/2)| {Bz_stagger_r,Bz_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j,k)        | {phi}
   * ==========================
   ******************************************/
  // Next compute phi at (i+1/2,j,k+1/2):
#pragma omp parallel for
  for(int k=1;k<cctk_lsh[2]-2;k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=1;i<cctk_lsh[0]-2;i++) {
        temporary[CCTK_GFINDEX3D(cctkGH,i,j,k)]=
          IPH(IPH(phi_bssn[CCTK_GFINDEX3D(cctkGH,i-1,j,k-1)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j,k-1)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i+1,j,k-1)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i+2,j,k-1)]),
              IPH(phi_bssn[CCTK_GFINDEX3D(cctkGH,i-1,j,k  )],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j,k  )],phi_bssn[CCTK_GFINDEX3D(cctkGH,i+1,j,k  )],phi_bssn[CCTK_GFINDEX3D(cctkGH,i+2,j,k  )]),
              IPH(phi_bssn[CCTK_GFINDEX3D(cctkGH,i-1,j,k+1)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j,k+1)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i+1,j,k+1)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i+2,j,k+1)]),
              IPH(phi_bssn[CCTK_GFINDEX3D(cctkGH,i-1,j,k+2)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i,j,k+2)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i+1,j,k+2)],phi_bssn[CCTK_GFINDEX3D(cctkGH,i+2,j,k+2)]));
      }


  int A_directiony=2;
  A_i_rhs_no_gauge_terms(A_directiony,cctkGH,cctk_lsh,cctk_nghostzones,out_prims_r,out_prims_l,temporary,cmax_z,cmin_z,cmax_x,cmin_x, Ay_rhs);


  // Next compute psi6phi_rhs, and add gauge terms to A_i_rhs terms!
  //   Note that in the following function, we don't bother with reconstruction, instead interpolating.
  // We need A^i, but only have A_i. So we add gtupij to the list of input variables.
  CCTK_REAL *interp_vars[MAXNUMINTERP];
  ww=0;
  interp_vars[ww]=betax;   ww++;
  interp_vars[ww]=betay;   ww++;
  interp_vars[ww]=betaz;   ww++;
  interp_vars[ww]=gtupxx;  ww++;
  interp_vars[ww]=gtupxy;  ww++;
  interp_vars[ww]=gtupxz;  ww++;
  interp_vars[ww]=gtupyy;  ww++;
  interp_vars[ww]=gtupyz;  ww++;
  interp_vars[ww]=gtupzz;  ww++;
  interp_vars[ww]=psi_bssn;ww++;
  interp_vars[ww]=lapm1;   ww++;
  interp_vars[ww]=Ax;      ww++;
  interp_vars[ww]=Ay;      ww++;
  interp_vars[ww]=Az;      ww++;
  int max_num_interp_variables=ww;
  if(max_num_interp_variables>MAXNUMINTERP) {CCTK_VError(VERR_DEF_PARAMS,"Error: Didn't allocate enough space for interp_vars[]."); }
  // We are FINISHED with v{x,y,z}{r,l} and P{r,l} so we use these 8 gridfunctions' worth of space as temp storage.
  Lorenz_psi6phi_rhs__add_gauge_terms_to_A_i_rhs(cctkGH,cctk_lsh,cctk_nghostzones,dX,interp_vars,psi6phi,
                                                 vxr,vyr,vzr,vxl,vyl,vzl,Pr,Pl,
                                                 psi6phi_rhs,Ax_rhs,Ay_rhs,Az_rhs);

  if( CCTK_IsThornActive("NRPyLeakageET") ) {
    // Convert rho, Y_e, T, and velocities to HydroBase
    // because they are needed by NRPyLeakage
#pragma omp parallel for
    for(int k=0;k<cctk_lsh[2];k++) {
      for(int j=0;j<cctk_lsh[1];j++) {
        for(int i=0;i<cctk_lsh[0];i++) {
          const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

          // Read from main memory
          const CCTK_REAL invalpL      = 1.0/alp[index];
          const CCTK_REAL betaxL       = betax[index];
          const CCTK_REAL betayL       = betay[index];
          const CCTK_REAL betazL       = betaz[index];
          const CCTK_REAL rhoL         = rho_b[index];
          const CCTK_REAL Y_eL         = igm_Ye[index];
          const CCTK_REAL temperatureL = igm_temperature[index];
          const CCTK_REAL vxL          = vx[index];
          const CCTK_REAL vyL          = vy[index];
          const CCTK_REAL vzL          = vz[index];

          // Write to main memory, converting to HydroBase
          rho[index]         = rhoL;
          Y_e[index]         = Y_eL;
          temperature[index] = temperatureL;
          velx[index]        = (vxL + betaxL)*invalpL;
          vely[index]        = (vyL + betayL)*invalpL;
          velz[index]        = (vzL + betazL)*invalpL;
        }
      }
    }
  }

  return;
  /*
  // FUN DEBUGGING TOOL (trust me!):
  #pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
  int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
  //st_x_rhs[index]=0.0;
  //st_y_rhs[index]=0.0;
  //st_z_rhs[index]=0.0;
  //rho_star_rhs[index]=0.0;
  //tau_rhs[index]=0.0;

  psi6phi_rhs[index] = 0.0;
  Ax_rhs[index] = 0.0;
  Ay_rhs[index] = 0.0;
  Az_rhs[index] = 0.0;
  }
  */
}

// We add #include's here instead of compiling these separately to help ensure that functions are properly inlined.
//    These files only include about 800 lines of code in total (~1200 lines in total), but it's arguably more
//    convenient to edit a 600 line file than an 1800 line file, so I'd prefer to leave this unconventional structure
//    alone.
#include "reconstruct_set_of_prims_PPM.C"
#include "compute_tau_rhs_extrinsic_curvature_terms_and_TUPmunu.C"
#include "add_fluxes_and_source_terms_to_hydro_rhss.C"
#include "mhdflux.C"
#include "A_i_rhs_no_gauge_terms.C"
#include "Lorenz_psi6phi_rhs__add_gauge_terms_to_A_i_rhs.C"
static void GRHayLMHD_tabulated_entropy_evaluate_sources_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int imin = cctk_nghostzones[0];
  const int jmin = cctk_nghostzones[1];
  const int kmin = cctk_nghostzones[2];
  const int imax = cctk_lsh[0] - cctk_nghostzones[0];
  const int jmax = cctk_lsh[1] - cctk_nghostzones[1];
  const int kmax = cctk_lsh[2] - cctk_nghostzones[2];
  const CCTK_REAL dxi = 1.0/CCTK_DELTA_SPACE(0);
  const CCTK_REAL dyi = 1.0/CCTK_DELTA_SPACE(1);
  const CCTK_REAL dzi = 1.0/CCTK_DELTA_SPACE(2);

#pragma omp parallel for
  for(int k=kmin; k<kmax; k++) {
    for(int j=jmin; j<jmax; j++) {
      for(int i=imin; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j ,k);

        // These variables have no source terms
        rho_star_rhs[index] = 0.0;
        S_star_rhs[index] = 0.0;
        Ye_star_rhs[index]  = 0.0;

        ghl_metric_quantities ADM_metric;
        ghl_initialize_metric(
              alp[index],
              betax[index], betay[index], betaz[index],
              gxx[index], gxy[index], gxz[index],
              gyy[index], gyz[index], gzz[index],
              &ADM_metric);

        ghl_extrinsic_curvature curv;
        ghl_initialize_extrinsic_curvature(
              kxx[index], kxy[index], kxz[index],
              kyy[index], kyz[index], kzz[index],
              &curv);

        ghl_primitive_quantities prims;
        prims.rho         = rho_b[index];
        prims.press       = P[index];
        prims.vU[0]       = vx[index];
        prims.vU[1]       = vy[index];
        prims.vU[2]       = vz[index];
        prims.BU[0]       = Bx[index] * ONE_OVER_SQRT_4PI;
        prims.BU[1]       = By[index] * ONE_OVER_SQRT_4PI;
        prims.BU[2]       = Bz[index] * ONE_OVER_SQRT_4PI;
        prims.entropy     = igm_entropy[index];
        prims.Y_e         = igm_Ye[index];
        prims.temperature = igm_temperature[index];
        const int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_limit_v_and_compute_u0(ghl_params, &ADM_metric, &prims);

        ghl_metric_quantities ADM_metric_derivs_x;
        GRHayLMHD_compute_metric_derivs(
              cctkGH, i, j, k,
              0, dxi, alp,
              betax, betay, betaz,
              gxx, gxy, gxz,
              gyy, gyz, gzz,
              &ADM_metric_derivs_x);

        ghl_metric_quantities ADM_metric_derivs_y;
        GRHayLMHD_compute_metric_derivs(
              cctkGH, i, j, k,
              1, dyi, alp,
              betax, betay, betaz,
              gxx, gxy, gxz,
              gyy, gyz, gzz,
              &ADM_metric_derivs_y);

        ghl_metric_quantities ADM_metric_derivs_z;
        GRHayLMHD_compute_metric_derivs(
              cctkGH, i, j, k,
              2, dzi, alp,
              betax, betay, betaz,
              gxx, gxy, gxz,
              gyy, gyz, gzz,
              &ADM_metric_derivs_z);

        ghl_conservative_quantities cons_source;
        ghl_calculate_source_terms(
              ghl_eos, &prims, &ADM_metric,
              &ADM_metric_derivs_x,
              &ADM_metric_derivs_y,
              &ADM_metric_derivs_z,
              &curv, &cons_source);

        tau_rhs [index] = cons_source.tau;
        st_x_rhs[index] = cons_source.SD[0];
        st_y_rhs[index] = cons_source.SD[1];
        st_z_rhs[index] = cons_source.SD[2];
      }
    }
  }
}

void GRHayLMHD_tabulated_entropy_calculate_flux_dir_rhs(
      const cGH *restrict cctkGH,
      const int flux_dir,
      const CCTK_REAL **B_center,
      const CCTK_REAL *restrict B_stagger,
      CCTK_REAL **vel_r,
      CCTK_REAL **vel_l,
      CCTK_REAL *restrict cmin,
      CCTK_REAL *restrict cmax) {
  DECLARE_CCTK_ARGUMENTS;

  const int imin = cctkGH->cctk_nghostzones[0];
  const int jmin = cctkGH->cctk_nghostzones[1];
  const int kmin = cctkGH->cctk_nghostzones[2];
  const int imax = cctkGH->cctk_lsh[0] - cctkGH->cctk_nghostzones[0];
  const int jmax = cctkGH->cctk_lsh[1] - cctkGH->cctk_nghostzones[1];
  const int kmax = cctkGH->cctk_lsh[2] - cctkGH->cctk_nghostzones[2];

  // Function pointer to allow for loop over fluxes and sources
  void (*calculate_characteristic_speed)(
        ghl_primitive_quantities *restrict prims_r,
        ghl_primitive_quantities *restrict prims_l,
        const ghl_eos_parameters *restrict eos,
        const ghl_metric_quantities *restrict ADM_metric_face,
        double *cmin, double *cmax);

  void (*calculate_HLLE_fluxes)(
        ghl_primitive_quantities *restrict prims_r,
        ghl_primitive_quantities *restrict prims_l,
        const ghl_eos_parameters *restrict eos,
        const ghl_metric_quantities *restrict ADM_metric_face,
        const double cmin,
        const double cmax,
        ghl_conservative_quantities *restrict cons_fluxes);

  const int xdir = (flux_dir == 0);
  const int ydir = (flux_dir == 1);
  const int zdir = (flux_dir == 2);

  const CCTK_REAL *v_flux;
  int B_recon[3];
  // Set function pointer to specific function for a given direction
  switch(flux_dir) {
    case 0:
      v_flux = vx;
      B_recon[0] = 0;
      B_recon[1] = 1;
      B_recon[2] = 2;
      calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn0;
      calculate_HLLE_fluxes = ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy;
      break;
    case 1:
      v_flux = vy;
      B_recon[0] = 1;
      B_recon[1] = 2;
      B_recon[2] = 0;
      calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn1;
      calculate_HLLE_fluxes = ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy;
      break;
    case 2:
      v_flux = vz;
      B_recon[0] = 2;
      B_recon[1] = 0;
      B_recon[2] = 1;
      calculate_characteristic_speed = &ghl_calculate_characteristic_speed_dirn2;
      calculate_HLLE_fluxes = ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy;
      break;
    default:
      CCTK_ERROR("Invalid flux_dir value (not 0, 1, or 2) has been passed to calculate_MHD_rhs.");
  }

  // This loop fills in all the data for reconstructed velocities. This loop is larger
  // than the others because we will need to reconstruct a second time.
  const int vimin = xdir*cctkGH->cctk_nghostzones[0];
  const int vjmin = ydir*cctkGH->cctk_nghostzones[1];
  const int vkmin = zdir*cctkGH->cctk_nghostzones[2];
  const int vimax = cctkGH->cctk_lsh[0] - xdir*(cctkGH->cctk_nghostzones[0] - 1);
  const int vjmax = cctkGH->cctk_lsh[1] - ydir*(cctkGH->cctk_nghostzones[1] - 1);
  const int vkmax = cctkGH->cctk_lsh[2] - zdir*(cctkGH->cctk_nghostzones[2] - 1);

#pragma omp parallel for
  for(int k=vkmin; k<vkmax; k++) {
    for(int j=vjmin; j<vjmax; j++) {
      for(int i=vimin; i<vimax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        CCTK_REAL press_stencil[6], v_flux_dir[6];
        CCTK_REAL vx_data[6], vy_data[6], vz_data[6];
        CCTK_REAL vxr, vxl, vyr, vyl, vzr, vzl;

        for(int ind=0; ind<6; ind++) {
          // Stencil from -3 to +2 reconstructs to e.g. i-1/2
          const int stencil = CCTK_GFINDEX3D(cctkGH, i+xdir*(ind-3), j+ydir*(ind-3), k+zdir*(ind-3));
          v_flux_dir[ind] = v_flux[stencil]; // Could be smaller; doesn't use full stencil
          press_stencil[ind] = P[stencil];
          vx_data[ind] = vx[stencil];
          vy_data[ind] = vy[stencil];
          vz_data[ind] = vz[stencil];
        }

        double ftilde[2];
        ghl_compute_ftilde(ghl_params, press_stencil, v_flux_dir, ftilde);

        ghl_ppm_reconstruction(ftilde, vx_data, &vxr, &vxl);
        ghl_ppm_reconstruction(ftilde, vy_data, &vyr, &vyl);
        ghl_ppm_reconstruction(ftilde, vz_data, &vzr, &vzl);

        vel_r[0][index] = vxr;
        vel_r[1][index] = vyr;
        vel_r[2][index] = vzr;

        vel_l[0][index] = vxl;
        vel_l[1][index] = vyl;
        vel_l[2][index] = vzl;
      }
    }
  }

  // This loop includes 1 ghostzone because the RHS calculation for e.g. the x direction
  // requires (i,j,k) and (i+1,j,k); if cmin/max weren't also needed for A_i, we could
  // technically have the loop only go 1 extra point in the flux_dir direction
#pragma omp parallel for
  for(int k=kmin; k<kmax+1; k++) {
    for(int j=jmin; j<jmax+1; j++) {
      for(int i=imin; i<imax+1; i++) {
        const int indm1 = CCTK_GFINDEX3D(cctkGH, i-xdir, j-ydir, k-zdir);
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        ghl_metric_quantities ADM_metric_face;
        GRHayLMHD_interpolate_metric_to_face(
              cctkGH, i, j, k,
              flux_dir, alp,
              betax, betay, betaz,
              gxx, gxy, gxz,
              gyy, gyz, gzz,
              &ADM_metric_face);

        CCTK_REAL rho_stencil[6], press_stencil[6], v_flux_dir[6];
        CCTK_REAL B1_stencil[6], B2_stencil[6], ent_stencil[6], Ye_stencil[6];
        ghl_primitive_quantities prims_r, prims_l;

        for(int ind=0; ind<6; ind++) {
          // Stencil from -3 to +2 reconstructs to e.g. i-1/2
          const int stencil  = CCTK_GFINDEX3D(cctkGH, i+xdir*(ind-3), j+ydir*(ind-3), k+zdir*(ind-3));
          v_flux_dir[ind]    = v_flux[stencil]; // Could be smaller; doesn't use full stencil
          rho_stencil[ind]   = rho_b[stencil];
          press_stencil[ind] = P[stencil];
          B1_stencil[ind]    = B_center[B_recon[1]][stencil] * ONE_OVER_SQRT_4PI;
          B2_stencil[ind]    = B_center[B_recon[2]][stencil] * ONE_OVER_SQRT_4PI;
          ent_stencil[ind]   = igm_entropy[stencil];
          Ye_stencil[ind]    = igm_Ye[stencil];
        }

        double ftilde[2];
        ghl_compute_ftilde(ghl_params, press_stencil, v_flux_dir, ftilde);

        const CCTK_REAL Gamma = get_Gamma_eff_tabulated(rho_b[index], P[index]);
        ghl_ppm_reconstruction_with_steepening(ghl_params, press_stencil, Gamma, ftilde, rho_stencil, &prims_r.rho, &prims_l.rho);

        ghl_ppm_reconstruction(ftilde, press_stencil, &prims_r.press, &prims_l.press);
        ghl_ppm_reconstruction(ftilde, B1_stencil, &prims_r.BU[B_recon[1]], &prims_l.BU[B_recon[1]]);
        ghl_ppm_reconstruction(ftilde, B2_stencil, &prims_r.BU[B_recon[2]], &prims_l.BU[B_recon[2]]);
        ghl_ppm_reconstruction(ftilde, ent_stencil, &prims_r.entropy, &prims_l.entropy);
        ghl_ppm_reconstruction(ftilde, Ye_stencil, &prims_r.Y_e, &prims_l.Y_e);

        // B_stagger is densitized, but B_center is not.
        prims_r.BU[B_recon[0]] = prims_l.BU[B_recon[0]] = B_stagger[indm1]/ADM_metric_face.sqrt_detgamma * ONE_OVER_SQRT_4PI;

        prims_r.vU[0] = vel_r[0][index];
        prims_r.vU[1] = vel_r[1][index];
        prims_r.vU[2] = vel_r[2][index];

        prims_l.vU[0] = vel_l[0][index];
        prims_l.vU[1] = vel_l[1][index];
        prims_l.vU[2] = vel_l[2][index];

        prims_r.temperature = prims_l.temperature = temperature[index];

        int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_limit_v_and_compute_u0(ghl_params, &ADM_metric_face, &prims_r);
        speed_limited = ghl_limit_v_and_compute_u0(ghl_params, &ADM_metric_face, &prims_l);

        // We must now compute eps and T
        ghl_tabulated_enforce_bounds_rho_Ye_P(ghl_eos, &prims_r.rho, &prims_r.Y_e, &prims_r.press);
        ghl_tabulated_compute_eps_T_from_P(ghl_eos, prims_r.rho, prims_r.Y_e, prims_r.press,
                                           &prims_r.eps, &prims_r.temperature);

        ghl_tabulated_enforce_bounds_rho_Ye_P(ghl_eos, &prims_l.rho, &prims_l.Y_e, &prims_l.press);
        ghl_tabulated_compute_eps_T_from_P(ghl_eos, prims_l.rho, prims_l.Y_e, prims_l.press,
                                           &prims_l.eps, &prims_l.temperature);

        ghl_conservative_quantities cons_fluxes;
        calculate_characteristic_speed(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, &cmin[index], &cmax[index]);
        calculate_HLLE_fluxes(&prims_r, &prims_l, ghl_eos, &ADM_metric_face, cmin[index], cmax[index], &cons_fluxes);

        rho_star_flux[index] = cons_fluxes.rho;
        tau_flux     [index] = cons_fluxes.tau;
        st_x_flux [index] = cons_fluxes.SD[0];
        st_y_flux [index] = cons_fluxes.SD[1];
        st_z_flux [index] = cons_fluxes.SD[2];
        S_star_flux[index] = cons_fluxes.entropy;
        Ye_star_flux [index] = cons_fluxes.Y_e;
      }
    }
  }

  const CCTK_REAL dxi = 1.0/CCTK_DELTA_SPACE(flux_dir);

#pragma omp parallel for
  for(int k=kmin; k<kmax; k++) {
    for(int j=jmin; j<jmax; j++) {
      for(int i=imin; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j ,k);
        const int indp1 = CCTK_GFINDEX3D(cctkGH, i+xdir, j+ydir, k+zdir);

        rho_star_rhs[index] += dxi*(rho_star_flux[index] - rho_star_flux[indp1]);
        tau_rhs     [index] += dxi*(tau_flux     [index] - tau_flux     [indp1]);
        st_x_rhs [index] += dxi*(st_x_flux [index] - st_x_flux [indp1]);
        st_y_rhs [index] += dxi*(st_y_flux [index] - st_y_flux [indp1]);
        st_z_rhs [index] += dxi*(st_z_flux [index] - st_z_flux [indp1]);
        S_star_rhs[index] += dxi*(S_star_flux[index] - S_star_flux[indp1]);
        Ye_star_rhs [index] += dxi*(Ye_star_flux [index] - Ye_star_flux [indp1]);
      }
    }
  }
}

void GRHayLMHD_evaluate_psi6phi_and_A_gauge_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Note that in this function, we don't bother with reconstruction, instead interpolating.
  // We need A^i, but only have A_i. So we use the BSSN metric gtupij.
  // The reconstruction variables are temporary variables and the data in them can be safely overwritten,
  // saving some memory.
  double *betax_interp = vxr;
  double *betay_interp = vyr;
  double *betaz_interp = vzr;
  double *alpha_interp = vxrr;
  double *alpha_Phi_minus_betaj_A_j_interp = vxll;
  double *sqrtg_Ax_interp = vxl;
  double *sqrtg_Ay_interp = vyl;
  double *sqrtg_Az_interp = vzl;

  /* Compute \partial_t psi6phi = -\partial_i (  \alpha psi^6 A^i - psi6phi \beta^i)
   *    (Eq 13 of http://arxiv.org/pdf/1110.4633.pdf), using Lorenz gauge.
   * Note that the RHS consists of a shift advection term on psi6phi and
   *    a term depending on the vector potential.
   * psi6phi is defined at (i+1/2,j+1/2,k+1/2), but instead of reconstructing
   *    to compute the RHS of \partial_t psi6phi, we instead use standard
   *    interpolations.
   */

  // We declare these values to be over the interior so the setting of ghostzone points is
  // more transparent.
  const int imin = cctkGH->cctk_nghostzones[0];
  const int jmin = cctkGH->cctk_nghostzones[1];
  const int kmin = cctkGH->cctk_nghostzones[2];
  const int imax = cctkGH->cctk_lsh[0]-cctkGH->cctk_nghostzones[0];
  const int jmax = cctkGH->cctk_lsh[1]-cctkGH->cctk_nghostzones[1];
  const int kmax = cctkGH->cctk_lsh[2]-cctkGH->cctk_nghostzones[2];

  // The RHS loop requires 2 ghostzones, so this loop must set those values, hence this loop
  // goes into the ghostzones. This loop requires a stencil of {-1,1},{-1,1},{-1,1},
  // which uses the last of the 3 ghostzones required by the simulation.
  // Note that ALL input variables are defined at ALL gridpoints, so no
  // worries about ghostzones.
#pragma omp parallel for
  for(int k=kmin-2; k<kmax+2; k++) {
    for(int j=jmin-2; j<jmax+2; j++) {
      for(int i=imin-2; i<imax+2; i++) {
        const int index=CCTK_GFINDEX3D(cctkGH,i,j,k);

        // First compute \partial_j \alpha \sqrt{\gamma} A^j (RHS of \partial_i psi6phi)
        // FIXME: Would be much cheaper & easier to unstagger A_i, raise, then interpolate A^i.
        //        However, we keep it this way to be completely compatible with the original
        //        Illinois GRMHD thorn, called mhd_evolve.
        //
        //Step 1) j=x: Need to raise A_i, but to do that, we must have all variables at the same gridpoints:
        // The goal is to compute \partial_j (\alpha \sqrt{\gamma} A^j) at (i+1/2,j+1/2,k+1/2)
        //    We do this by first interpolating (RHS1x) = (\alpha \sqrt{\gamma} A^x) at
        //    (i,j+1/2,k+1/2)and (i+1,j+1/2,k+1/2), then taking \partial_x (RHS1x) =
        //    [ RHS1x(i+1,j+1/2,k+1/2) - RHS1x(i,j+1/2,k+1/2) ]/dX.
        // First bring gtup's, psi, and alpha to (i,j+1/2,k+1/2):
        ghl_metric_quantities metric_stencil[2][2][2];
        double Ax_stencil[3][3][3];
        double Ay_stencil[3][3][3];
        double Az_stencil[3][3][3];
        induction_interp_vars interp_vars;

        // Read in variable at interpolation stencil points from main memory.
        for(int iterz=0; iterz<2; iterz++)
          for(int itery=0; itery<2; itery++)
            for(int iterx=0; iterx<2; iterx++) {
              const int ind = CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz);

              ghl_initialize_metric(
                    alp[ind], betax[ind], betay[ind], betaz[ind],
                    gxx[ind], gxy[ind], gxz[ind],
                    gyy[ind], gyz[ind], gzz[ind],
                    &metric_stencil[iterz][itery][iterx]);
        }
        // A_x needs a stencil s.t. interp_limits={ 0,1,-1,1,-1,1}.
        // A_y needs a stencil s.t. interp_limits={-1,1, 0,1,-1,1}.
        // A_z needs a stencil s.t. interp_limits={-1,1,-1,1, 0,1}.
        // We could fill only the needed elements, but it is cleaner
        // to fill in the whole 3x3x3 array.
        // TODO: the old code explicitly only filled in the necessary
        // elements. If we want to remove ~15 memcopies, do that here.
        for(int iterz=-1; iterz<2; iterz++)
          for(int itery=-1; itery<2; itery++)
            for(int iterx=-1; iterx<2; iterx++) {
              const int ind = CCTK_GFINDEX3D(cctkGH,i+iterx,j+itery,k+iterz);
              Ax_stencil[iterz+1][itery+1][iterx+1] = Ax[ind];
              Ay_stencil[iterz+1][itery+1][iterx+1] = Ay[ind];
              Az_stencil[iterz+1][itery+1][iterx+1] = Az[ind];
        }
// This code should only copy the needed data that isn't copied in the loop for other variables, but it is untested.
//        for(int iter2=0; iter2<2; iter2++)
//        for(int iter1=0; iter1<2; iter1++) {
//          gauge_vars.A_x[iter2+1][0][iter1+1] = in_vars[A_XI][CCTK_GFINDEX3D(cctkGH, i+iter1,     j-1, k+iter2)]; // { (0,1),    -1, (0,1)}
//          gauge_vars.A_x[0][iter2+1][iter1+1] = in_vars[A_XI][CCTK_GFINDEX3D(cctkGH, i+iter1, j+iter2, k-1    )]; // { (0,1), (0,1), -1}
//          gauge_vars.A_y[iter2+1][iter1+1][0] = in_vars[A_YI][CCTK_GFINDEX3D(cctkGH,     i-1, j+iter1, k+iter2)]; // { -1,    (0,1), (0,1)}
//          gauge_vars.A_y[0][iter1+1][iter2+1] = in_vars[A_YI][CCTK_GFINDEX3D(cctkGH, i+iter2, j+iter1, k-1    )]; // { (0,1), (0,1), -1}
//          gauge_vars.A_z[iter1+1][iter2+1][0] = in_vars[A_ZI][CCTK_GFINDEX3D(cctkGH,     i-1, j+iter2, k+iter1)]; // { -1,    (0,1), (0,1)}
//          gauge_vars.A_z[iter1+1][0][iter2+1] = in_vars[A_ZI][CCTK_GFINDEX3D(cctkGH, i+iter2,     j-1, k+iter1)]; // { (0,1),    -1, (0,1)}
//        }

        ghl_interpolate_with_cell_centered_ADM(metric_stencil, Ax_stencil, Ay_stencil, Az_stencil, psi6phi[index], &interp_vars);

        alpha_interp[index] = interp_vars.alpha;
        sqrtg_Ax_interp[index] = interp_vars.sqrtg_Ai[0];
        sqrtg_Ay_interp[index] = interp_vars.sqrtg_Ai[1];
        sqrtg_Az_interp[index] = interp_vars.sqrtg_Ai[2];
        alpha_Phi_minus_betaj_A_j_interp[index] = interp_vars.alpha_Phi_minus_betaj_A_j;
        betax_interp[index] = interp_vars.betai[0];
        betay_interp[index] = interp_vars.betai[1];
        betaz_interp[index] = interp_vars.betai[2];
      }
    }
  }

  const CCTK_REAL dxi[3] = { 1.0/CCTK_DELTA_SPACE(0), 1.0/CCTK_DELTA_SPACE(1), 1.0/CCTK_DELTA_SPACE(2) };

#pragma omp parallel for
  for(int k=kmin; k<kmax; k++) {
    for(int j=jmin; j<jmax; j++) {
      for(int i=imin; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // \partial_t A_i = [reconstructed stuff] + [gauge stuff],
        //    where [gauge stuff] = -\partial_i (\alpha \Phi - \beta^j A_j)
        Ax_rhs[index] += dxi[0]*(alpha_Phi_minus_betaj_A_j_interp[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] - alpha_Phi_minus_betaj_A_j_interp[index]);
        Ay_rhs[index] += dxi[1]*(alpha_Phi_minus_betaj_A_j_interp[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] - alpha_Phi_minus_betaj_A_j_interp[index]);
        Az_rhs[index] += dxi[2]*(alpha_Phi_minus_betaj_A_j_interp[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] - alpha_Phi_minus_betaj_A_j_interp[index]);

        double betax_stencil[5], betay_stencil[5], betaz_stencil[5];
        double psi6phi_stencil[3][5], sqrtg_Ai_stencil[3][2];

        sqrtg_Ai_stencil[0][0] = sqrtg_Ax_interp[index];
        sqrtg_Ai_stencil[1][0] = sqrtg_Ay_interp[index];
        sqrtg_Ai_stencil[2][0] = sqrtg_Az_interp[index];

        sqrtg_Ai_stencil[0][1] = sqrtg_Ax_interp[CCTK_GFINDEX3D(cctkGH,i+1,j,k)];
        sqrtg_Ai_stencil[1][1] = sqrtg_Ay_interp[CCTK_GFINDEX3D(cctkGH,i,j+1,k)];
        sqrtg_Ai_stencil[2][1] = sqrtg_Az_interp[CCTK_GFINDEX3D(cctkGH,i,j,k+1)];

        for(int iter=-2; iter<3; iter++) {
          const int indexx = CCTK_GFINDEX3D(cctkGH,i+iter,j,     k     );
          const int indexy = CCTK_GFINDEX3D(cctkGH,i,     j+iter,k     );
          const int indexz = CCTK_GFINDEX3D(cctkGH,i,     j,     k+iter);
          betax_stencil[iter+2] = betax_interp[indexx];
          betay_stencil[iter+2] = betay_interp[indexy];
          betaz_stencil[iter+2] = betaz_interp[indexz];
          psi6phi_stencil[0][iter+2] = psi6phi[indexx];
          psi6phi_stencil[1][iter+2] = psi6phi[indexy];
          psi6phi_stencil[2][iter+2] = psi6phi[indexz];
        }
        psi6phi_rhs[index] = ghl_calculate_phitilde_rhs(dxi, ghl_params->Lorenz_damping_factor, alpha_interp[index], betax_stencil, betay_stencil, betaz_stencil, sqrtg_Ai_stencil, psi6phi_stencil);
      }
    }
  }
}

static void GRHayLMHD_reconstruction_loop(const cGH *restrict cctkGH, const int flux_dir, const int num_vars,
                         const int *restrict var_indices,
                         const double *rho_b,
                         const double *pressure,
                         const double *v_flux,
                         const double **in_prims,
                         double **out_prims_r,
                         double **out_prims_l) {

  const int xdir = (flux_dir == 0);
  const int ydir = (flux_dir == 1);
  const int zdir = (flux_dir == 2);

  const int imin = cctkGH->cctk_nghostzones[0];
  const int imax = (cctkGH->cctk_lsh[0] - cctkGH->cctk_nghostzones[0]) + 1;
  const int jmin = cctkGH->cctk_nghostzones[1];
  const int jmax = (cctkGH->cctk_lsh[1] - cctkGH->cctk_nghostzones[1]) + 1;
  const int kmin = cctkGH->cctk_nghostzones[2];
  const int kmax = (cctkGH->cctk_lsh[2] - cctkGH->cctk_nghostzones[2]) + 1;

#pragma omp parallel for
  for(int k=kmin; k<kmax; k++) {
    for(int j=jmin; j<jmax; j++) {
      for(int i=imin; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
        double press_stencil[6], v_flux_stencil[6];
        double var_data[num_vars][6], vars_r[num_vars], vars_l[num_vars];

        for(int ind=0; ind<6; ind++) {
          // Stencil from -3 to +2 reconstructs to e.g. i-1/2
          const int stencil = CCTK_GFINDEX3D(cctkGH, i+xdir*(ind-3), j+ydir*(ind-3), k+zdir*(ind-3));
          v_flux_stencil[ind] = v_flux[stencil]; // Could be smaller; doesn't use full stencil
          press_stencil[ind] = pressure[stencil];
          for(int var=0; var<num_vars; var++) {
            var_data[var][ind] = in_prims[var_indices[var]][stencil];
          }
        }

        double ftilde[2];
        ghl_compute_ftilde(ghl_params, press_stencil, v_flux_stencil, ftilde);

        for(int var=0; var<num_vars; var++) {
          ghl_ppm_reconstruction(ftilde, var_data[var], &vars_r[var], &vars_l[var]);
          out_prims_r[var_indices[var]][index] = vars_r[var];
          out_prims_l[var_indices[var]][index] = vars_l[var];
        }
      }
    }
  }
}

void GRHayLMHD_A_flux_rhs(
      const cGH *restrict cctkGH,
      const int flux_dir,
      /*const*/ double **in_prims_r,
      /*const*/ double **in_prims_l,
      /*const*/ double **cmin,
      /*const*/ double **cmax,
      double *restrict A_rhs) {

  const int xdir = (flux_dir==0);
  const int ydir = (flux_dir==1);
  const int zdir = (flux_dir==2);

  // These are used to determine which components of v,
  // B, cmin, and cmax are used in the computation.
  const int dir1_offset = (flux_dir+1)%3, dir2_offset = (flux_dir+2)%3;

  // This offsets the index by +1 in the perpendicular directions
  const int v_offset[3] = { !xdir, !ydir, !zdir };

  // This offsets the index by +1 in the permuted direction (x<-y<-z)
  const int B1_offset[3] = { ydir, zdir, xdir };

  // This offsets the index by +1 in the permuted direction (x->y->z)
  const int B2_offset[3] = { zdir, xdir, ydir };

  const int imin = cctkGH->cctk_nghostzones[0];
  const int imax = cctkGH->cctk_lsh[0]-cctkGH->cctk_nghostzones[0];
  const int jmin = cctkGH->cctk_nghostzones[1];
  const int jmax = cctkGH->cctk_lsh[1]-cctkGH->cctk_nghostzones[1];
  const int kmin = cctkGH->cctk_nghostzones[2];
  const int kmax = cctkGH->cctk_lsh[2]-cctkGH->cctk_nghostzones[2];

#pragma omp parallel for
  for(int k=kmin; k<kmax; k++) {
    for(int j=jmin; j<jmax; j++) {
      for(int i=imin; i<imax; i++) {
        const int index    = CCTK_GFINDEX3D(cctkGH,i,j,k);
        const int index_v  = CCTK_GFINDEX3D(cctkGH,i+v_offset[0], j+v_offset[1], k+v_offset[2]);
        const int index_B1 = CCTK_GFINDEX3D(cctkGH,i+B1_offset[0],j+B1_offset[1],k+B1_offset[2]);
        const int index_B2 = CCTK_GFINDEX3D(cctkGH,i+B2_offset[0],j+B2_offset[1],k+B2_offset[2]);

        HLL_2D_vars vars;

        vars.v1rr = in_prims_r[VXR+dir1_offset][index_v];
        vars.v1rl = in_prims_l[VXR+dir1_offset][index_v];
        vars.v1lr = in_prims_r[VXL+dir1_offset][index_v];
        vars.v1ll = in_prims_l[VXL+dir1_offset][index_v];

        vars.v2rr = in_prims_r[VXR+dir2_offset][index_v];
        vars.v2rl = in_prims_l[VXR+dir2_offset][index_v];
        vars.v2lr = in_prims_r[VXL+dir2_offset][index_v];
        vars.v2ll = in_prims_l[VXL+dir2_offset][index_v];

        vars.B1r = in_prims_r[BX_STAGGER+dir1_offset][index_B1];
        vars.B1l = in_prims_l[BX_STAGGER+dir1_offset][index_B1];

        vars.B2r = in_prims_r[BX_STAGGER+dir2_offset][index_B2];
        vars.B2l = in_prims_l[BX_STAGGER+dir2_offset][index_B2];

        /*
          Note that cmax/cmin (\alpha^{\pm}  as defined in Del Zanna et al) is at a slightly DIFFERENT
          point (e.g., (i+1/2,j,k) instead of (i+1/2,j+1/2,k) for F3). Yuk Tung Liu discussed this point
          with M. Shibata, who found that the effect is negligible.
        */
        vars.c1_min = cmin[dir1_offset][index_B2];
        vars.c1_max = cmax[dir1_offset][index_B2];
        vars.c2_min = cmin[dir2_offset][index_B1];
        vars.c2_max = cmax[dir2_offset][index_B1];

        A_rhs[index] = ghl_HLL_2D_flux_with_Btilde(&vars);
      }
    }
  }
}

static void GRHayLMHD_compute_metric_derivs(
      const cGH *cctkGH,
      const int i, const int j, const int k,
      const int flux_dir,
      const CCTK_REAL dxi,
      const CCTK_REAL *restrict lapse,
      const CCTK_REAL *restrict betax,
      const CCTK_REAL *restrict betay,
      const CCTK_REAL *restrict betaz,
      const CCTK_REAL *restrict gxx,
      const CCTK_REAL *restrict gxy,
      const CCTK_REAL *restrict gxz,
      const CCTK_REAL *restrict gyy,
      const CCTK_REAL *restrict gyz,
      const CCTK_REAL *restrict gzz,
      ghl_metric_quantities *restrict metric_derivs) {

  const int xdir = (flux_dir == 0);
  const int ydir = (flux_dir == 1);
  const int zdir = (flux_dir == 2);

  const int indm2  = CCTK_GFINDEX3D(cctkGH, i-2*xdir, j-2*ydir, k-2*zdir);
  const int indm1  = CCTK_GFINDEX3D(cctkGH, i-xdir,   j-ydir,   k-zdir);
  const int indp1  = CCTK_GFINDEX3D(cctkGH, i+xdir,   j+ydir,   k+zdir);
  const int indp2  = CCTK_GFINDEX3D(cctkGH, i+2*xdir, j+2*ydir, k+2*zdir);

  const CCTK_REAL d_lapse = dxi*COMPUTE_DERIV(lapse[indm2], lapse[indm1], lapse[indp1], lapse[indp2]);
  const CCTK_REAL d_betax = dxi*COMPUTE_DERIV(betax[indm2], betax[indm1], betax[indp1], betax[indp2]);
  const CCTK_REAL d_betay = dxi*COMPUTE_DERIV(betay[indm2], betay[indm1], betay[indp1], betay[indp2]);
  const CCTK_REAL d_betaz = dxi*COMPUTE_DERIV(betaz[indm2], betaz[indm1], betaz[indp1], betaz[indp2]);

  const CCTK_REAL d_gxx = dxi*COMPUTE_DERIV(gxx[indm2], gxx[indm1], gxx[indp1], gxx[indp2]);
  const CCTK_REAL d_gxy = dxi*COMPUTE_DERIV(gxy[indm2], gxy[indm1], gxy[indp1], gxy[indp2]);
  const CCTK_REAL d_gxz = dxi*COMPUTE_DERIV(gxz[indm2], gxz[indm1], gxz[indp1], gxz[indp2]);
  const CCTK_REAL d_gyy = dxi*COMPUTE_DERIV(gyy[indm2], gyy[indm1], gyy[indp1], gyy[indp2]);
  const CCTK_REAL d_gyz = dxi*COMPUTE_DERIV(gyz[indm2], gyz[indm1], gyz[indp1], gyz[indp2]);
  const CCTK_REAL d_gzz = dxi*COMPUTE_DERIV(gzz[indm2], gzz[indm1], gzz[indp1], gzz[indp2]);

  ghl_initialize_metric(
        d_lapse,
        d_betax, d_betay, d_betaz,
        d_gxx, d_gxy, d_gxz,
        d_gyy, d_gyz, d_gzz,
        metric_derivs);
}

static void GRHayLMHD_interpolate_metric_to_face(
      const cGH *cctkGH,
      const int i, const int j, const int k,
      const int flux_dir,
      const CCTK_REAL *restrict lapse,
      const CCTK_REAL *restrict betax,
      const CCTK_REAL *restrict betay,
      const CCTK_REAL *restrict betaz,
      const CCTK_REAL *restrict gxx,
      const CCTK_REAL *restrict gxy,
      const CCTK_REAL *restrict gxz,
      const CCTK_REAL *restrict gyy,
      const CCTK_REAL *restrict gyz,
      const CCTK_REAL *restrict gzz,
      ghl_metric_quantities *restrict metric) {

  const int xdir = (flux_dir == 0);
  const int ydir = (flux_dir == 1);
  const int zdir = (flux_dir == 2);

  const int indm2  = CCTK_GFINDEX3D(cctkGH, i-2*xdir, j-2*ydir, k-2*zdir);
  const int indm1  = CCTK_GFINDEX3D(cctkGH, i-xdir,   j-ydir,   k-zdir);
  const int index  = CCTK_GFINDEX3D(cctkGH, i,        j ,       k);
  const int indp1  = CCTK_GFINDEX3D(cctkGH, i+xdir,   j+ydir,   k+zdir);

  const CCTK_REAL face_lapse = COMPUTE_FCVAL(lapse[indm2], lapse[indm1], lapse[index], lapse[indp1]);
  const CCTK_REAL face_betax = COMPUTE_FCVAL(betax[indm2], betax[indm1], betax[index], betax[indp1]);
  const CCTK_REAL face_betay = COMPUTE_FCVAL(betay[indm2], betay[indm1], betay[index], betay[indp1]);
  const CCTK_REAL face_betaz = COMPUTE_FCVAL(betaz[indm2], betaz[indm1], betaz[index], betaz[indp1]);

  const CCTK_REAL face_gxx = COMPUTE_FCVAL(gxx[indm2], gxx[indm1], gxx[index], gxx[indp1]);
  const CCTK_REAL face_gxy = COMPUTE_FCVAL(gxy[indm2], gxy[indm1], gxy[index], gxy[indp1]);
  const CCTK_REAL face_gxz = COMPUTE_FCVAL(gxz[indm2], gxz[indm1], gxz[index], gxz[indp1]);
  const CCTK_REAL face_gyy = COMPUTE_FCVAL(gyy[indm2], gyy[indm1], gyy[index], gyy[indp1]);
  const CCTK_REAL face_gyz = COMPUTE_FCVAL(gyz[indm2], gyz[indm1], gyz[index], gyz[indp1]);
  const CCTK_REAL face_gzz = COMPUTE_FCVAL(gzz[indm2], gzz[indm1], gzz[index], gzz[indp1]);

  ghl_initialize_metric(
        face_lapse,
        face_betax, face_betay, face_betaz,
        face_gxx, face_gxy, face_gxz,
        face_gyy, face_gyz, face_gzz,
        metric);
}
