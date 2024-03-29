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

#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

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
  compute_tau_rhs_extrinsic_curvature_terms_and_TUPmunu(eos,cctkGH,cctk_lsh,cctk_nghostzones,dX,metric,in_prims,TUPmunu,
                                                        gtupxy,gtupxz,gtupyz,
                                                        kxx,kxy,kxz,kyy,kyz,kzz,
                                                        tau_rhs, s_tau);

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
