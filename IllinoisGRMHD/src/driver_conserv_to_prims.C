#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include "IllinoisGRMHD.h"

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#include "IllinoisGRMHD_headers.h"
#include "con2prim_headers.h"
#include "inlined_functions.h"
#include "apply_tau_floor__enforce_limits_on_primitives_and_recompute_conservs.C"

extern "C" void IllinoisGRMHD_conserv_to_prims(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // We use proper C++ here, for file I/O later.
  using namespace std;

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

  // These BSSN-based variables are not evolved, and so are not defined anywhere that the grid has moved.
  // Here we convert ADM variables (from ADMBase) to the BSSN-based variables expected by this routine.
  IllinoisGRMHD_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij(cctkGH,cctk_lsh,  gxx,gxy,gxz,gyy,gyz,gzz,alp,
                                                                gtxx,gtxy,gtxz,gtyy,gtyz,gtzz,
                                                                gtupxx,gtupxy,gtupxz,gtupyy,gtupyz,gtupzz,
                                                                phi_bssn,psi_bssn,lapm1);


#ifndef ENABLE_STANDALONE_IGM_C2P_SOLVER
  if(CCTK_EQUALS(Symmetry,"equatorial")) {
    // SET SYMMETRY GHOSTZONES ON ALL CONSERVATIVE VARIABLES!
    int ierr=0;
    ierr+=CartSymGN(cctkGH,"IllinoisGRMHD::grmhd_conservatives");
    // FIXME: UGLY. Filling metric ghostzones is needed for, e.g., Cowling runs.
    ierr+=CartSymGN(cctkGH,"lapse::lapse_vars");
    ierr+=CartSymGN(cctkGH,"bssn::BSSN_vars");
    ierr+=CartSymGN(cctkGH,"bssn::BSSN_AH");
    ierr+=CartSymGN(cctkGH,"shift::shift_vars");
    if(ierr!=0) CCTK_VError(VERR_DEF_PARAMS,"IllinoisGRMHD ERROR (grep for it, foo!)  :(");
  }
#endif


  //Start the timer, so we can benchmark the primitives solver during evolution.
  //  Slower solver -> harder to find roots -> things may be going crazy!
  //FIXME: Replace this timing benchmark with something more meaningful, like the avg # of Newton-Raphson iterations per gridpoint!
  /*
    struct timeval start, end;
    long mtime, seconds, useconds;
    gettimeofday(&start, NULL);
  */

  int failures=0,font_fixes=0,vel_limited_ptcount=0,atm_resets=0,rho_star_fix_applied=0;
  int pointcount=0;
  int failures_inhoriz=0;
  int pointcount_inhoriz=0;
  int backup1=0,backup2=0,backup3=0;

  CCTK_REAL error_int_numer=0,error_int_denom=0;

  int imin=0,imax=cctk_lsh[0];
  int jmin=0,jmax=cctk_lsh[1];
  int kmin=0,kmax=cctk_lsh[2];

  // Whenever we get a conservative-to-primitive major failure, i.e. all
  // the routines and backups failed to recover the primitives from the
  // input conservatives, we will introduce a new fix, in which we will
  // reset the conservative variables at the given point by a weighted
  // average of the conservative variables at the neighboring points.
  // After that, the con2prim attempt will be retried. This mask allows
  // us to flag points in which the averaging procedure must be performed.
  int npoints = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];

  // We now add an integer to count the number of
  // points in which the averaging fix is required.
  // We initialize it to a nonzero value so that
  // the while condition below is triggered at least
  // once.
  int cons_avgs = 0;
  int loop_count = 0;
  int nan_found = 0;

#pragma omp parallel for reduction(+:failures,vel_limited_ptcount,font_fixes,pointcount,failures_inhoriz,pointcount_inhoriz,error_int_numer,error_int_denom,rho_star_fix_applied,atm_resets,backup1,backup2,backup3,nan_found) schedule(static)
  for(int k=kmin;k<kmax;k++) {
    for(int j=jmin;j<jmax;j++) {
      for(int i=imin;i<imax;i++) {

        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        ghl_metric_quantities ADM_metric;
        ghl_enforce_detgtij_and_initialize_ADM_metric(
              alp[index],
              betax[index], betay[index], betaz[index],
              gxx[index], gxy[index], gxz[index],
              gyy[index], gyz[index], gzz[index],
              &ADM_metric);

        ghl_ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        // Read in BSSN metric quantities from gridfunctions
        CCTK_REAL METRIC[NUMVARS_FOR_METRIC];
        METRIC[PHI   ] = phi_bssn[index];
        METRIC[GXX   ] = gtxx[index];
        METRIC[GXY   ] = gtxy[index];
        METRIC[GXZ   ] = gtxz[index];
        METRIC[GYY   ] = gtyy[index];
        METRIC[GYZ   ] = gtyz[index];
        METRIC[GZZ   ] = gtzz[index];
        METRIC[LAPM1 ] = lapm1[index];
        METRIC[SHIFTX] = betax[index];
        METRIC[SHIFTY] = betay[index];
        METRIC[SHIFTZ] = betaz[index];
        METRIC[GUPXX ] = gtupxx[index];
        METRIC[GUPYY ] = gtupyy[index];
        METRIC[GUPZZ ] = gtupzz[index];
        METRIC[GUPXY ] = gtupxy[index];
        METRIC[GUPXZ ] = gtupxz[index];
        METRIC[GUPYZ ] = gtupyz[index];

        // Set auxiliary metric quantities
        CCTK_REAL METRIC_LAP_PSI4[NUMVARS_METRIC_AUX];
        SET_LAPSE_PSI4(METRIC_LAP_PSI4,METRIC);

        // Now set the ADM metric quantities
        CCTK_REAL METRIC_PHYS[NUMVARS_FOR_METRIC];
        METRIC_PHYS[GXX  ] = METRIC[GXX  ]*METRIC_LAP_PSI4[PSI4 ];
        METRIC_PHYS[GXY  ] = METRIC[GXY  ]*METRIC_LAP_PSI4[PSI4 ];
        METRIC_PHYS[GXZ  ] = METRIC[GXZ  ]*METRIC_LAP_PSI4[PSI4 ];
        METRIC_PHYS[GYY  ] = METRIC[GYY  ]*METRIC_LAP_PSI4[PSI4 ];
        METRIC_PHYS[GYZ  ] = METRIC[GYZ  ]*METRIC_LAP_PSI4[PSI4 ];
        METRIC_PHYS[GZZ  ] = METRIC[GZZ  ]*METRIC_LAP_PSI4[PSI4 ];
        METRIC_PHYS[GUPXX] = METRIC[GUPXX]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPXY] = METRIC[GUPXY]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPXZ] = METRIC[GUPXZ]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPYY] = METRIC[GUPYY]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPYZ] = METRIC[GUPYZ]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPZZ] = METRIC[GUPZZ]*METRIC_LAP_PSI4[PSIM4];

        // Read in primitive variables from gridfunctions
        // FIXME: this seems wasteful as we won't use these values anyway
        CCTK_REAL PRIMS[old_MAXNUMVARS];
        PRIMS[RHOB         ] = rho[index];
        PRIMS[PRESSURE     ] = press[index];
        PRIMS[VX           ] = vx[index];
        PRIMS[VY           ] = vy[index];
        PRIMS[VZ           ] = vz[index];
        PRIMS[BX_CENTER    ] = Bx_center[index];
        PRIMS[BY_CENTER    ] = By_center[index];
        PRIMS[BZ_CENTER    ] = Bz_center[index];
        PRIMS[EPSILON      ] = eps[index];
        PRIMS[ENTROPY      ] = entropy[index];

        // Read in conservative variables from gridfunctions
        CCTK_REAL CONSERVS[NUM_CONSERVS],CONSERVS_avg_neighbors[NUM_CONSERVS];
        CONSERVS[RHOSTAR  ] = rho_star[index];
        CONSERVS[STILDEX  ] = Stildex[index];
        CONSERVS[STILDEY  ] = Stildey[index];
        CONSERVS[STILDEZ  ] = Stildez[index];
        CONSERVS[TAUENERGY] = tau     [index];
        CONSERVS[YESTAR   ] = Ye_star [index];
        CONSERVS[ENTSTAR  ] = ent_star  [index];

        // Tabulated EOS quantities
        if( eos.is_Tabulated ) {
          // Primitives
          PRIMS[YEPRIM     ] = Y_e[index];
          PRIMS[TEMPERATURE] = temperature[index];
        }

        CCTK_REAL shift_xL = METRIC_PHYS[GXX]*METRIC[SHIFTX] + METRIC_PHYS[GXY]*METRIC[SHIFTY] + METRIC_PHYS[GXZ]*METRIC[SHIFTZ];
        CCTK_REAL shift_yL = METRIC_PHYS[GXY]*METRIC[SHIFTX] + METRIC_PHYS[GYY]*METRIC[SHIFTY] + METRIC_PHYS[GYZ]*METRIC[SHIFTZ];
        CCTK_REAL shift_zL = METRIC_PHYS[GXZ]*METRIC[SHIFTX] + METRIC_PHYS[GYZ]*METRIC[SHIFTY] + METRIC_PHYS[GZZ]*METRIC[SHIFTZ];
        CCTK_REAL beta2L   = shift_xL*METRIC[SHIFTX] + shift_yL*METRIC[SHIFTY] + shift_zL*METRIC[SHIFTZ];


        // Compute 4-metric, both g_{\mu \nu} and g^{\mu \nu}.
        // This is for computing T_{\mu \nu} and T^{\mu \nu}. Also the HARM con2prim lowlevel function requires them.
        CCTK_REAL g4dn[4][4],g4up[4][4];
        g4dn[0][0] = -SQR(METRIC_LAP_PSI4[LAPSE]) + beta2L;
        g4dn[0][1] = g4dn[1][0] = shift_xL;
        g4dn[0][2] = g4dn[2][0] = shift_yL;
        g4dn[0][3] = g4dn[3][0] = shift_zL;
        g4dn[1][1]              = METRIC_PHYS[GXX];
        g4dn[1][2] = g4dn[2][1] = METRIC_PHYS[GXY];
        g4dn[1][3] = g4dn[3][1] = METRIC_PHYS[GXZ];
        g4dn[2][2]              = METRIC_PHYS[GYY];
        g4dn[2][3] = g4dn[3][2] = METRIC_PHYS[GYZ];
        g4dn[3][3]              = METRIC_PHYS[GZZ];

        CCTK_REAL alpha_inv_squared=SQR(METRIC_LAP_PSI4[LAPSEINV]);
        g4up[0][0] = -1.0*alpha_inv_squared;
        g4up[0][1] = g4up[1][0] = METRIC[SHIFTX]*alpha_inv_squared;
        g4up[0][2] = g4up[2][0] = METRIC[SHIFTY]*alpha_inv_squared;
        g4up[0][3] = g4up[3][0] = METRIC[SHIFTZ]*alpha_inv_squared;
        g4up[1][1]              = METRIC_PHYS[GUPXX] - METRIC[SHIFTX]*METRIC[SHIFTX]*alpha_inv_squared;
        g4up[1][2] = g4up[2][1] = METRIC_PHYS[GUPXY] - METRIC[SHIFTX]*METRIC[SHIFTY]*alpha_inv_squared;
        g4up[1][3] = g4up[3][1] = METRIC_PHYS[GUPXZ] - METRIC[SHIFTX]*METRIC[SHIFTZ]*alpha_inv_squared;
        g4up[2][2]              = METRIC_PHYS[GUPYY] - METRIC[SHIFTY]*METRIC[SHIFTY]*alpha_inv_squared;
        g4up[2][3] = g4up[3][2] = METRIC_PHYS[GUPYZ] - METRIC[SHIFTY]*METRIC[SHIFTZ]*alpha_inv_squared;
        g4up[3][3]              = METRIC_PHYS[GUPZZ] - METRIC[SHIFTZ]*METRIC[SHIFTZ]*alpha_inv_squared;


        //FIXME: might slow down the code.
        if(robust_isnan(CONSERVS[RHOSTAR]*CONSERVS[STILDEX]*CONSERVS[STILDEY]*CONSERVS[STILDEZ]*CONSERVS[TAUENERGY]*PRIMS[BX_CENTER]*PRIMS[BY_CENTER]*PRIMS[BZ_CENTER])) {
          CCTK_VWARN(CCTK_WARN_ALERT,"NAN FOUND: i,j,k = %d %d %d, x,y,z = %e %e %e , index=%d st_i = %e %e %e, rhostar = %e, tau = %e, Bi = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e",
                     i,j,k,x[index],y[index],z[index],index,
                     CONSERVS[STILDEX],CONSERVS[STILDEY],CONSERVS[STILDEZ],CONSERVS[RHOSTAR],CONSERVS[TAUENERGY],
                     PRIMS[BX_CENTER],PRIMS[BY_CENTER],PRIMS[BZ_CENTER],METRIC_PHYS[GXX],METRIC_PHYS[GXY],METRIC_PHYS[GXZ],METRIC_PHYS[GYY],METRIC_PHYS[GYZ],METRIC_PHYS[GZZ],METRIC_LAP_PSI4[PSI6]);
          nan_found++;
        }

        // Here we use _flux variables as temp storage for original values of conservative variables.. This is used for debugging purposes only.
        rho_star_flux[index]    = CONSERVS[RHOSTAR  ];
        Stildex_flux[index]        = CONSERVS[STILDEX  ];
        Stildey_flux[index]        = CONSERVS[STILDEY  ];
        Stildez_flux[index]        = CONSERVS[STILDEZ  ];
        tau_flux[index]         = CONSERVS[TAUENERGY];
        Ye_star_flux[index]     = CONSERVS[YESTAR   ];
        ent_star_flux[index]      = CONSERVS[ENTSTAR  ];

        CCTK_REAL rho_star_orig = CONSERVS[RHOSTAR  ];
        CCTK_REAL Stildex_orig = CONSERVS[STILDEX  ];
        CCTK_REAL Stildey_orig = CONSERVS[STILDEY  ];
        CCTK_REAL Stildez_orig = CONSERVS[STILDEZ  ];
        CCTK_REAL tau_orig      = CONSERVS[TAUENERGY];
        CCTK_REAL Ye_star_orig  = 0.0;
        CCTK_REAL ent_star_orig   = 0.0;
        if( eos.is_Tabulated) {
          Ye_star_orig          = CONSERVS[YESTAR   ];
        }
        if( eos.evolve_entropy ) {
          ent_star_orig           = CONSERVS[ENTSTAR  ];
        }

        int check=0;
        struct output_stats stats;
        stats.vel_limited    = 0;
        stats.failure_checker= 0;
        stats.font_fixed     = 0;
        stats.atm_reset      = 0;
        stats.backup[0]      = 0;
        stats.backup[1]      = 0;
        stats.backup[2]      = 0;
        stats.c2p_failed     = 0;
        stats.which_routine  = None;
        stats.dx[0]          = CCTK_DELTA_SPACE(0);
        stats.dx[1]          = CCTK_DELTA_SPACE(1);
        stats.dx[2]          = CCTK_DELTA_SPACE(2);
        stats.nan_found      = 0;
        if(CONSERVS[RHOSTAR]>0.0) {
          // Apply the tau floor
          if( eos.is_Hybrid ) {
            apply_tau_floor(index,Psi6threshold,PRIMS,METRIC,METRIC_PHYS,METRIC_LAP_PSI4,stats,eos,  CONSERVS);
          }

          for(int ii=0;ii<3;ii++) {
            check = con2prim(eos,
                             index,i,j,k,x,y,z,
                             METRIC,METRIC_PHYS,METRIC_LAP_PSI4,g4dn,g4up,
                             CONSERVS,PRIMS,
                             stats);
            if(check==0) ii=4;
            else stats.failure_checker+=100000;
          }
        } else {
          stats.failure_checker+=1;
          reset_prims_to_atmosphere( eos, PRIMS );
          rho_star_fix_applied++;
        }

        if( check != 0 ) {
          //--------------------------------------------------
          //----------- Primitive recovery failed ------------
          //--------------------------------------------------
          // Increment the failure flag
          // Sigh, reset to atmosphere
          reset_prims_to_atmosphere( eos, PRIMS );
          atm_resets++;
          // Then flag this point as a "success"
          check = 0;
          con2prim_failed_flag[index] = 0;
          if( eos.is_Hybrid ) {
            CCTK_VInfo(CCTK_THORNSTRING,"Couldn't find root from: %e %e %e %e %e, rhob approx=%e, rho_b_atm=%e, Bx=%e, By=%e, Bz=%e, gij_phys=%e %e %e %e %e %e, alpha=%e",
                       tau_orig,rho_star_orig,Stildex_orig,Stildey_orig,Stildez_orig,rho_star_orig/METRIC_LAP_PSI4[PSI6],eos.rho_atm,PRIMS[BX_CENTER],PRIMS[BY_CENTER],PRIMS[BZ_CENTER],METRIC_PHYS[GXX],METRIC_PHYS[GXY],METRIC_PHYS[GXZ],METRIC_PHYS[GYY],METRIC_PHYS[GYZ],METRIC_PHYS[GZZ],METRIC_LAP_PSI4[LAPSE]);
          }
          else if( eos.is_Tabulated ) {
            CCTK_VInfo(CCTK_THORNSTRING,"Couldn't find root from: %e %e %e %e %e %e %e, rhob approx=%e, rho_b_atm=%e, Bx=%e, By=%e, Bz=%e, gij_phys=%e %e %e %e %e %e, alpha=%e",
                       tau_orig,rho_star_orig,Stildex_orig,Stildey_orig,Stildez_orig,Ye_star_orig,ent_star_orig,rho_star_orig/METRIC_LAP_PSI4[PSI6],eos.rho_atm,PRIMS[BX_CENTER],PRIMS[BY_CENTER],PRIMS[BZ_CENTER],METRIC_PHYS[GXX],METRIC_PHYS[GXY],METRIC_PHYS[GXZ],METRIC_PHYS[GYY],METRIC_PHYS[GYZ],METRIC_PHYS[GZZ],METRIC_LAP_PSI4[LAPSE]);
          }
        }

        if( check == 0 ) {
          //--------------------------------------------------
          //---------- Primitive recovery succeeded ----------
          //--------------------------------------------------

          ghl_primitive_quantities prims;
          prims.rho         = PRIMS[RHOB        ];
          prims.press       = PRIMS[PRESSURE    ];
          prims.BU[0]       = Bx_center[index];
          prims.BU[1]       = By_center[index];
          prims.BU[2]       = Bz_center[index];
          prims.vU[0]       = PRIMS[VX          ];
          prims.vU[1]       = PRIMS[VY          ];
          prims.vU[2]       = PRIMS[VZ          ];
          prims.entropy     = PRIMS[ENTROPY     ];
          prims.Y_e         = PRIMS[YEPRIM     ];
          prims.temperature = PRIMS[TEMPERATURE];
          const int speed_limited = ghl_enforce_primitive_limits_and_compute_u0(
                ghl_params, ghl_eos, &ADM_metric, &prims);

          rho[index]         = prims.rho;
          press[index]       = prims.press;
          eps[index]         = prims.eps;
          u0[index]          = prims.u0;
          vx[index]          = prims.vU[0];
          vy[index]          = prims.vU[1];
          vz[index]          = prims.vU[2];
          entropy[index]     = prims.entropy;
          Y_e[index]         = prims.Y_e;
          temperature[index] = prims.temperature;

          //Now we compute the difference between original & new conservatives, for diagnostic purposes:
          error_int_numer += fabs(tau[index] - tau_orig) + fabs(rho_star[index] - rho_star_orig) +
            fabs(Stildex[index] - Stildex_orig) + fabs(Stildey[index] - Stildey_orig) + fabs(Stildez[index] - Stildez_orig);
          error_int_denom += tau_orig + rho_star_orig + fabs(Stildex_orig) + fabs(Stildey_orig) + fabs(Stildez_orig);

          if( eos.is_Tabulated ) {
            error_int_numer += fabs(Ye_star[index] - Ye_star_orig);
            error_int_denom += Ye_star_orig;
          }

          if(stats.atm_reset==1) {
            atm_resets++;
            stats.which_routine = -1;
          }
          igm_c2p_mask[index] = stats.which_routine;
          if(stats.backup[0]==1) backup1++;
          if(stats.backup[1]==1) backup2++;
          if(stats.backup[2]==1) backup3++;
          if(stats.font_fixed==1) font_fixes++;
          if(stats.nan_found==1) { CCTK_VWARN(CCTK_WARN_ALERT,"Found NAN while imposing speed limit"); nan_found++; }
          vel_limited_ptcount+=stats.vel_limited;
          if(check!=0) {
            failures++;
            if(exp(METRIC[PHI]*6.0)>Psi6threshold) {
              failures_inhoriz++;
              pointcount_inhoriz++;
            }
          }
          pointcount++;
          /***************************************************************************************************************************/
          failure_checker[index] = stats.failure_checker;
        }
      } // for(int i=imin;i<imax;i++)
    } // for(int j=jmin;j<jmax;j++)
  } // for(int k=kmin;k<kmax;k++)

  CCTK_REAL error_rho_numer = 0;
  CCTK_REAL error_tau_numer = 0;
  CCTK_REAL error_Sx_numer = 0;
  CCTK_REAL error_Sy_numer = 0;
  CCTK_REAL error_Sz_numer = 0;
  CCTK_REAL error_ent_numer = 0;
  CCTK_REAL error_Ye_numer = 0;

  CCTK_REAL error_rho_denom = 0;
  CCTK_REAL error_tau_denom = 0;
  CCTK_REAL error_Sx_denom = 0;
  CCTK_REAL error_Sy_denom = 0;
  CCTK_REAL error_Sz_denom = 0;
  CCTK_REAL error_ent_denom = 0;
  CCTK_REAL error_Ye_denom = 0;

#pragma omp parallel for reduction(+:                                   \
      error_rho_numer, error_tau_numer, error_Sx_numer, error_Sy_numer, \
      error_Sz_numer, error_rho_denom, error_tau_denom, error_Sx_denom, \
      error_Sy_denom, error_Sz_denom, error_ent_numer, error_ent_denom, \
      error_Ye_numer, error_Ye_denom) \
      schedule(static)
  for(int k=0; k<kmax; k++) {
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        ghl_metric_quantities ADM_metric;
        ghl_enforce_detgtij_and_initialize_ADM_metric(
              alp[index],
              betax[index], betay[index], betaz[index],
              gxx[index], gxy[index], gxz[index],
              gyy[index], gyz[index], gzz[index],
              &ADM_metric);

        ghl_ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        ghl_primitive_quantities prims;
        prims.rho         = rho[index];
        prims.press       = press[index];
        prims.eps         = eps[index];
        prims.u0          = u0[index];
        prims.vU[0]       = vx[index];
        prims.vU[1]       = vy[index];
        prims.vU[2]       = vz[index];
        prims.BU[0]       = Bx_center[index];
        prims.BU[1]       = By_center[index];
        prims.BU[2]       = Bz_center[index];
        prims.entropy     = entropy[index];
        prims.Y_e         = Y_e[index];
        prims.temperature = temperature[index];

        ghl_conservative_quantities cons, cons_orig;
        cons_orig.rho     = rho_star[index];
        cons_orig.tau     = tau[index];
        cons_orig.SD[0]   = Stildex[index];
        cons_orig.SD[1]   = Stildey[index];
        cons_orig.SD[2]   = Stildez[index];
        cons_orig.entropy = ent_star[index];
        cons_orig.Y_e     = Ye_star[index];

        ghl_compute_conservs(&ADM_metric, &metric_aux, &prims, &cons);

        rho_star[index] = cons.rho;
        tau[index]      = cons.tau;
        Stildex[index]  = cons.SD[0];
        Stildey[index]  = cons.SD[1];
        Stildez[index]  = cons.SD[2];
        ent_star[index] = cons.entropy;
        Ye_star[index]  = cons.Y_e;

        // for diagnostic purposes:
        error_rho_numer += fabs(cons.rho - cons_orig.rho);
        error_tau_numer += fabs(cons.tau - cons_orig.tau);
        error_Sx_numer  += fabs(cons.SD[0] - cons_orig.SD[0]);
        error_Sy_numer  += fabs(cons.SD[1] - cons_orig.SD[1]);
        error_Sz_numer  += fabs(cons.SD[2] - cons_orig.SD[2]);
        error_ent_numer += fabs(cons.entropy - cons_orig.entropy);
        error_Ye_numer  += fabs(cons.Y_e - cons_orig.Y_e);
        error_rho_denom += cons_orig.rho;
        error_tau_denom += cons_orig.tau;
        error_Sx_denom  += fabs(cons_orig.SD[0]);
        error_Sy_denom  += fabs(cons_orig.SD[1]);
        error_Sz_denom  += fabs(cons_orig.SD[2]);
        error_ent_denom += cons_orig.entropy;
        error_Ye_denom  += cons_orig.Y_e;
      }
    }
  }

  /*
    Failure checker decoder:
       1: atmosphere reset when rho_star < 0
      10: Limiting velocity u~ after C2P/Font Fix or v in ghl_enforce_primitive_limits_and_compute_u0
     100: Both C2P and Font Fix failed
      1k: backups used
     10k: tau~ was reset in ghl_apply_conservative_limits
    100k: S~ was reset in ghl_apply_conservative_limits
  */
  if(CCTK_Equals(verbose, "essential") || CCTK_Equals(verbose, "essential+iteration output")) {
    CCTK_VInfo(CCTK_THORNSTRING,"C2P: Lev: %d NumPts= %d | Fixes: BU: %d %d %d Font= %d VL= %d rho*= %d AVG= %d ATM= %d | Failures: %d InHoriz= %d / %d | Error: %.3e, ErrDenom: %.3e",
               (int)GetRefinementLevel(cctkGH),pointcount,
               backup1,backup2,backup3,
               font_fixes,vel_limited_ptcount,rho_star_fix_applied,cons_avgs,atm_resets,
               failures,
               failures_inhoriz,pointcount_inhoriz,
               error_int_numer/error_int_denom,error_int_denom);
  }
  if( nan_found ) {
    if( GetRefinementLevel(cctkGH) > 6 ) {
      CCTK_ERROR("Found NAN during con2prim driver. See error messages above. ABORTING!");
    }
    else {
      CCTK_VWARN(CCTK_WARN_ALERT,"Found NAN during con2prim driver, but not at finest level. Proceeding with caution...");
    }
  }
}

#include "harm_u2p_util.h"
#include "con2prim_wrapper.h"
