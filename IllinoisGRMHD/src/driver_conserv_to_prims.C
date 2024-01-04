/* We evolve forward in time a set of functions called the
 *  "conservative variables", and any time the conserv's
 *  are updated, we must solve for the primitive variables
 *  (rho, pressure, velocities) using a Newton-Raphson
 *  technique, before reconstructing & evaluating the RHSs
 *  of the MHD equations again.
 *
 * This file contains the driver routine for this Newton-
 *  Raphson solver. Truncation errors in conservative
 *  variables can lead to no physical solutions in
 *  primitive variables. We correct for these errors here
 *  through a number of tricks described in the appendices
 *  of http://arxiv.org/pdf/1112.0568.pdf.
 *
 * This is a wrapper for the 2d solver of Noble et al. See
 *  harm_utoprim_2d.c for references and copyright notice
 *  for that solver. This wrapper was primarily written by
 *  Zachariah Etienne & Yuk Tung Liu, in 2011-2013.
 *
 * For optimal compatibility, this wrapper is licensed under
 *  the GPL v2 or any later version.
 *
 * Note that this code assumes a simple gamma law for the
 *  moment, though it would be easy to extend to a piecewise
 *  polytrope. */

// Standard #include's
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>


#ifdef ENABLE_STANDALONE_IGM_C2P_SOLVER
#include "standalone_conserv_to_prims_main_function.h"
#else
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
#endif

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

  // Initialization of the masks above. Flag meaning:
  //
  // 0 -> Primitive recovery succeeded
  // 1 -> Primitive recovery failed
#pragma omp parallel for
  for(int i=0;i<npoints;i++)
    con2prim_failed_flag[i] = 0;

  // We now add an integer to count the number of
  // points in which the averaging fix is required.
  // We initialize it to a nonzero value so that
  // the while condition below is triggered at least
  // once.
  int num_of_conservative_averagings_needed = 1;
  int cons_avgs = 0;
  int loop_count = 0;
  int nan_found = 0;

  while( num_of_conservative_averagings_needed > 0 ) {

    // Now set the number of conserved averages to zero
    num_of_conservative_averagings_needed = 0;

#pragma omp parallel for reduction(+:failures,vel_limited_ptcount,font_fixes,pointcount,failures_inhoriz,pointcount_inhoriz,error_int_numer,error_int_denom,rho_star_fix_applied,atm_resets,backup1,backup2,backup3,num_of_conservative_averagings_needed,nan_found) schedule(static)
    for(int k=kmin;k<kmax;k++) {
      for(int j=jmin;j<jmax;j++) {
        for(int i=imin;i<imax;i++) {

          int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
          int c2p_fail_flag = con2prim_failed_flag[index];

          // Only attempt a primitive recovery if this is the first
          // attempt at doing so or if the previous attempt failed.
          if( (loop_count == 0) || (c2p_fail_flag != 0) ) {

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
            CCTK_REAL PRIMS[MAXNUMVARS];
            PRIMS[RHOB         ] = rho_b[index];
            PRIMS[PRESSURE     ] = P[index];
            PRIMS[VX           ] = vx[index];
            PRIMS[VY           ] = vy[index];
            PRIMS[VZ           ] = vz[index];
            PRIMS[BX_CENTER    ] = Bx[index];
            PRIMS[BY_CENTER    ] = By[index];
            PRIMS[BZ_CENTER    ] = Bz[index];
            PRIMS[EPSILON      ] = igm_eps[index];
            PRIMS[ENTROPY      ] = igm_entropy[index];

            // Read in conservative variables from gridfunctions
            CCTK_REAL CONSERVS[NUM_CONSERVS],CONSERVS_avg_neighbors[NUM_CONSERVS];
            CONSERVS[RHOSTAR  ] = rho_star[index];
            CONSERVS[STILDEX  ] = mhd_st_x[index];
            CONSERVS[STILDEY  ] = mhd_st_y[index];
            CONSERVS[STILDEZ  ] = mhd_st_z[index];
            CONSERVS[TAUENERGY] = tau     [index];
            CONSERVS[YESTAR   ] = Ye_star [index];
            CONSERVS[ENTSTAR  ] = S_star  [index];

            // Check if we need to perform the new conservative averaging fix
            if( c2p_fail_flag != 0 ) {

              // Average the conserved variables using the neighboring values
              int __attribute__((unused)) neighbors = con2prim_average_neighbor_conservatives(cctkGH,i,j,k,index,cctk_lsh,con2prim_failed_flag,
                                                                                              rho_star,mhd_st_x,mhd_st_y,mhd_st_z,tau,Ye_star,S_star,
                                                                                              CONSERVS_avg_neighbors);
	      if( neighbors > 0 ) {
		// We will attempt 4 fixes after the first attempt fails. These
		// new attempts will use the following conservatives:
		//
		// C = 0.75*C + 0.25*<C_neighbors>
		// C = 0.50*C + 0.50*<C_neighbors>
		// C = 0.25*C + 0.75*<C_neighbors>
		// C = 0.00*C + 1.00*<C_neighbors>
		//
		// Set the weights
		CCTK_REAL C_weight           = 0.25*c2p_fail_flag;
		CCTK_REAL one_minus_C_weight = 1.0 - C_weight;

		// Then update the conservs
		for(int which_con=0;which_con<NUM_CONSERVS;which_con++)
		  CONSERVS[which_con] = CONSERVS[which_con]*one_minus_C_weight + CONSERVS_avg_neighbors[which_con] * C_weight;

		// Assume we won't need to fix this point again, so
		// decrement the number of points that require the
		// conservative averaging fix
		num_of_conservative_averagings_needed--;
              }
	      else {
		// Probably should terminate in this case, but let us reset to ATM for now.
		// Simplest way is to set rhostar to a negative value.
		CONSERVS[RHOSTAR] = -1.0;
		num_of_conservative_averagings_needed--;
		con2prim_failed_flag[index] = 0;
		atm_resets++;
              }

            }

            // Tabulated EOS quantities
            if( eos.is_Tabulated ) {
              // Primitives
              PRIMS[YEPRIM     ] = igm_Ye[index];
              PRIMS[TEMPERATURE] = igm_temperature[index];
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
            st_x_flux[index]        = CONSERVS[STILDEX  ];
            st_y_flux[index]        = CONSERVS[STILDEY  ];
            st_z_flux[index]        = CONSERVS[STILDEZ  ];
            tau_flux[index]         = CONSERVS[TAUENERGY];
            Ye_star_flux[index]     = CONSERVS[YESTAR   ];
            S_star_flux[index]      = CONSERVS[ENTSTAR  ];

            CCTK_REAL rho_star_orig = CONSERVS[RHOSTAR  ];
            CCTK_REAL mhd_st_x_orig = CONSERVS[STILDEX  ];
            CCTK_REAL mhd_st_y_orig = CONSERVS[STILDEY  ];
            CCTK_REAL mhd_st_z_orig = CONSERVS[STILDEZ  ];
            CCTK_REAL tau_orig      = CONSERVS[TAUENERGY];
            CCTK_REAL Ye_star_orig  = 0.0;
            CCTK_REAL S_star_orig   = 0.0;
            if( eos.is_Tabulated) {
              Ye_star_orig          = CONSERVS[YESTAR   ];
            }
            if( eos.evolve_entropy ) {
              S_star_orig           = CONSERVS[ENTSTAR  ];
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
              con2prim_failed_flag[index] += 1;
              if( con2prim_failed_flag[index] > 4 ) {
                // Sigh, reset to atmosphere
                reset_prims_to_atmosphere( eos, PRIMS );
                atm_resets++;
                // Then flag this point as a "success"
                check = 0;
                con2prim_failed_flag[index] = 0;
                if( eos.is_Hybrid ) {
                  CCTK_VInfo(CCTK_THORNSTRING,"Couldn't find root from: %e %e %e %e %e, rhob approx=%e, rho_b_atm=%e, Bx=%e, By=%e, Bz=%e, gij_phys=%e %e %e %e %e %e, alpha=%e",
                             tau_orig,rho_star_orig,mhd_st_x_orig,mhd_st_y_orig,mhd_st_z_orig,rho_star_orig/METRIC_LAP_PSI4[PSI6],eos.rho_atm,PRIMS[BX_CENTER],PRIMS[BY_CENTER],PRIMS[BZ_CENTER],METRIC_PHYS[GXX],METRIC_PHYS[GXY],METRIC_PHYS[GXZ],METRIC_PHYS[GYY],METRIC_PHYS[GYZ],METRIC_PHYS[GZZ],METRIC_LAP_PSI4[LAPSE]);
                }
                else if( eos.is_Tabulated ) {
                  CCTK_VInfo(CCTK_THORNSTRING,"Couldn't find root from: %e %e %e %e %e %e %e, rhob approx=%e, rho_b_atm=%e, Bx=%e, By=%e, Bz=%e, gij_phys=%e %e %e %e %e %e, alpha=%e",
                             tau_orig,rho_star_orig,mhd_st_x_orig,mhd_st_y_orig,mhd_st_z_orig,Ye_star_orig,S_star_orig,rho_star_orig/METRIC_LAP_PSI4[PSI6],eos.rho_atm,PRIMS[BX_CENTER],PRIMS[BY_CENTER],PRIMS[BZ_CENTER],METRIC_PHYS[GXX],METRIC_PHYS[GXY],METRIC_PHYS[GXZ],METRIC_PHYS[GYY],METRIC_PHYS[GYZ],METRIC_PHYS[GZZ],METRIC_LAP_PSI4[LAPSE]);
                }
              }
              else {
                // Increment the number of gridpoints which will need the average fix
                num_of_conservative_averagings_needed++;
              }
            }

            if( check == 0 ) {
              //--------------------------------------------------
              //---------- Primitive recovery succeeded ----------
              //--------------------------------------------------
              // Enforce limits on primitive variables and recompute conservatives.
              static const int already_computed_physical_metric_and_inverse=1;
              CCTK_REAL TUPMUNU[10],TDNMUNU[10];
              IllinoisGRMHD_enforce_limits_on_primitives_and_recompute_conservs(already_computed_physical_metric_and_inverse,PRIMS,stats,eos,METRIC,g4dn,g4up, TUPMUNU,TDNMUNU,CONSERVS);

              rho_star   [index] = CONSERVS[RHOSTAR  ];
              mhd_st_x   [index] = CONSERVS[STILDEX  ];
              mhd_st_y   [index] = CONSERVS[STILDEY  ];
              mhd_st_z   [index] = CONSERVS[STILDEZ  ];
              tau        [index] = CONSERVS[TAUENERGY];

              // Set primitives, and/or provide a better guess.
              rho_b      [index] = PRIMS[RHOB        ];
              P          [index] = PRIMS[PRESSURE    ];
              vx         [index] = PRIMS[VX          ];
              vy         [index] = PRIMS[VY          ];
              vz         [index] = PRIMS[VZ          ];
              igm_eps    [index] = PRIMS[EPSILON     ];
              igm_entropy[index] = PRIMS[ENTROPY     ];

              // Tabulated EOS quantities
              if( eos.is_Tabulated ) {
                // Primitives
                igm_Ye[index]          = PRIMS[YEPRIM     ];
                igm_temperature[index] = PRIMS[TEMPERATURE];
                // Conservatives
                Ye_star[index]         = CONSERVS[YESTAR  ];
              }

              // Entropy evolution quantities
              if( eos.evolve_entropy ) {
                S_star[index]          = CONSERVS[ENTSTAR ];
              }

              if(update_Tmunu) {
                int ww=0;
                eTtt[index] = TDNMUNU[ww++];
                eTtx[index] = TDNMUNU[ww++];
                eTty[index] = TDNMUNU[ww++];
                eTtz[index] = TDNMUNU[ww++];
                eTxx[index] = TDNMUNU[ww++];
                eTxy[index] = TDNMUNU[ww++];
                eTxz[index] = TDNMUNU[ww++];
                eTyy[index] = TDNMUNU[ww++];
                eTyz[index] = TDNMUNU[ww++];
                eTzz[index] = TDNMUNU[ww  ];
              }

              //Now we compute the difference between original & new conservatives, for diagnostic purposes:
              error_int_numer += fabs(tau[index] - tau_orig) + fabs(rho_star[index] - rho_star_orig) +
                fabs(mhd_st_x[index] - mhd_st_x_orig) + fabs(mhd_st_y[index] - mhd_st_y_orig) + fabs(mhd_st_z[index] - mhd_st_z_orig);
              error_int_denom += tau_orig + rho_star_orig + fabs(mhd_st_x_orig) + fabs(mhd_st_y_orig) + fabs(mhd_st_z_orig);

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
          } // if( (loop_count == 0) || (c2p_fail_flag != 0) )
        } // for(int i=imin;i<imax;i++)
      } // for(int j=jmin;j<jmax;j++)
    } // for(int k=kmin;k<kmax;k++)

    if( num_of_conservative_averagings_needed > cons_avgs ) cons_avgs = num_of_conservative_averagings_needed;

    // Increment the loop counter. This counter is used to avoid
    // attempting to recover primitives at points in which we
    // have already done so successfully.
    loop_count++;

  } // while( num_of_conservative_averagings_needed > 0 )

  // fclose(c2pmaskfile);

  if(CCTK_Equals(verbose, "essential") || CCTK_Equals(verbose, "essential+iteration output")) {
    CCTK_VInfo(CCTK_THORNSTRING,"C2P: Lev: %d NumPts= %d | Fixes: BU: %d %d %d Font= %d VL= %d rho*= %d AVG= %d ATM= %d | Failures: %d InHoriz= %d / %d | Error: %.3e, ErrDenom: %.3e",
               (int)GetRefinementLevel(cctkGH),pointcount,
               backup1,backup2,backup3,
               font_fixes,vel_limited_ptcount,rho_star_fix_applied,cons_avgs,atm_resets,
               failures,
               failures_inhoriz,pointcount_inhoriz,
               error_int_numer/error_int_denom,error_int_denom);
  }
  // if( nan_found ) {
  //   if( GetRefinementLevel(cctkGH) > 6 ) {
  //     CCTK_ERROR("Found NAN during con2prim driver. See error messages above. ABORTING!");
  //   }
  //   else {
  //     CCTK_VWARN(CCTK_WARN_ALERT,"Found NAN during con2prim driver, but not at finest level. Proceeding with caution...");
  //   }
  // }

  // Very useful con2prim debugger. If the primitives (con2prim) solver fails, this will output all data needed to
  //     debug where and why the solver failed. Strongly suggested for experimenting with new fixes.
  if( ( (conserv_to_prims_debug==1) && (error_int_numer/error_int_denom > 0.05) ) ||
      ( atm_resets != 0 ) ) {

    ofstream myfile;
    char filename[100];
    srand(time(NULL));
    sprintf(filename,"primitives_debug-%e.dat",error_int_numer/error_int_denom);
    //Alternative, for debugging purposes as well:
    //srand(time(NULL));
    //sprintf(filename,"primitives_debug-%d.dat",rand());
    myfile.open (filename, ios::out | ios::binary);
    //myfile.open ("data.bin", ios::out | ios::binary);

    // This checker value will be printed last, and will
    // allow us to make sure we have read the dump file
    // correctly when debugging it.
    int checker=1063;

    // Grid information
    int fullsize=cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];
    myfile.write((char*)cctk_lsh,                             3*sizeof(int));
    myfile.write((char*)x,                           (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)y,                           (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)z,                           (fullsize)*sizeof(CCTK_REAL));

    // IllinoisGRMHD parameters
    myfile.write((char*)&GAMMA_SPEED_LIMIT,                   1*sizeof(CCTK_REAL));
    myfile.write((char*)&rho_b_max,                           1*sizeof(CCTK_REAL));
    myfile.write((char*)&rho_b_atm,                           1*sizeof(CCTK_REAL));
    myfile.write((char*)&tau_atm,                             1*sizeof(CCTK_REAL));
    myfile.write((char*)&Psi6threshold,                       1*sizeof(CCTK_REAL));
    myfile.write((char*)&update_Tmunu,                        1*sizeof(bool));
    myfile.write((char*)&neos,                                1*sizeof(int));
    myfile.write((char*)&Gamma_th,                            1*sizeof(CCTK_REAL));
    myfile.write((char*)&K_ppoly_tab0,                        1*sizeof(CCTK_REAL));
    myfile.write((char*)Gamma_ppoly_tab_in,                neos*sizeof(CCTK_REAL));
    myfile.write((char*)rho_ppoly_tab_in,              (neos-1)*sizeof(CCTK_REAL));
    myfile.write((char*)&igm_Ye_atm,                          1*sizeof(CCTK_REAL));
    myfile.write((char*)&igm_T_atm,                           1*sizeof(CCTK_REAL));
    myfile.write((char*)&igm_eos_table_ceiling_safety_factor, 1*sizeof(CCTK_REAL));
    myfile.write((char*)&igm_eos_table_floor_safety_factor,   1*sizeof(CCTK_REAL));
    myfile.write((char*)&igm_eos_root_finding_precision,      1*sizeof(CCTK_REAL));
    myfile.write((char*)&igm_eos_root_finding_precision,      1*sizeof(CCTK_REAL));
    myfile.write((char*)&igm_evolve_temperature,              1*sizeof(bool));
    myfile.write((char*)&igm_evolve_entropy,                  1*sizeof(bool));

    // Failure checker gridfunction
    myfile.write((char*)failure_checker,             (fullsize)*sizeof(CCTK_REAL));

    // Energy momentum tensor
    myfile.write((char*)eTtt,                        (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)eTtx,                        (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)eTty,                        (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)eTtz,                        (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)eTxx,                        (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)eTxy,                        (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)eTxz,                        (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)eTyy,                        (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)eTyz,                        (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)eTzz,                        (fullsize)*sizeof(CCTK_REAL));

    // Metric quantities
    myfile.write((char*)alp,                         (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gxx,                         (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gxy,                         (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gxz,                         (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gyy,                         (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gyz,                         (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gzz,                         (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)psi_bssn,                    (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)phi_bssn,                    (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtxx,                        (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtxy,                        (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtxz,                        (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtyy,                        (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtyz,                        (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtzz,                        (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtupxx,                      (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtupxy,                      (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtupxz,                      (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtupyy,                      (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtupyz,                      (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtupzz,                      (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)lapm1,                       (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)betax,                       (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)betay,                       (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)betaz,                       (fullsize)*sizeof(CCTK_REAL));

    // Original conservative variables (we used the flux gridfunctions to store these)
    myfile.write((char*)rho_star_flux,               (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)tau_flux,                    (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)st_x_flux,                   (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)st_y_flux,                   (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)st_z_flux,                   (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)Ye_star_flux,                (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)S_star_flux,                 (fullsize)*sizeof(CCTK_REAL));

    // Magnetic fields
    myfile.write((char*)Bx,                          (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)By,                          (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)Bz,                          (fullsize)*sizeof(CCTK_REAL));

    // Primitive variables
    myfile.write((char*)rho_b,                       (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)P,                           (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)vx,                          (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)vy,                          (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)vz,                          (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)igm_Ye,                      (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)igm_temperature,             (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)igm_eps,                     (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)igm_entropy,                 (fullsize)*sizeof(CCTK_REAL));

    // Checker value
    myfile.write((char*)&checker,                             1*sizeof(int));

    // All done! Close the file.
    myfile.close();
    CCTK_VInfo(CCTK_THORNSTRING,"Finished writing %s",filename);
  }

#ifdef ENABLE_STANDALONE_IGM_C2P_SOLVER
  return 0; // int main() requires an integer be returned
#endif

}

#include "harm_u2p_util.h"
#include "con2prim_wrapper.h"
