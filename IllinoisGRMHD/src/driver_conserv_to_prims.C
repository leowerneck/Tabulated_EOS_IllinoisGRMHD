#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "GRHayLib.h"
#include "Symmetry.h"

#include "IllinoisGRMHD_headers.h"
//#include "con2prim_headers.h"
#include "inlined_functions.h"
#include "apply_tau_floor__enforce_limits_on_primitives_and_recompute_conservs.C"

static inline int GRHayLMHD_local_avg(
      const cGH *restrict cctkGH,
      const int i,
      const int j,
      const int k,
      const int weight,
      const CCTK_INT  *restrict needs_average,
      const CCTK_REAL *restrict rho_star,
      const CCTK_REAL *restrict tau,
      const CCTK_REAL *restrict mhd_st_x,
      const CCTK_REAL *restrict mhd_st_y,
      const CCTK_REAL *restrict mhd_st_z,
      const CCTK_REAL *restrict S_star,
      const CCTK_REAL *restrict Ye_star,
      ghl_conservative_quantities *restrict cons) {

  // This sets how many neighbors must exist to consider
  // the averaging scheme a success. 1 means we only
  // require at least 1 valid neighboring point to continue
  // with the averaging method.
  const int min_neighbors = 1;

  const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
  const double wfac = weight/4.0;

  const int iavg_min = MAX(0, i-1);
  const int javg_min = MAX(0, j-1);
  const int kavg_min = MAX(0, k-1);
  const int iavg_max = MIN(cctkGH->cctk_lsh[0], i+2);
  const int javg_max = MIN(cctkGH->cctk_lsh[1], j+2);
  const int kavg_max = MIN(cctkGH->cctk_lsh[2], k+2);
  int num_avg = 0;

  const double cfac = 1.0 - wfac;
  cons->rho     = cfac*rho_star[index];
  cons->tau     = cfac*tau[index];
  cons->SD[0]   = cfac*mhd_st_x[index];
  cons->SD[1]   = cfac*mhd_st_y[index];
  cons->SD[2]   = cfac*mhd_st_z[index];
  cons->entropy = cfac*S_star[index];
  cons->Y_e     = cfac*Ye_star[index];

  double rhotmp, tautmp, Sxtmp, Sytmp, Sztmp, enttmp, yetmp;
  rhotmp = tautmp = Sxtmp = Sytmp = Sztmp = enttmp = yetmp = 0;

  for(int kavg=kavg_min; kavg<kavg_max; kavg++) {
    for(int javg=javg_min; javg<javg_max; javg++) {
      for(int iavg=iavg_min; iavg<iavg_max; iavg++) {
        // We're only averaging over neighbors, so skip this point.
        // Also, we skip over any points that also failed.
        const int inavg = CCTK_GFINDEX3D(cctkGH,iavg,javg,kavg);
        if((i==iavg && j==javg && k==kavg) || needs_average[inavg])
          continue;
        rhotmp += rho_star[inavg];
        tautmp += tau[inavg];
        Sxtmp  += mhd_st_x[inavg];
        Sytmp  += mhd_st_y[inavg];
        Sztmp  += mhd_st_z[inavg];
        enttmp += S_star[inavg];
        yetmp  += Ye_star[inavg];
        num_avg++;
      }
    }
  }

  if(num_avg < min_neighbors)
    return 1;

  cons->rho     += wfac*rhotmp;
  cons->tau     += wfac*tautmp;
  cons->SD[0]   += wfac*Sxtmp;
  cons->SD[1]   += wfac*Sytmp;
  cons->SD[2]   += wfac*Sztmp;
  cons->entropy += wfac*enttmp;
  cons->Y_e     += wfac*yetmp;

  // If weight=4, then the central point isn't part of the average
  num_avg += (weight != 4);

  cons->rho     /= num_avg;
  cons->tau     /= num_avg;
  cons->SD[0]   /= num_avg;
  cons->SD[1]   /= num_avg;
  cons->SD[2]   /= num_avg;
  cons->entropy /= num_avg;
  cons->Y_e     /= num_avg;

  return 0;
}

void IllinoisGRMHD_conserv_to_prims(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int imax = cctk_lsh[0];
  const int jmax = cctk_lsh[1];
  const int kmax = cctk_lsh[2];

  if(CCTK_EQUALS(Symmetry,"equatorial")) {
    // SET SYMMETRY GHOSTZONES ON ALL CONSERVATIVE VARIABLES!
    int ierr=0;
    ierr+=CartSymGN(cctkGH,"GRHayLMHD::grmhd_conservatives");
    // FIXME: UGLY. Filling metric ghostzones is needed for, e.g., Cowling runs.
    ierr+=CartSymGN(cctkGH,"lapse::lapse_vars");
    ierr+=CartSymGN(cctkGH,"bssn::BSSN_vars");
    ierr+=CartSymGN(cctkGH,"bssn::BSSN_AH");
    ierr+=CartSymGN(cctkGH,"shift::shift_vars");
    if(ierr!=0) CCTK_VERROR("Error with setting equatorial symmetries in con2prim.");
  }

  // Diagnostic variables.
  int failures = 0;
  int vel_limited_ptcount = 0;
  int rho_star_fix_applied = 0;
  int pointcount = 0;
  int navg_total = 0;
  int failures_inhoriz = 0;
  int pointcount_inhoriz = 0;
  int backup0 = 0;
  int backup1 = 0;
  int backup2 = 0;
  double error_rho_numer = 0;
  double error_tau_numer = 0;
  double error_Sx_numer = 0;
  double error_Sy_numer = 0;
  double error_Sz_numer = 0;
  double error_ent_numer = 0;
  double error_Ye_numer = 0;

  double error_rho_denom = 0;
  double error_tau_denom = 0;
  double error_Sx_denom = 0;
  double error_Sy_denom = 0;
  double error_Sz_denom = 0;
  double error_ent_denom = 0;
  double error_Ye_denom = 0;
  int n_iter = 0;
  int n_avg = 0;

#pragma omp parallel for reduction(+: \
      pointcount, backup0, backup1, backup2, vel_limited_ptcount, rho_star_fix_applied, failures, failures_inhoriz, pointcount_inhoriz, n_iter, \
      error_rho_numer, error_tau_numer, error_Sx_numer, error_Sy_numer, error_Sz_numer, error_ent_numer, error_Ye_numer, \
      error_rho_denom, error_tau_denom, error_Sx_denom, error_Sy_denom, error_Sz_denom, error_ent_denom, error_Ye_denom) schedule(static)
  for(int k=0; k<kmax; k++) {
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        double local_failure_checker = 0;

        ghl_con2prim_diagnostics diagnostics;
        ghl_initialize_diagnostics(&diagnostics);

        // Read in ADM metric quantities from gridfunctions and
        // set auxiliary and ADM metric quantities
        ghl_metric_quantities ADM_metric;
        ghl_enforce_detgtij_and_initialize_ADM_metric(
              alp[index],
              betax[index], betay[index], betaz[index],
              gxx[index], gxy[index], gxz[index],
              gyy[index], gyz[index], gzz[index],
              &ADM_metric);

        ghl_ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        // Read in primitive variables from gridfunctions
        // The code has only ever been tested using the default GRHayL guess,
        // so using the previous timelevel as an initial guess would need to
        // be implemented here.
        ghl_primitive_quantities prims;
        prims.BU[0] = Bx[index] * ONE_OVER_SQRT_4PI;
        prims.BU[1] = By[index] * ONE_OVER_SQRT_4PI;
        prims.BU[2] = Bz[index] * ONE_OVER_SQRT_4PI;

        // Read in conservative variables from gridfunctions
        ghl_conservative_quantities cons, cons_orig;
        cons.rho     = rho_star[index];
        cons.tau     = tau[index];
        cons.SD[0]   = mhd_st_x[index];
        cons.SD[1]   = mhd_st_y[index];
        cons.SD[2]   = mhd_st_z[index];
        cons.entropy = S_star[index];
        cons.Y_e     = Ye_star[index];

        // Here we save the original values of conservative variables in cons_orig for debugging purposes.
        cons_orig = cons;

        if(isnan(cons.rho*cons.tau*cons.SD[0]*cons.SD[1]*cons.SD[2]*cons.entropy*cons.Y_e*prims.BU[0]*prims.BU[1]*prims.BU[2])) {
          needs_average[index] = 1;
          n_avg++;
          continue;
        }

        /************* Main conservative-to-primitive logic ************/
        if(cons.rho>0.0) {
          // declare some variables for the C2P routine.
          ghl_conservative_quantities cons_undens;

          // Set the conserved variables required by the con2prim routine
          ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);

          prims.Y_e         = cons.Y_e / cons.rho;
          prims.rho         = rho_b[index];
          prims.temperature = igm_temperature[index];
          prims.vU[0]       = vx[index];
          prims.vU[1]       = vy[index];
          prims.vU[2]       = vz[index];
          if(isnan(prims.rho*prims.temperature*prims.vU[0]*prims.vU[1]*prims.vU[2])) {
            ghl_guess_primitives(ghl_eos, &ADM_metric, &cons, &prims);
          }


          /************* Conservative-to-primitive recovery ************/
          const int check = ghl_con2prim_multi_method(
                ghl_params, ghl_eos, &ADM_metric, &metric_aux,
                &cons_undens, &prims, &diagnostics);

          if(check || isnan(prims.rho*prims.press*prims.eps*prims.vU[0]*prims.vU[1]*prims.vU[2]*prims.entropy*prims.Y_e*prims.temperature)) {
            needs_average[index] = 1;
            n_avg++;
            continue;
          }
        } else {
          local_failure_checker += 1;
          ghl_set_prims_to_constant_atm(ghl_eos, &prims);
          rho_star_fix_applied++;
        } // if rho_star>0
        /***************************************************************/

        //--------------------------------------------------
        //---------- Primitive recovery succeeded ----------
        //--------------------------------------------------
        // Enforce limits on primitive variables and recompute conservatives.
        diagnostics.speed_limited += ghl_enforce_primitive_limits_and_compute_u0(
              ghl_params, ghl_eos, &ADM_metric, &prims);

        ghl_stress_energy Tmunu;
        ghl_compute_conservs_and_Tmunu(
                                       &ADM_metric, &metric_aux, &prims, &cons, &Tmunu);

        rho_b[index]           = prims.rho;
        P[index]               = prims.press;
        igm_eps[index]         = prims.eps;
        vx[index]              = prims.vU[0];
        vy[index]              = prims.vU[1];
        vz[index]              = prims.vU[2];
        igm_entropy[index]     = prims.entropy;
        igm_Ye[index]          = prims.Y_e;
        igm_temperature[index] = prims.temperature;

        rho_star[index] = cons.rho;
        tau[index]      = cons.tau;
        mhd_st_x[index]  = cons.SD[0];
        mhd_st_y[index]  = cons.SD[1];
        mhd_st_z[index]  = cons.SD[2];
        S_star[index] = cons.entropy;
        Ye_star[index]  = cons.Y_e;

        eTtt[index] = Tmunu.T4[0][0];
        eTtx[index] = Tmunu.T4[0][1];
        eTty[index] = Tmunu.T4[0][2];
        eTtz[index] = Tmunu.T4[0][3];
        eTxx[index] = Tmunu.T4[1][1];
        eTxy[index] = Tmunu.T4[1][2];
        eTxz[index] = Tmunu.T4[1][3];
        eTyy[index] = Tmunu.T4[2][2];
        eTyz[index] = Tmunu.T4[2][3];
        eTzz[index] = Tmunu.T4[3][3];

        needs_average[index] = 0;

        //Now we compute the difference between original & new conservatives, for diagnostic purposes:
        error_rho_numer += fabs(cons.rho - cons_orig.rho);
        error_tau_numer += fabs(cons.tau - cons_orig.tau);
        error_Sx_numer  += fabs(cons.SD[0] - cons_orig.SD[0]);
        error_Sy_numer  += fabs(cons.SD[1] - cons_orig.SD[1]);
        error_Sz_numer  += fabs(cons.SD[2] - cons_orig.SD[2]);
        error_ent_numer += fabs(cons.entropy - cons_orig.entropy);
        error_Ye_numer  += fabs(cons.Y_e - cons_orig.Y_e);
        error_rho_denom += cons_orig.tau;
        error_tau_denom += cons_orig.rho;
        error_Sx_denom  += fabs(cons_orig.SD[0]);
        error_Sy_denom  += fabs(cons_orig.SD[1]);
        error_Sz_denom  += fabs(cons_orig.SD[2]);
        error_ent_denom += cons_orig.entropy;
        error_Ye_denom  += cons_orig.Y_e;

        pointcount++;
        if(diagnostics.speed_limited) {
          local_failure_checker += 10;
          vel_limited_ptcount++;
        }
        backup0 += diagnostics.backup[0];
        backup1 += diagnostics.backup[1];
        backup2 += diagnostics.backup[2];
        n_iter += diagnostics.n_iter;
        failure_checker[index] = local_failure_checker
                               + 1000*diagnostics.backup[0]
                               + 10000*diagnostics.tau_fix
                               + 100000*diagnostics.Stilde_fix;
      }
    }
  }

  if(n_avg>0) {
    CCTK_VINFO("Con2Prim failed for %d points. Beginning averaging method...", n_avg);
    int ind_vals[n_avg], i_vals[n_avg], j_vals[n_avg], k_vals[n_avg];
    int counter = 0;
    int avg_weight = 1;
    for(int k=0; k<kmax; k++) {
      for(int j=0; j<jmax; j++) {
        for(int i=0; i<imax; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
          if(needs_average[index]) {
            ind_vals[counter] = index;
            i_vals[counter] = i;
            j_vals[counter] = j;
            k_vals[counter] = k;
            counter++;
          }
        }
      }
    }
    while(n_avg > 0 && avg_weight < 5) {
      int new_avg = 0;
      for(int iter=0; iter<n_avg; iter++) {
        const int i = i_vals[iter];
        const int j = j_vals[iter];
        const int k = k_vals[iter];
        const int index = ind_vals[iter];

        ghl_con2prim_diagnostics diagnostics;
        ghl_initialize_diagnostics(&diagnostics);

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
        prims.BU[0] = Bx[index];
        prims.BU[1] = By[index];
        prims.BU[2] = Bz[index];

        ghl_conservative_quantities cons, cons_orig;
        const int avg_fail = GRHayLMHD_local_avg(
              cctkGH, i, j, k, avg_weight,
              needs_average, rho_star, tau,
              mhd_st_x, mhd_st_y, mhd_st_z,
              S_star, Ye_star, &cons);

        if(avg_fail ||
           isnan(cons.rho*cons.tau*cons.SD[0]*cons.SD[1]*cons.SD[2]*cons.entropy*cons.Y_e*prims.BU[0]*prims.BU[1]*prims.BU[2])) {
          ind_vals[new_avg] = index;
          i_vals[new_avg] = i;
          j_vals[new_avg] = j;
          k_vals[new_avg] = k;
          new_avg++;
          continue;
        }

        cons_orig.rho     = rho_star[index];
        cons_orig.tau     = tau[index];
        cons_orig.SD[0]   = mhd_st_x[index];
        cons_orig.SD[1]   = mhd_st_y[index];
        cons_orig.SD[2]   = mhd_st_z[index];
        cons_orig.entropy = S_star[index];
        cons_orig.Y_e     = Ye_star[index];

        /************* Main conservative-to-primitive logic ************/
        if(cons.rho>0.0) {
          // declare some variables for the C2P routine.
          ghl_conservative_quantities cons_undens;

          // Set the conserved variables required by the con2prim routine
          ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);

          /************* Conservative-to-primitive recovery ************/
          const int check = ghl_con2prim_multi_method(
                ghl_params, ghl_eos, &ADM_metric, &metric_aux,
                &cons_undens, &prims, &diagnostics);

          //Check for NAN!
          if(check || isnan(prims.rho*prims.press*prims.eps*prims.vU[0]*prims.vU[1]*prims.vU[2]*prims.entropy*prims.Y_e*prims.temperature) ) {
            ind_vals[new_avg] = index;
            i_vals[new_avg] = i;
            j_vals[new_avg] = j;
            k_vals[new_avg] = k;
            new_avg++;
            continue;
          }
        } else {
          failure_checker[index] += 1;
          ghl_set_prims_to_constant_atm(ghl_eos, &prims);
          rho_star_fix_applied++;
        } // if rho_star>0

        diagnostics.speed_limited += ghl_enforce_primitive_limits_and_compute_u0(
              ghl_params, ghl_eos, &ADM_metric, &prims);
        ghl_compute_conservs(
              &ADM_metric, &metric_aux, &prims, &cons);

        rho_b[index]       = prims.rho;
        P[index]           = prims.press;
        eps[index]         = prims.eps;
        vx[index]          = prims.vU[0];
        vy[index]          = prims.vU[1];
        vz[index]          = prims.vU[2];
        igm_entropy[index]     = prims.entropy;
        igm_Ye[index]         = prims.Y_e;
        igm_temperature[index] = prims.temperature;

        rho_star[index] = cons.rho;
        tau[index]      = cons.tau;
        mhd_st_x[index]  = cons.SD[0];
        mhd_st_y[index]  = cons.SD[1];
        mhd_st_z[index]  = cons.SD[2];
        S_star[index] = cons.entropy;
        Ye_star[index]  = cons.Y_e;

        needs_average[index] = 0;
        navg_total++;

        //Now we compute the difference between original & new conservatives, for diagnostic purposes:
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

        pointcount++;
        if(diagnostics.speed_limited) {
          failure_checker[index] += 10;
          vel_limited_ptcount++;
        }
        backup0 += diagnostics.backup[0];
        backup1 += diagnostics.backup[1];
        backup2 += diagnostics.backup[2];
        n_iter += diagnostics.n_iter;
        failure_checker[index] = 1000*diagnostics.backup[0]
                               + 10000*diagnostics.tau_fix
                               + 100000*diagnostics.Stilde_fix;
      } // iter
      n_avg = new_avg;
      avg_weight++;
    } // while n_avg
    if(n_avg != 0) {
      for(int iter=0; iter<n_avg; iter++) {
        const int index = ind_vals[iter];
        //--------------------------------------------------
        //----------- Primitive recovery failed ------------
        //--------------------------------------------------
        failure_checker[index] += 100;

        ghl_con2prim_diagnostics diagnostics;
        ghl_initialize_diagnostics(&diagnostics);

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
        prims.BU[0] = Bx[index];
        prims.BU[1] = By[index];
        prims.BU[2] = Bz[index];

        ghl_conservative_quantities cons, cons_orig;
        cons_orig.rho     = rho_star[index];
        cons_orig.tau     = tau[index];
        cons_orig.SD[0]   = mhd_st_x[index];
        cons_orig.SD[1]   = mhd_st_y[index];
        cons_orig.SD[2]   = mhd_st_z[index];
        cons_orig.entropy = S_star[index];
        cons_orig.Y_e     = Ye_star[index];

        ghl_set_prims_to_constant_atm(ghl_eos, &prims);

        failures++;
        if(ADM_metric.sqrt_detgamma > ghl_params->psi6threshold) {
          failures_inhoriz++;
          pointcount_inhoriz++;
        }
        CCTK_VINFO("***********************************************************"
                   "Con2Prim and averaging backups failed! Resetting to atmosphere...\n"
                   "position = %e %e %e\n"
                   "lapse, shift = %e, %e, %e, %e\n"
                   "gij = %e, %e, %e, %e, %e, %e\n"
                   "B^i = %e, %e, %e\n"
                   "rho_*, ~tau, ~S_{i}: %e, %e, %e, %e, %e\n"
                   "~DS, ~DYe: %e, %e\n"
                   "***********************************************************",
                   x[index], y[index], z[index],
                   ADM_metric.lapse, ADM_metric.betaU[0], ADM_metric.betaU[1], ADM_metric.betaU[2],
                   ADM_metric.gammaDD[0][0], ADM_metric.gammaDD[0][1], ADM_metric.gammaDD[0][2],
                   ADM_metric.gammaDD[1][1], ADM_metric.gammaDD[1][2], ADM_metric.gammaDD[2][2],
                   prims.BU[0], prims.BU[1], prims.BU[2],
                   cons.rho, cons.tau, cons.SD[0], cons.SD[1], cons.SD[2],
                   cons.entropy, cons.Y_e);

        diagnostics.speed_limited += ghl_enforce_primitive_limits_and_compute_u0(
              ghl_params, ghl_eos, &ADM_metric, &prims);
        ghl_compute_conservs(
              &ADM_metric, &metric_aux, &prims, &cons);

        rho_b[index]       = prims.rho;
        P[index]           = prims.press;
        eps[index]         = prims.eps;
        vx[index]          = prims.vU[0];
        vy[index]          = prims.vU[1];
        vz[index]          = prims.vU[2];
        igm_entropy[index]     = prims.entropy;
        igm_Ye[index]         = prims.Y_e;
        igm_temperature[index] = prims.temperature;

        rho_star[index] = cons.rho;
        tau[index]      = cons.tau;
        mhd_st_x[index]  = cons.SD[0];
        mhd_st_y[index]  = cons.SD[1];
        mhd_st_z[index]  = cons.SD[2];
        S_star[index] = cons.entropy;
        Ye_star[index]  = cons.Y_e;

        //Now we compute the difference between original & new conservatives, for diagnostic purposes:
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

        pointcount++;
        if(diagnostics.speed_limited) {
          failure_checker[index] += 10;
          vel_limited_ptcount++;
        }
        backup0 += diagnostics.backup[0];
        backup1 += diagnostics.backup[1];
        backup2 += diagnostics.backup[2];
        n_iter += diagnostics.n_iter;
        failure_checker[index] = 1000*diagnostics.backup[0]
                               + 10000*diagnostics.tau_fix
                               + 100000*diagnostics.Stilde_fix;
      }
    }
  } // if n_avg

  const double rho_error     = (error_rho_denom==0) ? error_rho_numer : error_rho_numer/error_rho_denom;
  const double tau_error     = (error_tau_denom==0) ? error_tau_numer : error_tau_numer/error_tau_denom;
  const double Sx_error      = (error_Sx_denom==0)  ? error_Sx_numer  : error_Sx_numer/error_Sx_denom;
  const double Sy_error      = (error_Sy_denom==0)  ? error_Sy_numer  : error_Sy_numer/error_Sy_denom;
  const double Sz_error      = (error_Sz_denom==0)  ? error_Sz_numer  : error_Sz_numer/error_Sz_denom;
  const double entropy_error = (error_ent_denom==0) ? error_ent_numer : error_ent_numer/error_ent_denom;
  const double Ye_error      = (error_Ye_denom==0)  ? error_Ye_numer  : error_Ye_numer/error_Ye_denom;
  /*
    Failure checker decoder:
       1: atmosphere reset when rho_star < 0
      10: Limiting velocity u~ after C2P/Font Fix or v in ghl_enforce_primitive_limits_and_compute_u0
     100: Both C2P and Font Fix failed
      1k: backups used
     10k: tau~ was reset in ghl_apply_conservative_limits
    100k: S~ was reset in ghl_apply_conservative_limits
  */
  if(CCTK_Equals(verbose, "yes") || CCTK_Equals(verbose, "essential") || CCTK_Equals(verbose, "essential+iteration output")) {
    CCTK_VINFO("C2P: Iter. # %d, Lev: %d NumPts= %d | Backups: %d %d %d | Fixes: VL= %d rho*= %d | Failures: %d InHoriz= %d / %d | %.2f iters/gridpt\n"
               "   Error, Sum: rho %.3e, %.3e | tau %.3e, %.3e | entropy %.3e, %.3e | Y_e %.3e, %.3e\n"
               "               Sx %.3e, %.3e | Sy %.3e, %.3e | Sz %.3e, %.3e\n",
               cctk_iteration, (int)GetRefinementLevel(cctkGH), pointcount,
               backup0, backup1, backup2,
               vel_limited_ptcount, rho_star_fix_applied,
               failures, failures_inhoriz, pointcount_inhoriz,
               (double)n_iter/( (double)(cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]) ),
               rho_error, error_rho_denom,
               tau_error, error_tau_denom,
               entropy_error, error_ent_denom,
               Ye_error, error_Ye_denom,
               Sx_error, error_Sx_denom,
               Sy_error, error_Sy_denom,
               Sz_error, error_Sz_denom);
  }
}
