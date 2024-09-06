// Thorn      : IllinoisGRMHD
// File       : con2prim_test_suite.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: In this file we provide an extensive test suite of
//              the con2prim routines available in the code.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "IllinoisGRMHD_headers.h"
#include "con2prim_headers.h"
#include "GRHayLib.h"

#include <fstream>
using namespace std;

static inline CCTK_REAL relative_error(const CCTK_REAL a, const CCTK_REAL b) {
  if(a != 0) {
    return (fabs(1.0 - b / a));
  }
  else if(b != 0) {
    return (fabs(1.0 - a / b));
  }
  else {
    return (0.0);
  }
}

static inline void
set_BSSN_metric_from_ADM_struct(const ghl_metric_quantities *ADM_metric, CCTK_REAL *METRIC) {
  // Recall that gamma_{ij} = \psi^{4} \tilde\gamma_{ij}. Therefore
  // \gamma = \psi^{12} \tilde\gamma \implies phi = 1/12
  // log(\gamma/\tilde\gamma), and that \tilde\gamma=1.
  const double phi = (1.0 / 12.0) * log(ADM_metric->detgamma);
  const double psi = exp(phi);
  const double psi4 = pow(psi, 4.0);
  const double psim4 = 1.0 / psi4;
  METRIC[PHI] = phi;
  METRIC[LAPM1] = ADM_metric->lapse - 1.0;
  METRIC[SHIFTX] = ADM_metric->betaU[0];
  METRIC[SHIFTY] = ADM_metric->betaU[1];
  METRIC[SHIFTZ] = ADM_metric->betaU[2];
  METRIC[GXX] = ADM_metric->gammaDD[0][0] * psim4;
  METRIC[GXY] = ADM_metric->gammaDD[0][1] * psim4;
  METRIC[GXZ] = ADM_metric->gammaDD[0][2] * psim4;
  METRIC[GYY] = ADM_metric->gammaDD[1][1] * psim4;
  METRIC[GYZ] = ADM_metric->gammaDD[1][2] * psim4;
  METRIC[GZZ] = ADM_metric->gammaDD[2][2] * psim4;
  METRIC[GUPXX] = ADM_metric->gammaUU[0][0] * psi4;
  METRIC[GUPYY] = ADM_metric->gammaUU[0][1] * psi4;
  METRIC[GUPZZ] = ADM_metric->gammaUU[0][2] * psi4;
  METRIC[GUPXY] = ADM_metric->gammaUU[1][1] * psi4;
  METRIC[GUPXZ] = ADM_metric->gammaUU[1][2] * psi4;
  METRIC[GUPYZ] = ADM_metric->gammaUU[2][2] * psi4;
}

static inline void
set_PRIMS_from_prims_struct(const ghl_primitive_quantities *prims, CCTK_REAL *PRIMS) {
  // Now set the primitive variables array, following IGM's
  // standards
  PRIMS[RHOB] = prims->rho;
  PRIMS[YEPRIM] = prims->Y_e;
  PRIMS[TEMPERATURE] = prims->temperature;
  PRIMS[PRESSURE] = prims->press;
  PRIMS[EPSILON] = prims->eps;
  PRIMS[VX] = prims->vU[0];
  PRIMS[VY] = prims->vU[1];
  PRIMS[VZ] = prims->vU[2];
  PRIMS[BX_CENTER] = prims->BU[0];
  PRIMS[BY_CENTER] = prims->BU[1];
  PRIMS[BZ_CENTER] = prims->BU[2];
}

static inline void
set_prims_struct_from_PRIMS(const CCTK_REAL *PRIMS, ghl_primitive_quantities *prims) {
  // Now set the primitive variables array, following IGM's
  // standards
  prims->rho = PRIMS[RHOB];
  prims->Y_e = PRIMS[YEPRIM];
  prims->temperature = PRIMS[TEMPERATURE];
  prims->press = PRIMS[PRESSURE];
  prims->eps = PRIMS[EPSILON];
  prims->vU[0] = PRIMS[VX];
  prims->vU[1] = PRIMS[VY];
  prims->vU[2] = PRIMS[VZ];
  prims->BU[0] = PRIMS[BX_CENTER];
  prims->BU[1] = PRIMS[BY_CENTER];
  prims->BU[2] = PRIMS[BZ_CENTER];
}

static void run_unit_test(
      const char *routine,
      igm_eos_parameters &eos) {

  CCTK_VINFO("Beginning unit test for %s", routine);

  const char *vars_string = "rho_vs_T";
  CCTK_VINFO("  Running test %s", vars_string);
  char filename[256];
  sprintf(filename, "con2prim_tabulated_%s_%s_unperturbed.bin", routine, vars_string);
  FILE *fp_unpert = fopen(filename, "rb");
  sprintf(filename, "con2prim_tabulated_%s_%s_perturbed.bin", routine, vars_string);
  FILE *fp_pert = fopen(filename, "rb");

  int n1, n2;
  if(fread(&n1, sizeof(int), 1, fp_unpert) != 1) {
    CCTK_ERROR("Failed to read from file");
  };
  if(fread(&n2, sizeof(int), 1, fp_pert) != 1) {
    CCTK_ERROR("Failed to read from file");
  };
  if(n1 != n2) {
    CCTK_VERROR("Problem reading input data files (%d != %d)", n1, n2);
  }

  const int npoints = n1;
  for(int ir = 0; ir < npoints; ir++) {
    for(int it = 0; it < npoints; it++) {
      for(int iy = 0; iy < npoints; iy++) {
        ghl_con2prim_diagnostics diagnostics;
        ghl_initialize_diagnostics(&diagnostics);

        // Read input metric from unperturbed data file
        ghl_metric_quantities ADM_metric;
        if(fread(&ADM_metric, sizeof(ghl_metric_quantities), 1, fp_unpert) != 1) {
          CCTK_ERROR("Failed to read input metric from file");
        }

        ghl_ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        // Read input primitives from unperturbed data file
        ghl_primitive_quantities prims;
        if(fread(&prims, sizeof(ghl_primitive_quantities), 1, fp_unpert) != 1) {
          CCTK_ERROR("Failed to read input primitives from file");
        }

        // Compute conserved variables and Tmunu
        ghl_conservative_quantities cons;
        __attribute__((unused)) ghl_stress_energy dummy;
        ghl_compute_conservs_and_Tmunu(&ADM_metric, &metric_aux, &prims, &cons, &dummy);

        // Undensitize the conserved variables
        ghl_conservative_quantities cons_undens;
        ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);

        // Now perform the con2prim
        if(ghl_con2prim_tabulated_multi_method(
                 ghl_params, ghl_eos, &ADM_metric, &metric_aux, &cons_undens, &prims, &diagnostics)) {
          CCTK_INFO("Con2Prim failed");
        }

        prims.vU[0] = prims.vU[0] / prims.u0;
        prims.vU[1] = prims.vU[1] / prims.u0;
        prims.vU[2] = prims.vU[2] / prims.u0;

        // Read unperturbed and perturbed results from file
        ghl_primitive_quantities prims_trusted, prims_pert;
        if(fread(&prims_trusted, sizeof(ghl_primitive_quantities), 1, fp_unpert) != 1) {
          CCTK_ERROR("Failed to read trusted primitives from file");
        }
        if(fread(&prims_pert, sizeof(ghl_primitive_quantities), 1, fp_pert) != 1) {
          CCTK_ERROR("Failed to read perturbed primitives from file");
        }

        // Validate results
        // ghl_pert_test_fail(prims_trusted.rho, prims.rho, prims_pert.rho);
        // ghl_pert_test_fail(prims_trusted.Y_e, prims.Y_e, prims_pert.Y_e);
        // ghl_pert_test_fail(prims_trusted.temperature, prims.temperature, prims_pert.temperature);
        // ghl_pert_test_fail(prims_trusted.press, prims.press, prims_pert.press);
        // ghl_pert_test_fail(prims_trusted.eps, prims.eps, prims_pert.eps);
        // ghl_pert_test_fail(prims_trusted.vU[0], prims.vU[0], prims_pert.vU[0]);
        // ghl_pert_test_fail(prims_trusted.vU[1], prims.vU[1], prims_pert.vU[1]);
        // ghl_pert_test_fail(prims_trusted.vU[2], prims.vU[2], prims_pert.vU[2]);
      }
    }
  }
  fclose(fp_unpert);
  fclose(fp_pert);
}

extern "C" void IllinoisGRMHD_con2prim_test_suit(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Initialize the EOS parameters
  igm_eos_parameters eos;
  initialize_igm_eos_parameters_from_input(igm_eos_key, cctk_time, eos);

  // Print information
  if(CCTK_EQUALS(igm_eos_type, "Tabulated") || CCTK_EQUALS(igm_eos_type, "nuc_eos")) {
    CCTK_INFO("EOS type        : Tabulated");
  }

  // Count number of routines tested
  CCTK_INT num_routines_tested = 1;
  CCTK_INT con2prim_test_keys[4];
  char con2prim_test_names[4][100];
  con2prim_test_keys[0] = eos.c2p_routine;
  sprintf(con2prim_test_names[0], "%s", igm_con2prim_routine);
  if(eos.c2p_backup[0] != igm_None) {
    num_routines_tested++;
    con2prim_test_keys[1] = eos.c2p_backup[0];
    sprintf(con2prim_test_names[1], "%s", igm_con2prim_backup_routine[0]);
    if(eos.c2p_backup[1] != igm_None) {
      num_routines_tested++;
      con2prim_test_keys[2] = eos.c2p_backup[1];
      sprintf(con2prim_test_names[2], "%s", igm_con2prim_backup_routine[1]);
      if(eos.c2p_backup[2] != igm_None) {
        num_routines_tested++;
        con2prim_test_keys[3] = eos.c2p_backup[2];
        sprintf(con2prim_test_names[3], "%s", igm_con2prim_backup_routine[2]);
      }
    }
  }

  const double W_test = 2;

  // Now perform one test for each of the selected routines
  for(int which_routine = 0; which_routine < num_routines_tested; which_routine++) {

    const char *vars_string = "rho_vs_T";
    CCTK_VINFO("  Running test %s", vars_string);
    char filename[256];
    char *routine = con2prim_test_names[which_routine];
    sprintf(filename, "con2prim_tabulated_%s_%s_unperturbed.bin", routine, vars_string);
    FILE *fp = fopen(filename, "rb");


    int npoints;
    if(fread(&npoints, sizeof(int), 1, fp) != 1) {
      CCTK_ERROR("Failed to read from file");
    };

    int failures = 0;

    for(int ir = 0; ir < npoints; ir++) {
      for(int it = 0; it < npoints; it++) {
        for(int iy = 0; iy < npoints; iy++) {

          ghl_metric_quantities ADM_metric;
          if(fread(&ADM_metric, sizeof(ghl_metric_quantities), 1, fp) != 1) {
            CCTK_ERROR("Failed to read input metric from file");
          }

          CCTK_REAL METRIC[NUMVARS_FOR_METRIC];
          set_BSSN_metric_from_ADM_struct(&ADM_metric, METRIC);

          // Read input primitives from unperturbed data file
          ghl_primitive_quantities prims;
          if(fread(&prims, sizeof(ghl_primitive_quantities), 1, fp) != 1) {
            CCTK_ERROR("Failed to read input primitives from file");
          }
          CCTK_REAL PRIMS[MAXNUMVARS];
          set_PRIMS_from_prims_struct(&prims, PRIMS);

          // We'll also need the auxilary metric variables array
          CCTK_REAL METRIC_PHYS[NUMVARS_FOR_METRIC];
          CCTK_REAL METRIC_LAP_PSI4[NUMVARS_METRIC_AUX];
          SET_LAPSE_PSI4(METRIC_LAP_PSI4, METRIC);

          // Now the physical metric
          METRIC_PHYS[GXX] = METRIC[GXX] * METRIC_LAP_PSI4[PSI4];
          METRIC_PHYS[GXY] = METRIC[GXY] * METRIC_LAP_PSI4[PSI4];
          METRIC_PHYS[GXZ] = METRIC[GXZ] * METRIC_LAP_PSI4[PSI4];
          METRIC_PHYS[GYY] = METRIC[GYY] * METRIC_LAP_PSI4[PSI4];
          METRIC_PHYS[GYZ] = METRIC[GYZ] * METRIC_LAP_PSI4[PSI4];
          METRIC_PHYS[GZZ] = METRIC[GZZ] * METRIC_LAP_PSI4[PSI4];
          METRIC_PHYS[GUPXX] = METRIC[GUPXX] * METRIC_LAP_PSI4[PSIM4];
          METRIC_PHYS[GUPXY] = METRIC[GUPXY] * METRIC_LAP_PSI4[PSIM4];
          METRIC_PHYS[GUPXZ] = METRIC[GUPXZ] * METRIC_LAP_PSI4[PSIM4];
          METRIC_PHYS[GUPYY] = METRIC[GUPYY] * METRIC_LAP_PSI4[PSIM4];
          METRIC_PHYS[GUPYZ] = METRIC[GUPYZ] * METRIC_LAP_PSI4[PSIM4];
          METRIC_PHYS[GUPZZ] = METRIC[GUPZZ] * METRIC_LAP_PSI4[PSIM4];

          // Then set the conservative variables array
          const int already_computed_physical_metric_and_inverse = 0;
          CCTK_REAL TUPMUNU[10], TDNMUNU[10];
          CCTK_REAL CONSERVS[NUM_CONSERVS];
          CCTK_REAL g4dn[4][4]
                = { { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 } };
          CCTK_REAL g4up[4][4]
                = { { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 } };
          struct output_stats stats;
          IllinoisGRMHD_enforce_limits_on_primitives_and_recompute_conservs(
                already_computed_physical_metric_and_inverse, PRIMS, stats, eos, METRIC, g4dn,
                g4up, TUPMUNU, TDNMUNU, CONSERVS);

          // The con2prim routines require different conservative
          // variables than those that IllinoisGRMHD evolve. Therefore, we
          // must first convert them into a new set of variables, suitable
          // for primitive recovery.
          CCTK_REAL cons[numcons];
          set_cons_from_PRIMS_and_CONSERVS(
                eos, con2prim_test_keys[which_routine], METRIC, METRIC_LAP_PSI4, PRIMS, CONSERVS,
                cons);

          // The con2prim routines require primitive guesses in order to
          // perform the recovery. In IllinoisGRMHD, we do not keep track
          // of the primitives in between time steps, and therefore our
          // guesses are *not* the values of the primitives in the
          // previous time level. Instead, we provide guesses based on the
          // conservative variables.
          CCTK_INT check = 0;
          CCTK_REAL prim[numprims];
          for(int which_guess = 1; which_guess <= 2; which_guess++) {
            set_prim_from_PRIMS_and_CONSERVS(
                  eos, con2prim_test_keys[which_routine], which_guess, METRIC, METRIC_LAP_PSI4,
                  PRIMS, CONSERVS, cons, prim);
            for(int i = 0; i < numprims; i++) {
              prim[i] = 1e300;
            }
            prim[TEMP] = PRIMS[TEMPERATURE] * 0.95;
            prim[WLORENTZ] = W_test * 0.95;
            check = con2prim_select(
                  eos, con2prim_test_keys[which_routine], METRIC_PHYS, g4dn, g4up, cons, prim,
                  stats);
            if(check == 0) {
              break;
            }
          }

          CCTK_REAL ERRORS[MAXNUMVARS], accumulated_error = 0.0;
          if(check != 0) {
            failures++;
            CCTK_VINFO("Recovery FAILED!");
            accumulated_error = 1e300;
          }
          else {

            // Now that we have found some solution, we first limit
            // velocity:
            // FIXME: Probably want to use exactly the same velocity
            // limiter function here as in mhdflux.C
            CCTK_REAL utx_new = prim[UTCON1];
            CCTK_REAL uty_new = prim[UTCON2];
            CCTK_REAL utz_new = prim[UTCON3];

            // //Velocity limiter:
            CCTK_REAL gijuiuj
                  = METRIC_PHYS[GXX] * SQR(utx_new) + 2.0 * METRIC_PHYS[GXY] * utx_new * uty_new
                    + 2.0 * METRIC_PHYS[GXZ] * utx_new * utz_new + METRIC_PHYS[GYY] * SQR(uty_new)
                    + 2.0 * METRIC_PHYS[GYZ] * uty_new * utz_new + METRIC_PHYS[GZZ] * SQR(utz_new);
            CCTK_REAL au0m1 = gijuiuj / (1.0 + sqrt(1.0 + gijuiuj));
            CCTK_REAL u0L = (au0m1 + 1.0) * METRIC_LAP_PSI4[LAPSEINV];
            PRIMS[YEPRIM] = prim[YE];
            PRIMS[TEMPERATURE] = prim[TEMP];
            PRIMS[PRESSURE] = prim[PRESS];
            PRIMS[EPSILON] = prim[EPS];
            PRIMS[VX] = utx_new / u0L - METRIC[SHIFTX];
            PRIMS[VY] = uty_new / u0L - METRIC[SHIFTY];
            PRIMS[VZ] = utz_new / u0L - METRIC[SHIFTZ];

            // CCTK_VINFO("Recovery SUCCEEDED!");
            // for (int which_prim = 0; which_prim < num_prims_in_error;
            //      which_prim++) {
            //   int primL = which_prims_in_error[which_prim];
            //   ERRORS[primL] = relative_error(PRIMS[primL], PRIMS_ORIG[primL]);
            //   CCTK_VINFO("Relative error for prim %s: %.3e (%e -> %e)",
            //              primnames[primL], ERRORS[primL], PRIMS_ORIG[primL],
            //              PRIMS[primL]);
            //   accumulated_error += ERRORS[primL];
            // }
            // CCTK_VINFO("Total accumulated error    : %e\n", accumulated_error);
          }
          // fprintf(
          //       outfile, "%e %e %e\n", log10(PRIMS_ORIG[RHOB]), log10(PRIMS_ORIG[TEMPERATURE]),
          //       log10(MAX(accumulated_error, 1e-16)));
        }
      }
      // fprintf(outfile, "\n");
    }

    // fclose(outfile);

    const int ntotal = npoints * npoints * npoints;

    CCTK_VINFO("Completed test for routine %s", routine);
    CCTK_VINFO("Final report:");
    CCTK_VINFO("    Number of recovery attempts: %d", ntotal);
    CCTK_VINFO("    Number of failed recoveries: %d", failures);
    CCTK_VINFO(
          "    Recovery failure rate      : %.2lf%%",
          ((CCTK_REAL)failures) / ((CCTK_REAL)ntotal) * 100.0);
  }

  CCTK_VINFO("All done! Terminating the run.");
  exit(1);
}
