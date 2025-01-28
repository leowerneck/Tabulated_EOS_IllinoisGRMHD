// Thorn      : IllinoisGRMHD
// File       : con2prim_test_suite.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: In this file we provide an extensive test suite of
//              the con2prim routines available in the code.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "GRHayLib.h"
#include "IllinoisGRMHD_headers.h"
#include "con2prim_headers.h"

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

static inline void ghl_adm_to_igm_bssn(const ghl_metric_quantities *ADM_metric, CCTK_REAL *METRIC) {
  const double phi   = (1.0 / 12.0) * log(ADM_metric->detgamma);
  const double psi   = exp(phi);
  const double psi4  = pow(psi, 4.0);
  const double psim4 = 1.0 / psi4;
  METRIC[PHI]        = phi;
  METRIC[LAPM1]      = ADM_metric->lapse - 1.0;
  METRIC[SHIFTX]     = ADM_metric->betaU[0];
  METRIC[SHIFTY]     = ADM_metric->betaU[1];
  METRIC[SHIFTZ]     = ADM_metric->betaU[2];
  METRIC[GXX]        = ADM_metric->gammaDD[0][0] * psim4;
  METRIC[GXY]        = ADM_metric->gammaDD[0][1] * psim4;
  METRIC[GXZ]        = ADM_metric->gammaDD[0][2] * psim4;
  METRIC[GYY]        = ADM_metric->gammaDD[1][1] * psim4;
  METRIC[GYZ]        = ADM_metric->gammaDD[1][2] * psim4;
  METRIC[GZZ]        = ADM_metric->gammaDD[2][2] * psim4;
  METRIC[GUPXX]      = ADM_metric->gammaUU[0][0] * psi4;
  METRIC[GUPYY]      = ADM_metric->gammaUU[0][1] * psi4;
  METRIC[GUPZZ]      = ADM_metric->gammaUU[0][2] * psi4;
  METRIC[GUPXY]      = ADM_metric->gammaUU[1][1] * psi4;
  METRIC[GUPXZ]      = ADM_metric->gammaUU[1][2] * psi4;
  METRIC[GUPYZ]      = ADM_metric->gammaUU[2][2] * psi4;
}

static inline void ghl_prims_to_igm(const ghl_primitive_quantities *prims, CCTK_REAL *PRIMS) {
  PRIMS[RHOB]        = prims->rho;
  PRIMS[YEPRIM]      = prims->Y_e;
  PRIMS[TEMPERATURE] = prims->temperature;
  PRIMS[PRESSURE]    = prims->press;
  PRIMS[EPSILON]     = prims->eps;
  PRIMS[VX]          = prims->vU[0];
  PRIMS[VY]          = prims->vU[1];
  PRIMS[VZ]          = prims->vU[2];
  PRIMS[BX_CENTER]   = prims->BU[0];
  PRIMS[BY_CENTER]   = prims->BU[1];
  PRIMS[BZ_CENTER]   = prims->BU[2];
}

static inline void igm_prims_to_ghl(const CCTK_REAL *PRIMS, ghl_primitive_quantities *prims) {
  prims->rho         = PRIMS[RHOB];
  prims->Y_e         = PRIMS[YEPRIM];
  prims->temperature = PRIMS[TEMPERATURE];
  prims->press       = PRIMS[PRESSURE];
  prims->eps         = PRIMS[EPSILON];
  prims->vU[0]       = PRIMS[VX];
  prims->vU[1]       = PRIMS[VY];
  prims->vU[2]       = PRIMS[VZ];
  prims->BU[0]       = PRIMS[BX_CENTER];
  prims->BU[1]       = PRIMS[BY_CENTER];
  prims->BU[2]       = PRIMS[BZ_CENTER];
}

static inline void ghl_cons_to_igm(const ghl_conservative_quantities *ghl_cons, CCTK_REAL *igm_cons) {
  igm_cons[RHOSTAR]   = ghl_cons->rho;
  igm_cons[STILDEX]   = ghl_cons->SD[0];
  igm_cons[STILDEY]   = ghl_cons->SD[1];
  igm_cons[STILDEZ]   = ghl_cons->SD[2];
  igm_cons[TAUENERGY] = ghl_cons->tau;
  igm_cons[YESTAR]    = ghl_cons->Y_e;
  igm_cons[ENTSTAR]   = ghl_cons->entropy;
}

static inline void igm_cons_to_ghl(const CCTK_REAL *igm_cons, ghl_conservative_quantities *ghl_cons) {
  ghl_cons->rho     = igm_cons[RHOSTAR];
  ghl_cons->SD[0]   = igm_cons[STILDEX];
  ghl_cons->SD[1]   = igm_cons[STILDEY];
  ghl_cons->SD[2]   = igm_cons[STILDEZ];
  ghl_cons->tau     = igm_cons[TAUENERGY];
  ghl_cons->Y_e     = igm_cons[YESTAR];
  ghl_cons->entropy = igm_cons[ENTSTAR];
}

extern "C" void IllinoisGRMHD_con2prim_test_suit(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  igm_eos_parameters eos;
  initialize_igm_eos_parameters_from_input(igm_eos_key, cctk_time, eos);

  FILE *fp = fopen(con2prim_test_input_file_filename, "rb");
  if(fp == NULL) {
    CCTK_VERROR("Could not open input file %s", con2prim_test_input_file_filename);
  }

  int nrho, nt, nye;
  if(fread(&nrho, sizeof(int), 1, fp) != 1) {
    fclose(fp);
    CCTK_ERROR("Failed to read nrho");
  }
  if(fread(&nt, sizeof(int), 1, fp) != 1) {
    fclose(fp);
    CCTK_ERROR("Failed to read nt");
  }
  if(fread(&nye, sizeof(int), 1, fp) != 1) {
    fclose(fp);
    CCTK_ERROR("Failed to read nye");
  }
  assert(nrho == 50);
  assert(nt == 50);
  assert(nye == 50);

  int failures = 0;
  for(int k = 0; k < nye; k++) {
    for(int j = 0; j < nt; j++) {
      for(int i = 0; i < nrho; i++) {
        ghl_metric_quantities ghl_ADM_metric;
        if(fread(&ADM_metric, sizeof(ghl_ADM_metric), 1, fp) != 1) {
          fclose(fp);
          CCTK_VERROR("Failed to read ADM metric at index %d, %d, %d", i, j, k);
        }

        ghl_ADM_aux_quantities ghl_AUX_metric;
        if(fread(&AUX_metric, sizeof(ghl_AUX_metric), 1, fp) != 1) {
          fclose(fp);
          CCTK_VERROR("Failed to read AUX metric at index %d, %d, %d", i, j, k);
        }

        ghl_primitive_quantities ghl_prims_orig;
        if(fread(&prims_orig, sizeof(ghl_prims_orig), 1, fp) != 1) {
          fclose(fp);
          CCTK_VERROR("Failed to read prims at index %d, %d, %d", i, j, k);
        }

        ghl_conservative_quantities ghl_cons;
        if(fread(&cons, sizeof(ghl_cons), 1, fp) != 1) {
          fclose(fp);
          CCTK_VERROR("Failed to read cons at index %d, %d, %d", i, j, k);
        }

        CCTK_REAL metric[NUMVARS_FOR_METRIC];
        ghl_adm_to_igm_bssn(ghl_ADM_metric, metric);

        CCTK_REAL metric_aux[NUMVARS_METRIC_AUX];
        SET_LAPSE_PSI4(metric_aux, metric);

        CCTK_REAL METRIC_PHYS[NUMVARS_FOR_METRIC];
        METRIC_PHYS[GXX]   = METRIC[GXX] * METRIC_LAP_PSI4[PSI4];
        METRIC_PHYS[GXY]   = METRIC[GXY] * METRIC_LAP_PSI4[PSI4];
        METRIC_PHYS[GXZ]   = METRIC[GXZ] * METRIC_LAP_PSI4[PSI4];
        METRIC_PHYS[GYY]   = METRIC[GYY] * METRIC_LAP_PSI4[PSI4];
        METRIC_PHYS[GYZ]   = METRIC[GYZ] * METRIC_LAP_PSI4[PSI4];
        METRIC_PHYS[GZZ]   = METRIC[GZZ] * METRIC_LAP_PSI4[PSI4];
        METRIC_PHYS[GUPXX] = METRIC[GUPXX] * METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPXY] = METRIC[GUPXY] * METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPXZ] = METRIC[GUPXZ] * METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPYY] = METRIC[GUPYY] * METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPYZ] = METRIC[GUPYZ] * METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPZZ] = METRIC[GUPZZ] * METRIC_LAP_PSI4[PSIM4];

        CCTK_REAL cons[NUM_CONSERVS];
        ghl_cons_to_igm(ghl_cons, cons);

        CCTK_REAL g4dn[4][4], g4up[4][4];
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

        CCTK_REAL alpha_inv_squared = SQR(METRIC_LAP_PSI4[LAPSEINV]);
        g4up[0][0]                  = -1.0 * alpha_inv_squared;
        g4up[0][1] = g4up[1][0] = METRIC[SHIFTX] * alpha_inv_squared;
        g4up[0][2] = g4up[2][0] = METRIC[SHIFTY] * alpha_inv_squared;
        g4up[0][3] = g4up[3][0] = METRIC[SHIFTZ] * alpha_inv_squared;
        g4up[1][1]              = METRIC_PHYS[GUPXX] - METRIC[SHIFTX] * METRIC[SHIFTX] * alpha_inv_squared;
        g4up[1][2] = g4up[2][1] = METRIC_PHYS[GUPXY] - METRIC[SHIFTX] * METRIC[SHIFTY] * alpha_inv_squared;
        g4up[1][3] = g4up[3][1] = METRIC_PHYS[GUPXZ] - METRIC[SHIFTX] * METRIC[SHIFTZ] * alpha_inv_squared;
        g4up[2][2]              = METRIC_PHYS[GUPYY] - METRIC[SHIFTY] * METRIC[SHIFTY] * alpha_inv_squared;
        g4up[2][3] = g4up[3][2] = METRIC_PHYS[GUPYZ] - METRIC[SHIFTY] * METRIC[SHIFTZ] * alpha_inv_squared;
        g4up[3][3]              = METRIC_PHYS[GUPZZ] - METRIC[SHIFTZ] * METRIC[SHIFTZ] * alpha_inv_squared;

        // for(int n = 0; n < 4; n++) {
        //   params.main_routine = methods[n];
        //   CCTK_REAL prims[MAXNUMVARS];
        //
        //   ghl_primitive_quantities prims;
        //   if(ghl_con2prim_tabulated_multi_method(
        //            &params, &eos, &ADM_metric, &AUX_metric, &cons_undens, &prims, &diagnostics)) {
        //     fails[n]++;
        //   }
        // }
      }
    }
  }
  fclose(fp);

  const int ntotal = npoints * npoints * npoints;

  // CCTK_VINFO("Completed test for routine %s", routine);
  CCTK_VINFO("Final report:");
  CCTK_VINFO("    Number of recovery attempts: %d", ntotal);
  CCTK_VINFO("    Number of failed recoveries: %d", failures);
  CCTK_VINFO("    Recovery failure rate      : %.2lf%%", ((CCTK_REAL)failures) / ((CCTK_REAL)ntotal) * 100.0);

  CCTK_VINFO("All done! Terminating the run.");
  exit(1);
}
