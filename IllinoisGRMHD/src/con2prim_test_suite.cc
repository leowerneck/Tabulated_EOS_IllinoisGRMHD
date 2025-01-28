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

static inline void ghl_adm_to_igm_bssn(const ghl_metric_quantities *ghl_adm, CCTK_REAL *igm_bssn) {
  const double phi   = (1.0 / 12.0) * log(ghl_adm->detgamma);
  const double psi   = exp(phi);
  const double psi4  = pow(psi, 4.0);
  const double psim4 = 1.0 / psi4;
  igm_bssn[PHI]      = phi;
  igm_bssn[LAPM1]    = ghl_adm->lapse - 1.0;
  igm_bssn[SHIFTX]   = ghl_adm->betaU[0];
  igm_bssn[SHIFTY]   = ghl_adm->betaU[1];
  igm_bssn[SHIFTZ]   = ghl_adm->betaU[2];
  igm_bssn[GXX]      = ghl_adm->gammaDD[0][0] * psim4;
  igm_bssn[GXY]      = ghl_adm->gammaDD[0][1] * psim4;
  igm_bssn[GXZ]      = ghl_adm->gammaDD[0][2] * psim4;
  igm_bssn[GYY]      = ghl_adm->gammaDD[1][1] * psim4;
  igm_bssn[GYZ]      = ghl_adm->gammaDD[1][2] * psim4;
  igm_bssn[GZZ]      = ghl_adm->gammaDD[2][2] * psim4;
  igm_bssn[GUPXX]    = ghl_adm->gammaUU[0][0] * psi4;
  igm_bssn[GUPYY]    = ghl_adm->gammaUU[0][1] * psi4;
  igm_bssn[GUPZZ]    = ghl_adm->gammaUU[0][2] * psi4;
  igm_bssn[GUPXY]    = ghl_adm->gammaUU[1][1] * psi4;
  igm_bssn[GUPXZ]    = ghl_adm->gammaUU[1][2] * psi4;
  igm_bssn[GUPYZ]    = ghl_adm->gammaUU[2][2] * psi4;
}

static inline void ghl_prims_to_igm(const ghl_primitive_quantities *ghl_prims, CCTK_REAL *igm_prims) {
  igm_prims[RHOB]        = ghl_prims->rho;
  igm_prims[YEPRIM]      = ghl_prims->Y_e;
  igm_prims[TEMPERATURE] = ghl_prims->temperature;
  igm_prims[PRESSURE]    = ghl_prims->press;
  igm_prims[EPSILON]     = ghl_prims->eps;
  igm_prims[VX]          = ghl_prims->vU[0];
  igm_prims[VY]          = ghl_prims->vU[1];
  igm_prims[VZ]          = ghl_prims->vU[2];
  igm_prims[BX_CENTER]   = ghl_prims->BU[0];
  igm_prims[BY_CENTER]   = ghl_prims->BU[1];
  igm_prims[BZ_CENTER]   = ghl_prims->BU[2];
}

static inline void igm_prims_to_ghl(const CCTK_REAL *igm_prims, ghl_primitive_quantities *ghl_prims) {
  ghl_prims->rho         = igm_prims[RHOB];
  ghl_prims->Y_e         = igm_prims[YEPRIM];
  ghl_prims->temperature = igm_prims[TEMPERATURE];
  ghl_prims->press       = igm_prims[PRESSURE];
  ghl_prims->eps         = igm_prims[EPSILON];
  ghl_prims->vU[0]       = igm_prims[VX];
  ghl_prims->vU[1]       = igm_prims[VY];
  ghl_prims->vU[2]       = igm_prims[VZ];
  ghl_prims->BU[0]       = igm_prims[BX_CENTER];
  ghl_prims->BU[1]       = igm_prims[BY_CENTER];
  ghl_prims->BU[2]       = igm_prims[BZ_CENTER];
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

  const char *filename = "ghl_unit_test_con2prim_tabulated.bin";
  FILE *fp             = fopen(filename, "rb");
  if(fp == NULL) {
    CCTK_VERROR("Could not open input file %s", filename);
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
        if(fread(&ghl_ADM_metric, sizeof(ghl_ADM_metric), 1, fp) != 1) {
          fclose(fp);
          CCTK_VERROR("Failed to read ADM metric at index %d, %d, %d", i, j, k);
        }

        ghl_ADM_aux_quantities ghl_AUX_metric;
        if(fread(&ghl_AUX_metric, sizeof(ghl_AUX_metric), 1, fp) != 1) {
          fclose(fp);
          CCTK_VERROR("Failed to read AUX metric at index %d, %d, %d", i, j, k);
        }

        ghl_primitive_quantities ghl_prims_orig;
        if(fread(&ghl_prims_orig, sizeof(ghl_prims_orig), 1, fp) != 1) {
          fclose(fp);
          CCTK_VERROR("Failed to read prims at index %d, %d, %d", i, j, k);
        }

        ghl_conservative_quantities ghl_cons;
        if(fread(&ghl_cons, sizeof(ghl_cons), 1, fp) != 1) {
          fclose(fp);
          CCTK_VERROR("Failed to read cons at index %d, %d, %d", i, j, k);
        }

        CCTK_REAL metric[NUMVARS_FOR_METRIC];
        ghl_adm_to_igm_bssn(&ghl_ADM_metric, metric);

        CCTK_REAL metric_aux[NUMVARS_METRIC_AUX];
        SET_LAPSE_PSI4(metric_aux, metric);

        CCTK_REAL metric_phys[NUMVARS_FOR_METRIC];
        metric_phys[GXX]   = metric[GXX] * metric_aux[PSI4];
        metric_phys[GXY]   = metric[GXY] * metric_aux[PSI4];
        metric_phys[GXZ]   = metric[GXZ] * metric_aux[PSI4];
        metric_phys[GYY]   = metric[GYY] * metric_aux[PSI4];
        metric_phys[GYZ]   = metric[GYZ] * metric_aux[PSI4];
        metric_phys[GZZ]   = metric[GZZ] * metric_aux[PSI4];
        metric_phys[GUPXX] = metric[GUPXX] * metric_aux[PSIM4];
        metric_phys[GUPXY] = metric[GUPXY] * metric_aux[PSIM4];
        metric_phys[GUPXZ] = metric[GUPXZ] * metric_aux[PSIM4];
        metric_phys[GUPYY] = metric[GUPYY] * metric_aux[PSIM4];
        metric_phys[GUPYZ] = metric[GUPYZ] * metric_aux[PSIM4];
        metric_phys[GUPZZ] = metric[GUPZZ] * metric_aux[PSIM4];

        CCTK_REAL igm_cons[NUM_CONSERVS];
        ghl_cons_to_igm(&ghl_cons, igm_cons);

        CCTK_REAL shift_xL = metric_phys[GXX] * metric[SHIFTX] + metric_phys[GXY] * metric[SHIFTY]
                             + metric_phys[GXZ] * metric[SHIFTZ];
        CCTK_REAL shift_yL = metric_phys[GXY] * metric[SHIFTX] + metric_phys[GYY] * metric[SHIFTY]
                             + metric_phys[GYZ] * metric[SHIFTZ];
        CCTK_REAL shift_zL = metric_phys[GXZ] * metric[SHIFTX] + metric_phys[GYZ] * metric[SHIFTY]
                             + metric_phys[GZZ] * metric[SHIFTZ];
        CCTK_REAL beta2L = shift_xL * metric[SHIFTX] + shift_yL * metric[SHIFTY] + shift_zL * metric[SHIFTZ];

        CCTK_REAL g4dn[4][4], g4up[4][4];
        g4dn[0][0] = -SQR(metric_aux[LAPSE]) + beta2L;
        g4dn[0][1] = g4dn[1][0] = metric[SHIFTX];
        g4dn[0][2] = g4dn[2][0] = metric[SHIFTY];
        g4dn[0][3] = g4dn[3][0] = metric[SHIFTZ];
        g4dn[1][1]              = metric_phys[GXX];
        g4dn[1][2] = g4dn[2][1] = metric_phys[GXY];
        g4dn[1][3] = g4dn[3][1] = metric_phys[GXZ];
        g4dn[2][2]              = metric_phys[GYY];
        g4dn[2][3] = g4dn[3][2] = metric_phys[GYZ];
        g4dn[3][3]              = metric_phys[GZZ];

        CCTK_REAL alpha_inv_squared = SQR(metric_aux[LAPSEINV]);
        g4up[0][0]                  = -1.0 * alpha_inv_squared;
        g4up[0][1] = g4up[1][0] = metric[SHIFTX] * alpha_inv_squared;
        g4up[0][2] = g4up[2][0] = metric[SHIFTY] * alpha_inv_squared;
        g4up[0][3] = g4up[3][0] = metric[SHIFTZ] * alpha_inv_squared;
        g4up[1][1]              = metric_phys[GUPXX] - metric[SHIFTX] * metric[SHIFTX] * alpha_inv_squared;
        g4up[1][2] = g4up[2][1] = metric_phys[GUPXY] - metric[SHIFTX] * metric[SHIFTY] * alpha_inv_squared;
        g4up[1][3] = g4up[3][1] = metric_phys[GUPXZ] - metric[SHIFTX] * metric[SHIFTZ] * alpha_inv_squared;
        g4up[2][2]              = metric_phys[GUPYY] - metric[SHIFTY] * metric[SHIFTY] * alpha_inv_squared;
        g4up[2][3] = g4up[3][2] = metric_phys[GUPYZ] - metric[SHIFTY] * metric[SHIFTZ] * alpha_inv_squared;
        g4up[3][3]              = metric_phys[GUPZZ] - metric[SHIFTZ] * metric[SHIFTZ] * alpha_inv_squared;
      }
    }
  }
  fclose(fp);

  const int ntotal = nrho * nye * nt;

  // CCTK_VINFO("Completed test for routine %s", routine);
  CCTK_VINFO("Final report:");
  CCTK_VINFO("    Number of recovery attempts: %d", ntotal);
  CCTK_VINFO("    Number of failed recoveries: %d", failures);
  CCTK_VINFO("    Recovery failure rate      : %.2lf%%", ((CCTK_REAL)failures) / ((CCTK_REAL)ntotal) * 100.0);

  CCTK_VINFO("All done! Terminating the run.");
  exit(0);
}
