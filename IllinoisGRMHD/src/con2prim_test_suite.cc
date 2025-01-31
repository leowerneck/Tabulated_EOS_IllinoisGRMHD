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

typedef struct {
  CCTK_REAL g4DD[4][4], g4UU[4][4];
} igm_aux_metric;

static inline void ghl_metric_aux_to_igm(const ghl_ADM_aux_quantities *ghl_aux, igm_aux_metric *igm_aux) {
  for(int mu = 0; mu < 4; mu++) {
    for(int nu = 0; nu < 4; nu++) {
      igm_aux->g4DD[mu][nu] = ghl_aux->g4DD[mu][nu];
      igm_aux->g4UU[mu][nu] = ghl_aux->g4UU[mu][nu];
    }
  }
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
  eos.T_max = 100;

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

  int num_routines           = 4;
  int routines[4]            = { igm_Palenzuela1D, igm_Palenzuela1D_entropy, igm_Newman1D, igm_Newman1D_entropy };
  const char *methodnames[4] = { "Palenzuela1D", "Palenzuela1D_entropy", "Newman1D", "Newman1D_entropy" };
  int fails[4]               = { 0, 0, 0, 0 };
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

        ghl_primitive_quantities ghl_prims;
        if(fread(&ghl_prims, sizeof(ghl_prims), 1, fp) != 1) {
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

        CCTK_REAL igm_prims[MAXNUMVARS];
        ghl_prims_to_igm(&ghl_prims, igm_prims);

        CCTK_REAL igm_cons[NUM_CONSERVS];
        ghl_cons_to_igm(&ghl_cons, igm_cons);

        igm_aux_metric igm_aux;
        ghl_metric_aux_to_igm(&ghl_AUX_metric, &igm_aux);

        struct output_stats stats = {};
        stats.which_routine       = igm_None;
        stats.dx[0]               = CCTK_DELTA_SPACE(0);
        stats.dx[1]               = CCTK_DELTA_SPACE(1);
        stats.dx[2]               = CCTK_DELTA_SPACE(2);

        for(int n = 0; n < num_routines; n++) {
          CCTK_REAL c2p_cons[numcons];
          set_cons_from_PRIMS_and_CONSERVS(eos, eos.c2p_routine, metric, metric_aux, igm_prims, igm_cons, c2p_cons);

          CCTK_REAL c2p_prims[numprims];
          set_prim_from_PRIMS_and_CONSERVS(
                eos, eos.c2p_routine, 0, metric, metric_aux, igm_prims, igm_cons, c2p_cons, c2p_prims);

          if(con2prim_select(eos, routines[n], metric_phys, igm_aux.g4DD, igm_aux.g4UU, c2p_cons, c2p_prims, stats)) {
            fails[n]++;
          }
        }
      }
    }
  }
  fclose(fp);

  CCTK_VINFO("T_max = %g", eos.T_max);

  const int npts = nrho * nt * nye;
  CCTK_VINFO("Failure rates:");
  for(int n = 0; n < num_routines; n++) {
    const char *name    = methodnames[n];
    const int failcount = fails[n];
    const double pct    = ((double)failcount) / ((double)npts) * 100;
    CCTK_VINFO("    %-20s : %06d/%06d : %5.1lf%%", name, failcount, npts, pct);
  }

  CCTK_VINFO("All done!");
  exit(0);
}
