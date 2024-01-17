#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "GRHayLib.h"
#include "IllinoisGRMHD_headers.h"
#include "inlined_functions.h"
#include "apply_tau_floor__enforce_limits_on_primitives_and_recompute_conservs.C"

/*********************************************
 * Apply outer boundary conditions on A_{\mu}
 ********************************************/
extern "C" void IllinoisGRMHD_outer_boundaries_on_A_mu(CCTK_ARGUMENTS) {
}

extern "C" void IllinoisGRMHD_outer_boundaries_on_P_rho_b_vx_vy_vz(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_EQUALS(Matter_BC,"frozen")) return;

  bool Symmetry_none=false; if(CCTK_EQUALS(Symmetry,"none")) Symmetry_none=true;

  int levelnumber = GetRefinementLevel(cctkGH);

  // Don't apply approximate outer boundary conditions on initial data, which should be defined everywhere, or on levels != [coarsest level].
  if(cctk_iteration==0 || levelnumber!=0) return;

  if(cctk_nghostzones[0]!=cctk_nghostzones[1] || cctk_nghostzones[0]!=cctk_nghostzones[2])
    CCTK_VError(VERR_DEF_PARAMS,"ERROR: IllinoisGRMHD outer BC driver does not support unequal number of ghostzones in different directions!");
  for(int which_bdry_pt=0;which_bdry_pt<cctk_nghostzones[0];which_bdry_pt++) {

    // Order here is for compatibility with old version of this code.
    /* XMIN & XMAX */
    // i=imax=outer boundary
    if(cctk_bbox[1]) {
      int imax = cctk_lsh[0] - cctk_nghostzones[0] + which_bdry_pt;
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int j=0; j<cctk_lsh[1]; j++) {
          const int index = CCTK_GFINDEX3D(cctkGH, imax,   j, k);
          const int indm1 = CCTK_GFINDEX3D(cctkGH, imax-1, j, k);

          // removed ENABLE from IGM and do_outflow from GRHayLHD for all inflow logic
          ghl_primitive_quantities prims;
          prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;
          prims.rho         = rho[indm1];
          prims.press       = press[indm1];
          prims.vU[0]       = (vx[indm1] < 0.0) ? 0.0 : vx[indm1];
          prims.vU[1]       = vy[indm1];
          prims.vU[2]       = vz[indm1];
          prims.Y_e         = Y_e[indm1];
          prims.temperature = temperature[indm1];

          rho[index]       = prims.rho;
          press[index]           = prims.press;
          vx[index]          = prims.vU[0];
          vy[index]          = prims.vU[1];
          vz[index]          = prims.vU[2];
          Y_e[index]         = prims.Y_e;
          temperature[index] = prims.temperature;
        }
      }
    }

    // i=imin=outer boundary
    if(cctk_bbox[0]) {
      int imin = cctk_nghostzones[0] - which_bdry_pt - 1;
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int j=0; j<cctk_lsh[1]; j++) { 
          const int index = CCTK_GFINDEX3D(cctkGH, imin, j, k);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, imin+1, j, k);

          ghl_primitive_quantities prims;
          prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;
          prims.rho         = rho[indp1];
          prims.press       = press[indp1];
          prims.vU[0]       = (vx[indp1] > 0.0) ? 0.0 : vx[indp1];
          prims.vU[1]       = vy[indp1];
          prims.vU[2]       = vz[indp1];
          prims.Y_e         = Y_e[indp1];
          prims.temperature = temperature[indp1];

          rho[index]       = prims.rho;
          press[index]           = prims.press;
          vx[index]          = prims.vU[0];
          vy[index]          = prims.vU[1];
          vz[index]          = prims.vU[2];
          Y_e[index]         = prims.Y_e;
          temperature[index] = prims.temperature;
        }
      }
    }

    /* YMIN & YMAX */
    // j=jmax=outer boundary
    if(cctk_bbox[3]) {
      int jmax = cctk_lsh[1] - cctk_nghostzones[1] + which_bdry_pt;
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, jmax, k);
          const int indm1 = CCTK_GFINDEX3D(cctkGH, i, jmax-1, k);

          ghl_primitive_quantities prims;
          prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;
          prims.rho         = rho[indm1];
          prims.press       = press[indm1];
          prims.vU[0]       = vx[indm1];
          prims.vU[1]       = (vy[indm1] < 0.0) ? 0.0 : vy[indm1];
          prims.vU[2]       = vz[indm1];
          prims.Y_e         = Y_e[indm1];
          prims.temperature = temperature[indm1];

          rho[index]       = prims.rho;
          press[index]           = prims.press;
          vx[index]          = prims.vU[0];
          vy[index]          = prims.vU[1];
          vz[index]          = prims.vU[2];
          Y_e[index]         = prims.Y_e;
          temperature[index] = prims.temperature;
        }
      }
    }

    // j=jmin=outer boundary
    if(cctk_bbox[2]) {
      int jmin = cctk_nghostzones[1] - which_bdry_pt - 1;
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, jmin, k);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, i, jmin+1, k);

          ghl_primitive_quantities prims;
          prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;
          prims.rho         = rho[indp1];
          prims.press       = press[indp1];
          prims.vU[0]       = vx[indp1];
          prims.vU[1]       = (vy[indp1] > 0.0) ? 0.0 : vy[indp1];
          prims.vU[2]       = vz[indp1];
          prims.Y_e         = Y_e[indp1];
          prims.temperature = temperature[indp1];

          rho[index]       = prims.rho;
          press[index]           = prims.press;
          vx[index]          = prims.vU[0];
          vy[index]          = prims.vU[1];
          vz[index]          = prims.vU[2];
          Y_e[index]         = prims.Y_e;
          temperature[index] = prims.temperature;
        }
      }
    }

    /* ZMIN & ZMAX */
    // k=kmax=outer boundary
    if(cctk_bbox[5]) {
      int kmax = cctk_lsh[2] - cctk_nghostzones[2] + which_bdry_pt;
      for(int j=0; j<cctk_lsh[1]; j++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, j, kmax);
          const int indm1 = CCTK_GFINDEX3D(cctkGH, i, j, kmax-1);

          ghl_primitive_quantities prims;
          prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;
          prims.rho         = rho[indm1];
          prims.press       = press[indm1];
          prims.vU[0]       = vx[indm1];
          prims.vU[1]       = vy[indm1];
          prims.vU[2]       = (vz[indm1] < 0.0) ? 0.0 : vz[indm1];
          prims.Y_e         = Y_e[indm1];
          prims.temperature = temperature[indm1];

          rho[index]       = prims.rho;
          press[index]           = prims.press;
          vx[index]          = prims.vU[0];
          vy[index]          = prims.vU[1];
          vz[index]          = prims.vU[2];
          Y_e[index]         = prims.Y_e;
          temperature[index] = prims.temperature;
        }
      }
    }

    // k=kmin=outer boundary
    if((cctk_bbox[4]) && Symmetry_none) {
      int kmin = cctk_nghostzones[2] - which_bdry_pt - 1;
      for(int j=0; j<cctk_lsh[1]; j++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, j, kmin);
          const int indp1 = CCTK_GFINDEX3D(cctkGH, i, j, kmin+1);

          ghl_primitive_quantities prims;
          prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;
          prims.rho         = rho[indp1];
          prims.press       = press[indp1];
          prims.vU[0]       = vx[indp1];
          prims.vU[1]       = vy[indp1];
          prims.vU[2]       = (vz[indp1] > 0.0) ? 0.0 : vz[indp1];
          prims.Y_e         = Y_e[indp1];
          prims.temperature = temperature[indp1];

          rho[index] = rho[indp1];
          press[index] = press[indp1];
          vx[index] = vx[indp1];
          vy[index] = vy[indp1];
          vz[index] = vz[indp1];
          Y_e[index] = Y_e[indp1];
          temperature[index] = temperature[indp1];
        }
      }
    }
  }

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        if(((cctk_bbox[0]) && i<cctk_nghostzones[0]) ||
           ((cctk_bbox[1]) && i>=cctk_lsh[0]-cctk_nghostzones[0]) ||
           ((cctk_bbox[2]) && j<cctk_nghostzones[1]) ||
           ((cctk_bbox[3]) && j>=cctk_lsh[1]-cctk_nghostzones[1]) ||
           ((cctk_bbox[4]) && k<cctk_nghostzones[2] && CCTK_EQUALS(Symmetry,"none")) ||
           ((cctk_bbox[5]) && k>=cctk_lsh[2]-cctk_nghostzones[2])) {
          int index = CCTK_GFINDEX3D(cctkGH,i, j, k);

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
          prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;
          prims.rho         = rho[index];
          prims.press       = press[index];
          prims.vU[0]       = vx[index];
          prims.vU[1]       = vy[index];
          prims.vU[2]       = vz[index];
          prims.Y_e         = Y_e[index];
          prims.temperature = temperature[index];

          ghl_conservative_quantities cons;
          const int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_enforce_primitive_limits_and_compute_u0(
                ghl_params, ghl_eos, &ADM_metric, &prims);

          ghl_compute_conservs(
                &ADM_metric, &metric_aux, &prims, &cons);

          rho[index]       = prims.rho;
          press[index]           = prims.press;
          vx[index]          = prims.vU[0];
          vy[index]          = prims.vU[1];
          vz[index]          = prims.vU[2];
          Y_e[index]         = prims.Y_e;
          temperature[index] = prims.temperature;

          rho_star[index] = cons.rho;
          tau[index]      = cons.tau;
          mhd_st_x[index]  = cons.SD[0];
          mhd_st_y[index]  = cons.SD[1];
          mhd_st_z[index]  = cons.SD[2];
          Ye_star[index]  = cons.Y_e;

          // removed setTmunu from here
        }
      }
    }
  }
}
