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
    CCTK_VERROR("ERROR: IllinoisGRMHD outer BC driver does not support unequal number of ghostzones in different directions!");
  for(int which_bdry_pt=0;which_bdry_pt<cctk_nghostzones[0];which_bdry_pt++) {

    // Order here is for compatibility with old version of this code.
    /* XMIN & XMAX */
    // i=imax=outer boundary
    if(cctk_bbox[1]) {
      int imax = cctk_lsh[0] - cctk_nghostzones[0] + which_bdry_pt;
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int j=0; j<cctk_lsh[1]; j++) {
          rho[CCTK_GFINDEX3D(cctkGH, imax, j, k)] = rho[CCTK_GFINDEX3D(cctkGH, imax-1, j, k)];
          press[CCTK_GFINDEX3D(cctkGH, imax, j, k)] = press[CCTK_GFINDEX3D(cctkGH, imax-1, j, k)];
          vx[CCTK_GFINDEX3D(cctkGH, imax, j, k)] = vx[CCTK_GFINDEX3D(cctkGH, imax-1, j, k)];
          vy[CCTK_GFINDEX3D(cctkGH, imax, j, k)] = vy[CCTK_GFINDEX3D(cctkGH, imax-1, j, k)];
          vz[CCTK_GFINDEX3D(cctkGH, imax, j, k)] = vz[CCTK_GFINDEX3D(cctkGH, imax-1, j, k)];
          Y_e[CCTK_GFINDEX3D(cctkGH, imax, j, k)] = Y_e[CCTK_GFINDEX3D(cctkGH, imax-1, j, k)];
          temperature[CCTK_GFINDEX3D(cctkGH, imax, j, k)] = temperature[CCTK_GFINDEX3D(cctkGH, imax-1, j, k)];
          if(vx[CCTK_GFINDEX3D(cctkGH, imax, j, k)]<0.) vx[CCTK_GFINDEX3D(cctkGH, imax, j, k)]=0.;
        }
      }
    }

    // i=imin=outer boundary
    if(cctk_bbox[0]) {
      int imin = cctk_nghostzones[0] - which_bdry_pt - 1;
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int j=0; j<cctk_lsh[1]; j++) { 
          rho[CCTK_GFINDEX3D(cctkGH, imin, j, k)] = rho[CCTK_GFINDEX3D(cctkGH, imin+1, j, k)];
          press[CCTK_GFINDEX3D(cctkGH, imin, j, k)] = press[CCTK_GFINDEX3D(cctkGH, imin+1, j, k)];
          vx[CCTK_GFINDEX3D(cctkGH, imin, j, k)] = vx[CCTK_GFINDEX3D(cctkGH, imin+1, j, k)];
          vy[CCTK_GFINDEX3D(cctkGH, imin, j, k)] = vy[CCTK_GFINDEX3D(cctkGH, imin+1, j, k)];
          vz[CCTK_GFINDEX3D(cctkGH, imin, j, k)] = vz[CCTK_GFINDEX3D(cctkGH, imin+1, j, k)];
          Y_e[CCTK_GFINDEX3D(cctkGH, imin, j, k)] = Y_e[CCTK_GFINDEX3D(cctkGH, imin+1, j, k)];
          temperature[CCTK_GFINDEX3D(cctkGH, imin, j, k)] = temperature[CCTK_GFINDEX3D(cctkGH, imin+1, j, k)];
          if(vx[CCTK_GFINDEX3D(cctkGH, imin, j, k)]>0.) vx[CCTK_GFINDEX3D(cctkGH, imin, j, k)]=0.;
        }
      }
    }

    /* YMIN & YMAX */
    // j=jmax=outer boundary
    if(cctk_bbox[3]) {
      int jmax = cctk_lsh[1] - cctk_nghostzones[1] + which_bdry_pt;
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          rho[CCTK_GFINDEX3D(cctkGH, i, jmax, k)] = rho[CCTK_GFINDEX3D(cctkGH, i, jmax-1, k)];
          press[CCTK_GFINDEX3D(cctkGH, i, jmax, k)] = press[CCTK_GFINDEX3D(cctkGH, i, jmax-1, k)];
          vx[CCTK_GFINDEX3D(cctkGH, i, jmax, k)] = vx[CCTK_GFINDEX3D(cctkGH, i, jmax-1, k)];
          vy[CCTK_GFINDEX3D(cctkGH, i, jmax, k)] = vy[CCTK_GFINDEX3D(cctkGH, i, jmax-1, k)];
          vz[CCTK_GFINDEX3D(cctkGH, i, jmax, k)] = vz[CCTK_GFINDEX3D(cctkGH, i, jmax-1, k)];
          Y_e[CCTK_GFINDEX3D(cctkGH, i, jmax, k)] = Y_e[CCTK_GFINDEX3D(cctkGH, i, jmax-1, k)];
          temperature[CCTK_GFINDEX3D(cctkGH, i, jmax, k)] = temperature[CCTK_GFINDEX3D(cctkGH, i, jmax-1, k)];
          if(vy[CCTK_GFINDEX3D(cctkGH, i, jmax, k)]<0.) vy[CCTK_GFINDEX3D(cctkGH, i, jmax, k)]=0.;
        }
      }
    }

    // j=jmin=outer boundary
    if(cctk_bbox[2]) {
      int jmin = cctk_nghostzones[1] - which_bdry_pt - 1;
      for(int k=0; k<cctk_lsh[2]; k++) {
        for(int i=0; i<cctk_lsh[0]; i++) { 
          rho[CCTK_GFINDEX3D(cctkGH, i, jmin, k)] = rho[CCTK_GFINDEX3D(cctkGH, i, jmin+1, k)];
          press[CCTK_GFINDEX3D(cctkGH, i, jmin, k)] = press[CCTK_GFINDEX3D(cctkGH, i, jmin+1, k)];
          vx[CCTK_GFINDEX3D(cctkGH, i, jmin, k)] = vx[CCTK_GFINDEX3D(cctkGH, i, jmin+1, k)];
          vy[CCTK_GFINDEX3D(cctkGH, i, jmin, k)] = vy[CCTK_GFINDEX3D(cctkGH, i, jmin+1, k)];
          vz[CCTK_GFINDEX3D(cctkGH, i, jmin, k)] = vz[CCTK_GFINDEX3D(cctkGH, i, jmin+1, k)];
          Y_e[CCTK_GFINDEX3D(cctkGH, i, jmin, k)] = Y_e[CCTK_GFINDEX3D(cctkGH, i, jmin+1, k)];
          temperature[CCTK_GFINDEX3D(cctkGH, i, jmin, k)] = temperature[CCTK_GFINDEX3D(cctkGH, i, jmin+1, k)];
          if(vy[CCTK_GFINDEX3D(cctkGH, i, jmin, k)]>0.) vy[CCTK_GFINDEX3D(cctkGH, i, jmin, k)]=0.;
        }
      }
    }

    /* ZMIN & ZMAX */
    // k=kmax=outer boundary
    if(cctk_bbox[5]) {
      int kmax = cctk_lsh[2] - cctk_nghostzones[2] + which_bdry_pt;
      for(int j=0; j<cctk_lsh[1]; j++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          rho[CCTK_GFINDEX3D(cctkGH, i, j, kmax)] = rho[CCTK_GFINDEX3D(cctkGH, i, j, kmax-1)];
          press[CCTK_GFINDEX3D(cctkGH, i, j, kmax)] = press[CCTK_GFINDEX3D(cctkGH, i, j, kmax-1)];
          vx[CCTK_GFINDEX3D(cctkGH, i, j, kmax)] = vx[CCTK_GFINDEX3D(cctkGH, i, j, kmax-1)];
          vy[CCTK_GFINDEX3D(cctkGH, i, j, kmax)] = vy[CCTK_GFINDEX3D(cctkGH, i, j, kmax-1)];
          Y_e[CCTK_GFINDEX3D(cctkGH, i, j, kmax)] = Y_e[CCTK_GFINDEX3D(cctkGH, i, j, kmax-1)];
          temperature[CCTK_GFINDEX3D(cctkGH, i, j, kmax)] = temperature[CCTK_GFINDEX3D(cctkGH, i, j, kmax-1)];
          if(vz[CCTK_GFINDEX3D(cctkGH, i, j, kmax)]<0.) vz[CCTK_GFINDEX3D(cctkGH, i, j, kmax)]=0.;
        }
      }
    }

    // k=kmin=outer boundary
    if((cctk_bbox[4]) && Symmetry_none) {
      int kmin = cctk_nghostzones[2] - which_bdry_pt - 1;
      for(int j=0; j<cctk_lsh[1]; j++) {
        for(int i=0; i<cctk_lsh[0]; i++) {
          rho[CCTK_GFINDEX3D(cctkGH, i, j, kmin)] = rho[CCTK_GFINDEX3D(cctkGH, i, j, kmin+1)];
          press[CCTK_GFINDEX3D(cctkGH, i, j, kmin)] = press[CCTK_GFINDEX3D(cctkGH, i, j, kmin+1)];
          vx[CCTK_GFINDEX3D(cctkGH, i, j, kmin)] = vx[CCTK_GFINDEX3D(cctkGH, i, j, kmin+1)];
          vy[CCTK_GFINDEX3D(cctkGH, i, j, kmin)] = vy[CCTK_GFINDEX3D(cctkGH, i, j, kmin+1)];
          vz[CCTK_GFINDEX3D(cctkGH, i, j, kmin)] = vz[CCTK_GFINDEX3D(cctkGH, i, j, kmin+1)];
          Y_e[CCTK_GFINDEX3D(cctkGH, i, j, kmin)] = Y_e[CCTK_GFINDEX3D(cctkGH, i, j, kmin+1)];
          temperature[CCTK_GFINDEX3D(cctkGH, i, j, kmin)] = temperature[CCTK_GFINDEX3D(cctkGH, i, j, kmin+1)];
          if(vz[CCTK_GFINDEX3D(cctkGH, i, j, kmin)]>0.) vz[CCTK_GFINDEX3D(cctkGH, i, j, kmin)]=0.;
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

          rho[index]         = prims.rho;
          press[index]       = prims.press;
          vx[index]          = prims.vU[0];
          vy[index]          = prims.vU[1];
          vz[index]          = prims.vU[2];
          Y_e[index]         = prims.Y_e;
          temperature[index] = prims.temperature;

          rho_star[index] = cons.rho;
          tau[index]      = cons.tau;
          mhd_st_x[index] = cons.SD[0];
          mhd_st_y[index] = cons.SD[1];
          mhd_st_z[index] = cons.SD[2];
          Ye_star[index]  = cons.Y_e;

          // removed setTmunu from here
        }
      }
    }
  }
}
