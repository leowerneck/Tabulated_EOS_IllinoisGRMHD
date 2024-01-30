
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "GRHayLib.h"

void IllinoisGRMHD_Prim2Con(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int imax = cctk_lsh[0];
  const int jmax = cctk_lsh[1];
  const int kmax = cctk_lsh[2];

#pragma omp parallel for
  for(int k=0; k<kmax; k++) {
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        ghl_metric_quantities ADM_metric;
        ghl_initialize_metric(
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
        prims.entropy     = entropy[index];
        prims.Y_e         = Y_e[index];
        prims.temperature = temperature[index];

        ghl_conservative_quantities cons;
        const int speed_limited CCTK_ATTRIBUTE_UNUSED = ghl_enforce_primitive_limits_and_compute_u0(
              ghl_params, ghl_eos, &ADM_metric, &prims);
        ghl_compute_conservs(
              &ADM_metric, &metric_aux, &prims, &cons);

        rho[index]         = prims.rho;
        press[index]       = prims.press;
        eps[index]         = prims.eps;
        vx[index]          = prims.vU[0];
        vy[index]          = prims.vU[1];
        vz[index]          = prims.vU[2];
        entropy[index]     = prims.entropy;
        Y_e[index]         = prims.Y_e;
        temperature[index] = prims.temperature;

        rho_star[index] = cons.rho;
        tau[index]      = cons.tau;
        mhd_st_x[index]  = cons.SD[0];
        mhd_st_y[index]  = cons.SD[1];
        mhd_st_z[index]  = cons.SD[2];
        S_star[index] = cons.entropy;
        Ye_star[index]  = cons.Y_e;
      }
    }
  }
}
