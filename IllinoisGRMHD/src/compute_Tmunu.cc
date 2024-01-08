#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "GRHayLib.h"

void IllinoisGRMHD_compute_Tmunu(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel for
  for (int k = 0; k < cctk_lsh[2]; k++) {
    for (int j = 0; j < cctk_lsh[1]; j++) {
      for (int i = 0; i < cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        // Read in ADM metric quantities from gridfunctions and
        // set auxiliary and ADM metric quantities
        ghl_metric_quantities ADM_metric;
        ghl_enforce_detgtij_and_initialize_ADM_metric(
            alp[index], betax[index], betay[index], betaz[index], gxx[index],
            gxy[index], gxz[index], gyy[index], gyz[index], gzz[index],
            &ADM_metric);

        ghl_ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        // Read in primitive variables from gridfunctions
        ghl_primitive_quantities prims;
        prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0;
        prims.rho = rho[index];
        prims.press = press[index];
        prims.eps = eps[index];
        prims.vU[0] = vx[index];
        prims.vU[1] = vy[index];
        prims.vU[2] = vz[index];

        ghl_enforce_primitive_limits_and_compute_u0(ghl_params, ghl_eos,
                                                    &ADM_metric, &prims);

        ghl_stress_energy Tmunu;
        ghl_compute_TDNmunu(&ADM_metric, &metric_aux, &prims, &Tmunu);

        eTtt[index] += Tmunu.T4[0][0];
        eTtx[index] += Tmunu.T4[0][1];
        eTty[index] += Tmunu.T4[0][2];
        eTtz[index] += Tmunu.T4[0][3];
        eTxx[index] += Tmunu.T4[1][1];
        eTxy[index] += Tmunu.T4[1][2];
        eTxz[index] += Tmunu.T4[1][3];
        eTyy[index] += Tmunu.T4[2][2];
        eTyz[index] += Tmunu.T4[2][3];
        eTzz[index] += Tmunu.T4[3][3];
      }
    }
  }
}
