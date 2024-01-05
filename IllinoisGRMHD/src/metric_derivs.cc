#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void compute_four_metric_from_adm_quantities(
    const int index, const double *alpha, const double *betax,
    const double *betay, const double *betaz, const double *gxx,
    const double *gxy, const double *gxz, const double *gyy, const double *gyz,
    const double *gzz, double g4DD[4][4]) {

  const double alpha_sqr = alpha[index] * alpha[index];
  const double betaU[3] = {betax[index], betay[index], betaz[index]};
  const double gammaDD[3][3] = {
      gxx[index], gxy[index], gxz[index], gxy[index], gyy[index],
      gyz[index], gxz[index], gyz[index], gzz[index],
  };

  double betaD[3] = {0, 0, 0};
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      betaD[i] += gammaDD[i][j] * betaU[j];
    }
  }

  double beta_sqr = 0;
  for (int i = 0; i < 3; i++) {
    beta_sqr += betaU[i] * betaD[i];
  }

  g4DD[0][0] = -alpha_sqr + beta_sqr;
  for (int i = 1; i < 4; i++) {
    g4DD[0][i] = g4DD[i][0] = betaD[i - 1];
    for (int j = 1; j < 4; j++) {
      g4DD[i][j] = gammaDD[i - 1][j - 1];
    }
  }
}

void IllinoisGRMHD_compute_four_metric_time_derivatives(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int imin = cctk_nghostzones[0];
  const int imax = cctk_lsh[0] - imin;
  const int jmin = cctk_nghostzones[1];
  const int jmax = cctk_lsh[1] - jmin;
  const int kmin = cctk_nghostzones[2];
  const int kmax = cctk_lsh[2];

  for (int k = kmin; k < kmax; k++) {
    for (int j = jmin; j < jmax; j++) {
      for (int i = imin; i < imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
        const double inv_dt = 1.0 / CCTK_DELTA_TIME;

        double g4DD[4][4];
        compute_four_metric_from_adm_quantities(index, alp, betax, betay, betaz,
                                                gxx, gxy, gxz, gyy, gyz, gzz,
                                                g4DD);

        double g4DD_p[4][4];
        compute_four_metric_from_adm_quantities(index, alp_p, betax_p, betay_p,
                                                betaz_p, gxx_p, gxy_p, gxz_p,
                                                gyy_p, gyz_p, gzz_p, g4DD_p);

        gtt_dt[index] = (g4DD[0][0] - g4DD_p[0][0]) * inv_dt;
        gtx_dt[index] = (g4DD[0][1] - g4DD_p[0][1]) * inv_dt;
        gty_dt[index] = (g4DD[0][2] - g4DD_p[0][2]) * inv_dt;
        gtz_dt[index] = (g4DD[0][3] - g4DD_p[0][3]) * inv_dt;
        gxx_dt[index] = (g4DD[1][1] - g4DD_p[1][1]) * inv_dt;
        gxy_dt[index] = (g4DD[1][2] - g4DD_p[1][2]) * inv_dt;
        gxz_dt[index] = (g4DD[1][3] - g4DD_p[1][3]) * inv_dt;
        gyy_dt[index] = (g4DD[2][2] - g4DD_p[2][2]) * inv_dt;
        gyz_dt[index] = (g4DD[2][3] - g4DD_p[2][3]) * inv_dt;
        gzz_dt[index] = (g4DD[3][3] - g4DD_p[3][3]) * inv_dt;
      }
    }
  }
}
