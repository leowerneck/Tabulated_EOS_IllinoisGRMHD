#ifndef INTEGRANDS__VOLUMEINTEGRALS_VACUUM__
#define INTEGRANDS__VOLUMEINTEGRALS_VACUUM__

/* Integrand for L2 norms */
void VolumeIntegrals_vacuum_L2_integrand(double *VolIntegrand1, int index,double *f,double *x,double *y,double *z) {
  double fL = f[index];
  VolIntegrand1[index] = fL*fL;
}


/* Center of Lapse: */
void VolumeIntegrals_vacuum_CoL_integrand(double *VolIntegrand1,double *VolIntegrand2,double *VolIntegrand3,double *VolIntegrand4, int index,double *lapse,double *x,double *y,double *z) {
  double one_minus_lapseL = pow(1.0 - lapse[index],80); // <- Yields pretty consistent results with CoM integrand.
  VolIntegrand1[index] = one_minus_lapseL*x[index];
  VolIntegrand2[index] = one_minus_lapseL*y[index];
  VolIntegrand3[index] = one_minus_lapseL*z[index];
  VolIntegrand4[index] = one_minus_lapseL;
}

/* ADM Mass */
void VolumeIntegrals_vacuum_ADM_Mass_integrand_eval_derivs(double *ADM_M_integrand_x, double *ADM_M_integrand_y, double *ADM_M_integrand_z, const cGH *cctkGH, int i,int j,int k, double idx,double idy,double idz, double *alp,double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz) {
  const CCTK_REAL cm1 = -0.5;

  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

  // Read in \gamma_{i j}
  CCTK_REAL g11L = gxx[index];
  CCTK_REAL g12L = gxy[index];
  CCTK_REAL g13L = gxz[index];
  CCTK_REAL g22L = gyy[index];
  CCTK_REAL g23L = gyz[index];
  CCTK_REAL g33L = gzz[index];

  // Metric determinant
  CCTK_REAL detgL = -g13L * g13L * g22L + 2 * g12L * g13L * g23L - g11L * g23L * g23L - g12L * g12L * g33L + g11L * g22L * g33L;

  CCTK_REAL ginv[3][3];

  // Calculate inverse metric \gamma^{i j}
  ginv[0][0] = (g22L * g33L - g23L * g23L) / detgL;
  ginv[0][1] = (g13L * g23L - g12L * g33L) / detgL;
  ginv[0][2] = (g12L * g23L - g13L * g22L) / detgL;
  ginv[1][1] = (g11L * g33L - g13L * g13L) / detgL;
  ginv[1][2] = (g12L * g13L - g11L * g23L) / detgL;
  ginv[2][2] = (g11L * g22L - g12L * g12L) / detgL;

  ginv[1][0] = ginv[0][1];
  ginv[2][0] = ginv[0][2];
  ginv[2][1] = ginv[1][2];

  CCTK_REAL g_d1[3][3][3]; // g_d1[i][j][k] = d_i g_{j k}

  g_d1[0][0][0] = (cm1 * (gxx[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] - gxx[CCTK_GFINDEX3D(cctkGH,i+1,j,k)])) * idx;
  g_d1[0][0][1] = (cm1 * (gxy[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] - gxy[CCTK_GFINDEX3D(cctkGH,i+1,j,k)])) * idx;
  g_d1[0][0][2] = (cm1 * (gxz[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] - gxz[CCTK_GFINDEX3D(cctkGH,i+1,j,k)])) * idx;
  g_d1[0][1][1] = (cm1 * (gyy[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] - gyy[CCTK_GFINDEX3D(cctkGH,i+1,j,k)])) * idx;
  g_d1[0][1][2] = (cm1 * (gyz[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] - gyz[CCTK_GFINDEX3D(cctkGH,i+1,j,k)])) * idx;
  g_d1[0][2][2] = (cm1 * (gzz[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] - gzz[CCTK_GFINDEX3D(cctkGH,i+1,j,k)])) * idx;

  g_d1[0][1][0] = g_d1[0][0][1];
  g_d1[0][2][0] = g_d1[0][0][2];
  g_d1[0][2][1] = g_d1[0][1][2];

  g_d1[1][0][0] = (cm1 * (gxx[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] - gxx[CCTK_GFINDEX3D(cctkGH,i,j+1,k)])) * idy;
  g_d1[1][0][1] = (cm1 * (gxy[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] - gxy[CCTK_GFINDEX3D(cctkGH,i,j+1,k)])) * idy;
  g_d1[1][0][2] = (cm1 * (gxz[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] - gxz[CCTK_GFINDEX3D(cctkGH,i,j+1,k)])) * idy;
  g_d1[1][1][1] = (cm1 * (gyy[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] - gyy[CCTK_GFINDEX3D(cctkGH,i,j+1,k)])) * idy;
  g_d1[1][1][2] = (cm1 * (gyz[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] - gyz[CCTK_GFINDEX3D(cctkGH,i,j+1,k)])) * idy;
  g_d1[1][2][2] = (cm1 * (gzz[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] - gzz[CCTK_GFINDEX3D(cctkGH,i,j+1,k)])) * idy;

  g_d1[1][1][0] = g_d1[1][0][1];
  g_d1[1][2][0] = g_d1[1][0][2];
  g_d1[1][2][1] = g_d1[1][1][2];

  g_d1[2][0][0] = (cm1 * (gxx[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] - gxx[CCTK_GFINDEX3D(cctkGH,i,j,k+1)])) * idz;
  g_d1[2][0][1] = (cm1 * (gxy[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] - gxy[CCTK_GFINDEX3D(cctkGH,i,j,k+1)])) * idz;
  g_d1[2][0][2] = (cm1 * (gxz[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] - gxz[CCTK_GFINDEX3D(cctkGH,i,j,k+1)])) * idz;
  g_d1[2][1][1] = (cm1 * (gyy[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] - gyy[CCTK_GFINDEX3D(cctkGH,i,j,k+1)])) * idz;
  g_d1[2][1][2] = (cm1 * (gyz[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] - gyz[CCTK_GFINDEX3D(cctkGH,i,j,k+1)])) * idz;
  g_d1[2][2][2] = (cm1 * (gzz[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] - gzz[CCTK_GFINDEX3D(cctkGH,i,j,k+1)])) * idz;

  g_d1[2][1][0] = g_d1[2][0][1];
  g_d1[2][2][0] = g_d1[2][0][2];
  g_d1[2][2][1] = g_d1[2][1][2];

  ADM_M_integrand_x[index] = 0;
  ADM_M_integrand_y[index] = 0;
  ADM_M_integrand_z[index] = 0;

  for(int i2 = 0; i2 < 3; i2++)
    {
      for(int j2 = 0; j2 < 3; j2++)
        {
          for(int k2 = 0; k2 < 3; k2++)
            {
              ADM_M_integrand_x[index] += ginv[i2][j2] * ginv[k2][0] * (g_d1[j2][i2][k2] - g_d1[k2][i2][j2]);
              ADM_M_integrand_y[index] += ginv[i2][j2] * ginv[k2][1] * (g_d1[j2][i2][k2] - g_d1[k2][i2][j2]);
              ADM_M_integrand_z[index] += ginv[i2][j2] * ginv[k2][2] * (g_d1[j2][i2][k2] - g_d1[k2][i2][j2]);
            }
        }
    }

  ADM_M_integrand_x[index] *= alp[index] * sqrt(detgL);
  ADM_M_integrand_y[index] *= alp[index] * sqrt(detgL);
  ADM_M_integrand_z[index] *= alp[index] * sqrt(detgL);
}

void VolumeIntegrals_vacuum_ADM_Mass_integrand(double *ADM_M_integrand, const cGH *cctkGH, int i,int j,int k, double idx,double idy,double idz, double *ADM_M_integrand_x, double *ADM_M_integrand_y, double *ADM_M_integrand_z) {
  const CCTK_REAL cm1 = -0.5;

  ADM_M_integrand[CCTK_GFINDEX3D(cctkGH,i,j,k)] = 0.0625 / M_PI *
    (
     (cm1 * (ADM_M_integrand_x[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] - ADM_M_integrand_x[CCTK_GFINDEX3D(cctkGH,i+1,j,k)])) * idx +
     (cm1 * (ADM_M_integrand_y[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] - ADM_M_integrand_y[CCTK_GFINDEX3D(cctkGH,i,j+1,k)])) * idy +
     (cm1 * (ADM_M_integrand_z[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] - ADM_M_integrand_z[CCTK_GFINDEX3D(cctkGH,i,j,k+1)])) * idz
     );
}

/* ADM Momentum */
void VolumeIntegrals_vacuum_ADM_Momentum_integrand_eval_derivs(double *ADM_Py_integrand_x, double *ADM_Py_integrand_y, double *ADM_Py_integrand_z, const cGH *cctkGH, int i,int j,int k, double idx,double idy,double idz, double *alp,double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,double *kzz) {

  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

  // Read in \gamma_{i j}
  CCTK_REAL g11L = gxx[index];
  CCTK_REAL g12L = gxy[index];
  CCTK_REAL g13L = gxz[index];
  CCTK_REAL g22L = gyy[index];
  CCTK_REAL g23L = gyz[index];
  CCTK_REAL g33L = gzz[index];

  // Metric determinant
  CCTK_REAL detgL = -g13L * g13L * g22L + 2 * g12L * g13L * g23L - g11L * g23L * g23L - g12L * g12L * g33L + g11L * g22L * g33L;

  CCTK_REAL ginv[3][3];

  // Calculate inverse metric \gamma^{i j}
  ginv[0][0] = (g22L * g33L - g23L * g23L) / detgL;
  ginv[0][1] = (g13L * g23L - g12L * g33L) / detgL;
  ginv[0][2] = (g12L * g23L - g13L * g22L) / detgL;
  ginv[1][1] = (g11L * g33L - g13L * g13L) / detgL;
  ginv[1][2] = (g12L * g13L - g11L * g23L) / detgL;
  ginv[2][2] = (g11L * g22L - g12L * g12L) / detgL;

  ginv[1][0] = ginv[0][1];
  ginv[2][0] = ginv[0][2];
  ginv[2][1] = ginv[1][2];

  CCTK_REAL Kdown[3][3];

  // Read in covariant extrinsic curvature K_{i j}
  Kdown[0][0] = kxx[index];
  Kdown[0][1] = kxy[index];
  Kdown[0][2] = kxz[index];
  Kdown[1][1] = kyy[index];
  Kdown[1][2] = kyz[index];
  Kdown[2][2] = kzz[index];

  Kdown[1][0] = Kdown[0][1];
  Kdown[2][0] = Kdown[0][2];
  Kdown[2][1] = Kdown[1][2];

  CCTK_REAL Kup[3][3];

  // Calculate contravariant extrinsic curvature K^{i j}
  for(int i2 = 0; i2 < 3; i2++)
    {
      for(int j2 = 0; j2 < 3; j2++)
        {
	  Kup[i2][j2] = 0;

          for(int k2 = 0; k2 < 3; k2++)
            {
	      for(int l2 = 0; l2 < 3; l2++)
		{
		  Kup[i2][j2] += ginv[i2][k2] * ginv[j2][l2] * Kdown[k2][l2];
		}
	    }
	}
    }

  CCTK_REAL K = 0;

  // Calculate mean curvature K = g^{i j} K_{i j}
  for(int i2 = 0; i2 < 3; i2++)
    {
      for(int j2 = 0; j2 < 3; j2++)
        {
	  K += ginv[i2][j2] * Kdown[i2][j2];
	}
    }

  CCTK_REAL surface_integrand[3][3];

  for(int i2 = 0; i2 < 3; i2++)
    {
      for(int j2 = 0; j2 < 3; j2++)
        {
	  surface_integrand[i2][j2] = alp[index] * sqrt(detgL) * (Kup[i2][j2] - ginv[i2][j2] * K);
	}
    }

  ADM_Py_integrand_x[index] = surface_integrand[1][0];
  ADM_Py_integrand_y[index] = surface_integrand[1][1];
  ADM_Py_integrand_z[index] = surface_integrand[1][2];
}

void VolumeIntegrals_vacuum_ADM_Momentum_integrand(double *ADM_Py_integrand, const cGH *cctkGH, int i,int j,int k, double idx,double idy,double idz, double *ADM_Py_integrand_x, double *ADM_Py_integrand_y, double *ADM_Py_integrand_z) {
  const CCTK_REAL cm1 = -0.5;

  ADM_Py_integrand[CCTK_GFINDEX3D(cctkGH,i,j,k)] = 0.125 / M_PI *
    (
     (cm1 * (ADM_Py_integrand_x[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] - ADM_Py_integrand_x[CCTK_GFINDEX3D(cctkGH,i+1,j,k)])) * idx +
     (cm1 * (ADM_Py_integrand_y[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] - ADM_Py_integrand_y[CCTK_GFINDEX3D(cctkGH,i,j+1,k)])) * idy +
     (cm1 * (ADM_Py_integrand_z[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] - ADM_Py_integrand_z[CCTK_GFINDEX3D(cctkGH,i,j,k+1)])) * idz
     );
}

/* ADM Angular Momentum */
void VolumeIntegrals_vacuum_ADM_Angular_Momentum_integrand_eval_derivs(double *ADM_Jz_integrand_x, double *ADM_Jz_integrand_y, double *ADM_Jz_integrand_z, const cGH *cctkGH, int i,int j,int k, double idx,double idy,double idz, double *x,double *y,double *z,double *alp,double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,double *kzz) {

  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

  // Read in \gamma_{i j}
  CCTK_REAL g11L = gxx[index];
  CCTK_REAL g12L = gxy[index];
  CCTK_REAL g13L = gxz[index];
  CCTK_REAL g22L = gyy[index];
  CCTK_REAL g23L = gyz[index];
  CCTK_REAL g33L = gzz[index];

  // Metric determinant
  CCTK_REAL detgL = -g13L * g13L * g22L + 2 * g12L * g13L * g23L - g11L * g23L * g23L - g12L * g12L * g33L + g11L * g22L * g33L;

  CCTK_REAL ginv[3][3];

  // Calculate inverse metric \gamma^{i j}
  ginv[0][0] = (g22L * g33L - g23L * g23L) / detgL;
  ginv[0][1] = (g13L * g23L - g12L * g33L) / detgL;
  ginv[0][2] = (g12L * g23L - g13L * g22L) / detgL;
  ginv[1][1] = (g11L * g33L - g13L * g13L) / detgL;
  ginv[1][2] = (g12L * g13L - g11L * g23L) / detgL;
  ginv[2][2] = (g11L * g22L - g12L * g12L) / detgL;

  ginv[1][0] = ginv[0][1];
  ginv[2][0] = ginv[0][2];
  ginv[2][1] = ginv[1][2];

  CCTK_REAL Kdown[3][3];

  // Read in covariant extrinsic curvature K_{i j}
  Kdown[0][0] = kxx[index];
  Kdown[0][1] = kxy[index];
  Kdown[0][2] = kxz[index];
  Kdown[1][1] = kyy[index];
  Kdown[1][2] = kyz[index];
  Kdown[2][2] = kzz[index];

  Kdown[1][0] = Kdown[0][1];
  Kdown[2][0] = Kdown[0][2];
  Kdown[2][1] = Kdown[1][2];

  CCTK_REAL Kup[3][3];

  // Calculate contravariant extrinsic curvature K^{i j}
  for(int i2 = 0; i2 < 3; i2++)
    {
      for(int j2 = 0; j2 < 3; j2++)
        {
	  Kup[i2][j2] = 0;

          for(int k2 = 0; k2 < 3; k2++)
            {
	      for(int l2 = 0; l2 < 3; l2++)
		{
		  Kup[i2][j2] += ginv[i2][k2] * ginv[j2][l2] * Kdown[k2][l2];
		}
	    }
	}
    }

  CCTK_REAL K = 0;

  // Calculate mean curvature K = g^{i j} K_{i j}
  for(int i2 = 0; i2 < 3; i2++)
    {
      for(int j2 = 0; j2 < 3; j2++)
        {
	  K += ginv[i2][j2] * Kdown[i2][j2];
	}
    }

  // Levi-Civita tensor
  CCTK_REAL LCT[3][3][3] = {{{0, 0, 0}, {0, 0, 1}, {0, -1, 0}}, {{0, 0, -1}, {0, 0, 0}, {1, 0, 0}}, {{0, 1, 0}, {-1, 0, 0}, {0, 0, 0}}};

  // Coordinate vector
  CCTK_REAL xx[3] = {x[index], y[index], z[index]};

  CCTK_REAL surface_integrand[3][3];

  for(int i2 = 0; i2 < 3; i2++)
    {
      for(int j2 = 0; j2 < 3; j2++)
        {
	  surface_integrand[i2][j2] = 0;

          for(int k2 = 0; k2 < 3; k2++)
            {
	      for(int l2 = 0; l2 < 3; l2++)
		{
		  surface_integrand[i2][j2] += LCT[i2][k2][l2] * xx[k2] * (Kup[l2][j2] - ginv[l2][j2] * K);
		}
	    }

	  surface_integrand[i2][j2] *= alp[index] * sqrt(detgL);
	}
    }

  ADM_Jz_integrand_x[index] = surface_integrand[2][0];
  ADM_Jz_integrand_y[index] = surface_integrand[2][1];
  ADM_Jz_integrand_z[index] = surface_integrand[2][2];
}

void VolumeIntegrals_vacuum_ADM_Angular_Momentum_integrand(double *ADM_Jz_integrand, const cGH *cctkGH, int i,int j,int k, double idx,double idy,double idz, double *ADM_Jz_integrand_x, double *ADM_Jz_integrand_y, double *ADM_Jz_integrand_z) {
  const CCTK_REAL cm1 = -0.5;

  ADM_Jz_integrand[CCTK_GFINDEX3D(cctkGH,i,j,k)] = 0.125 / M_PI *
    (
     (cm1 * (ADM_Jz_integrand_x[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] - ADM_Jz_integrand_x[CCTK_GFINDEX3D(cctkGH,i+1,j,k)])) * idx +
     (cm1 * (ADM_Jz_integrand_y[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] - ADM_Jz_integrand_y[CCTK_GFINDEX3D(cctkGH,i,j+1,k)])) * idy +
     (cm1 * (ADM_Jz_integrand_z[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] - ADM_Jz_integrand_z[CCTK_GFINDEX3D(cctkGH,i,j,k+1)])) * idz
     );
}

#endif // INTEGRANDS__VOLUMEINTEGRALS_VACUUM__
