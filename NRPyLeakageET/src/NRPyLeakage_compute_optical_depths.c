#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "NRPyLeakage.h"


/*
 * (c) Leo Werneck
 * Compute GRMHD source terms following Ruffert et al. (1996)
 * https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
void NRPyLeakage_compute_optical_depths(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(0);

  // Step 0: Loop over the grid computing the optical depth
#pragma omp parallel for
  for(int k=cctk_nghostzones[2];k<cctk_lsh[2]-cctk_nghostzones[2];k++) {
    for(int j=cctk_nghostzones[1];j<cctk_lsh[1]-cctk_nghostzones[1];j++) {
      for(int i=cctk_nghostzones[0];i<cctk_lsh[0]-cctk_nghostzones[0];i++) {

        // Step 1: Set gridpoint indices
        const int i_j_k   = CCTK_GFINDEX3D(cctkGH,i  ,j,k  );
        const int ip1_j_k = CCTK_GFINDEX3D(cctkGH,i+1,j,k  );
        const int im1_j_k = CCTK_GFINDEX3D(cctkGH,i-1,j,k  );
        const int i_jp1_k = CCTK_GFINDEX3D(cctkGH,i,j+1,k  );
        const int i_jm1_k = CCTK_GFINDEX3D(cctkGH,i,j-1,k  );
        const int i_j_kp1 = CCTK_GFINDEX3D(cctkGH,i,j  ,k+1);
        const int i_j_km1 = CCTK_GFINDEX3D(cctkGH,i,j  ,k-1);

        // Step 2: Read in metric gfs from main memory
        const CCTK_REAL gxx_i_j_k = gxx[i_j_k];
        const CCTK_REAL gxx_ip1_j_k = gxx[ip1_j_k];
        const CCTK_REAL gxx_im1_j_k = gxx[im1_j_k];

        const CCTK_REAL gyy_i_j_k = gyy[i_j_k];
        const CCTK_REAL gyy_i_jp1_k = gyy[i_jp1_k];
        const CCTK_REAL gyy_i_jm1_k = gyy[i_jm1_k];

        const CCTK_REAL gzz_i_j_k = gzz[i_j_k];
        const CCTK_REAL gzz_i_j_kp1 = gzz[i_j_kp1];
        const CCTK_REAL gzz_i_j_km1 = gzz[i_j_km1];


        // Step 3: Read in opacity gfs from main memory
        const CCTK_REAL kappa_0_nue_i_j_k = kappa_0_nue[i_j_k];
        const CCTK_REAL kappa_0_nue_ip1_j_k = kappa_0_nue[ip1_j_k];
        const CCTK_REAL kappa_0_nue_im1_j_k = kappa_0_nue[im1_j_k];
        const CCTK_REAL kappa_0_nue_i_jp1_k = kappa_0_nue[i_jp1_k];
        const CCTK_REAL kappa_0_nue_i_jm1_k = kappa_0_nue[i_jm1_k];
        const CCTK_REAL kappa_0_nue_i_j_kp1 = kappa_0_nue[i_j_kp1];
        const CCTK_REAL kappa_0_nue_i_j_km1 = kappa_0_nue[i_j_km1];

        const CCTK_REAL kappa_1_nue_i_j_k = kappa_1_nue[i_j_k];
        const CCTK_REAL kappa_1_nue_ip1_j_k = kappa_1_nue[ip1_j_k];
        const CCTK_REAL kappa_1_nue_im1_j_k = kappa_1_nue[im1_j_k];
        const CCTK_REAL kappa_1_nue_i_jp1_k = kappa_1_nue[i_jp1_k];
        const CCTK_REAL kappa_1_nue_i_jm1_k = kappa_1_nue[i_jm1_k];
        const CCTK_REAL kappa_1_nue_i_j_kp1 = kappa_1_nue[i_j_kp1];
        const CCTK_REAL kappa_1_nue_i_j_km1 = kappa_1_nue[i_j_km1];

        const CCTK_REAL kappa_0_anue_i_j_k = kappa_0_anue[i_j_k];
        const CCTK_REAL kappa_0_anue_ip1_j_k = kappa_0_anue[ip1_j_k];
        const CCTK_REAL kappa_0_anue_im1_j_k = kappa_0_anue[im1_j_k];
        const CCTK_REAL kappa_0_anue_i_jp1_k = kappa_0_anue[i_jp1_k];
        const CCTK_REAL kappa_0_anue_i_jm1_k = kappa_0_anue[i_jm1_k];
        const CCTK_REAL kappa_0_anue_i_j_kp1 = kappa_0_anue[i_j_kp1];
        const CCTK_REAL kappa_0_anue_i_j_km1 = kappa_0_anue[i_j_km1];

        const CCTK_REAL kappa_1_anue_i_j_k = kappa_1_anue[i_j_k];
        const CCTK_REAL kappa_1_anue_ip1_j_k = kappa_1_anue[ip1_j_k];
        const CCTK_REAL kappa_1_anue_im1_j_k = kappa_1_anue[im1_j_k];
        const CCTK_REAL kappa_1_anue_i_jp1_k = kappa_1_anue[i_jp1_k];
        const CCTK_REAL kappa_1_anue_i_jm1_k = kappa_1_anue[i_jm1_k];
        const CCTK_REAL kappa_1_anue_i_j_kp1 = kappa_1_anue[i_j_kp1];
        const CCTK_REAL kappa_1_anue_i_j_km1 = kappa_1_anue[i_j_km1];

        const CCTK_REAL kappa_0_nux_i_j_k = kappa_0_nux[i_j_k];
        const CCTK_REAL kappa_0_nux_ip1_j_k = kappa_0_nux[ip1_j_k];
        const CCTK_REAL kappa_0_nux_im1_j_k = kappa_0_nux[im1_j_k];
        const CCTK_REAL kappa_0_nux_i_jp1_k = kappa_0_nux[i_jp1_k];
        const CCTK_REAL kappa_0_nux_i_jm1_k = kappa_0_nux[i_jm1_k];
        const CCTK_REAL kappa_0_nux_i_j_kp1 = kappa_0_nux[i_j_kp1];
        const CCTK_REAL kappa_0_nux_i_j_km1 = kappa_0_nux[i_j_km1];

        const CCTK_REAL kappa_1_nux_i_j_k = kappa_1_nux[i_j_k];
        const CCTK_REAL kappa_1_nux_ip1_j_k = kappa_1_nux[ip1_j_k];
        const CCTK_REAL kappa_1_nux_im1_j_k = kappa_1_nux[im1_j_k];
        const CCTK_REAL kappa_1_nux_i_jp1_k = kappa_1_nux[i_jp1_k];
        const CCTK_REAL kappa_1_nux_i_jm1_k = kappa_1_nux[i_jm1_k];
        const CCTK_REAL kappa_1_nux_i_j_kp1 = kappa_1_nux[i_j_kp1];
        const CCTK_REAL kappa_1_nux_i_j_km1 = kappa_1_nux[i_j_km1];


        // Step 4: Compute metric at cell faces
        const CCTK_REAL gxx_iphalf_j_k = 0.5*(gxx_i_j_k + gxx_ip1_j_k);
        const CCTK_REAL gxx_imhalf_j_k = 0.5*(gxx_i_j_k + gxx_im1_j_k);

        const CCTK_REAL gyy_i_jphalf_k = 0.5*(gyy_i_j_k + gyy_i_jp1_k);
        const CCTK_REAL gyy_i_jmhalf_k = 0.5*(gyy_i_j_k + gyy_i_jm1_k);

        const CCTK_REAL gzz_i_j_kphalf = 0.5*(gzz_i_j_k + gzz_i_j_kp1);
        const CCTK_REAL gzz_i_j_kmhalf = 0.5*(gzz_i_j_k + gzz_i_j_km1);


        // Step 5: Compute ds^{i} = sqrt(gamma_{ii}dx^{i}dx^{i})
        const CCTK_REAL ds_iphalf_j_k = sqrt(dx*dx*gxx_iphalf_j_k);
        const CCTK_REAL ds_imhalf_j_k = sqrt(dx*dx*gxx_imhalf_j_k);
        const CCTK_REAL ds_i_jphalf_k = sqrt(dy*dy*gyy_i_jphalf_k);
        const CCTK_REAL ds_i_jmhalf_k = sqrt(dy*dy*gyy_i_jmhalf_k);
        const CCTK_REAL ds_i_j_kphalf = sqrt(dz*dz*gzz_i_j_kphalf);
        const CCTK_REAL ds_i_j_kmhalf = sqrt(dz*dz*gzz_i_j_kmhalf);

        // Step 6: Compute opacities at cell faces
        const CCTK_REAL kappa_0_nue_iphalf_j_k = 0.5*(kappa_0_nue_i_j_k + kappa_0_nue_ip1_j_k);
        const CCTK_REAL kappa_0_nue_imhalf_j_k = 0.5*(kappa_0_nue_i_j_k + kappa_0_nue_im1_j_k);
        const CCTK_REAL kappa_0_nue_i_jphalf_k = 0.5*(kappa_0_nue_i_j_k + kappa_0_nue_i_jp1_k);
        const CCTK_REAL kappa_0_nue_i_jmhalf_k = 0.5*(kappa_0_nue_i_j_k + kappa_0_nue_i_jm1_k);
        const CCTK_REAL kappa_0_nue_i_j_kphalf = 0.5*(kappa_0_nue_i_j_k + kappa_0_nue_i_j_kp1);
        const CCTK_REAL kappa_0_nue_i_j_kmhalf = 0.5*(kappa_0_nue_i_j_k + kappa_0_nue_i_j_km1);

        const CCTK_REAL kappa_1_nue_iphalf_j_k = 0.5*(kappa_1_nue_i_j_k + kappa_1_nue_ip1_j_k);
        const CCTK_REAL kappa_1_nue_imhalf_j_k = 0.5*(kappa_1_nue_i_j_k + kappa_1_nue_im1_j_k);
        const CCTK_REAL kappa_1_nue_i_jphalf_k = 0.5*(kappa_1_nue_i_j_k + kappa_1_nue_i_jp1_k);
        const CCTK_REAL kappa_1_nue_i_jmhalf_k = 0.5*(kappa_1_nue_i_j_k + kappa_1_nue_i_jm1_k);
        const CCTK_REAL kappa_1_nue_i_j_kphalf = 0.5*(kappa_1_nue_i_j_k + kappa_1_nue_i_j_kp1);
        const CCTK_REAL kappa_1_nue_i_j_kmhalf = 0.5*(kappa_1_nue_i_j_k + kappa_1_nue_i_j_km1);

        const CCTK_REAL kappa_0_anue_iphalf_j_k = 0.5*(kappa_0_anue_i_j_k + kappa_0_anue_ip1_j_k);
        const CCTK_REAL kappa_0_anue_imhalf_j_k = 0.5*(kappa_0_anue_i_j_k + kappa_0_anue_im1_j_k);
        const CCTK_REAL kappa_0_anue_i_jphalf_k = 0.5*(kappa_0_anue_i_j_k + kappa_0_anue_i_jp1_k);
        const CCTK_REAL kappa_0_anue_i_jmhalf_k = 0.5*(kappa_0_anue_i_j_k + kappa_0_anue_i_jm1_k);
        const CCTK_REAL kappa_0_anue_i_j_kphalf = 0.5*(kappa_0_anue_i_j_k + kappa_0_anue_i_j_kp1);
        const CCTK_REAL kappa_0_anue_i_j_kmhalf = 0.5*(kappa_0_anue_i_j_k + kappa_0_anue_i_j_km1);

        const CCTK_REAL kappa_1_anue_iphalf_j_k = 0.5*(kappa_1_anue_i_j_k + kappa_1_anue_ip1_j_k);
        const CCTK_REAL kappa_1_anue_imhalf_j_k = 0.5*(kappa_1_anue_i_j_k + kappa_1_anue_im1_j_k);
        const CCTK_REAL kappa_1_anue_i_jphalf_k = 0.5*(kappa_1_anue_i_j_k + kappa_1_anue_i_jp1_k);
        const CCTK_REAL kappa_1_anue_i_jmhalf_k = 0.5*(kappa_1_anue_i_j_k + kappa_1_anue_i_jm1_k);
        const CCTK_REAL kappa_1_anue_i_j_kphalf = 0.5*(kappa_1_anue_i_j_k + kappa_1_anue_i_j_kp1);
        const CCTK_REAL kappa_1_anue_i_j_kmhalf = 0.5*(kappa_1_anue_i_j_k + kappa_1_anue_i_j_km1);

        const CCTK_REAL kappa_0_nux_iphalf_j_k = 0.5*(kappa_0_nux_i_j_k + kappa_0_nux_ip1_j_k);
        const CCTK_REAL kappa_0_nux_imhalf_j_k = 0.5*(kappa_0_nux_i_j_k + kappa_0_nux_im1_j_k);
        const CCTK_REAL kappa_0_nux_i_jphalf_k = 0.5*(kappa_0_nux_i_j_k + kappa_0_nux_i_jp1_k);
        const CCTK_REAL kappa_0_nux_i_jmhalf_k = 0.5*(kappa_0_nux_i_j_k + kappa_0_nux_i_jm1_k);
        const CCTK_REAL kappa_0_nux_i_j_kphalf = 0.5*(kappa_0_nux_i_j_k + kappa_0_nux_i_j_kp1);
        const CCTK_REAL kappa_0_nux_i_j_kmhalf = 0.5*(kappa_0_nux_i_j_k + kappa_0_nux_i_j_km1);

        const CCTK_REAL kappa_1_nux_iphalf_j_k = 0.5*(kappa_1_nux_i_j_k + kappa_1_nux_ip1_j_k);
        const CCTK_REAL kappa_1_nux_imhalf_j_k = 0.5*(kappa_1_nux_i_j_k + kappa_1_nux_im1_j_k);
        const CCTK_REAL kappa_1_nux_i_jphalf_k = 0.5*(kappa_1_nux_i_j_k + kappa_1_nux_i_jp1_k);
        const CCTK_REAL kappa_1_nux_i_jmhalf_k = 0.5*(kappa_1_nux_i_j_k + kappa_1_nux_i_jm1_k);
        const CCTK_REAL kappa_1_nux_i_j_kphalf = 0.5*(kappa_1_nux_i_j_k + kappa_1_nux_i_j_kp1);
        const CCTK_REAL kappa_1_nux_i_j_kmhalf = 0.5*(kappa_1_nux_i_j_k + kappa_1_nux_i_j_km1);


        // Step 7: Compute optical depth at neighboring points
        const CCTK_REAL tau_0_nue_ip1_j_k = tau_0_nue[ip1_j_k] + ds_iphalf_j_k*kappa_0_nue_iphalf_j_k;
        const CCTK_REAL tau_0_nue_im1_j_k = tau_0_nue[im1_j_k] + ds_imhalf_j_k*kappa_0_nue_imhalf_j_k;
        const CCTK_REAL tau_0_nue_i_jp1_k = tau_0_nue[i_jp1_k] + ds_i_jphalf_k*kappa_0_nue_i_jphalf_k;
        const CCTK_REAL tau_0_nue_i_jm1_k = tau_0_nue[i_jm1_k] + ds_i_jmhalf_k*kappa_0_nue_i_jmhalf_k;
        const CCTK_REAL tau_0_nue_i_j_kp1 = tau_0_nue[i_j_kp1] + ds_i_j_kphalf*kappa_0_nue_i_j_kphalf;
        const CCTK_REAL tau_0_nue_i_j_km1 = tau_0_nue[i_j_km1] + ds_i_j_kmhalf*kappa_0_nue_i_j_kmhalf;

        const CCTK_REAL tau_1_nue_ip1_j_k = tau_1_nue[ip1_j_k] + ds_iphalf_j_k*kappa_1_nue_iphalf_j_k;
        const CCTK_REAL tau_1_nue_im1_j_k = tau_1_nue[im1_j_k] + ds_imhalf_j_k*kappa_1_nue_imhalf_j_k;
        const CCTK_REAL tau_1_nue_i_jp1_k = tau_1_nue[i_jp1_k] + ds_i_jphalf_k*kappa_1_nue_i_jphalf_k;
        const CCTK_REAL tau_1_nue_i_jm1_k = tau_1_nue[i_jm1_k] + ds_i_jmhalf_k*kappa_1_nue_i_jmhalf_k;
        const CCTK_REAL tau_1_nue_i_j_kp1 = tau_1_nue[i_j_kp1] + ds_i_j_kphalf*kappa_1_nue_i_j_kphalf;
        const CCTK_REAL tau_1_nue_i_j_km1 = tau_1_nue[i_j_km1] + ds_i_j_kmhalf*kappa_1_nue_i_j_kmhalf;

        const CCTK_REAL tau_0_anue_ip1_j_k = tau_0_anue[ip1_j_k] + ds_iphalf_j_k*kappa_0_anue_iphalf_j_k;
        const CCTK_REAL tau_0_anue_im1_j_k = tau_0_anue[im1_j_k] + ds_imhalf_j_k*kappa_0_anue_imhalf_j_k;
        const CCTK_REAL tau_0_anue_i_jp1_k = tau_0_anue[i_jp1_k] + ds_i_jphalf_k*kappa_0_anue_i_jphalf_k;
        const CCTK_REAL tau_0_anue_i_jm1_k = tau_0_anue[i_jm1_k] + ds_i_jmhalf_k*kappa_0_anue_i_jmhalf_k;
        const CCTK_REAL tau_0_anue_i_j_kp1 = tau_0_anue[i_j_kp1] + ds_i_j_kphalf*kappa_0_anue_i_j_kphalf;
        const CCTK_REAL tau_0_anue_i_j_km1 = tau_0_anue[i_j_km1] + ds_i_j_kmhalf*kappa_0_anue_i_j_kmhalf;

        const CCTK_REAL tau_1_anue_ip1_j_k = tau_1_anue[ip1_j_k] + ds_iphalf_j_k*kappa_1_anue_iphalf_j_k;
        const CCTK_REAL tau_1_anue_im1_j_k = tau_1_anue[im1_j_k] + ds_imhalf_j_k*kappa_1_anue_imhalf_j_k;
        const CCTK_REAL tau_1_anue_i_jp1_k = tau_1_anue[i_jp1_k] + ds_i_jphalf_k*kappa_1_anue_i_jphalf_k;
        const CCTK_REAL tau_1_anue_i_jm1_k = tau_1_anue[i_jm1_k] + ds_i_jmhalf_k*kappa_1_anue_i_jmhalf_k;
        const CCTK_REAL tau_1_anue_i_j_kp1 = tau_1_anue[i_j_kp1] + ds_i_j_kphalf*kappa_1_anue_i_j_kphalf;
        const CCTK_REAL tau_1_anue_i_j_km1 = tau_1_anue[i_j_km1] + ds_i_j_kmhalf*kappa_1_anue_i_j_kmhalf;

        const CCTK_REAL tau_0_nux_ip1_j_k = tau_0_nux[ip1_j_k] + ds_iphalf_j_k*kappa_0_nux_iphalf_j_k;
        const CCTK_REAL tau_0_nux_im1_j_k = tau_0_nux[im1_j_k] + ds_imhalf_j_k*kappa_0_nux_imhalf_j_k;
        const CCTK_REAL tau_0_nux_i_jp1_k = tau_0_nux[i_jp1_k] + ds_i_jphalf_k*kappa_0_nux_i_jphalf_k;
        const CCTK_REAL tau_0_nux_i_jm1_k = tau_0_nux[i_jm1_k] + ds_i_jmhalf_k*kappa_0_nux_i_jmhalf_k;
        const CCTK_REAL tau_0_nux_i_j_kp1 = tau_0_nux[i_j_kp1] + ds_i_j_kphalf*kappa_0_nux_i_j_kphalf;
        const CCTK_REAL tau_0_nux_i_j_km1 = tau_0_nux[i_j_km1] + ds_i_j_kmhalf*kappa_0_nux_i_j_kmhalf;

        const CCTK_REAL tau_1_nux_ip1_j_k = tau_1_nux[ip1_j_k] + ds_iphalf_j_k*kappa_1_nux_iphalf_j_k;
        const CCTK_REAL tau_1_nux_im1_j_k = tau_1_nux[im1_j_k] + ds_imhalf_j_k*kappa_1_nux_imhalf_j_k;
        const CCTK_REAL tau_1_nux_i_jp1_k = tau_1_nux[i_jp1_k] + ds_i_jphalf_k*kappa_1_nux_i_jphalf_k;
        const CCTK_REAL tau_1_nux_i_jm1_k = tau_1_nux[i_jm1_k] + ds_i_jmhalf_k*kappa_1_nux_i_jmhalf_k;
        const CCTK_REAL tau_1_nux_i_j_kp1 = tau_1_nux[i_j_kp1] + ds_i_j_kphalf*kappa_1_nux_i_j_kphalf;
        const CCTK_REAL tau_1_nux_i_j_km1 = tau_1_nux[i_j_km1] + ds_i_j_kmhalf*kappa_1_nux_i_j_kmhalf;


        // Step 8: Select path of least resistance
        const CCTK_REAL new_tau_0_nue_i_j_k = MIN(MIN(MIN(MIN(MIN(tau_0_nue_ip1_j_k,tau_0_nue_im1_j_k),tau_0_nue_i_jp1_k),tau_0_nue_i_jm1_k),tau_0_nue_i_j_kp1),tau_0_nue_i_j_km1);
        const CCTK_REAL new_tau_1_nue_i_j_k = MIN(MIN(MIN(MIN(MIN(tau_1_nue_ip1_j_k,tau_1_nue_im1_j_k),tau_1_nue_i_jp1_k),tau_1_nue_i_jm1_k),tau_1_nue_i_j_kp1),tau_1_nue_i_j_km1);
        const CCTK_REAL new_tau_0_anue_i_j_k = MIN(MIN(MIN(MIN(MIN(tau_0_anue_ip1_j_k,tau_0_anue_im1_j_k),tau_0_anue_i_jp1_k),tau_0_anue_i_jm1_k),tau_0_anue_i_j_kp1),tau_0_anue_i_j_km1);
        const CCTK_REAL new_tau_1_anue_i_j_k = MIN(MIN(MIN(MIN(MIN(tau_1_anue_ip1_j_k,tau_1_anue_im1_j_k),tau_1_anue_i_jp1_k),tau_1_anue_i_jm1_k),tau_1_anue_i_j_kp1),tau_1_anue_i_j_km1);
        const CCTK_REAL new_tau_0_nux_i_j_k = MIN(MIN(MIN(MIN(MIN(tau_0_nux_ip1_j_k,tau_0_nux_im1_j_k),tau_0_nux_i_jp1_k),tau_0_nux_i_jm1_k),tau_0_nux_i_j_kp1),tau_0_nux_i_j_km1);
        const CCTK_REAL new_tau_1_nux_i_j_k = MIN(MIN(MIN(MIN(MIN(tau_1_nux_ip1_j_k,tau_1_nux_im1_j_k),tau_1_nux_i_jp1_k),tau_1_nux_i_jm1_k),tau_1_nux_i_j_kp1),tau_1_nux_i_j_km1);

         // Step 9: Write results to main memory
        tau_0_nue[i_j_k] = new_tau_0_nue_i_j_k;
        tau_1_nue[i_j_k] = new_tau_1_nue_i_j_k;
        tau_0_anue[i_j_k] = new_tau_0_anue_i_j_k;
        tau_1_anue[i_j_k] = new_tau_1_anue_i_j_k;
        tau_0_nux[i_j_k] = new_tau_0_nux_i_j_k;
        tau_1_nux[i_j_k] = new_tau_1_nux_i_j_k;

      } // for(int i=cctk_nghostzones[0];i<cctk_lsh[0]-cctk_nghostzones[0];i++)
    } // for(int j=cctk_nghostzones[1];j<cctk_lsh[1]-cctk_nghostzones[1];j++)
  } // for(int k=cctk_nghostzones[2];k<cctk_lsh[2]-cctk_nghostzones[2];k++)
}
