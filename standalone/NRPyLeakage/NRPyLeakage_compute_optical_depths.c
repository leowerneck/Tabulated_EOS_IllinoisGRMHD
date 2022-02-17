#include "Basic_defines.h"


#ifndef IDX3D
#define IDX3D(i0,i1,i2) ( (i0) + (N0)*( (i1) + (N1)*(i2) ) )
#endif

/*
 * (c) Leo Werneck
 * Compute GRMHD source terms following Ruffert et al. (1996)
 * https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
void NRPyLeakage_compute_optical_depths(const int N0,
                                        const int N1,
                                        const int N2,
                                        const int Ng0,
                                        const int Ng1,
                                        const int Ng2,
                                        const int dxx0,
                                        const int dxx1,
                                        const int dxx2,
                                        const REAL *restrict gammaDD00,
                                        const REAL *restrict gammaDD11,
                                        const REAL *restrict gammaDD22,
                                        const REAL *restrict kappa_0_nue,
                                        const REAL *restrict kappa_1_nue,
                                        const REAL *restrict kappa_0_anue,
                                        const REAL *restrict kappa_1_anue,
                                        const REAL *restrict kappa_0_nux,
                                        const REAL *restrict kappa_1_nux,
                                        REAL *restrict tau_0_nue,
                                        REAL *restrict tau_1_nue,
                                        REAL *restrict tau_0_anue,
                                        REAL *restrict tau_1_anue,
                                        REAL *restrict tau_0_nux,
                                        REAL *restrict tau_1_nux) {


  // Step 0: Loop over the grid computing the optical depth
#pragma omp parallel for
  for(int i2=Ng2;i2<N2-Ng2;i2++) {
    for(int i1=Ng1;i1<N1-Ng1;i1++) {
      for(int i0=Ng0;i0<N0-Ng0;i0++) {

        // Step 1: Set gridpoint indices
        const int i0_i1_i2   = IDX3D(i0  ,i1,i2  );
        const int i0p1_i1_i2 = IDX3D(i0+1,i1,i2  );
        const int i0m1_i1_i2 = IDX3D(i0-1,i1,i2  );
        const int i0_i1p1_i2 = IDX3D(i0,i1+1,i2  );
        const int i0_i1m1_i2 = IDX3D(i0,i1-1,i2  );
        const int i0_i1_i2p1 = IDX3D(i0,i1  ,i2+1);
        const int i0_i1_i2m1 = IDX3D(i0,i1  ,i2-1);

        // Step 2: Read in metric gfs from main memory
        const REAL gammaDD00_i0_i1_i2 = gammaDD00[i0_i1_i2];
        const REAL gammaDD00_i0p1_i1_i2 = gammaDD00[i0p1_i1_i2];
        const REAL gammaDD00_i0m1_i1_i2 = gammaDD00[i0m1_i1_i2];

        const REAL gammaDD11_i0_i1_i2 = gammaDD11[i0_i1_i2];
        const REAL gammaDD11_i0_i1p1_i2 = gammaDD11[i0_i1p1_i2];
        const REAL gammaDD11_i0_i1m1_i2 = gammaDD11[i0_i1m1_i2];

        const REAL gammaDD22_i0_i1_i2 = gammaDD22[i0_i1_i2];
        const REAL gammaDD22_i0_i1_i2p1 = gammaDD22[i0_i1_i2p1];
        const REAL gammaDD22_i0_i1_i2m1 = gammaDD22[i0_i1_i2m1];


        // Step 3: Read in opacity gfs from main memory
        const REAL kappa_0_nue_i0_i1_i2 = kappa_0_nue[i0_i1_i2];
        const REAL kappa_0_nue_i0p1_i1_i2 = kappa_0_nue[i0p1_i1_i2];
        const REAL kappa_0_nue_i0m1_i1_i2 = kappa_0_nue[i0m1_i1_i2];
        const REAL kappa_0_nue_i0_i1p1_i2 = kappa_0_nue[i0_i1p1_i2];
        const REAL kappa_0_nue_i0_i1m1_i2 = kappa_0_nue[i0_i1m1_i2];
        const REAL kappa_0_nue_i0_i1_i2p1 = kappa_0_nue[i0_i1_i2p1];
        const REAL kappa_0_nue_i0_i1_i2m1 = kappa_0_nue[i0_i1_i2m1];

        const REAL kappa_1_nue_i0_i1_i2 = kappa_1_nue[i0_i1_i2];
        const REAL kappa_1_nue_i0p1_i1_i2 = kappa_1_nue[i0p1_i1_i2];
        const REAL kappa_1_nue_i0m1_i1_i2 = kappa_1_nue[i0m1_i1_i2];
        const REAL kappa_1_nue_i0_i1p1_i2 = kappa_1_nue[i0_i1p1_i2];
        const REAL kappa_1_nue_i0_i1m1_i2 = kappa_1_nue[i0_i1m1_i2];
        const REAL kappa_1_nue_i0_i1_i2p1 = kappa_1_nue[i0_i1_i2p1];
        const REAL kappa_1_nue_i0_i1_i2m1 = kappa_1_nue[i0_i1_i2m1];

        const REAL kappa_0_anue_i0_i1_i2 = kappa_0_anue[i0_i1_i2];
        const REAL kappa_0_anue_i0p1_i1_i2 = kappa_0_anue[i0p1_i1_i2];
        const REAL kappa_0_anue_i0m1_i1_i2 = kappa_0_anue[i0m1_i1_i2];
        const REAL kappa_0_anue_i0_i1p1_i2 = kappa_0_anue[i0_i1p1_i2];
        const REAL kappa_0_anue_i0_i1m1_i2 = kappa_0_anue[i0_i1m1_i2];
        const REAL kappa_0_anue_i0_i1_i2p1 = kappa_0_anue[i0_i1_i2p1];
        const REAL kappa_0_anue_i0_i1_i2m1 = kappa_0_anue[i0_i1_i2m1];

        const REAL kappa_1_anue_i0_i1_i2 = kappa_1_anue[i0_i1_i2];
        const REAL kappa_1_anue_i0p1_i1_i2 = kappa_1_anue[i0p1_i1_i2];
        const REAL kappa_1_anue_i0m1_i1_i2 = kappa_1_anue[i0m1_i1_i2];
        const REAL kappa_1_anue_i0_i1p1_i2 = kappa_1_anue[i0_i1p1_i2];
        const REAL kappa_1_anue_i0_i1m1_i2 = kappa_1_anue[i0_i1m1_i2];
        const REAL kappa_1_anue_i0_i1_i2p1 = kappa_1_anue[i0_i1_i2p1];
        const REAL kappa_1_anue_i0_i1_i2m1 = kappa_1_anue[i0_i1_i2m1];

        const REAL kappa_0_nux_i0_i1_i2 = kappa_0_nux[i0_i1_i2];
        const REAL kappa_0_nux_i0p1_i1_i2 = kappa_0_nux[i0p1_i1_i2];
        const REAL kappa_0_nux_i0m1_i1_i2 = kappa_0_nux[i0m1_i1_i2];
        const REAL kappa_0_nux_i0_i1p1_i2 = kappa_0_nux[i0_i1p1_i2];
        const REAL kappa_0_nux_i0_i1m1_i2 = kappa_0_nux[i0_i1m1_i2];
        const REAL kappa_0_nux_i0_i1_i2p1 = kappa_0_nux[i0_i1_i2p1];
        const REAL kappa_0_nux_i0_i1_i2m1 = kappa_0_nux[i0_i1_i2m1];

        const REAL kappa_1_nux_i0_i1_i2 = kappa_1_nux[i0_i1_i2];
        const REAL kappa_1_nux_i0p1_i1_i2 = kappa_1_nux[i0p1_i1_i2];
        const REAL kappa_1_nux_i0m1_i1_i2 = kappa_1_nux[i0m1_i1_i2];
        const REAL kappa_1_nux_i0_i1p1_i2 = kappa_1_nux[i0_i1p1_i2];
        const REAL kappa_1_nux_i0_i1m1_i2 = kappa_1_nux[i0_i1m1_i2];
        const REAL kappa_1_nux_i0_i1_i2p1 = kappa_1_nux[i0_i1_i2p1];
        const REAL kappa_1_nux_i0_i1_i2m1 = kappa_1_nux[i0_i1_i2m1];


        // Step 4: Compute metric at cell faces
        const REAL gammaDD00_i0phalf_i1_i2 = 0.5*(gammaDD00_i0_i1_i2 + gammaDD00_i0p1_i1_i2);
        const REAL gammaDD00_i0mhalf_i1_i2 = 0.5*(gammaDD00_i0_i1_i2 + gammaDD00_i0m1_i1_i2);

        const REAL gammaDD11_i0_i1phalf_i2 = 0.5*(gammaDD11_i0_i1_i2 + gammaDD11_i0_i1p1_i2);
        const REAL gammaDD11_i0_i1mhalf_i2 = 0.5*(gammaDD11_i0_i1_i2 + gammaDD11_i0_i1m1_i2);

        const REAL gammaDD22_i0_i1_i2phalf = 0.5*(gammaDD22_i0_i1_i2 + gammaDD22_i0_i1_i2p1);
        const REAL gammaDD22_i0_i1_i2mhalf = 0.5*(gammaDD22_i0_i1_i2 + gammaDD22_i0_i1_i2m1);


        // Step 5: Compute ds^{i} = sqrt(gamma_{ii}dx^{i}dx^{i})
        const REAL ds_i0phalf_i1_i2 = sqrt(dxx0*dxx0*gammaDD00_i0phalf_i1_i2);
        const REAL ds_i0mhalf_i1_i2 = sqrt(dxx0*dxx0*gammaDD00_i0mhalf_i1_i2);
        const REAL ds_i0_i1phalf_i2 = sqrt(dxx1*dxx1*gammaDD11_i0_i1phalf_i2);
        const REAL ds_i0_i1mhalf_i2 = sqrt(dxx1*dxx1*gammaDD11_i0_i1mhalf_i2);
        const REAL ds_i0_i1_i2phalf = sqrt(dxx2*dxx2*gammaDD22_i0_i1_i2phalf);
        const REAL ds_i0_i1_i2mhalf = sqrt(dxx2*dxx2*gammaDD22_i0_i1_i2mhalf);

        // Step 6: Compute opacities at cell faces
        const REAL kappa_0_nue_i0phalf_i1_i2 = 0.5*(kappa_0_nue_i0_i1_i2 + kappa_0_nue_i0p1_i1_i2);
        const REAL kappa_0_nue_i0mhalf_i1_i2 = 0.5*(kappa_0_nue_i0_i1_i2 + kappa_0_nue_i0m1_i1_i2);
        const REAL kappa_0_nue_i0_i1phalf_i2 = 0.5*(kappa_0_nue_i0_i1_i2 + kappa_0_nue_i0_i1p1_i2);
        const REAL kappa_0_nue_i0_i1mhalf_i2 = 0.5*(kappa_0_nue_i0_i1_i2 + kappa_0_nue_i0_i1m1_i2);
        const REAL kappa_0_nue_i0_i1_i2phalf = 0.5*(kappa_0_nue_i0_i1_i2 + kappa_0_nue_i0_i1_i2p1);
        const REAL kappa_0_nue_i0_i1_i2mhalf = 0.5*(kappa_0_nue_i0_i1_i2 + kappa_0_nue_i0_i1_i2m1);

        const REAL kappa_1_nue_i0phalf_i1_i2 = 0.5*(kappa_1_nue_i0_i1_i2 + kappa_1_nue_i0p1_i1_i2);
        const REAL kappa_1_nue_i0mhalf_i1_i2 = 0.5*(kappa_1_nue_i0_i1_i2 + kappa_1_nue_i0m1_i1_i2);
        const REAL kappa_1_nue_i0_i1phalf_i2 = 0.5*(kappa_1_nue_i0_i1_i2 + kappa_1_nue_i0_i1p1_i2);
        const REAL kappa_1_nue_i0_i1mhalf_i2 = 0.5*(kappa_1_nue_i0_i1_i2 + kappa_1_nue_i0_i1m1_i2);
        const REAL kappa_1_nue_i0_i1_i2phalf = 0.5*(kappa_1_nue_i0_i1_i2 + kappa_1_nue_i0_i1_i2p1);
        const REAL kappa_1_nue_i0_i1_i2mhalf = 0.5*(kappa_1_nue_i0_i1_i2 + kappa_1_nue_i0_i1_i2m1);

        const REAL kappa_0_anue_i0phalf_i1_i2 = 0.5*(kappa_0_anue_i0_i1_i2 + kappa_0_anue_i0p1_i1_i2);
        const REAL kappa_0_anue_i0mhalf_i1_i2 = 0.5*(kappa_0_anue_i0_i1_i2 + kappa_0_anue_i0m1_i1_i2);
        const REAL kappa_0_anue_i0_i1phalf_i2 = 0.5*(kappa_0_anue_i0_i1_i2 + kappa_0_anue_i0_i1p1_i2);
        const REAL kappa_0_anue_i0_i1mhalf_i2 = 0.5*(kappa_0_anue_i0_i1_i2 + kappa_0_anue_i0_i1m1_i2);
        const REAL kappa_0_anue_i0_i1_i2phalf = 0.5*(kappa_0_anue_i0_i1_i2 + kappa_0_anue_i0_i1_i2p1);
        const REAL kappa_0_anue_i0_i1_i2mhalf = 0.5*(kappa_0_anue_i0_i1_i2 + kappa_0_anue_i0_i1_i2m1);

        const REAL kappa_1_anue_i0phalf_i1_i2 = 0.5*(kappa_1_anue_i0_i1_i2 + kappa_1_anue_i0p1_i1_i2);
        const REAL kappa_1_anue_i0mhalf_i1_i2 = 0.5*(kappa_1_anue_i0_i1_i2 + kappa_1_anue_i0m1_i1_i2);
        const REAL kappa_1_anue_i0_i1phalf_i2 = 0.5*(kappa_1_anue_i0_i1_i2 + kappa_1_anue_i0_i1p1_i2);
        const REAL kappa_1_anue_i0_i1mhalf_i2 = 0.5*(kappa_1_anue_i0_i1_i2 + kappa_1_anue_i0_i1m1_i2);
        const REAL kappa_1_anue_i0_i1_i2phalf = 0.5*(kappa_1_anue_i0_i1_i2 + kappa_1_anue_i0_i1_i2p1);
        const REAL kappa_1_anue_i0_i1_i2mhalf = 0.5*(kappa_1_anue_i0_i1_i2 + kappa_1_anue_i0_i1_i2m1);

        const REAL kappa_0_nux_i0phalf_i1_i2 = 0.5*(kappa_0_nux_i0_i1_i2 + kappa_0_nux_i0p1_i1_i2);
        const REAL kappa_0_nux_i0mhalf_i1_i2 = 0.5*(kappa_0_nux_i0_i1_i2 + kappa_0_nux_i0m1_i1_i2);
        const REAL kappa_0_nux_i0_i1phalf_i2 = 0.5*(kappa_0_nux_i0_i1_i2 + kappa_0_nux_i0_i1p1_i2);
        const REAL kappa_0_nux_i0_i1mhalf_i2 = 0.5*(kappa_0_nux_i0_i1_i2 + kappa_0_nux_i0_i1m1_i2);
        const REAL kappa_0_nux_i0_i1_i2phalf = 0.5*(kappa_0_nux_i0_i1_i2 + kappa_0_nux_i0_i1_i2p1);
        const REAL kappa_0_nux_i0_i1_i2mhalf = 0.5*(kappa_0_nux_i0_i1_i2 + kappa_0_nux_i0_i1_i2m1);

        const REAL kappa_1_nux_i0phalf_i1_i2 = 0.5*(kappa_1_nux_i0_i1_i2 + kappa_1_nux_i0p1_i1_i2);
        const REAL kappa_1_nux_i0mhalf_i1_i2 = 0.5*(kappa_1_nux_i0_i1_i2 + kappa_1_nux_i0m1_i1_i2);
        const REAL kappa_1_nux_i0_i1phalf_i2 = 0.5*(kappa_1_nux_i0_i1_i2 + kappa_1_nux_i0_i1p1_i2);
        const REAL kappa_1_nux_i0_i1mhalf_i2 = 0.5*(kappa_1_nux_i0_i1_i2 + kappa_1_nux_i0_i1m1_i2);
        const REAL kappa_1_nux_i0_i1_i2phalf = 0.5*(kappa_1_nux_i0_i1_i2 + kappa_1_nux_i0_i1_i2p1);
        const REAL kappa_1_nux_i0_i1_i2mhalf = 0.5*(kappa_1_nux_i0_i1_i2 + kappa_1_nux_i0_i1_i2m1);


        // Step 7: Compute optical depth at neighboring points
        const REAL tau_0_nue_i0p1_i1_i2 = tau_0_nue[i0p1_i1_i2] + ds_i0phalf_i1_i2*kappa_0_nue_i0phalf_i1_i2;
        const REAL tau_0_nue_i0m1_i1_i2 = tau_0_nue[i0m1_i1_i2] + ds_i0mhalf_i1_i2*kappa_0_nue_i0mhalf_i1_i2;
        const REAL tau_0_nue_i0_i1p1_i2 = tau_0_nue[i0_i1p1_i2] + ds_i0_i1phalf_i2*kappa_0_nue_i0_i1phalf_i2;
        const REAL tau_0_nue_i0_i1m1_i2 = tau_0_nue[i0_i1m1_i2] + ds_i0_i1mhalf_i2*kappa_0_nue_i0_i1mhalf_i2;
        const REAL tau_0_nue_i0_i1_i2p1 = tau_0_nue[i0_i1_i2p1] + ds_i0_i1_i2phalf*kappa_0_nue_i0_i1_i2phalf;
        const REAL tau_0_nue_i0_i1_i2m1 = tau_0_nue[i0_i1_i2m1] + ds_i0_i1_i2mhalf*kappa_0_nue_i0_i1_i2mhalf;

        const REAL tau_1_nue_i0p1_i1_i2 = tau_1_nue[i0p1_i1_i2] + ds_i0phalf_i1_i2*kappa_1_nue_i0phalf_i1_i2;
        const REAL tau_1_nue_i0m1_i1_i2 = tau_1_nue[i0m1_i1_i2] + ds_i0mhalf_i1_i2*kappa_1_nue_i0mhalf_i1_i2;
        const REAL tau_1_nue_i0_i1p1_i2 = tau_1_nue[i0_i1p1_i2] + ds_i0_i1phalf_i2*kappa_1_nue_i0_i1phalf_i2;
        const REAL tau_1_nue_i0_i1m1_i2 = tau_1_nue[i0_i1m1_i2] + ds_i0_i1mhalf_i2*kappa_1_nue_i0_i1mhalf_i2;
        const REAL tau_1_nue_i0_i1_i2p1 = tau_1_nue[i0_i1_i2p1] + ds_i0_i1_i2phalf*kappa_1_nue_i0_i1_i2phalf;
        const REAL tau_1_nue_i0_i1_i2m1 = tau_1_nue[i0_i1_i2m1] + ds_i0_i1_i2mhalf*kappa_1_nue_i0_i1_i2mhalf;

        const REAL tau_0_anue_i0p1_i1_i2 = tau_0_anue[i0p1_i1_i2] + ds_i0phalf_i1_i2*kappa_0_anue_i0phalf_i1_i2;
        const REAL tau_0_anue_i0m1_i1_i2 = tau_0_anue[i0m1_i1_i2] + ds_i0mhalf_i1_i2*kappa_0_anue_i0mhalf_i1_i2;
        const REAL tau_0_anue_i0_i1p1_i2 = tau_0_anue[i0_i1p1_i2] + ds_i0_i1phalf_i2*kappa_0_anue_i0_i1phalf_i2;
        const REAL tau_0_anue_i0_i1m1_i2 = tau_0_anue[i0_i1m1_i2] + ds_i0_i1mhalf_i2*kappa_0_anue_i0_i1mhalf_i2;
        const REAL tau_0_anue_i0_i1_i2p1 = tau_0_anue[i0_i1_i2p1] + ds_i0_i1_i2phalf*kappa_0_anue_i0_i1_i2phalf;
        const REAL tau_0_anue_i0_i1_i2m1 = tau_0_anue[i0_i1_i2m1] + ds_i0_i1_i2mhalf*kappa_0_anue_i0_i1_i2mhalf;

        const REAL tau_1_anue_i0p1_i1_i2 = tau_1_anue[i0p1_i1_i2] + ds_i0phalf_i1_i2*kappa_1_anue_i0phalf_i1_i2;
        const REAL tau_1_anue_i0m1_i1_i2 = tau_1_anue[i0m1_i1_i2] + ds_i0mhalf_i1_i2*kappa_1_anue_i0mhalf_i1_i2;
        const REAL tau_1_anue_i0_i1p1_i2 = tau_1_anue[i0_i1p1_i2] + ds_i0_i1phalf_i2*kappa_1_anue_i0_i1phalf_i2;
        const REAL tau_1_anue_i0_i1m1_i2 = tau_1_anue[i0_i1m1_i2] + ds_i0_i1mhalf_i2*kappa_1_anue_i0_i1mhalf_i2;
        const REAL tau_1_anue_i0_i1_i2p1 = tau_1_anue[i0_i1_i2p1] + ds_i0_i1_i2phalf*kappa_1_anue_i0_i1_i2phalf;
        const REAL tau_1_anue_i0_i1_i2m1 = tau_1_anue[i0_i1_i2m1] + ds_i0_i1_i2mhalf*kappa_1_anue_i0_i1_i2mhalf;

        const REAL tau_0_nux_i0p1_i1_i2 = tau_0_nux[i0p1_i1_i2] + ds_i0phalf_i1_i2*kappa_0_nux_i0phalf_i1_i2;
        const REAL tau_0_nux_i0m1_i1_i2 = tau_0_nux[i0m1_i1_i2] + ds_i0mhalf_i1_i2*kappa_0_nux_i0mhalf_i1_i2;
        const REAL tau_0_nux_i0_i1p1_i2 = tau_0_nux[i0_i1p1_i2] + ds_i0_i1phalf_i2*kappa_0_nux_i0_i1phalf_i2;
        const REAL tau_0_nux_i0_i1m1_i2 = tau_0_nux[i0_i1m1_i2] + ds_i0_i1mhalf_i2*kappa_0_nux_i0_i1mhalf_i2;
        const REAL tau_0_nux_i0_i1_i2p1 = tau_0_nux[i0_i1_i2p1] + ds_i0_i1_i2phalf*kappa_0_nux_i0_i1_i2phalf;
        const REAL tau_0_nux_i0_i1_i2m1 = tau_0_nux[i0_i1_i2m1] + ds_i0_i1_i2mhalf*kappa_0_nux_i0_i1_i2mhalf;

        const REAL tau_1_nux_i0p1_i1_i2 = tau_1_nux[i0p1_i1_i2] + ds_i0phalf_i1_i2*kappa_1_nux_i0phalf_i1_i2;
        const REAL tau_1_nux_i0m1_i1_i2 = tau_1_nux[i0m1_i1_i2] + ds_i0mhalf_i1_i2*kappa_1_nux_i0mhalf_i1_i2;
        const REAL tau_1_nux_i0_i1p1_i2 = tau_1_nux[i0_i1p1_i2] + ds_i0_i1phalf_i2*kappa_1_nux_i0_i1phalf_i2;
        const REAL tau_1_nux_i0_i1m1_i2 = tau_1_nux[i0_i1m1_i2] + ds_i0_i1mhalf_i2*kappa_1_nux_i0_i1mhalf_i2;
        const REAL tau_1_nux_i0_i1_i2p1 = tau_1_nux[i0_i1_i2p1] + ds_i0_i1_i2phalf*kappa_1_nux_i0_i1_i2phalf;
        const REAL tau_1_nux_i0_i1_i2m1 = tau_1_nux[i0_i1_i2m1] + ds_i0_i1_i2mhalf*kappa_1_nux_i0_i1_i2mhalf;


        // Step 8: Select path of least resistance
        const REAL new_tau_0_nue_i0_i1_i2 = MIN(MIN(MIN(MIN(MIN(tau_0_nue_i0p1_i1_i2,tau_0_nue_i0m1_i1_i2),tau_0_nue_i0_i1p1_i2),tau_0_nue_i0_i1m1_i2),tau_0_nue_i0_i1_i2p1),tau_0_nue_i0_i1_i2m1);
        const REAL new_tau_1_nue_i0_i1_i2 = MIN(MIN(MIN(MIN(MIN(tau_1_nue_i0p1_i1_i2,tau_1_nue_i0m1_i1_i2),tau_1_nue_i0_i1p1_i2),tau_1_nue_i0_i1m1_i2),tau_1_nue_i0_i1_i2p1),tau_1_nue_i0_i1_i2m1);
        const REAL new_tau_0_anue_i0_i1_i2 = MIN(MIN(MIN(MIN(MIN(tau_0_anue_i0p1_i1_i2,tau_0_anue_i0m1_i1_i2),tau_0_anue_i0_i1p1_i2),tau_0_anue_i0_i1m1_i2),tau_0_anue_i0_i1_i2p1),tau_0_anue_i0_i1_i2m1);
        const REAL new_tau_1_anue_i0_i1_i2 = MIN(MIN(MIN(MIN(MIN(tau_1_anue_i0p1_i1_i2,tau_1_anue_i0m1_i1_i2),tau_1_anue_i0_i1p1_i2),tau_1_anue_i0_i1m1_i2),tau_1_anue_i0_i1_i2p1),tau_1_anue_i0_i1_i2m1);
        const REAL new_tau_0_nux_i0_i1_i2 = MIN(MIN(MIN(MIN(MIN(tau_0_nux_i0p1_i1_i2,tau_0_nux_i0m1_i1_i2),tau_0_nux_i0_i1p1_i2),tau_0_nux_i0_i1m1_i2),tau_0_nux_i0_i1_i2p1),tau_0_nux_i0_i1_i2m1);
        const REAL new_tau_1_nux_i0_i1_i2 = MIN(MIN(MIN(MIN(MIN(tau_1_nux_i0p1_i1_i2,tau_1_nux_i0m1_i1_i2),tau_1_nux_i0_i1p1_i2),tau_1_nux_i0_i1m1_i2),tau_1_nux_i0_i1_i2p1),tau_1_nux_i0_i1_i2m1);

         // Step 9: Write results to main memory
        tau_0_nue[i0_i1_i2] = new_tau_0_nue_i0_i1_i2;
        tau_1_nue[i0_i1_i2] = new_tau_1_nue_i0_i1_i2;
        tau_0_anue[i0_i1_i2] = new_tau_0_anue_i0_i1_i2;
        tau_1_anue[i0_i1_i2] = new_tau_1_anue_i0_i1_i2;
        tau_0_nux[i0_i1_i2] = new_tau_0_nux_i0_i1_i2;
        tau_1_nux[i0_i1_i2] = new_tau_1_nux_i0_i1_i2;

      } // for(int i0=Ng0;i0<N0-Ng0;i0++)
    } // for(int i1=Ng1;i1<N1-Ng1;i1++)
  } // for(int i2=Ng2;i2<N2-Ng2;i2++)
}
