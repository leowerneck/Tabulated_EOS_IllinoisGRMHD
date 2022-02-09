#include "Basic_defines.h"


#ifdef IDX3D
#undef IDX3D
#endif
#define IDX3D(i0,i1,i2) ( (i0) + (Nt0)*( (i1) + (Nt1)*(i2) ) )

/*
 * (c) Leo Werneck
 * Compute GRMHD source terms following Ruffert et al. (1996)
 * https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
void NRPyLeakage_compute_optical_depths(const int Nt0,
                                        const int Nt1,
                                        const int Nt2,
                                        const int Ng0,
                                        const int Ng1,
                                        const int Ng2,
                                        const REAL dxx0,
                                        const REAL dxx1,
                                        const REAL dxx2,
                                        REAL *gammaDD00,
                                        REAL *gammaDD11,
                                        REAL *gammaDD22,
                                        REAL *kappa_0_nue,
                                        REAL *kappa_1_nue,
                                        REAL *kappa_0_anue,
                                        REAL *kappa_1_anue,
                                        REAL *kappa_0_nux,
                                        REAL *kappa_1_nux,
                                        REAL *tau_0_nue,
                                        REAL *tau_1_nue,
                                        REAL *tau_0_anue,
                                        REAL *tau_1_anue,
                                        REAL *tau_0_nux,
                                        REAL *tau_1_nux) {

  printf("Inside NRPyLeakage_compute_optical_depths\n");
  printf("%e %e %e\n",gammaDD00[0],gammaDD11[0],gammaDD22[0]);
  getchar();

  // Step 0: Loop over the grid computing the optical depth
// #pragma omp parallel for
  for(int i2=Ng2;i2<Nt2-Ng2;i2++) {
    for(int i1=Ng1;i1<Nt1-Ng1;i1++) {
      for(int i0=Ng0;i0<Nt0-Ng0;i0++) {

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

        printf("Read metric\n");
        getchar();


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

        // if( i0==Nt0/2 && i1==Nt1/2 && i2==Nt2/2 ) {
        //   printf("dx: %e %e %e\n",dxx0,dxx1,dxx2);
        //   printf("sg: %e %e %e %e %e %e\n",
        //          gammaDD00_i0phalf_i1_i2,
        //          gammaDD00_i0mhalf_i1_i2,
        //          gammaDD11_i0_i1phalf_i2,
        //          gammaDD11_i0_i1mhalf_i2,
        //          gammaDD22_i0_i1_i2phalf,
        //          gammaDD22_i0_i1_i2mhalf);
        //   printf("ds: %e %e %e %e %e %e\n",
        //          ds_i0phalf_i1_i2,ds_i0mhalf_i1_i2,
        //          ds_i0_i1phalf_i2,ds_i0_i1mhalf_i2,
        //          ds_i0_i1_i2phalf,ds_i0_i1_i2mhalf);
        // }

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

        // if( i0==Nt0/2 && i1==Nt1/2 && i2==Nt2/2 ) {
        //   printf("%e %e %e %e\n",tau_0_nue[i0p1_i1_i2],ds_i0phalf_i1_i2,kappa_0_nue_i0phalf_i1_i2,tau_0_nue_i0p1_i1_i2);
        //   printf("%e %e %e %e\n",tau_0_nue[i0m1_i1_i2],ds_i0mhalf_i1_i2,kappa_0_nue_i0mhalf_i1_i2,tau_0_nue_i0m1_i1_i2);
        //   printf("%e %e %e %e\n",tau_0_nue[i0_i1p1_i2],ds_i0_i1phalf_i2,kappa_0_nue_i0_i1phalf_i2,tau_0_nue_i0_i1p1_i2);
        //   printf("%e %e %e %e\n",tau_0_nue[i0_i1m1_i2],ds_i0_i1mhalf_i2,kappa_0_nue_i0_i1mhalf_i2,tau_0_nue_i0_i1m1_i2);
        //   printf("%e %e %e %e\n",tau_0_nue[i0_i1_i2p1],ds_i0_i1_i2phalf,kappa_0_nue_i0_i1_i2phalf,tau_0_nue_i0_i1_i2p1);
        //   printf("%e %e %e %e\n",tau_0_nue[i0_i1_i2m1],ds_i0_i1_i2mhalf,kappa_0_nue_i0_i1_i2mhalf,tau_0_nue_i0_i1_i2m1);
        //   getchar();
        // }

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
