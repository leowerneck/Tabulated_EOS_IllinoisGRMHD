#include <cmath>
#include <cstdio>
#include <cstdlib>
#ifndef ENABLE_STANDALONE_IGM_C2P_SOLVER
#include "cctk.h"
#include "cctk_Parameters.h"
#endif

void IllinoisGRMHD_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij(
    const cGH *restrict cctkGH,
    const int *restrict cctk_lsh,
    CCTK_REAL *restrict gxx, CCTK_REAL *gxy,
    CCTK_REAL *restrict gxz,
    CCTK_REAL *restrict gyy,
    CCTK_REAL *restrict gyz,
    CCTK_REAL *restrict gzz,
    const CCTK_REAL *restrict alp,
    CCTK_REAL *restrict gtxx,
    CCTK_REAL *restrict gtxy,
    CCTK_REAL *restrict gtxz,
    CCTK_REAL *restrict gtyy,
    CCTK_REAL *restrict gtyz,
    CCTK_REAL *restrict gtzz,
    CCTK_REAL *restrict gtupxx,
    CCTK_REAL *restrict gtupxy,
    CCTK_REAL *restrict gtupxz,
    CCTK_REAL *restrict gtupyy,
    CCTK_REAL *restrict gtupyz,
    CCTK_REAL *restrict gtupzz,
    CCTK_REAL *restrict phi,
    CCTK_REAL *restrict psi,
    CCTK_REAL *restrict lapm1 ) {

#ifndef ENABLE_STANDALONE_IGM_C2P_SOLVER
 DECLARE_CCTK_PARAMETERS;
#endif


#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
        int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
        CCTK_REAL gxx_physL=gxx[index];
        CCTK_REAL gxy_physL=gxy[index];
        CCTK_REAL gxz_physL=gxz[index];
        CCTK_REAL gyy_physL=gyy[index];
        CCTK_REAL gyz_physL=gyz[index];
        CCTK_REAL gzz_physL=gzz[index];

        /**********************************************************************
         * Compute \tilde{\gamma_{ij}}, phi, and psi (BSSN) from g_{ij} (ADM) *
         **********************************************************************/
        CCTK_REAL gijdet = gxx_physL * gyy_physL * gzz_physL + gxy_physL * gyz_physL * gxz_physL + gxz_physL * gxy_physL * gyz_physL
          - gxz_physL * gyy_physL * gxz_physL - gxy_physL * gxy_physL * gzz_physL - gxx_physL * gyz_physL * gyz_physL;

        gijdet = fabs(gijdet);


        CCTK_REAL phiL = (1.0/12.0) * log(gijdet);
        CCTK_REAL psiL = exp(phiL);

        CCTK_REAL Psim4 = 1.0/(psiL*psiL*psiL*psiL);
        CCTK_REAL gtxxL = gxx_physL*Psim4;
        CCTK_REAL gtxyL = gxy_physL*Psim4;
        CCTK_REAL gtxzL = gxz_physL*Psim4;
        CCTK_REAL gtyyL = gyy_physL*Psim4;
        CCTK_REAL gtyzL = gyz_physL*Psim4;
        CCTK_REAL gtzzL = gzz_physL*Psim4;


        /*********************************
         * Apply det gtij = 1 constraint *
         *********************************/
        CCTK_REAL gtijdet = gtxxL * gtyyL * gtzzL + gtxyL * gtyzL * gtxzL + gtxzL * gtxyL * gtyzL -
          gtxzL * gtyyL * gtxzL - gtxyL * gtxyL * gtzzL - gtxxL * gtyzL * gtyzL;

        CCTK_REAL gtijdet_Fm1o3 = fabs(1.0/cbrt(gtijdet));

        gtxxL = gtxxL * gtijdet_Fm1o3;
        gtxyL = gtxyL * gtijdet_Fm1o3;
        gtxzL = gtxzL * gtijdet_Fm1o3;
        gtyyL = gtyyL * gtijdet_Fm1o3;
        gtyzL = gtyzL * gtijdet_Fm1o3;
        gtzzL = gtzzL * gtijdet_Fm1o3;

        if(gtijdet<0.0) {
#ifndef ENABLE_STANDALONE_IGM_C2P_SOLVER
          CCTK_VWarn(CCTK_WARN_ALERT,__LINE__, __FILE__, CCTK_THORNSTRING,
#else
          fprintf(stderr,
#endif
"WARNING: det[3-metric]<0.0 at point  %d %d %d | cctk_lsh: %d %d %d. Hopefully this is occurring in gz's! gtij_phys = %.2e %.2e %.2e %.2e %.2e %.2e gtij_new = %.2e %.2e %.2e %.2e %.2e %.2e | gijdet = %.2e | gtijdet = %.2e\n",
i,j,k,cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],gxx_physL,gxy_physL,gxz_physL,gyy_physL,gyz_physL,gzz_physL,gtxxL,gtxyL,gtxzL,gtyyL,gtyzL,gtzzL,-gijdet,gtijdet);
}


        CCTK_REAL Psi4 = psiL*psiL*psiL*psiL;
        /*****************************************
         * Set all the needed BSSN gridfunctions *
         *****************************************/
        phi[index] = phiL;
        psi[index] = psiL;

        lapm1[index] = alp[index] - 1.0;

        gtxx[index] = gtxxL;
        gtxy[index] = gtxyL;
        gtxz[index] = gtxzL;
        gtyy[index] = gtyyL;
        gtyz[index] = gtyzL;
        gtzz[index] = gtzzL;

        gxx[index] = gtxxL*Psi4;
        gxy[index] = gtxyL*Psi4;
        gxz[index] = gtxzL*Psi4;
        gyy[index] = gtyyL*Psi4;
        gyz[index] = gtyzL*Psi4;
        gzz[index] = gtzzL*Psi4;

        gtupxx[index] =   ( gtyyL * gtzzL - gtyzL * gtyzL );
        gtupxy[index] = - ( gtxyL * gtzzL - gtyzL * gtxzL );
        gtupxz[index] =   ( gtxyL * gtyzL - gtyyL * gtxzL );
        gtupyy[index] =   ( gtxxL * gtzzL - gtxzL * gtxzL );
        gtupyz[index] = - ( gtxxL * gtyzL - gtxyL * gtxzL );
        gtupzz[index] =   ( gtxxL * gtyyL - gtxyL * gtxyL );
      }
}
