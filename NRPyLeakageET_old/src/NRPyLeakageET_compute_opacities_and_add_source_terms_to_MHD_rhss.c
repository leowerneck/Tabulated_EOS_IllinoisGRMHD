#include "cctk.h"
#include "cctk_Parameters.h"
#include "NRPyLeakageET.h"

#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

#define CHECK_POINTER(pointer,name)                                     \
  if( !pointer ) CCTK_VError(__LINE__,__FILE__,CCTK_THORNSTRING,"Failed to get pointer for gridfunction '%s'",name);

void NRPyLeakageET_compute_opacities_and_add_source_terms_to_MHD_rhss(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(verbosity_level>1) CCTK_VINFO("Inside NRPyLeakageET_compute_opacities_and_add_source_terms_to_MHD_rhss");

  const int timelevel = 0;

  // Step 1: Get pointers to opacity and optical depth gridfunctions
  CCTK_REAL *Ye_star_rhs = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,timelevel,GFstring_Ye_star_rhs));
  CCTK_REAL *tau_rhs     = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,timelevel,GFstring_tau_rhs));
  CCTK_REAL *st_x_rhs    = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,timelevel,GFstring_st_x_rhs));
  CCTK_REAL *st_y_rhs    = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,timelevel,GFstring_st_y_rhs));
  CCTK_REAL *st_z_rhs    = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,timelevel,GFstring_st_z_rhs));

  // Step 2: Check pointers are ok
  CHECK_POINTER(Ye_star_rhs,GFstring_Ye_star_rhs);
  CHECK_POINTER(tau_rhs    ,GFstring_tau_rhs    );
  CHECK_POINTER(st_x_rhs   ,GFstring_st_x_rhs   );
  CHECK_POINTER(st_y_rhs   ,GFstring_st_y_rhs   );
  CHECK_POINTER(st_z_rhs   ,GFstring_st_z_rhs   );

  // Step 3: Ghostzones begin and end index
  const int imin = cctk_nghostzones[0];
  const int imax = cctk_lsh[0] - cctk_nghostzones[0];
  const int jmin = cctk_nghostzones[1];
  const int jmax = cctk_lsh[1] - cctk_nghostzones[1];
  const int kmin = cctk_nghostzones[2];
  const int kmax = cctk_lsh[2] - cctk_nghostzones[2];


  // Step 3: Compute opacities and leakage source terms
  CCTK_REAL R_avg=0,Q_avg=0;
#pragma omp parallel for reduction(+:R_avg,Q_avg)
  for(int k=kmin;k<kmax;k++) {
    for(int j=jmin;j<jmax;j++) {
      for(int i=imin;i<imax;i++) {

        // Step 3.a: Set the index of the current gridpoint
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // Step 3.b: Check if we are within the threshold
        CCTK_REAL rhoL = rho[index];

        if( rhoL < rho_threshold ) {
          // Step 3.b.i: Below density threshold.
          //             Set opacities to zero; don't add anything to the RHSs
          kappa_0_nue [index] = 0.0;
          kappa_1_nue [index] = 0.0;
          kappa_0_anue[index] = 0.0;
          kappa_1_anue[index] = 0.0;
          kappa_0_nux [index] = 0.0;
          kappa_1_nux [index] = 0.0;
        }
        else {
          // Step 3.c: Read from main memory
          const CCTK_REAL alpL         = alp[index];
          const CCTK_REAL alpinvsqrdL  = 1.0/(alpL*alpL);
          const CCTK_REAL betaxL       = betax[index];
          const CCTK_REAL betayL       = betay[index];
          const CCTK_REAL betazL       = betaz[index];
          CCTK_REAL gxxL               = gxx[index];
          CCTK_REAL gxyL               = gxy[index];
          CCTK_REAL gxzL               = gxz[index];
          CCTK_REAL gyyL               = gyy[index];
          CCTK_REAL gyzL               = gyz[index];
          CCTK_REAL gzzL               = gzz[index];
          CCTK_REAL vxL                = alpL*velx[index] - betaxL;
          CCTK_REAL vyL                = alpL*vely[index] - betayL;
          CCTK_REAL vzL                = alpL*velz[index] - betazL;
          const CCTK_REAL Y_eL         = Y_e[index];
          const CCTK_REAL temperatureL = temperature[index];
          const CCTK_REAL tau_nueL [2] = {tau_0_nue_p [index],tau_1_nue_p [index]};
          const CCTK_REAL tau_anueL[2] = {tau_0_anue_p[index],tau_1_anue_p[index]};
          const CCTK_REAL tau_nuxL [2] = {tau_0_nux_p [index],tau_1_nux_p [index]};

          // Step 3.d: Compute BSSN quantities; enforce det(gammabar_{ij}) = 1
          // Step 3.d.i: Compute the determinant of the physical metric
          CCTK_REAL gdet = gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + gxzL * gxyL * gyzL
            - gxzL * gyyL * gxzL - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL;
          gdet = fabs(gdet);
          // Step 3.d.ii: Compute BSSN quantities (gb = "gbar" = conformal metric)
          const CCTK_REAL phiL   = (1.0/12.0) * log(gdet);
          const CCTK_REAL psiL   = exp(phiL);
          const CCTK_REAL psi2L  = psiL *psiL;
          const CCTK_REAL psi4L  = psi2L*psi2L;
          const CCTK_REAL psi6L  = psi4L*psi2L;
          const CCTK_REAL psim4L = 1.0/(psi4L);
          CCTK_REAL gbxxL        = gxxL*psim4L;
          CCTK_REAL gbxyL        = gxyL*psim4L;
          CCTK_REAL gbxzL        = gxzL*psim4L;
          CCTK_REAL gbyyL        = gyyL*psim4L;
          CCTK_REAL gbyzL        = gyzL*psim4L;
          CCTK_REAL gbzzL        = gzzL*psim4L;
          // Step 3.d.iii: Compute the determinant of the conformal metric
          CCTK_REAL gbdet = gbxxL * gbyyL * gbzzL + gbxyL * gbyzL * gbxzL + gbxzL * gbxyL * gbyzL
            - gbxzL * gbyyL * gbxzL - gbxyL * gbxyL * gbzzL - gbxxL * gbyzL * gbyzL;
          gbdet = fabs(gbdet);
          // Step 3.d.iv: Enforce det(gammabar) = 1 constraint
          CCTK_REAL gbdet_Fm1o3 = fabs(1.0/cbrt(gbdet));
          gbxxL *= gbdet_Fm1o3;
          gbxxL *= gbdet_Fm1o3;
          gbxxL *= gbdet_Fm1o3;
          gbxxL *= gbdet_Fm1o3;
          gbxxL *= gbdet_Fm1o3;
          gbxxL *= gbdet_Fm1o3;
          // Step 3.d.v: Recompute physical metric
          gxxL = gbxxL*psi4L;
          gxyL = gbxyL*psi4L;
          gxzL = gbxzL*psi4L;
          gyyL = gbyyL*psi4L;
          gyzL = gbyzL*psi4L;
          gzzL = gbzzL*psi4L;

          // Step 3.e: Compute u^{mu}
          // Step 3.e.i: Compute Lorentz factor W
          const CCTK_REAL one_minus_one_over_Wsqrd = (    gxxL*(vxL + betaxL)*(vxL + betaxL) +
                                                      2.0*gxyL*(vxL + betaxL)*(vyL + betayL) +
                                                      2.0*gxzL*(vxL + betaxL)*(vzL + betazL) +
                                                          gyyL*(vyL + betayL)*(vyL + betayL) +
                                                      2.0*gyzL*(vyL + betayL)*(vzL + betazL) +
                                                          gzzL*(vzL + betazL)*(vzL + betazL) )*alpinvsqrdL;
          CCTK_REAL W = 1.0/sqrt( 1 - one_minus_one_over_Wsqrd );

          // Step 3.e.ii: Impose a speed limit to the velocities
          if( W > W_max ) {
            const CCTK_REAL W_scale = W_max/W;
            vxL = (vxL + betaxL)*W_scale - betaxL;
            vyL = (vyL + betayL)*W_scale - betayL;
            vzL = (vzL + betazL)*W_scale - betazL;
            W   = W_max;
          }

          // Step 3.f: Compute u^{mu} using:
          //  - W = alpha u^{0}     => u^{0} = W / alpha
          //  - v^{i} = u^{i}/u^{i} => u^{i} = v^{i}u^{0}
          const CCTK_REAL u0L = W / alpL;
          const CCTK_REAL uxL = vxL * u0L;
          const CCTK_REAL uyL = vyL * u0L;
          const CCTK_REAL uzL = vzL * u0L;

          // Step 3.g: Compute R, Q, and the neutrino opacities
          CCTK_REAL R_sourceL, Q_sourceL, kappa_nueL[2],kappa_anueL[2],kappa_nuxL[2];
          NRPyLeakageET_compute_GRMHD_source_terms_and_opacities(constants_key,
                                                                 rhoL,Y_eL,temperatureL,
                                                                 tau_nueL,tau_anueL,tau_nuxL,
                                                                 &R_sourceL,&Q_sourceL,kappa_nueL,kappa_anueL,kappa_nuxL);

          // Step 3.h: Compute MHD right-hand sides
          const CCTK_REAL sqrtmgL      = alpL * psi6L;
          const CCTK_REAL Ye_star_rhsL = sqrtmgL * R_sourceL;
          const CCTK_REAL tau_rhsL     = sqrtmgL * Q_sourceL * u0L;
          const CCTK_REAL st_x_rhsL    = sqrtmgL * Q_sourceL * uxL;
          const CCTK_REAL st_y_rhsL    = sqrtmgL * Q_sourceL * uyL;
          const CCTK_REAL st_z_rhsL    = sqrtmgL * Q_sourceL * uzL;

          // Step 3.i: Write to main memory
          kappa_0_nue [index]  = kappa_nueL [0];
          kappa_1_nue [index]  = kappa_nueL [1];
          kappa_0_anue[index]  = kappa_anueL[0];
          kappa_1_anue[index]  = kappa_anueL[1];
          kappa_0_nux [index]  = kappa_nuxL [0];
          kappa_1_nux [index]  = kappa_nuxL [1];

          // Step 3.j: Update right-hand sides only in the grid interior
          Ye_star_rhs [index] += Ye_star_rhsL;
          tau_rhs     [index] += tau_rhsL;
          st_x_rhs    [index] += st_x_rhsL;
          st_y_rhs    [index] += st_y_rhsL;
          st_z_rhs    [index] += st_z_rhsL;

          R_avg += R_sourceL/rhoL;
          Q_avg += Q_sourceL/rhoL;
        }
      }
    }
  }
  if(verbosity_level>0) {
    CCTK_VInfo(CCTK_THORNSTRING,"***** Iter. # %d, Lev: %d, Averages -- R/rho: %e | Q/rho: %e *****",cctk_iteration,GetRefinementLevel(cctkGH),R_avg,Q_avg);
    if(verbosity_level>1) CCTK_VINFO("Finished NRPyLeakageET_compute_opacities_and_add_source_terms_to_MHD_rhss");
  }
}
