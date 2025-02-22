#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <sys/time.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "IllinoisGRMHD_headers.h"

void Convert_to_HydroBase(CCTK_ARGUMENTS)
{

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // Generally, we only need the HydroBase variables for diagnostic purposes, so we run the below
    // loop only at iterations in which diagnostics are run.
    if(Convert_to_HydroBase_every == 0 || cctk_iteration % Convert_to_HydroBase_every != 0) {
        return;
    }

    /***************
     * PPEOS Patch *
     ***************
     * We will need to set up our EOS in
     * order to be able to compute eps below
     */
    igm_eos_parameters eos;
    initialize_igm_eos_parameters_from_input(igm_eos_key, cctk_time, eos);

#pragma omp parallel for
    for(int k = 0; k < cctk_lsh[2]; k++) {
        for(int j = 0; j < cctk_lsh[1]; j++) {
            for(int i = 0; i < cctk_lsh[0]; i++) {
                int index = CCTK_GFINDEX3D(cctkGH, i, j, k);

                // IllinoisGRMHD defines v^i = u^i/u^0.

                // Meanwhile, the ET/HydroBase formalism, called the Valencia
                // formalism, splits the 4 velocity into a purely spatial part
                // and a part that is normal to the spatial hypersurface:
                // u^a = G (n^a + U^a), (Eq. 14 of arXiv:1304.5544; G=W, U^a=v^a)
                // where n^a is the unit normal vector to the spatial hypersurface,
                // n_a = {-\alpha,0,0,0}, and U^a is the purely spatial part, which
                // is defined in HydroBase as the vel[] vector gridfunction.
                // Then u^a n_a = - \alpha u^0 = G n^a n_a = -G, and
                // of course \alpha u^0 = 1/sqrt(1+γ^ij u_i u_j) = \Gamma,
                // the standard Lorentz factor.

                // Note that n^i = - \beta^i / \alpha, so
                // u^a = \Gamma (n^a + U^a)
                // -> u^i = \Gamma ( U^i - \beta^i / \alpha )
                // which implies
                // v^i = u^i/u^0
                //     = \Gamma/u^0 ( U^i - \beta^i / \alpha ) <- \Gamma = \alpha u^0
                //     = \alpha ( U^i - \beta^i / \alpha )
                //     = \alpha U^i - \beta^i
                const CCTK_REAL vxL        = vx[index];
                const CCTK_REAL vyL        = vy[index];
                const CCTK_REAL vzL        = vz[index];
                const CCTK_REAL lapseL     = alp[index];
                const CCTK_REAL lapseL_inv = 1.0 / lapseL;

                vel[CCTK_GFINDEX4D(cctkGH, i, j, k, 0)] = (vxL + betax[index]) * lapseL_inv;
                vel[CCTK_GFINDEX4D(cctkGH, i, j, k, 1)] = (vyL + betay[index]) * lapseL_inv;
                vel[CCTK_GFINDEX4D(cctkGH, i, j, k, 2)] = (vzL + betaz[index]) * lapseL_inv;

                // \alpha u^0 = 1/sqrt(1+γ^ij u_i u_j) = \Gamma = w_lorentz
                // First compute u^0:
                // Derivation of first equation:
                // \gamma_{ij} (v^i + \beta^i)(v^j + \beta^j)/(\alpha)^2
                //   = \gamma_{ij} 1/(u^0)^2 ( \gamma^{ik} u_k \gamma^{jl} u_l /(\alpha)^2 <- Using
                //   Eq. 53 of arXiv:astro-ph/0503420 = 1/(u^0 \alpha)^2 u_j u_l \gamma^{jl}  <-
                //   Since \gamma_{ij} \gamma^{ik} = \delta^k_j = 1/(u^0 \alpha)^2 ( (u^0 \alpha)^2
                //   - 1 ) <- Using Eq. 56 of arXiv:astro-ph/0503420 = 1 - 1/(u^0 \alpha)^2 <= 1
                CCTK_REAL shiftxL = betax[index];
                CCTK_REAL shiftyL = betay[index];
                CCTK_REAL shiftzL = betaz[index];

                CCTK_REAL gxxL = gxx[index];
                CCTK_REAL gxyL = gxy[index];
                CCTK_REAL gxzL = gxz[index];
                CCTK_REAL gyyL = gyy[index];
                CCTK_REAL gyzL = gyz[index];
                CCTK_REAL gzzL = gzz[index];

                CCTK_REAL one_minus_one_over_alpha_u0_squared
                    = (gxxL * SQR(vxL + shiftxL) + 2.0 * gxyL * (vxL + shiftxL) * (vyL + shiftyL)
                       + 2.0 * gxzL * (vxL + shiftxL) * (vzL + shiftzL) + gyyL * SQR(vyL + shiftyL)
                       + 2.0 * gyzL * (vyL + shiftyL) * (vzL + shiftzL)
                       + gzzL * SQR(vzL + shiftzL))
                    * SQR(lapseL_inv);
                /*** Check for superluminal velocity ***/
                // FIXME: Instead of >1.0, should be one_minus_one_over_alpha_u0_squared >
                // ONE_MINUS_ONE_OVER_GAMMA_SPEED_LIMIT_SQUARED, for consistency with
                // conserv_to_prims routines

                if(one_minus_one_over_alpha_u0_squared > 1.0) {
                    CCTK_VInfo(CCTK_THORNSTRING,
                               "Convert_to_HydroBase WARNING: Found superluminal velocity. This "
                               "should have been caught by IllinoisGRMHD.");
                }

                // A = 1.0-one_minus_one_over_alpha_u0_squared = 1-(1-1/(al u0)^2) = 1/(al u0)^2
                // 1/sqrt(A) = al u0
                CCTK_REAL W = 1.0 / sqrt(1.0 - one_minus_one_over_alpha_u0_squared);
                if(std::isnan(W * lapseL_inv)) {
                    printf("BAD FOUND NAN ALPHAU0 CALC: %.15e %.15e %.15e\n",
                           W,
                           lapseL_inv,
                           one_minus_one_over_alpha_u0_squared);
                }

                w_lorentz[index] = W;

                Bvec[CCTK_GFINDEX4D(cctkGH, i, j, k, 0)] = Bx_center[index];
                Bvec[CCTK_GFINDEX4D(cctkGH, i, j, k, 1)] = By_center[index];
                Bvec[CCTK_GFINDEX4D(cctkGH, i, j, k, 2)] = Bz_center[index];
            }
        }
    }
}
