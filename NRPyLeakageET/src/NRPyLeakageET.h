#ifndef NRPYLEAKAGE_H_
#define NRPYLEAKAGE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cctk.h"

#ifdef MIN
#undef MIN
#endif

#ifdef MAX
#undef MAX
#endif

#define MAX(a,b) ( (a) > (b) ? (a) : (b) )
#define MIN(a,b) ( (a) < (b) ? (a) : (b) )

// "Primary" parameters
#define NRPyLeakageET_enable_beta_nue (1)
#define NRPyLeakageET_enable_beta_anue (1)
#define NRPyLeakageET_enable_pair_nue_anue (1)
#define NRPyLeakageET_enable_pair_nux_anux (1)
#define NRPyLeakageET_enable_plasmon_nue_anue (1)
#define NRPyLeakageET_enable_plasmon_nux_anux (1)
#define NRPyLeakageET_enable_brems_nui_anui (1)
#define USE_NRPY_CONSTANTS (0)
#define USE_HARM_CONSTANTS (1)
#define USE_ZELMANI_CONSTANTS (2)
#define NRPyLeakageET_Q_npmass (1.2935)
#define NRPyLeakageET_gamma_0 (5.565e-2)
#define NRPyLeakageET_sigma_0 (1.76e-44)
#define NRPyLeakageET_alpha (1.25)
#define NRPyLeakageET_C_A (0.5)
#define NRPyLeakageET_sinthw2 (0.23)
#define NRPyLeakageET_Brems_C1 (2.9988e7)
#define NRPyLeakageET_Brems_C2 (6.5428e7)
#define NRPyLeakageET_Brems_zeta (0.5)
#define NRPyLeakageET_eta_nue_0 (0.0)
#define NRPyLeakageET_eta_anue_0 (0.0)
#define NRPyLeakageET_c_light (2.997924580000000e+10)
#define NRPyLeakageET_N_A (6.022140760000000e+23)
#define NRPyLeakageET_alpha_fs (7.297352569300000e-03)
#define NRPyLeakageET_amu (1.660539066600000e-24)
#define NRPyLeakageET_hc3 (1.905895198207216e-30)
#define NRPyLeakageET_m_e_c2 (5.109989499961642e-01)
#define NRPyLeakageET_units_geom_to_cgs_D (6.175828479261933e+17)
#define NRPyLeakageET_units_cgs_to_geom_D (1.619215953548485e-18)
#define NRPyLeakageET_units_cgs_to_geom_R (7.975433521479384e-24)
#define NRPyLeakageET_units_cgs_to_geom_Q (1.421750164721599e-50)
#define NRPyLeakageET_units_geom_to_cgs_M (1.988409870698051e+33)
#define NRPyLeakageET_units_geom_to_cgs_L (1.476625038050125e+05)
#define NRPyLeakageET_units_geom_to_cgs_T (4.925490947641267e-06)
#define NRPyLeakageET_units_cgs_to_geom_M (5.029144215870041e-34)
#define NRPyLeakageET_units_cgs_to_geom_L (6.772199944005382e-06)
#define NRPyLeakageET_units_cgs_to_geom_T (2.030254467280836e+05)
#define NRPyLeakageET_harm_c_light (2.99792458e10)
#define NRPyLeakageET_harm_N_A (6.0221415e23)
#define NRPyLeakageET_harm_alpha_fs (0.00729735252051)
#define NRPyLeakageET_harm_amu (1.66053886e-24)
#define NRPyLeakageET_harm_hc3 (1.90589514992e-30)
#define NRPyLeakageET_harm_m_e_c2 (5.109989461e-01)
#define NRPyLeakageET_harm_units_geom_to_cgs_D (6.17714470405638e17)
#define NRPyLeakageET_harm_units_cgs_to_geom_D (1.61887093132742e-18)
#define NRPyLeakageET_harm_units_cgs_to_geom_R (7.97315453692e-24)
#define NRPyLeakageET_harm_units_cgs_to_geom_Q (1.42134688491e-50)
#define NRPyLeakageET_harm_units_geom_to_cgs_L (1.476625038050125e+05)
#define NRPyLeakageET_ZL_Q_npmass (1.293333)
#define NRPyLeakageET_ZL_alpha (1.23)
#define NRPyLeakageET_ZL_N_A (6.0221367e+23)
#define NRPyLeakageET_ZL_alpha_fs (7.297352520505561e-03)
#define NRPyLeakageET_ZL_amu (1.674927211e-24)
#define NRPyLeakageET_ZL_hc3 (1.905893979207552e-30)
#define NRPyLeakageET_ZL_m_e_c2 (5.1099891e-01)
// "Derived" parameters
#define NRPyLeakageET_C_V (NRPyLeakageET_C_A + 2*NRPyLeakageET_sinthw2)
#define NRPyLeakageET_beta (NRPyLeakageET_c_light*NRPyLeakageET_sigma_0/((NRPyLeakageET_m_e_c2)*(NRPyLeakageET_m_e_c2)))
#define NRPyLeakageET_harm_beta (NRPyLeakageET_harm_c_light*NRPyLeakageET_sigma_0/((NRPyLeakageET_harm_m_e_c2)*(NRPyLeakageET_harm_m_e_c2)))
#define NRPyLeakageET_C1pC2_nue_anue (((-NRPyLeakageET_C_A + NRPyLeakageET_C_V)*(-NRPyLeakageET_C_A + NRPyLeakageET_C_V)) + ((NRPyLeakageET_C_A + NRPyLeakageET_C_V)*(NRPyLeakageET_C_A + NRPyLeakageET_C_V)))
#define NRPyLeakageET_C1pC2_nux_anux (((-NRPyLeakageET_C_A + NRPyLeakageET_C_V)*(-NRPyLeakageET_C_A + NRPyLeakageET_C_V)) + ((NRPyLeakageET_C_A + NRPyLeakageET_C_V - 2)*(NRPyLeakageET_C_A + NRPyLeakageET_C_V - 2)))


#ifdef __cplusplus
extern "C" {
#endif
// Function prototypes
CCTK_REAL NRPyLeakageET_Fermi_Dirac_integrals(const int k, const CCTK_REAL z);

void NRPyLeakageET_compute_GRMHD_source_terms_and_opacities_nrpy_constants(const CCTK_REAL rho_b,
                                                                           const CCTK_REAL Y_e,
                                                                           const CCTK_REAL T,
                                                                           const CCTK_REAL *restrict tau_nue,
                                                                           const CCTK_REAL *restrict tau_anue,
                                                                           const CCTK_REAL *restrict tau_nux,
                                                                           CCTK_REAL *restrict R_source,
                                                                           CCTK_REAL *restrict Q_source,
                                                                           CCTK_REAL *restrict kappa_nue,
                                                                           CCTK_REAL *restrict kappa_anue,
                                                                           CCTK_REAL *restrict kappa_nux);

void NRPyLeakageET_compute_GRMHD_source_terms_and_opacities_harm_constants(const CCTK_REAL rho_b,
                                                                           const CCTK_REAL Y_e,
                                                                           const CCTK_REAL T,
                                                                           const CCTK_REAL *restrict tau_nue,
                                                                           const CCTK_REAL *restrict tau_anue,
                                                                           const CCTK_REAL *restrict tau_nux,
                                                                           CCTK_REAL *restrict R_source,
                                                                           CCTK_REAL *restrict Q_source,
                                                                           CCTK_REAL *restrict kappa_nue,
                                                                           CCTK_REAL *restrict kappa_anue,
                                                                           CCTK_REAL *restrict kappa_nux);

void NRPyLeakageET_compute_GRMHD_source_terms_and_opacities(const int which_constants_to_use,
                                                            const CCTK_REAL rho_b,
                                                            const CCTK_REAL Y_e,
                                                            const CCTK_REAL T,
                                                            const CCTK_REAL *restrict tau_nue,
                                                            const CCTK_REAL *restrict tau_anue,
                                                            const CCTK_REAL *restrict tau_nux,
                                                            CCTK_REAL *restrict R_source,
                                                            CCTK_REAL *restrict Q_source,
                                                            CCTK_REAL *restrict kappa_nue,
                                                            CCTK_REAL *restrict kappa_anue,
                                                            CCTK_REAL *restrict kappa_nux);

void NRPyLeakageET_compute_opacities_nrpy_constants(const CCTK_REAL rho_b,
                                                    const CCTK_REAL Y_e,
                                                    const CCTK_REAL T,
                                                    const CCTK_REAL *restrict tau_nue,
                                                    const CCTK_REAL *restrict tau_anue,
                                                    const CCTK_REAL *restrict tau_nux,
                                                    CCTK_REAL *restrict kappa_nue,
                                                    CCTK_REAL *restrict kappa_anue,
                                                    CCTK_REAL *restrict kappa_nux);

void NRPyLeakageET_compute_opacities_harm_constants(const CCTK_REAL rho_b,
                                                    const CCTK_REAL Y_e,
                                                    const CCTK_REAL T,
                                                    const CCTK_REAL *restrict tau_nue,
                                                    const CCTK_REAL *restrict tau_anue,
                                                    const CCTK_REAL *restrict tau_nux,
                                                    CCTK_REAL *restrict kappa_nue,
                                                    CCTK_REAL *restrict kappa_anue,
                                                    CCTK_REAL *restrict kappa_nux);

int  NRPyLeakageET_ProcessOwnsData();
void NRPyLeakageET_GetMaxSize(CCTK_ARGUMENTS,int *IterationCounter);
void NRPyLeakageET_optical_depths_initialize_to_zero(CCTK_ARGUMENTS);
void NRPyLeakageET_compute_opacities(CCTK_ARGUMENTS);
void NRPyLeakageET_optical_depths_PathOfLeastResistance(CCTK_ARGUMENTS);
void NRPyLeakageET_copy_opacities_and_optical_depths_to_previous_time_levels(CCTK_ARGUMENTS);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // NRPYLEAKAGE_H_
