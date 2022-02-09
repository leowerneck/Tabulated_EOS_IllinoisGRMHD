#ifndef NRPYLEAKAGE_H_
#define NRPYLEAKAGE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cctk.h"

#define REAL CCTK_REAL

#ifdef MIN
#undef MIN
#endif

#ifdef MAX
#undef MAX
#endif

#define MAX(a,b) ( (a) > (b) ? (a) : (b) )
#define MIN(a,b) ( (a) < (b) ? (a) : (b) )

// "Primary" parameters
#define NRPyLeakage_enable_beta_nue (1)
#define NRPyLeakage_enable_beta_anue (1)
#define NRPyLeakage_enable_pair_nue_anue (1)
#define NRPyLeakage_enable_pair_nux_anux (1)
#define NRPyLeakage_enable_plasmon_nue_anue (1)
#define NRPyLeakage_enable_plasmon_nux_anux (1)
#define NRPyLeakage_enable_brems_nui_anui (1)
#define USE_NRPY_CONSTANTS (0)
#define USE_HARM_CONSTANTS (1)
#define USE_ZELMANI_CONSTANTS (2)
#define NRPyLeakage_Q_npmass (1.2935)
#define NRPyLeakage_gamma_0 (5.565e-2)
#define NRPyLeakage_sigma_0 (1.76e-44)
#define NRPyLeakage_alpha (1.25)
#define NRPyLeakage_C_A (0.5)
#define NRPyLeakage_sinthw2 (0.23)
#define NRPyLeakage_Brems_C1 (2.9988e7)
#define NRPyLeakage_Brems_C2 (6.5428e7)
#define NRPyLeakage_Brems_zeta (0.5)
#define NRPyLeakage_eta_nue_0 (0.0)
#define NRPyLeakage_eta_anue_0 (0.0)
#define NRPyLeakage_c_light (2.997924580000000e+10)
#define NRPyLeakage_N_A (6.022140760000000e+23)
#define NRPyLeakage_alpha_fs (7.297352569300000e-03)
#define NRPyLeakage_amu (1.660539066600000e-24)
#define NRPyLeakage_hc3 (7.838450084421221e-48)
#define NRPyLeakage_m_e_c2 (5.109989499961642e-01)
#define NRPyLeakage_units_geom_to_cgs_D (6.175828479261933e+17)
#define NRPyLeakage_units_cgs_to_geom_D (1.619215953548485e-18)
#define NRPyLeakage_units_cgs_to_geom_R (7.975433521479384e-24)
#define NRPyLeakage_units_cgs_to_geom_Q (1.421750164721599e-50)
#define NRPyLeakage_units_geom_to_cgs_M (1.988409870698051e+33)
#define NRPyLeakage_units_geom_to_cgs_L (1.476625038050125e+05)
#define NRPyLeakage_units_geom_to_cgs_T (4.925490947641267e-06)
#define NRPyLeakage_units_cgs_to_geom_M (5.029144215870041e-34)
#define NRPyLeakage_units_cgs_to_geom_L (6.772199944005382e-06)
#define NRPyLeakage_units_cgs_to_geom_T (2.030254467280836e+05)
#define NRPyLeakage_harm_c_light (2.99792458e10)
#define NRPyLeakage_harm_N_A (6.0221415e23)
#define NRPyLeakage_harm_alpha_fs (0.00729735252051)
#define NRPyLeakage_harm_amu (1.66053886e-24)
#define NRPyLeakage_harm_hc3 (1.90589514992e-30)
#define NRPyLeakage_harm_m_e_c2 (5.109989461e-01)
#define NRPyLeakage_harm_units_geom_to_cgs_D (6.17714470405638e17)
#define NRPyLeakage_harm_units_cgs_to_geom_D (1.61887093132742e-18)
#define NRPyLeakage_harm_units_cgs_to_geom_R (7.97315453692e-24)
#define NRPyLeakage_harm_units_cgs_to_geom_Q (1.42134688491e-50)
#define NRPyLeakage_harm_units_geom_to_cgs_L (1.476625038050125e+05)
#define NRPyLeakage_ZL_Q_npmass (1.293333)
#define NRPyLeakage_ZL_alpha (1.23)
#define NRPyLeakage_ZL_N_A (6.0221367e+23)
#define NRPyLeakage_ZL_alpha_fs (7.297352520505561e-03)
#define NRPyLeakage_ZL_amu (1.674927211e-24)
#define NRPyLeakage_ZL_hc3 (1.905893979207552e-30)
#define NRPyLeakage_ZL_m_e_c2 (5.1099891e-01)
// "Derived" parameters
#define NRPyLeakage_C_V (NRPyLeakage_C_A + 2*NRPyLeakage_sinthw2)
#define NRPyLeakage_beta (NRPyLeakage_c_light*NRPyLeakage_sigma_0/((NRPyLeakage_m_e_c2)*(NRPyLeakage_m_e_c2)))
#define NRPyLeakage_harm_beta (NRPyLeakage_harm_c_light*NRPyLeakage_sigma_0/((NRPyLeakage_harm_m_e_c2)*(NRPyLeakage_harm_m_e_c2)))
#define NRPyLeakage_C1pC2_nue_anue (((-NRPyLeakage_C_A + NRPyLeakage_C_V)*(-NRPyLeakage_C_A + NRPyLeakage_C_V)) + ((NRPyLeakage_C_A + NRPyLeakage_C_V)*(NRPyLeakage_C_A + NRPyLeakage_C_V)))
#define NRPyLeakage_C1pC2_nux_anux (((-NRPyLeakage_C_A + NRPyLeakage_C_V)*(-NRPyLeakage_C_A + NRPyLeakage_C_V)) + ((NRPyLeakage_C_A + NRPyLeakage_C_V - 2)*(NRPyLeakage_C_A + NRPyLeakage_C_V - 2)))


// Function prototypes
REAL NRPyLeakage_Fermi_Dirac_integrals(const int k, const REAL z);

void NRPyLeakage_compute_GRMHD_source_terms_and_opacities_nrpy_constants(const REAL rho_b,
                                                                         const REAL Y_e,
                                                                         const REAL T,
                                                                         const REAL *restrict tau_nue,
                                                                         const REAL *restrict tau_anue,
                                                                         const REAL *restrict tau_nux,
                                                                         REAL *restrict R_source,
                                                                         REAL *restrict Q_source,
                                                                         REAL *restrict kappa_nue,
                                                                         REAL *restrict kappa_anue,
                                                                         REAL *restrict kappa_nux);

void NRPyLeakage_compute_GRMHD_source_terms_and_opacities_harm_constants(const REAL rho_b,
                                                                         const REAL Y_e,
                                                                         const REAL T,
                                                                         const REAL *restrict tau_nue,
                                                                         const REAL *restrict tau_anue,
                                                                         const REAL *restrict tau_nux,
                                                                         REAL *restrict R_source,
                                                                         REAL *restrict Q_source,
                                                                         REAL *restrict kappa_nue,
                                                                         REAL *restrict kappa_anue,
                                                                         REAL *restrict kappa_nux);

void NRPyLeakage_compute_GRMHD_source_terms_and_opacities(const int which_constants_to_use,
                                                          const REAL rho_b,
                                                          const REAL Y_e,
                                                          const REAL T,
                                                          const REAL *restrict tau_nue,
                                                          const REAL *restrict tau_anue,
                                                          const REAL *restrict tau_nux,
                                                          REAL *restrict R_source,
                                                          REAL *restrict Q_source,
                                                          REAL *restrict kappa_nue,
                                                          REAL *restrict kappa_anue,
                                                          REAL *restrict kappa_nux);

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
                                        REAL *restrict tau_1_nux);

void NRPyLeakage_compute_opacities_and_add_source_terms_to_MHD_rhss_impl( const CCTK_POINTER_TO_CONST cctkGH,
                                                                          const int *restrict cctk_lsh,
                                                                          const int *restrict cctk_nghostzones,
                                                                          const CCTK_REAL W_max,
                                                                          const CCTK_REAL *restrict alp,
                                                                          const CCTK_REAL *restrict betax,
                                                                          const CCTK_REAL *restrict betay,
                                                                          const CCTK_REAL *restrict betaz,
                                                                          const CCTK_REAL *restrict gxx,
                                                                          const CCTK_REAL *restrict gxy,
                                                                          const CCTK_REAL *restrict gxz,
                                                                          const CCTK_REAL *restrict gyy,
                                                                          const CCTK_REAL *restrict gyz,
                                                                          const CCTK_REAL *restrict gzz,
                                                                          const CCTK_REAL *restrict rho,
                                                                          const CCTK_REAL *restrict Y_e,
                                                                          const CCTK_REAL *restrict temperature,
                                                                          const CCTK_REAL *restrict vel,
                                                                          CCTK_REAL *restrict Ye_star_rhs,
                                                                          CCTK_REAL *restrict tau_rhs,
                                                                          CCTK_REAL *restrict st_x_rhs,
                                                                          CCTK_REAL *restrict st_y_rhs,
                                                                          CCTK_REAL *restrict st_z_rhs );

#endif // NRPYLEAKAGE_H_
