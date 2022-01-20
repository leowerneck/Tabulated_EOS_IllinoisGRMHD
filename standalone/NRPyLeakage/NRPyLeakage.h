#ifndef NRPYLEAKAGE_H_
#define NRPYLEAKAGE_H_

#include "NRPyEOS.h"

// "Primary" parameters
#define NRPyLeakage_Q_npmass (1.2935)
#define NRPyLeakage_gamma_0 (5.565e-2)
#define NRPyLeakage_c_light (2.997924580000000e+10)
#define NRPyLeakage_N_A (6.022140760000000e+23)
#define NRPyLeakage_sigma_0 (1.76e-44)
#define NRPyLeakage_alpha (1.25)
#define NRPyLeakage_alpha_fs (7.297352569300000e-03)
#define NRPyLeakage_C_A (0.5)
#define NRPyLeakage_sinthw2 (0.23)
#define NRPyLeakage_MeV_to_erg (1.602176634000000e-06)
#define NRPyLeakage_amu (1.660539066600000e-24)
#define NRPyLeakage_hc3 (1.905895198207216e-30)
#define NRPyLeakage_m_e_c2 (5.109989499961642e-01)
#define NRPyLeakage_Brems_C1 (2.9988e7)
#define NRPyLeakage_Brems_C2 (6.5428e7)
#define NRPyLeakage_Brems_zeta (0.5)
#define NRPyLeakage_eta_nue_0 (0.0)
#define NRPyLeakage_eta_anue_0 (0.0)
#define NRPyLeakage_units_geom_to_cgs_density (6.175828479261933e+17)
#define NRPyLeakage_units_cgs_to_geom_density (1.619215953548485e-18)
#define NRPyLeakage_units_MeV_to_geom (8.396457816868713e-16)
// "Derived" parameters
#define NRPyLeakage_C_V (NRPyLeakage_C_A + 2*NRPyLeakage_sinthw2)

#define NRPyLeakage_beta (NRPyLeakage_c_light*NRPyLeakage_sigma_0/((NRPyLeakage_m_e_c2)*(NRPyLeakage_m_e_c2)))

#define NRPyLeakage_C1pC2_nue_anue (((-NRPyLeakage_C_A + NRPyLeakage_C_V)*(-NRPyLeakage_C_A + NRPyLeakage_C_V)) + ((NRPyLeakage_C_A + NRPyLeakage_C_V)*(NRPyLeakage_C_A + NRPyLeakage_C_V)))

#define NRPyLeakage_C1pC2_nux_anux (((-NRPyLeakage_C_A + NRPyLeakage_C_V)*(-NRPyLeakage_C_A + NRPyLeakage_C_V)) + ((NRPyLeakage_C_A + NRPyLeakage_C_V - 2)*(NRPyLeakage_C_A + NRPyLeakage_C_V - 2)))



// Function prototypes
REAL NRPyLeakage_Fermi_Dirac_integrals(const int k, const REAL z);

void NRPyLeakage_compute_GRMHD_source_terms(const NRPyEOS_params *restrict eos_params,
                                            const REAL rho_b,
                                            const REAL Y_e,
                                            const REAL T,
                                            const REAL tau_nue,
                                            const REAL tau_anue,
                                            const REAL tau_nux,
                                            REAL *restrict R_source,
                                            REAL *restrict Q_source);

#endif // NRPYLEAKAGE_H_