#ifndef LEAKAGE_HH_
#define LEAKAGE_HH_

namespace Leakage {
  // "Primary" parameters
  constexpr REAL Q_npmass = 1.2935;
  constexpr REAL gamma_0 = 5.565e-2;
  constexpr REAL c_light = 2.997924580000000e+10;
  constexpr REAL N_A = 6.022140760000000e+23;
  constexpr REAL sigma_0 = 1.76e-44;
  constexpr REAL alpha = 1.25;
  constexpr REAL alpha_fs = 7.297352569300000e-03;
  constexpr REAL C_A = 0.5;
  constexpr REAL sinthw2 = 0.23;
  constexpr REAL MeV_to_erg = 1.602176634000000e-06;
  constexpr REAL amu = 1.660539066600000e-24;
  constexpr REAL hc3 = 1.905895198207216e-30;
  constexpr REAL m_e_c2 = 5.109989499961642e-01;
  constexpr REAL Brems_C1 = 2.9988e7;
  constexpr REAL Brems_C2 = 6.5428e7;
  constexpr REAL Brems_zeta = 0.5;
  constexpr REAL eta_nue_0 = 0.0;
  constexpr REAL eta_anue_0 = 0.0;
  constexpr REAL units_geom_to_cgs_density = 6.175828479261933e+17;
  // "Derived" parameters
  constexpr REAL C_V = C_A + 2*sinthw2;
  constexpr REAL beta = N_A*c_light*sigma_0/((m_e_c2)*(m_e_c2));
  constexpr REAL C1pC2_nue_anue = ((-C_A + C_V)*(-C_A + C_V)) + ((C_A + C_V)*(C_A + C_V));
  constexpr REAL C1pC2_nux_anux = ((-C_A + C_V)*(-C_A + C_V)) + ((C_A + C_V - 2)*(C_A + C_V - 2));
}

// Function prototypes
extern "C"
REAL Leakage_Fermi_Dirac_integrals(const int k, const REAL z);

extern "C"
void Leakage_compute_free_rates(const REAL rho_b,
                const REAL Y_e,
                const REAL T,
                const REAL tau_nue,
                const REAL tau_anue,
                REAL *restrict R_free_total_nue ,
                REAL *restrict R_free_total_anue,
                REAL *restrict R_free_total_nux ,
                REAL *restrict Q_free_total_nue ,
                REAL *restrict Q_free_total_anue,
                REAL *restrict Q_free_total_nux);

#endif // LEAKAGE_HH_