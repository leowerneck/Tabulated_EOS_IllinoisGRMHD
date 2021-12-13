// Compute emission and cooling rates of neutrinos
// following [1] Ruffert et al. (1996)
// https://adsabs.harvard.edu/pdf/1996A%26A...311..532R

#include "Leakage.hh"
#include "WVU_EOS_Tabulated_headers.hh"

void Leakage_emission_and_cooling_rates(const REAL rho,
                                        const REAL Y_e,
                                        const REAL T) {

  // Step 1: Compute needed Fermi-Dirac integrals
  const REAL F_4_nue  = Leakage_Fermi_Dirac_integral(4,eta_nue);
  const REAL F_5_nue  = Leakage_Fermi_Dirac_integral(5,eta_nue);
  const REAL F_4_anue = Leakage_Fermi_Dirac_integral(4,eta_anue);
  const REAL F_5_anue = Leakage_Fermi_Dirac_integral(5,eta_anue);

  // Step 2: Compute all needed blocking factors
  // Step 2.a: Electron capture by a proton; Eq. (A15) of [1]
  const REAL blocking_factor_beta_nue = 1.0/( 1.0 + exp(eta_nue - F_5_nue/F_4_nue) );

  // Step 2.b: Positron capture by a neutron; Eq. (A16) of [1]
  const REAL blocking_factor_beta_anue = 1.0/( 1.0 + exp(eta_anue - F_5_anue/F_4_anue) );

  // Step 2.c: Electron-positron pair annihilation; Eq. (B9) of [1]
  const REAL pair_annihilation_nue_anue_factor = 0.5 * ( F_4_nue/F_3_nue + F_4_anue/F_3_anue );
  const REAL blocking_factor_pair_nue          = 1.0/( 1.0 + exp( eta_nue  - pair_annihilation_nue_anue_factor ) );
  const REAL blocking_factor_pair_anue         = 1.0/( 1.0 + exp( eta_anue - pair_annihilation_nue_anue_factor ) );

  // Step 2.d: Plasmon-decay; Eq. (B13) of [1]
  const REAL gamma                        = LEAKAGE_PLASMON_GAMMA_0 * sqrt( (3.0*eta_e*eta_e + PISQR)/3.0 );
  const REAL plasmon_decay_gamma_factor   = 1.0 + 0.5*gamma*gamma/(1.0+gamma);
  const REAL blocking_factor_plasmon_nue  = 1.0/( eta_nue  - plasmon_decay_gamma_factor );
  const REAL blocking_factor_plasmon_anue = 1.0/( eta_anue - plasmon_decay_gamma_factor );
  const REAL blocking_factor_plasmon_nux  = 1.0/( eta_nux  - plasmon_decay_gamma_factor );

}
