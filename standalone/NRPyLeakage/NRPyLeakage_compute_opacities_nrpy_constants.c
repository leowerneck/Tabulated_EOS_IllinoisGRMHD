#include "Basic_defines.h"
#include "NRPyEOS.h"
#include "NRPyLeakage.h"

/*
 * (c) Leo Werneck
 * Compute GRMHD source terms following Ruffert et al. (1996)
 * https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
void NRPyLeakage_compute_opacities_nrpy_constants(const NRPyEOS_params *restrict eos_params,
                                                  const REAL rho_b,
                                                  const REAL Y_e,
                                                  const REAL T,
                                                  const REAL *restrict tau_nue,
                                                  const REAL *restrict tau_anue,
                                                  const REAL *restrict tau_nux,
                                                  REAL *restrict kappa_nue,
                                                  REAL *restrict kappa_anue,
                                                  REAL *restrict kappa_nux) {


  // Step 1: Get chemical potentials and mass
  //         fractions using the EOS
  REAL mu_e, mu_p, mu_n, muhat, X_p, X_n;
  NRPyEOS_mue_mup_mun_muhat_Xn_and_Xp_from_rho_Ye_T(eos_params, rho_b, Y_e, T,
                                                    &mu_e, &mu_p, &mu_n, &muhat, &X_p, &X_n);

  // Step 2: Compute rho_b in cgs units
  const REAL rho_b_cgs = rho_b * NRPyLeakage_units_geom_to_cgs_D;

  // Step 3: Compute Y_{pn} and Y_{np}
  const REAL Y_p = Y_e;
  const REAL Y_n = 1-Y_e;
  const REAL exp_metahat = exp(-muhat/T);
  // Step 3.a: Compute Y_{np}
  REAL Y_np = (Y_e < 0.5) ? (2.0*Y_e-1.0)/(exp_metahat-1.0) : Y_n;

  // Step 3.b: Compute Y_{pn}
  REAL Y_pn = (Y_e > 0.5) ? exp_metahat*(2.0*Y_e-1.0)/(exp_metahat-1.0) : Y_p;

  // Step 3.c: Make sure both Y_np and Y_pn are non-negative
  if( Y_np < 0.0 ) Y_np = Y_n;
  if( Y_pn < 0.0 ) Y_pn = Y_p;

  // Step 4: Compute the source terms
  //         Note: The code below is generated by NRPy+
  const REAL tmp_0 = exp(-tau_nue[0]);
  const REAL tmp_1 = (1.0/(T));
  const REAL tmp_2 = mu_e*tmp_1;
  const REAL tmp_4 = NRPyLeakage_eta_nue_0*tmp_0 + (1 - tmp_0)*(-muhat*tmp_1 + tmp_2);
  const REAL tmp_5 = NRPyLeakage_Fermi_Dirac_integrals(4, tmp_4);
  const REAL tmp_6 = NRPyLeakage_N_A*NRPyLeakage_sigma_0*((T)*(T))*rho_b_cgs/((NRPyLeakage_m_e_c2)*(NRPyLeakage_m_e_c2));
  const REAL tmp_7 = tmp_5*tmp_6/NRPyLeakage_Fermi_Dirac_integrals(2, tmp_4);
  const REAL tmp_9 = (5.0/24.0)*((NRPyLeakage_alpha)*(NRPyLeakage_alpha));
  const REAL tmp_10 = (1 - Y_e)*(tmp_9 + 1.0/24.0)/((2.0/3.0)*MAX(mu_n*tmp_1, 0) + 1);
  const REAL tmp_11 = Y_e*(tmp_9 + (1.0/6.0)*((NRPyLeakage_C_V - 1)*(NRPyLeakage_C_V - 1)))/((2.0/3.0)*MAX(mu_p*tmp_1, 0) + 1);
  const REAL tmp_12 = NRPyLeakage_Fermi_Dirac_integrals(5, tmp_4);
  const REAL tmp_13 = (3.0/4.0)*((NRPyLeakage_alpha)*(NRPyLeakage_alpha)) + 1.0/4.0;
  const REAL tmp_14 = Y_np*tmp_13/(exp(-tmp_12/tmp_5 + tmp_2) + 1);
  const REAL tmp_15 = tmp_12/NRPyLeakage_Fermi_Dirac_integrals(3, tmp_4);
  const REAL tmp_17 = tmp_11*tmp_6;
  const REAL tmp_18 = exp(-tau_anue[0]);
  const REAL tmp_20 = NRPyLeakage_eta_anue_0*tmp_18 + (1 - tmp_18)*(muhat*tmp_1 - tmp_2);
  const REAL tmp_21 = NRPyLeakage_Fermi_Dirac_integrals(4, tmp_20);
  const REAL tmp_22 = tmp_21/NRPyLeakage_Fermi_Dirac_integrals(2, tmp_20);
  const REAL tmp_24 = NRPyLeakage_Fermi_Dirac_integrals(5, tmp_20);
  const REAL tmp_25 = Y_pn*tmp_13/(exp(-tmp_2 - tmp_24/tmp_21) + 1);
  const REAL tmp_26 = tmp_24/NRPyLeakage_Fermi_Dirac_integrals(3, tmp_20);
  const REAL tmp_28 = NRPyLeakage_Fermi_Dirac_integrals(4, 0)/NRPyLeakage_Fermi_Dirac_integrals(2, 0);
  const REAL tmp_30 = NRPyLeakage_Fermi_Dirac_integrals(5, 0)/NRPyLeakage_Fermi_Dirac_integrals(3, 0);
  kappa_nue[0] = NRPyLeakage_units_geom_to_cgs_L*(EnsureFinite(tmp_10*tmp_7) + EnsureFinite(tmp_11*tmp_7) + EnsureFinite(tmp_14*tmp_7));
  kappa_nue[1] = NRPyLeakage_units_geom_to_cgs_L*(EnsureFinite(tmp_15*tmp_17) + EnsureFinite(tmp_10*tmp_15*tmp_6) + EnsureFinite(tmp_14*tmp_15*tmp_6));
  kappa_anue[0] = NRPyLeakage_units_geom_to_cgs_L*(EnsureFinite(tmp_17*tmp_22) + EnsureFinite(tmp_10*tmp_22*tmp_6) + EnsureFinite(tmp_22*tmp_25*tmp_6));
  kappa_anue[1] = NRPyLeakage_units_geom_to_cgs_L*(EnsureFinite(tmp_17*tmp_26) + EnsureFinite(tmp_10*tmp_26*tmp_6) + EnsureFinite(tmp_25*tmp_26*tmp_6));
  kappa_nux[0] = NRPyLeakage_units_geom_to_cgs_L*(EnsureFinite(tmp_17*tmp_28) + EnsureFinite(tmp_10*tmp_28*tmp_6));
  kappa_nux[1] = NRPyLeakage_units_geom_to_cgs_L*(EnsureFinite(tmp_17*tmp_30) + EnsureFinite(tmp_10*tmp_30*tmp_6));

  // Step 5: Make sure results are finite; if not reset to small value
  for(int i=0;i<2;i++) {
    if( !isfinite(kappa_nue [i]) ) kappa_nue [i] = 1e-15;
    if( !isfinite(kappa_anue[i]) ) kappa_anue[i] = 1e-15;
    if( !isfinite(kappa_nux [i]) ) kappa_nux [i] = 1e-15;
  }
}
