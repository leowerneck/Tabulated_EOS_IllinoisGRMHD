#include "cctk.h"
#include "NRPyLeakageET.h"

/*
 * (c) Leo Werneck
 * Compute GRMHD source terms following Ruffert et al. (1996)
 * https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
void NRPyLeakageET_compute_opacities_nrpy_constants(const CCTK_REAL rho_b,
                                                    const CCTK_REAL Y_e,
                                                    const CCTK_REAL T,
                                                    const CCTK_REAL *restrict tau_nue,
                                                    const CCTK_REAL *restrict tau_anue,
                                                    const CCTK_REAL *restrict tau_nux,
                                                    CCTK_REAL *restrict kappa_nue,
                                                    CCTK_REAL *restrict kappa_anue,
                                                    CCTK_REAL *restrict kappa_nux) {


  // Step 1: Get chemical potentials and mass
  //         fractions using the EOS
  CCTK_REAL mu_e, mu_p, mu_n, muhat, X_p, X_n;
  WVU_EOS_mue_mup_mun_muhat_Xn_and_Xp_from_rho_Ye_T(rho_b, Y_e, T, &mu_e, &mu_p, &mu_n, &muhat, &X_p, &X_n);

  // Step 2: Compute rho_b in cgs units
  const CCTK_REAL rho_b_cgs = rho_b * NRPyLeakageET_units_geom_to_cgs_D;

  // Step 3: Compute Y_{pn} and Y_{np}
  const CCTK_REAL Y_p = Y_e;
  const CCTK_REAL Y_n = 1-Y_e;
  const CCTK_REAL exp_metahat = exp(-muhat/T);
  // Step 3.a: Compute Y_{np}
  CCTK_REAL Y_np = (Y_e < 0.5) ? (2.0*Y_e-1.0)/(exp_metahat-1.0) : Y_n;

  // Step 3.b: Compute Y_{pn}
  CCTK_REAL Y_pn = (Y_e > 0.5) ? exp_metahat*(2.0*Y_e-1.0)/(exp_metahat-1.0) : Y_p;

  // Step 3.c: Make sure both Y_np and Y_pn are non-negative
  if( Y_np < 0.0 ) Y_np = Y_n;
  if( Y_pn < 0.0 ) Y_pn = Y_p;

  // Leo says: this is the way ZelmaniLeak computes Y_pn and Y_np
  // Step 3: Compute Y_{pn} and Y_{np}
  // Step 3.a: Compute Y_{pn} (See discussion below Eq. A8 in https://arxiv.org/pdf/1306.4953.pdf)
  // CCTK_REAL Y_pn;
  // if( rho_b_cgs > 2e12 ) {
  //   // Step 3.a.i: Use Eqs. (A13) and (A14) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  //   Y_pn = (2*Y_e - 1)/(1-exp(muhat/T));
  // }
  // else {
  //   Y_pn = 1-Y_e;
  // }
  // // Step 3.a.ii: Make sure Y_{pn} is nonzero
  // Y_pn = MAX(Y_pn,0.0);

  // // Step 3.b: Compute Y_{np} (Eq. A13 in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R)
  // const CCTK_REAL Y_np = exp(muhat/T) * Y_pn;

  // Step 4: Compute the source terms
  //         Note: The code below is generated by NRPy+
  const CCTK_REAL tmp_0 = (1.0/(T));
  const CCTK_REAL tmp_1 = mu_e*tmp_0;
  const CCTK_REAL tmp_3 = -muhat*tmp_0 + tmp_1;
  const CCTK_REAL tmp_4 = NRPyLeakageET_Fermi_Dirac_integrals(4, tmp_3);
  const CCTK_REAL tmp_5 = NRPyLeakageET_N_A*NRPyLeakageET_sigma_0*((T)*(T))*rho_b_cgs/((NRPyLeakageET_m_e_c2)*(NRPyLeakageET_m_e_c2));
  const CCTK_REAL tmp_6 = tmp_4*tmp_5/NRPyLeakageET_Fermi_Dirac_integrals(2, tmp_3);
  const CCTK_REAL tmp_8 = (5.0/24.0)*((NRPyLeakageET_alpha)*(NRPyLeakageET_alpha));
  const CCTK_REAL tmp_9 = (1 - Y_e)*(tmp_8 + 1.0/24.0)/((2.0/3.0)*MAX(mu_n*tmp_0, 0) + 1);
  const CCTK_REAL tmp_10 = Y_e*(tmp_8 + (1.0/6.0)*((NRPyLeakageET_C_V - 1)*(NRPyLeakageET_C_V - 1)))/((2.0/3.0)*MAX(mu_p*tmp_0, 0) + 1);
  const CCTK_REAL tmp_11 = NRPyLeakageET_Fermi_Dirac_integrals(5, tmp_3);
  const CCTK_REAL tmp_12 = (3.0/4.0)*((NRPyLeakageET_alpha)*(NRPyLeakageET_alpha)) + 1.0/4.0;
  const CCTK_REAL tmp_13 = Y_np*tmp_12/(exp(tmp_1 - tmp_11/tmp_4) + 1);
  const CCTK_REAL tmp_14 = tmp_11/NRPyLeakageET_Fermi_Dirac_integrals(3, tmp_3);
  const CCTK_REAL tmp_16 = tmp_10*tmp_5;
  const CCTK_REAL tmp_18 = muhat*tmp_0 - tmp_1;
  const CCTK_REAL tmp_19 = NRPyLeakageET_Fermi_Dirac_integrals(4, tmp_18);
  const CCTK_REAL tmp_20 = tmp_19/NRPyLeakageET_Fermi_Dirac_integrals(2, tmp_18);
  const CCTK_REAL tmp_22 = NRPyLeakageET_Fermi_Dirac_integrals(5, tmp_18);
  const CCTK_REAL tmp_23 = Y_pn*tmp_12/(exp(-tmp_1 - tmp_22/tmp_19) + 1);
  const CCTK_REAL tmp_24 = tmp_22/NRPyLeakageET_Fermi_Dirac_integrals(3, tmp_18);
  const CCTK_REAL tmp_26 = NRPyLeakageET_Fermi_Dirac_integrals(4, 0)/NRPyLeakageET_Fermi_Dirac_integrals(2, 0);
  const CCTK_REAL tmp_28 = NRPyLeakageET_Fermi_Dirac_integrals(5, 0)/NRPyLeakageET_Fermi_Dirac_integrals(3, 0);
  kappa_nue[0] = NRPyLeakageET_units_geom_to_cgs_L*(tmp_10*tmp_6 + tmp_13*tmp_6 + tmp_6*tmp_9);
  kappa_nue[1] = NRPyLeakageET_units_geom_to_cgs_L*(tmp_13*tmp_14*tmp_5 + tmp_14*tmp_16 + tmp_14*tmp_5*tmp_9);
  kappa_anue[0] = NRPyLeakageET_units_geom_to_cgs_L*(tmp_16*tmp_20 + tmp_20*tmp_23*tmp_5 + tmp_20*tmp_5*tmp_9);
  kappa_anue[1] = NRPyLeakageET_units_geom_to_cgs_L*(tmp_16*tmp_24 + tmp_23*tmp_24*tmp_5 + tmp_24*tmp_5*tmp_9);
  kappa_nux[0] = NRPyLeakageET_units_geom_to_cgs_L*(tmp_16*tmp_26 + tmp_26*tmp_5*tmp_9);
  kappa_nux[1] = NRPyLeakageET_units_geom_to_cgs_L*(tmp_16*tmp_28 + tmp_28*tmp_5*tmp_9);

  // Step 5: Make sure results are finite; if not reset to small value
  for(int i=0;i<2;i++) {
    if( !isfinite(kappa_nue [i]) ) kappa_nue [i] = 1e-15;
    if( !isfinite(kappa_anue[i]) ) kappa_anue[i] = 1e-15;
    if( !isfinite(kappa_nux [i]) ) kappa_nux [i] = 1e-15;
  }
}
