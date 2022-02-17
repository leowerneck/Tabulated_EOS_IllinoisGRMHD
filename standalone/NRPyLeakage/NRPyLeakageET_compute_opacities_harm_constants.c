#include "cctk.h"
#include "NRPyLeakageET.h"

/*
 * (c) Leo Werneck
 * Compute GRMHD source terms following Ruffert et al. (1996)
 * https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
void NRPyLeakageET_compute_opacities_harm_constants(const CCTK_REAL rho_b,
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
  const CCTK_REAL rho_b_cgs = rho_b * NRPyLeakage_harm_units_geom_to_cgs_D;

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

  // Step 4: Compute the source terms
  //         Note: The code below is generated by NRPy+

  CCTK_REAL R_beta_nue,R_beta_anue,R_pair_nue_anue,R_pair_nux_anux,R_plasmon_nue_anue,R_plasmon_nux_anux,R_Brems_nui_anui;
  CCTK_REAL Q_beta_nue,Q_beta_anue,Q_pair_nue_anue,Q_pair_nux_anux,Q_plasmon_nue_anue,Q_plasmon_nux_anux,Q_Brems_nui_anui;
  CCTK_REAL R_eff_nue,R_eff_anue,Q_eff_nue,Q_eff_anue,Q_eff_nux;
  CCTK_REAL kappa_s_n_0_nue,kappa_s_p_0_nue,kappa_a_0_nue,kappa_s_n_1_nue,kappa_s_p_1_nue,kappa_a_1_nue;
  CCTK_REAL kappa_s_n_0_anue,kappa_s_p_0_anue,kappa_a_0_anue,kappa_s_n_1_anue,kappa_s_p_1_anue,kappa_a_1_anue;
  CCTK_REAL kappa_s_n_0_nux,kappa_s_p_0_nux,kappa_s_n_1_nux,kappa_s_p_1_nux;
  const CCTK_REAL tmp_0 = (1.0/(T));
  const CCTK_REAL tmp_1 = mu_e*tmp_0;
  const CCTK_REAL tmp_3 = -muhat*tmp_0 + tmp_1;
  const CCTK_REAL tmp_4 = NRPyLeakage_Fermi_Dirac_integrals(4, tmp_3);
  const CCTK_REAL tmp_5 = NRPyLeakage_harm_N_A*NRPyLeakage_sigma_0*((T)*(T))*rho_b_cgs/((NRPyLeakage_harm_m_e_c2)*(NRPyLeakage_harm_m_e_c2));
  const CCTK_REAL tmp_6 = tmp_4*tmp_5/NRPyLeakage_Fermi_Dirac_integrals(2, tmp_3);
  const CCTK_REAL tmp_8 = (5.0/24.0)*((NRPyLeakage_alpha)*(NRPyLeakage_alpha));
  const CCTK_REAL tmp_9 = (1 - Y_e)*(tmp_8 + 1.0/24.0)/((2.0/3.0)*MAX(mu_n*tmp_0, 0) + 1);
  const CCTK_REAL tmp_10 = EnsureFinite(tmp_6*tmp_9);
  const CCTK_REAL tmp_11 = Y_e*(tmp_8 + (1.0/6.0)*((NRPyLeakage_C_V - 1)*(NRPyLeakage_C_V - 1)))/((2.0/3.0)*MAX(mu_p*tmp_0, 0) + 1);
  const CCTK_REAL tmp_12 = EnsureFinite(tmp_11*tmp_6);
  const CCTK_REAL tmp_13 = NRPyLeakage_Fermi_Dirac_integrals(5, tmp_3);
  const CCTK_REAL tmp_14 = (3.0/4.0)*((NRPyLeakage_alpha)*(NRPyLeakage_alpha)) + 1.0/4.0;
  const CCTK_REAL tmp_15 = Y_np*tmp_14/(exp(tmp_1 - tmp_13/tmp_4) + 1);
  const CCTK_REAL tmp_16 = EnsureFinite(tmp_15*tmp_6);
  const CCTK_REAL tmp_17 = tmp_5*tmp_9;
  const CCTK_REAL tmp_18 = tmp_13/NRPyLeakage_Fermi_Dirac_integrals(3, tmp_3);
  const CCTK_REAL tmp_19 = EnsureFinite(tmp_17*tmp_18);
  const CCTK_REAL tmp_21 = EnsureFinite(tmp_11*tmp_18*tmp_5);
  const CCTK_REAL tmp_22 = EnsureFinite(tmp_15*tmp_18*tmp_5);
  const CCTK_REAL tmp_24 = muhat*tmp_0 - tmp_1;
  const CCTK_REAL tmp_25 = NRPyLeakage_Fermi_Dirac_integrals(4, tmp_24);
  const CCTK_REAL tmp_26 = tmp_25/NRPyLeakage_Fermi_Dirac_integrals(2, tmp_24);
  const CCTK_REAL tmp_27 = EnsureFinite(tmp_17*tmp_26);
  const CCTK_REAL tmp_29 = EnsureFinite(tmp_11*tmp_26*tmp_5);
  const CCTK_REAL tmp_30 = NRPyLeakage_Fermi_Dirac_integrals(5, tmp_24);
  const CCTK_REAL tmp_31 = Y_pn*tmp_14/(exp(-tmp_1 - tmp_30/tmp_25) + 1);
  const CCTK_REAL tmp_32 = EnsureFinite(tmp_26*tmp_31*tmp_5);
  const CCTK_REAL tmp_33 = tmp_30/NRPyLeakage_Fermi_Dirac_integrals(3, tmp_24);
  const CCTK_REAL tmp_34 = EnsureFinite(tmp_17*tmp_33);
  const CCTK_REAL tmp_36 = EnsureFinite(tmp_11*tmp_33*tmp_5);
  const CCTK_REAL tmp_37 = EnsureFinite(tmp_31*tmp_33*tmp_5);
  const CCTK_REAL tmp_38 = NRPyLeakage_Fermi_Dirac_integrals(4, 0)/NRPyLeakage_Fermi_Dirac_integrals(2, 0);
  const CCTK_REAL tmp_39 = EnsureFinite(tmp_17*tmp_38);
  const CCTK_REAL tmp_41 = EnsureFinite(tmp_11*tmp_38*tmp_5);
  const CCTK_REAL tmp_42 = NRPyLeakage_Fermi_Dirac_integrals(5, 0)/NRPyLeakage_Fermi_Dirac_integrals(3, 0);
  const CCTK_REAL tmp_43 = EnsureFinite(tmp_17*tmp_42);
  const CCTK_REAL tmp_44 = EnsureFinite(tmp_11*tmp_42*tmp_5);
  kappa_s_n_0_nue = tmp_10;
  kappa_s_p_0_nue = tmp_12;
  kappa_a_0_nue = tmp_16;
  kappa_s_n_1_nue = tmp_19;
  kappa_s_p_1_nue = tmp_21;
  kappa_a_1_nue = tmp_22;
  kappa_s_n_0_anue = tmp_27;
  kappa_s_p_0_anue = tmp_29;
  kappa_a_0_anue = tmp_32;
  kappa_s_n_1_anue = tmp_34;
  kappa_s_p_1_anue = tmp_36;
  kappa_a_1_anue = tmp_37;
  kappa_s_n_0_nux = tmp_39;
  kappa_s_p_0_nux = tmp_41;
  kappa_s_n_1_nux = tmp_43;
  kappa_s_p_1_nux = tmp_44;
  kappa_nue[0] = NRPyLeakage_harm_units_geom_to_cgs_L*(tmp_10 + tmp_12 + tmp_16);
  kappa_nue[1] = NRPyLeakage_harm_units_geom_to_cgs_L*(tmp_19 + tmp_21 + tmp_22);
  kappa_anue[0] = NRPyLeakage_harm_units_geom_to_cgs_L*(tmp_27 + tmp_29 + tmp_32);
  kappa_anue[1] = NRPyLeakage_harm_units_geom_to_cgs_L*(tmp_34 + tmp_36 + tmp_37);
  kappa_nux[0] = NRPyLeakage_harm_units_geom_to_cgs_L*(tmp_39 + tmp_41);
  kappa_nux[1] = NRPyLeakage_harm_units_geom_to_cgs_L*(tmp_43 + tmp_44);

  // Step 5: Make sure results are finite; if not reset to small value
  for(int i=0;i<2;i++) {
    if( !isfinite(kappa_nue [i]) ) kappa_nue [i] = 1e-15;
    if( !isfinite(kappa_anue[i]) ) kappa_anue[i] = 1e-15;
    if( !isfinite(kappa_nux [i]) ) kappa_nux [i] = 1e-15;
  }
  fprintf(stderr,"(NRPyLeakage) kappa_s_n_0_nue  = %.15e\n",kappa_s_n_0_nue);
  fprintf(stderr,"(NRPyLeakage) kappa_s_p_0_nue  = %.15e\n",kappa_s_p_0_nue);
  fprintf(stderr,"(NRPyLeakage) kappa_a_0_nue    = %.15e\n",kappa_a_0_nue);
  fprintf(stderr,"(NRPyLeakage) kappa_s_n_1_nue  = %.15e\n",kappa_s_n_1_nue);
  fprintf(stderr,"(NRPyLeakage) kappa_s_p_1_nue  = %.15e\n",kappa_s_p_1_nue);
  fprintf(stderr,"(NRPyLeakage) kappa_a_1_nue    = %.15e\n",kappa_a_1_nue);
  fprintf(stderr,"(NRPyLeakage) kappa_s_n_0_anue = %.15e\n",kappa_s_n_0_anue);
  fprintf(stderr,"(NRPyLeakage) kappa_s_p_0_anue = %.15e\n",kappa_s_p_0_anue);
  fprintf(stderr,"(NRPyLeakage) kappa_a_0_anue   = %.15e\n",kappa_a_0_anue);
  fprintf(stderr,"(NRPyLeakage) kappa_s_n_1_anue = %.15e\n",kappa_s_n_1_anue);
  fprintf(stderr,"(NRPyLeakage) kappa_s_p_1_anue = %.15e\n",kappa_s_p_1_anue);
  fprintf(stderr,"(NRPyLeakage) kappa_a_1_anue   = %.15e\n",kappa_a_1_anue);
  fprintf(stderr,"(NRPyLeakage) kappa_s_n_0_nux  = %.15e\n",kappa_s_n_0_nux);
  fprintf(stderr,"(NRPyLeakage) kappa_s_p_0_nux  = %.15e\n",kappa_s_p_0_nux);
  fprintf(stderr,"(NRPyLeakage) kappa_s_n_1_nux  = %.15e\n",kappa_s_n_1_nux);
  fprintf(stderr,"(NRPyLeakage) kappa_s_p_1_nux  = %.15e\n",kappa_s_p_1_nux);
  fprintf(stderr,"(NRPyLeakage) kappa_nue[0]     = %.15e\n",kappa_nue[0]);
  fprintf(stderr,"(NRPyLeakage) kappa_nue[1]     = %.15e\n",kappa_nue[1]);
  fprintf(stderr,"(NRPyLeakage) kappa_anue[0]    = %.15e\n",kappa_anue[0]);
  fprintf(stderr,"(NRPyLeakage) kappa_anue[1]    = %.15e\n",kappa_anue[1]);
  fprintf(stderr,"(NRPyLeakage) kappa_nux[0]     = %.15e\n",kappa_nux[0]);
  fprintf(stderr,"(NRPyLeakage) kappa_nux[1]     = %.15e\n",kappa_nux[1]);
  fprintf(stderr,"(NRPyLeakage) ************************************************\n");
}