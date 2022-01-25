#include "Basic_defines.h"
#include "NRPyEOS.h"
#include "NRPyLeakage.h"

/*
 * (c) Leo Werneck
 * Compute GRMHD source terms following Ruffert et al. (1996)
 * https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
void NRPyLeakage_compute_GRMHD_source_terms_nrpy_constants(const NRPyEOS_params *restrict eos_params,
                                                           const REAL rho_b,
                                                           const REAL Y_e,
                                                           const REAL T,
                                                           const REAL tau_nue,
                                                           const REAL tau_anue,
                                                           const REAL tau_nux,
                                                           REAL *restrict R_source,
                                                           REAL *restrict Q_source) {


  // Step 1: Get chemical potentials and mass
  //         fractions using the EOS
  REAL mu_e, mu_p, mu_n, muhat, X_p, X_n;
  NRPyEOS_mue_mup_mun_muhat_Xn_and_Xp_from_rho_Ye_T(eos_params,rho_b, Y_e, T,
                                                    &mu_e, &mu_p, &mu_n, &muhat, &X_p, &X_n);

  // Step 2: Compute rho_b in cgs units
  const REAL rho_b_cgs = rho_b * NRPyLeakage_units_geom_to_cgs_D;

/* Leo says: this is the way ZelmaniLeak computes Y_pn and Y_np
  // Step 3: Compute Y_{pn} and Y_{np}
  // Step 3.a: Compute Y_{pn} (See discussion below Eq. A8 in https://arxiv.org/pdf/1306.4953.pdf)
  REAL Y_pn;
  if( rho_b_cgs > 2e12 ) {
    // Step 3.a.i: Use Eqs. (A13) and (A14) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
    Y_pn = (2*Y_e - 1)/(1-exp(muhat/T));
  }
  else {
    Y_pn = 1-Y_e;
  }
  // Step 3.a.ii: Make sure Y_{pn} is nonzero
  Y_pn = MAX(Y_pn,0.0);
  
  // Step 3.b: Compute Y_{np} (Eq. A13 in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R)
  const REAL Y_np = exp(muhat/T) * Y_pn;
*/

  // Leo says: This is the way HARM3D computes Y_pn and Y_np
  // Step 2: Compute Y_{pn} and Y_{np}
  const REAL exp_metahat = exp(-muhat/T);
  // Step 2.a: Compute Y_{np}
  REAL Y_np = (Y_e < 0.5) ? (2.0*Y_e-1.0)/(exp_metahat-1.0) : (1.0-Y_e);
    
  // Step 2.b: Compute Y_{pn}
  REAL Y_pn = (Y_e > 0.5) ? exp_metahat*(2.0*Y_e-1.0)/(exp_metahat-1.0) : Y_e;
    
  // Step 2.c: Make sure both Y_np and Y_pn are non-negative
  if( Y_np < 0.0 ) Y_np = 1.0-Y_e;
  if( Y_pn < 0.0 ) Y_pn = Y_e;
  
  // Step 3: Compute diffusion timescales (FIXME)
  const REAL t_diff_0_nue  = 0.0;
  const REAL t_diff_0_anue = 0.0;
  const REAL t_diff_1_nue  = 0.0;
  const REAL t_diff_1_anue = 0.0;
  const REAL t_diff_1_nux  = 0.0;

  // Step 4: Compute the source terms
  //         Note: The code below is generated by NRPy+
  const REAL tmp_0 = (1.0/(T));
  const REAL tmp_1 = mu_e*tmp_0;
  const REAL tmp_2 = NRPyLeakage_Fermi_Dirac_integrals(4, tmp_1);
  const REAL tmp_3 = 8*M_PI*NRPyLeakage_N_A*NRPyLeakage_beta*rho_b_cgs*((3.0/8.0)*((NRPyLeakage_alpha)*(NRPyLeakage_alpha)) + 1.0/8.0)/NRPyLeakage_hc3;
  const REAL tmp_4 = ((T)*(T)*(T)*(T)*(T))*tmp_3;
  const REAL tmp_5 = NRPyLeakage_Fermi_Dirac_integrals(5, tmp_1);
  const REAL tmp_7 = -muhat*tmp_0 + tmp_1;
  const REAL tmp_8 = NRPyLeakage_enable_beta_nue*Y_pn/(exp(tmp_7 - tmp_5/tmp_2) + 1);
  const REAL tmp_9 = NRPyLeakage_Brems_zeta*NRPyLeakage_enable_brems_nui_anui*rho_b_cgs*(((X_n)*(X_n)) + (28.0/3.0)*X_n*X_p + ((X_p)*(X_p)));
  const REAL tmp_11 = NRPyLeakage_Fermi_Dirac_integrals(3, -tmp_1);
  const REAL tmp_12 = ((M_PI)*(M_PI));
  const REAL tmp_13 = NRPyLeakage_Fermi_Dirac_integrals(3, tmp_1);
  const REAL tmp_15 = NRPyLeakage_Fermi_Dirac_integrals(4, -tmp_1);
  const REAL tmp_16 = -1.0/2.0*tmp_2/tmp_13 - 1.0/2.0*tmp_15/tmp_11;
  const REAL tmp_17 = muhat*tmp_0 - tmp_1;
  const REAL tmp_18 = NRPyLeakage_C1pC2_nue_anue*NRPyLeakage_enable_pair_nue_anue/((exp(tmp_16 + tmp_17) + 1)*(exp(tmp_16 + tmp_7) + 1));
  const REAL tmp_19 = (1.0/((NRPyLeakage_hc3)*(NRPyLeakage_hc3)));
  const REAL tmp_20 = NRPyLeakage_beta*pow(T, 8)*tmp_19;
  const REAL tmp_21 = (1.0/3.0)*tmp_12 + ((mu_e)*(mu_e))/((T)*(T));
  const REAL tmp_22 = NRPyLeakage_gamma_0*sqrt(tmp_21);
  const REAL tmp_24 = ((NRPyLeakage_gamma_0)*(NRPyLeakage_gamma_0))*tmp_21/(tmp_22 + 1);
  const REAL tmp_25 = -1.0/2.0*tmp_24 - 1;
  const REAL tmp_26 = ((M_PI)*(M_PI)*(M_PI))*pow(NRPyLeakage_gamma_0, 6)*((tmp_21)*(tmp_21)*(tmp_21))*(tmp_22 + 1)*exp(-tmp_22)/NRPyLeakage_alpha_fs;
  const REAL tmp_27 = ((NRPyLeakage_C_V)*(NRPyLeakage_C_V))*NRPyLeakage_enable_plasmon_nue_anue*tmp_26/((exp(tmp_17 + tmp_25) + 1)*(exp(tmp_25 + tmp_7) + 1));
  const REAL tmp_28 = NRPyLeakage_Brems_C1*pow(T, 4.5)*tmp_9 + (16.0/9.0)*tmp_11*tmp_12*tmp_13*tmp_18*tmp_20 + (1.0/3.0)*tmp_20*tmp_27;
  const REAL tmp_29 = tmp_2*tmp_4*tmp_8 + tmp_28;
  const REAL tmp_30 = (1.0/4.0)*NRPyLeakage_hc3/M_PI;
  const REAL tmp_31 = tmp_30/((T)*(T)*(T));
  const REAL tmp_32 = NRPyLeakage_Fermi_Dirac_integrals(5, -tmp_1);
  const REAL tmp_33 = NRPyLeakage_enable_beta_anue*Y_np/(exp(tmp_17 - tmp_32/tmp_15) + 1);
  const REAL tmp_34 = tmp_15*tmp_33*tmp_4 + tmp_28;
  const REAL tmp_35 = NRPyLeakage_Brems_C2*pow(T, 5.5)*tmp_9;
  const REAL tmp_36 = pow(T, 9)*tmp_19;
  const REAL tmp_38 = (1.0/36.0)*NRPyLeakage_beta*(32*tmp_11*tmp_12*tmp_2*tmp_36 + 32*tmp_12*tmp_13*tmp_15*tmp_36);
  const REAL tmp_39 = (1.0/6.0)*NRPyLeakage_beta*tmp_36*(tmp_24 + 2);
  const REAL tmp_40 = pow(T, 6)*tmp_3;
  const REAL tmp_41 = tmp_18*tmp_38 + tmp_27*tmp_39 + tmp_35;
  const REAL tmp_42 = tmp_40*tmp_5*tmp_8 + tmp_41;
  const REAL tmp_43 = tmp_30/((T)*(T)*(T)*(T));
  const REAL tmp_45 = tmp_32*tmp_33*tmp_40 + tmp_41;
  *R_source = NRPyLeakage_amu*NRPyLeakage_units_cgs_to_geom_R*(-tmp_29/(t_diff_0_nue*tmp_29*tmp_31/NRPyLeakage_Fermi_Dirac_integrals(2, tmp_7) + 1) + tmp_34/(t_diff_0_anue*tmp_31*tmp_34/NRPyLeakage_Fermi_Dirac_integrals(2, tmp_17) + 1));
  *Q_source = NRPyLeakage_units_cgs_to_geom_Q*(-tmp_42/(t_diff_1_nue*tmp_42*tmp_43/NRPyLeakage_Fermi_Dirac_integrals(3, tmp_7) + 1) - tmp_45/(t_diff_1_anue*tmp_43*tmp_45/NRPyLeakage_Fermi_Dirac_integrals(3, tmp_17) + 1) - 4*(NRPyLeakage_C1pC2_nux_anux*NRPyLeakage_enable_pair_nux_anux*tmp_38/((exp(tmp_16) + 1)*(exp(tmp_16) + 1)) + NRPyLeakage_enable_plasmon_nux_anux*tmp_26*tmp_39*((NRPyLeakage_C_V - 1)*(NRPyLeakage_C_V - 1))/((exp(tmp_25) + 1)*(exp(tmp_25) + 1)) + tmp_35)/(t_diff_1_nux*tmp_42*tmp_43/NRPyLeakage_Fermi_Dirac_integrals(3, 0) + 1));
}
