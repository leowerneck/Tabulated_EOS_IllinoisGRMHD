#include "Basic_defines.h"
#include "NRPyEOS.h"
#include "NRPyLeakage.h"

/*
 * (c) Leo Werneck
 * Compute GRMHD source terms following Ruffert et al. (1996)
 * https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
void NRPyLeakage_compute_GRMHD_source_terms_harm_constants(const NRPyEOS_params *restrict eos_params,
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
  const REAL rho_b_cgs = rho_b * NRPyLeakage_harm_units_geom_to_cgs_density;

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

  REAL R_beta_nue,R_beta_anue,R_pair_nue_anue,R_pair_nux_anux,R_plasmon_nue_anue,R_plasmon_nux_anux,R_Brems_nui_anui;
  REAL Q_beta_nue,Q_beta_anue,Q_pair_nue_anue,Q_pair_nux_anux,Q_plasmon_nue_anue,Q_plasmon_nux_anux,Q_Brems_nui_anui;
  REAL R_eff_nue,R_eff_anue,Q_eff_nue,Q_eff_anue,Q_eff_nux;
  const REAL tmp_0 = (1.0/(T));
  const REAL tmp_1 = mu_e*tmp_0;
  const REAL tmp_2 = NRPyLeakage_Fermi_Dirac_integrals(4, tmp_1);
  const REAL tmp_3 = 8*M_PI*NRPyLeakage_harm_N_A*NRPyLeakage_beta*rho_b_cgs*((3.0/8.0)*((NRPyLeakage_alpha)*(NRPyLeakage_alpha)) + 1.0/8.0)/NRPyLeakage_harm_hc3;
  const REAL tmp_4 = ((T)*(T)*(T)*(T)*(T))*tmp_3;
  const REAL tmp_5 = NRPyLeakage_Fermi_Dirac_integrals(5, tmp_1);
  const REAL tmp_7 = -muhat*tmp_0 + tmp_1;
  const REAL tmp_8 = NRPyLeakage_enable_beta_nue*Y_pn/(exp(tmp_7 - tmp_5/tmp_2) + 1);
  const REAL tmp_9 = tmp_2*tmp_4*tmp_8;
  const REAL tmp_11 = NRPyLeakage_Fermi_Dirac_integrals(4, -tmp_1);
  const REAL tmp_12 = NRPyLeakage_Fermi_Dirac_integrals(5, -tmp_1);
  const REAL tmp_13 = muhat*tmp_0 - tmp_1;
  const REAL tmp_14 = NRPyLeakage_enable_beta_anue*Y_np/(exp(tmp_13 - tmp_12/tmp_11) + 1);
  const REAL tmp_15 = tmp_11*tmp_14*tmp_4;
  const REAL tmp_16 = NRPyLeakage_Fermi_Dirac_integrals(3, tmp_1);
  const REAL tmp_17 = NRPyLeakage_Fermi_Dirac_integrals(3, -tmp_1);
  const REAL tmp_18 = -1.0/2.0*tmp_11/tmp_17 - 1.0/2.0*tmp_2/tmp_16;
  const REAL tmp_19 = NRPyLeakage_C1pC2_nue_anue*NRPyLeakage_enable_pair_nue_anue/((exp(tmp_13 + tmp_18) + 1)*(exp(tmp_18 + tmp_7) + 1));
  const REAL tmp_20 = ((M_PI)*(M_PI));
  const REAL tmp_22 = (1.0/((NRPyLeakage_harm_hc3)*(NRPyLeakage_harm_hc3)));
  const REAL tmp_23 = NRPyLeakage_beta*pow(T, 8)*tmp_22;
  const REAL tmp_24 = tmp_16*tmp_17*tmp_20*tmp_23;
  const REAL tmp_25 = (16.0/9.0)*tmp_19*tmp_24;
  const REAL tmp_26 = NRPyLeakage_C1pC2_nux_anux*NRPyLeakage_enable_pair_nux_anux/((exp(tmp_18) + 1)*(exp(tmp_18) + 1));
  const REAL tmp_27 = (1.0/3.0)*tmp_20 + ((mu_e)*(mu_e))/((T)*(T));
  const REAL tmp_28 = NRPyLeakage_gamma_0*sqrt(tmp_27);
  const REAL tmp_30 = ((M_PI)*(M_PI)*(M_PI))*pow(NRPyLeakage_gamma_0, 6)*((tmp_27)*(tmp_27)*(tmp_27))*(tmp_28 + 1)*exp(-tmp_28)/NRPyLeakage_harm_alpha_fs;
  const REAL tmp_31 = (1.0/3.0)*tmp_23*tmp_30;
  const REAL tmp_32 = ((NRPyLeakage_gamma_0)*(NRPyLeakage_gamma_0))*tmp_27/(tmp_28 + 1);
  const REAL tmp_33 = -1.0/2.0*tmp_32 - 1;
  const REAL tmp_34 = ((NRPyLeakage_C_V)*(NRPyLeakage_C_V))*NRPyLeakage_enable_plasmon_nue_anue/((exp(tmp_13 + tmp_33) + 1)*(exp(tmp_33 + tmp_7) + 1));
  const REAL tmp_36 = NRPyLeakage_enable_plasmon_nux_anux*((NRPyLeakage_C_V - 1)*(NRPyLeakage_C_V - 1))/((exp(tmp_33) + 1)*(exp(tmp_33) + 1));
  const REAL tmp_37 = NRPyLeakage_Brems_zeta*NRPyLeakage_enable_brems_nui_anui*rho_b_cgs*(((X_n)*(X_n)) + (28.0/3.0)*X_n*X_p + ((X_p)*(X_p)));
  const REAL tmp_38 = NRPyLeakage_Brems_C1*pow(T, 4.5)*tmp_37;
  const REAL tmp_39 = pow(T, 6)*tmp_3;
  const REAL tmp_40 = tmp_39*tmp_5*tmp_8;
  const REAL tmp_41 = tmp_12*tmp_14*tmp_39;
  const REAL tmp_42 = pow(T, 9)*tmp_22;
  const REAL tmp_44 = NRPyLeakage_beta*(32*tmp_11*tmp_16*tmp_20*tmp_42 + 32*tmp_17*tmp_2*tmp_20*tmp_42);
  const REAL tmp_45 = (1.0/36.0)*tmp_19*tmp_44;
  const REAL tmp_46 = (1.0/9.0)*tmp_26*tmp_44;
  const REAL tmp_47 = (1.0/6.0)*NRPyLeakage_beta*tmp_30*tmp_42*(tmp_32 + 2);
  const REAL tmp_50 = NRPyLeakage_Brems_C2*pow(T, 5.5)*tmp_37;
  const REAL tmp_51 = tmp_25 + tmp_31*tmp_34 + tmp_38;
  const REAL tmp_53 = (1.0/4.0)*NRPyLeakage_harm_hc3/M_PI;
  const REAL tmp_54 = tmp_53/((T)*(T)*(T));
  const REAL tmp_55 = (tmp_51 + tmp_9)/(t_diff_0_nue*tmp_54*(tmp_51 + tmp_9)/NRPyLeakage_Fermi_Dirac_integrals(2, tmp_7) + 1);
  const REAL tmp_57 = (tmp_15 + tmp_51)/(t_diff_0_anue*tmp_54*(tmp_15 + tmp_51)/NRPyLeakage_Fermi_Dirac_integrals(2, tmp_13) + 1);
  const REAL tmp_58 = tmp_34*tmp_47 + tmp_45 + tmp_50;
  const REAL tmp_60 = tmp_53/((T)*(T)*(T)*(T));
  const REAL tmp_61 = tmp_60*(tmp_40 + tmp_58);
  const REAL tmp_62 = (tmp_40 + tmp_58)/(t_diff_1_nue*tmp_61/NRPyLeakage_Fermi_Dirac_integrals(3, tmp_7) + 1);
  const REAL tmp_64 = (tmp_41 + tmp_58)/(t_diff_1_anue*tmp_60*(tmp_41 + tmp_58)/NRPyLeakage_Fermi_Dirac_integrals(3, tmp_13) + 1);
  const REAL tmp_65 = (tmp_36*tmp_47 + tmp_46 + tmp_50)/(t_diff_1_nux*tmp_61/NRPyLeakage_Fermi_Dirac_integrals(3, 0) + 1);
  R_beta_nue = tmp_9;
  R_beta_anue = tmp_15;
  R_pair_nue_anue = tmp_25;
  R_pair_nux_anux = (64.0/9.0)*tmp_24*tmp_26;
  R_plasmon_nue_anue = tmp_31*tmp_34;
  R_plasmon_nux_anux = tmp_31*tmp_36;
  R_Brems_nui_anui = tmp_38;
  Q_beta_nue = tmp_40;
  Q_beta_anue = tmp_41;
  Q_pair_nue_anue = tmp_45;
  Q_pair_nux_anux = tmp_46;
  Q_plasmon_nue_anue = tmp_34*tmp_47;
  Q_plasmon_nux_anux = tmp_36*tmp_47;
  Q_Brems_nui_anui = tmp_50;
  R_eff_nue = tmp_55;
  R_eff_anue = tmp_57;
  Q_eff_nue = tmp_62;
  Q_eff_anue = tmp_64;
  Q_eff_nux = tmp_65;
  *R_source = NRPyLeakage_amu*NRPyLeakage_harm_units_cgs_to_geom_R*(-tmp_55 + tmp_57);
  *Q_source = NRPyLeakage_harm_units_cgs_to_geom_Q*(-tmp_62 - tmp_64 - 4*tmp_65);

  fprintf(stderr,"(NRPyLeakage) ****************** debug mode ******************\n");
  // fprintf(stderr,"(NRPyLeakage) Y_pn               = %.15e\n",Y_pn);
  // fprintf(stderr,"(NRPyLeakage) mu_e               = %.15e\n",mu_e);
  // fprintf(stderr,"(NRPyLeakage) mu_p               = %.15e\n",mu_p);
  // fprintf(stderr,"(NRPyLeakage) mu_n               = %.15e\n",mu_n);
  // fprintf(stderr,"(NRPyLeakage) muhat              = %.15e\n",muhat);
  // fprintf(stderr,"(NRPyLeakage) X_p                = %.15e\n",X_p);
  // fprintf(stderr,"(NRPyLeakage) X_n                = %.15e\n",X_n);
  // fprintf(stderr,"(NRPyLeakage) eta_e              = %.15e\n",tmp_1);
  // fprintf(stderr,"(NRPyLeakage) eta_nue            = %.15e\n",tmp_7);
  // fprintf(stderr,"(NRPyLeakage) FD_5_nue           = %.15e\n",NRPyLeakage_Fermi_Dirac_integrals(5, tmp_7));
  // fprintf(stderr,"(NRPyLeakage) FD_4_nue           = %.15e\n",NRPyLeakage_Fermi_Dirac_integrals(4, tmp_7));
  // fprintf(stderr,"(NRPyLeakage) Ratio              = %.15e\n",tmp_5/tmp_2);
  // fprintf(stderr,"(NRPyLeakage) arg of exp         = %.15e\n",tmp_7 - tmp_5/tmp_2);
  // fprintf(stderr,"(NRPyLeakage) blocking_factor_ec = %.15e\n",1.0/(exp(tmp_7 - tmp_5/tmp_2) + 1));
  fprintf(stderr,"(NRPyLeakage) R_beta_nue         = %.15e\n",R_beta_nue);
  fprintf(stderr,"(NRPyLeakage) R_beta_anue        = %.15e\n",R_beta_anue);
  fprintf(stderr,"(NRPyLeakage) R_pair_nue_anue    = %.15e\n",R_pair_nue_anue);
  fprintf(stderr,"(NRPyLeakage) R_pair_nux_anux    = %.15e\n",R_pair_nux_anux);
  fprintf(stderr,"(NRPyLeakage) R_plasmon_nue_anue = %.15e\n",R_plasmon_nue_anue);
  fprintf(stderr,"(NRPyLeakage) R_plasmon_nux_anux = %.15e\n",R_plasmon_nux_anux);
  fprintf(stderr,"(NRPyLeakage) R_Brems_nui_anui   = %.15e\n",R_Brems_nui_anui);
  fprintf(stderr,"(NRPyLeakage) Q_beta_nue         = %.15e\n",Q_beta_nue);
  fprintf(stderr,"(NRPyLeakage) Q_beta_anue        = %.15e\n",Q_beta_anue);
  fprintf(stderr,"(NRPyLeakage) Q_pair_nue_anue    = %.15e\n",Q_pair_nue_anue);
  fprintf(stderr,"(NRPyLeakage) Q_pair_nux_anux    = %.15e\n",Q_pair_nux_anux);
  fprintf(stderr,"(NRPyLeakage) Q_plasmon_nue_anue = %.15e\n",Q_plasmon_nue_anue);
  fprintf(stderr,"(NRPyLeakage) Q_plasmon_nux_anux = %.15e\n",Q_plasmon_nux_anux);
  fprintf(stderr,"(NRPyLeakage) Q_Brems_nui_anui   = %.15e\n",Q_Brems_nui_anui);
  // fprintf(stderr,"(NRPyLeakage) R_eff_nue          = %.15e\n",R_eff_nue);
  // fprintf(stderr,"(NRPyLeakage) R_eff_anue         = %.15e\n",R_eff_anue);
  // fprintf(stderr,"(NRPyLeakage) Q_eff_nue          = %.15e\n",Q_eff_nue);
  // fprintf(stderr,"(NRPyLeakage) Q_eff_anue         = %.15e\n",Q_eff_anue);
  // fprintf(stderr,"(NRPyLeakage) Q_eff_nux          = %.15e\n",Q_eff_nux);
  fprintf(stderr,"(NRPyLeakage) *R_source          = %.15e\n",*R_source);
  fprintf(stderr,"(NRPyLeakage) *Q_source          = %.15e\n",*Q_source);
  fprintf(stderr,"(NRPyLeakage) ************************************************\n");
}
