#include "cctk.h"
#include "NRPyLeakage.h"
#include "WVU_EOS_Tabulated_headers.h"

/*
 * (c) Leo Werneck
 * Compute GRMHD source terms following Ruffert et al. (1996)
 * https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
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
                                                                         REAL *restrict kappa_nux) {


  // Step 1: Get chemical potentials and mass
  //         fractions using the EOS
  REAL mu_e, mu_p, mu_n, muhat, X_p, X_n;
  WVU_EOS_mue_mup_mun_muhat_Xn_and_Xp_from_rho_Ye_T(rho_b, Y_e, T, &mu_e, &mu_p, &mu_n, &muhat, &X_p, &X_n);

  // Step 2: Compute rho_b in cgs units
  const REAL rho_b_cgs = rho_b * NRPyLeakage_units_geom_to_cgs_D;

/* Leo says: this is the way ZelmaniLeak computes Y_pn and Y_np
  // Step 3: Compute Y_{pn} and Y_{np}
  // Step 3.a: Compute Y_{pn} (See discussion below Eq. A8 in https://arxiv.org/pdf/1306.4953.pdf)
  const REAL Y_p = Y_e;
  const REAL Y_n = 1-Y_e;
  REAL Y_pn;
  if( rho_b_cgs > 2e12 ) {
    // Step 3.a.i: Use Eqs. (A13) and (A14) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
    Y_pn = (2*Y_e - 1)/(1-exp(muhat/T));
  }
  else {
    Y_pn = Y_n;
  }
  // Step 3.a.ii: Make sure Y_{pn} is nonzero
  Y_pn = MAX(Y_pn,0.0);
  
  // Step 3.b: Compute Y_{np} (Eq. A13 in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R)
  const REAL Y_np = exp(muhat/T) * Y_pn;
*/

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

  printf("rho        = %.15e\n",rho_b);
  printf("Y_e        = %.15e\n",Y_e);
  printf(" T         = %.15e\n",T);
  printf("Y_np       = %.15e\n",Y_np);
  printf("Y_pn       = %.15e\n",Y_pn);
  printf("tau_0_nue  = %.15e\n",tau_nue [0]);
  printf("tau_1_nue  = %.15e\n",tau_nue [1]);
  printf("tau_0_anue = %.15e\n",tau_anue[0]);
  printf("tau_1_anue = %.15e\n",tau_anue[1]);
  printf("tau_0_nux  = %.15e\n",tau_nux [0]);
  printf("tau_1_nux  = %.15e\n",tau_nux [1]);

  // Step 4: Compute the source terms
  //         Note: The code below is generated by NRPy+
  const REAL tmp_0 = (1.0/(T));
  const REAL tmp_1 = mu_e*tmp_0;
  const REAL tmp_2 = NRPyLeakage_Fermi_Dirac_integrals(4, tmp_1);
  const REAL tmp_3 = ((NRPyLeakage_alpha)*(NRPyLeakage_alpha));
  const REAL tmp_5 = 8*M_PI*NRPyLeakage_N_A*NRPyLeakage_beta*rho_b_cgs*((3.0/8.0)*tmp_3 + 1.0/8.0)/NRPyLeakage_hc3;
  const REAL tmp_6 = ((T)*(T)*(T)*(T)*(T))*tmp_5;
  const REAL tmp_7 = NRPyLeakage_Fermi_Dirac_integrals(5, tmp_1);
  const REAL tmp_9 = -muhat*tmp_0 + tmp_1;
  const REAL tmp_10 = NRPyLeakage_enable_beta_nue*Y_pn/(exp(tmp_9 - tmp_7/tmp_2) + 1);
  const REAL tmp_11 = NRPyLeakage_Brems_zeta*NRPyLeakage_enable_brems_nui_anui*rho_b_cgs*(((X_n)*(X_n)) + (28.0/3.0)*X_n*X_p + ((X_p)*(X_p)));
  const REAL tmp_13 = NRPyLeakage_Fermi_Dirac_integrals(3, -tmp_1);
  const REAL tmp_14 = ((M_PI)*(M_PI));
  const REAL tmp_15 = NRPyLeakage_Fermi_Dirac_integrals(3, tmp_1);
  const REAL tmp_17 = NRPyLeakage_Fermi_Dirac_integrals(4, -tmp_1);
  const REAL tmp_18 = -1.0/2.0*tmp_2/tmp_15 - 1.0/2.0*tmp_17/tmp_13;
  const REAL tmp_19 = muhat*tmp_0 - tmp_1;
  const REAL tmp_20 = NRPyLeakage_C1pC2_nue_anue*NRPyLeakage_enable_pair_nue_anue/((exp(tmp_18 + tmp_19) + 1)*(exp(tmp_18 + tmp_9) + 1));
  const REAL tmp_21 = (1.0/((NRPyLeakage_hc3)*(NRPyLeakage_hc3)));
  const REAL tmp_22 = NRPyLeakage_beta*pow(T, 8)*tmp_21;
  const REAL tmp_24 = (1.0/3.0)*tmp_14 + ((mu_e)*(mu_e))/((T)*(T));
  const REAL tmp_25 = NRPyLeakage_gamma_0*sqrt(tmp_24);
  const REAL tmp_27 = ((NRPyLeakage_gamma_0)*(NRPyLeakage_gamma_0))*tmp_24/(tmp_25 + 1);
  const REAL tmp_28 = -1.0/2.0*tmp_27 - 1;
  const REAL tmp_29 = ((M_PI)*(M_PI)*(M_PI))*pow(NRPyLeakage_gamma_0, 6)*((tmp_24)*(tmp_24)*(tmp_24))*(tmp_25 + 1)*exp(-tmp_25)/NRPyLeakage_alpha_fs;
  const REAL tmp_30 = ((NRPyLeakage_C_V)*(NRPyLeakage_C_V))*NRPyLeakage_enable_plasmon_nue_anue*tmp_29/((exp(tmp_19 + tmp_28) + 1)*(exp(tmp_28 + tmp_9) + 1));
  const REAL tmp_31 = NRPyLeakage_Brems_C1*pow(T, 4.5)*tmp_11 + (16.0/9.0)*tmp_13*tmp_14*tmp_15*tmp_20*tmp_22 + (1.0/3.0)*tmp_22*tmp_30;
  const REAL tmp_32 = tmp_10*tmp_2*tmp_6 + tmp_31;
  const REAL // Not supported in C:
// MAX
tmp_34 = (1 - Y_e)*((5.0/24.0)*tmp_3 + 1.0/24.0)/((2.0/3.0)*MAX(mu_n*tmp_0, 0) + 1);
  const REAL tmp_35 = NRPyLeakage_N_A*NRPyLeakage_sigma_0*((T)*(T))*rho_b_cgs/((NRPyLeakage_m_e_c2)*(NRPyLeakage_m_e_c2));
  const REAL tmp_36 = (1.0/(NRPyLeakage_Fermi_Dirac_integrals(2, tmp_9)));
  const REAL tmp_37 = NRPyLeakage_Fermi_Dirac_integrals(4, tmp_9);
  const REAL tmp_39 = tmp_35*tmp_36*tmp_37;
  const REAL tmp_40 = (1.0/6.0)*((NRPyLeakage_C_V - 1)*(NRPyLeakage_C_V - 1));
  const REAL // Not supported in C:
// MAX
tmp_41 = Y_e*((5.0/24.0)*tmp_3 + tmp_40)/((2.0/3.0)*MAX(mu_p*tmp_0, 0) + 1);
  const REAL tmp_42 = tmp_35*tmp_41;
  const REAL tmp_43 = (3.0/4.0)*tmp_3 + 1.0/4.0;
  const REAL tmp_44 = NRPyLeakage_Fermi_Dirac_integrals(5, tmp_9);
  const REAL tmp_45 = Y_np*tmp_43/(exp(tmp_1 - tmp_44/tmp_37) + 1);
  const REAL tmp_46 = tmp_34*tmp_39 + tmp_36*tmp_37*tmp_42 + tmp_39*tmp_45;
  const REAL tmp_47 = (3.0/2.0)*NRPyLeakage_hc3/(M_PI*NRPyLeakage_units_geom_to_cgs_L);
  const REAL tmp_48 = tmp_47/((T)*(T)*(T));
  const REAL tmp_49 = NRPyLeakage_Fermi_Dirac_integrals(5, -tmp_1);
  const REAL tmp_50 = NRPyLeakage_enable_beta_anue*Y_np/(exp(tmp_19 - tmp_49/tmp_17) + 1);
  const REAL tmp_51 = tmp_17*tmp_50*tmp_6 + tmp_31;
  const REAL tmp_52 = (1.0/(NRPyLeakage_Fermi_Dirac_integrals(2, tmp_19)));
  const REAL tmp_53 = NRPyLeakage_Fermi_Dirac_integrals(4, tmp_19);
  const REAL tmp_54 = tmp_35*tmp_52*tmp_53;
  const REAL tmp_55 = NRPyLeakage_Fermi_Dirac_integrals(5, tmp_19);
  const REAL tmp_56 = Y_pn*tmp_43/(exp(-tmp_1 - tmp_55/tmp_53) + 1);
  const REAL tmp_57 = tmp_34*tmp_54 + tmp_41*tmp_54 + tmp_54*tmp_56;
  const REAL tmp_58 = NRPyLeakage_Brems_C2*pow(T, 5.5)*tmp_11;
  const REAL tmp_59 = pow(T, 9)*tmp_21;
  const REAL tmp_61 = (1.0/36.0)*NRPyLeakage_beta*(32*tmp_13*tmp_14*tmp_2*tmp_59 + 32*tmp_14*tmp_15*tmp_17*tmp_59);
  const REAL tmp_62 = NRPyLeakage_beta*tmp_59*(tmp_27 + 2);
  const REAL tmp_63 = (1.0/(NRPyLeakage_Fermi_Dirac_integrals(3, 0)));
  const REAL tmp_64 = tmp_63*NRPyLeakage_Fermi_Dirac_integrals(5, 0);
  const REAL tmp_66 = tmp_34*tmp_35*tmp_64 + tmp_42*tmp_64;
  const REAL tmp_67 = pow(T, 6)*tmp_5;
  const REAL tmp_68 = tmp_20*tmp_61 + (1.0/6.0)*tmp_30*tmp_62 + tmp_58;
  const REAL tmp_69 = tmp_10*tmp_67*tmp_7 + tmp_68;
  const REAL tmp_70 = tmp_47/((T)*(T)*(T)*(T));
  const REAL tmp_72 = (1.0/(NRPyLeakage_Fermi_Dirac_integrals(3, tmp_9)));
  const REAL tmp_74 = tmp_35*tmp_44*tmp_72;
  const REAL tmp_75 = tmp_34*tmp_74 + tmp_42*tmp_44*tmp_72 + tmp_45*tmp_74;
  const REAL tmp_76 = tmp_49*tmp_50*tmp_67 + tmp_68;
  const REAL tmp_77 = (1.0/(NRPyLeakage_Fermi_Dirac_integrals(3, tmp_19)));
  const REAL tmp_79 = tmp_35*tmp_55*tmp_77;
  const REAL tmp_80 = tmp_34*tmp_79 + tmp_42*tmp_55*tmp_77 + tmp_56*tmp_79;
  const REAL tmp_81 = NRPyLeakage_Fermi_Dirac_integrals(4, 0)/NRPyLeakage_Fermi_Dirac_integrals(2, 0);
  *R_source = NRPyLeakage_amu*NRPyLeakage_units_cgs_to_geom_R*(-tmp_32/(((tau_nue[0])*(tau_nue[0]))*tmp_32*tmp_36*tmp_48/tmp_46 + 1) + tmp_51/(((tau_anue[0])*(tau_anue[0]))*tmp_48*tmp_51*tmp_52/tmp_57 + 1));
  *Q_source = NRPyLeakage_units_cgs_to_geom_Q*(-tmp_69/(((tau_nue[1])*(tau_nue[1]))*tmp_69*tmp_70*tmp_72/tmp_75 + 1) - tmp_76/(((tau_anue[1])*(tau_anue[1]))*tmp_70*tmp_76*tmp_77/tmp_80 + 1) - 4*(NRPyLeakage_C1pC2_nux_anux*NRPyLeakage_enable_pair_nux_anux*tmp_61/((exp(tmp_18) + 1)*(exp(tmp_18) + 1)) + NRPyLeakage_enable_plasmon_nux_anux*tmp_29*tmp_40*tmp_62/((exp(tmp_28) + 1)*(exp(tmp_28) + 1)) + tmp_58)/(((tau_nux[1])*(tau_nux[1]))*tmp_63*tmp_69*tmp_70/tmp_66 + 1));
  kappa_nue[0] = NRPyLeakage_units_geom_to_cgs_L*tmp_46;
  kappa_nue[1] = NRPyLeakage_units_geom_to_cgs_L*tmp_75;
  kappa_anue[0] = NRPyLeakage_units_geom_to_cgs_L*tmp_57;
  kappa_anue[1] = NRPyLeakage_units_geom_to_cgs_L*tmp_80;
  kappa_nux[0] = NRPyLeakage_units_geom_to_cgs_L*(tmp_34*tmp_35*tmp_81 + tmp_42*tmp_81);
  kappa_nux[1] = NRPyLeakage_units_geom_to_cgs_L*tmp_66;

  printf("R = %.15e\n",*R_source);
  printf("Q = %.15e\n",*Q_source);
}
