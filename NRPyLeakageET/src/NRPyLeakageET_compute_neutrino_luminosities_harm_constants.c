#include "cctk.h"
#include "NRPyLeakageET.h"

static CCTK_REAL EnsureFinite(const CCTK_REAL x) {
  if(isfinite(x))
    return x;
  else
    return 1e-15;
}

/*
 * (c) Leo Werneck
 * Compute neutrino luminosities following Siegel & Metzger (2018)
 * https://arxiv.org/pdf/1711.00868.pdf
 * Neutrino rates: https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
void NRPyLeakageET_compute_neutrino_luminosities_harm_constants(const CCTK_REAL alpha,
                                                                const CCTK_REAL gammaxx,
                                                                const CCTK_REAL gammaxy,
                                                                const CCTK_REAL gammaxz,
                                                                const CCTK_REAL gammayy,
                                                                const CCTK_REAL gammayz,
                                                                const CCTK_REAL gammazz,
                                                                const CCTK_REAL rho_b,
                                                                const CCTK_REAL Y_e,
                                                                const CCTK_REAL T,
                                                                const CCTK_REAL W,
                                                                const CCTK_REAL *restrict tau_nue,
                                                                const CCTK_REAL *restrict tau_anue,
                                                                const CCTK_REAL *restrict tau_nux,
                                                                CCTK_REAL *restrict lum_nue,
                                                                CCTK_REAL *restrict lum_anue,
                                                                CCTK_REAL *restrict lum_nux) {


  // Step 1: Get chemical potentials and mass
  //         fractions using the EOS
  CCTK_REAL mu_e, mu_p, mu_n, muhat, X_p, X_n;
  WVU_EOS_mue_mup_mun_muhat_Xn_and_Xp_from_rho_Ye_T(rho_b, Y_e, T, &mu_e, &mu_p, &mu_n, &muhat, &X_p, &X_n);

  // Step 2: Compute rho_b in cgs units
  const CCTK_REAL rho_b_cgs = rho_b * NRPyLeakageET_harm_units_geom_to_cgs_D;

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
  const CCTK_REAL tmp_0 = (1.0/(T));
  const CCTK_REAL tmp_1 = mu_e*tmp_0;
  const CCTK_REAL tmp_2 = NRPyLeakageET_Fermi_Dirac_integrals(4, tmp_1);
  const CCTK_REAL tmp_3 = NRPyLeakageET_Fermi_Dirac_integrals(5, tmp_1)/tmp_2;
  const CCTK_REAL tmp_4 = exp(-tau_nue[0]);
  const CCTK_REAL tmp_6 = NRPyLeakageET_eta_nue_0*tmp_4 + (1 - tmp_4)*(-muhat*tmp_0 + tmp_1);
  const CCTK_REAL tmp_7 = ((NRPyLeakageET_alpha)*(NRPyLeakageET_alpha));
  const CCTK_REAL tmp_8 = M_PI/NRPyLeakageET_harm_hc3;
  const CCTK_REAL tmp_10 = 8*NRPyLeakageET_harm_N_A*NRPyLeakageET_beta*((T)*(T)*(T)*(T)*(T))*rho_b_cgs*tmp_8*((3.0/8.0)*tmp_7 + 1.0/8.0);
  const CCTK_REAL tmp_11 = EnsureFinite(NRPyLeakageET_Brems_C2*NRPyLeakageET_enable_brems_nui_anui*T*EnsureFinite(NRPyLeakageET_Brems_C1*NRPyLeakageET_Brems_zeta*pow(T, 4.5)*rho_b_cgs*(((X_n)*(X_n)) + (28.0/3.0)*X_n*X_p + ((X_p)*(X_p))))/NRPyLeakageET_Brems_C1);
  const CCTK_REAL tmp_12 = exp(-tau_anue[0]);
  const CCTK_REAL tmp_14 = NRPyLeakageET_eta_anue_0*tmp_12 + (1 - tmp_12)*(muhat*tmp_0 - tmp_1);
  const CCTK_REAL tmp_15 = ((M_PI)*(M_PI));
  const CCTK_REAL tmp_17 = (1.0/3.0)*tmp_15 + ((mu_e)*(mu_e))/((T)*(T));
  const CCTK_REAL tmp_18 = NRPyLeakageET_gamma_0*sqrt(tmp_17);
  const CCTK_REAL tmp_20 = ((NRPyLeakageET_gamma_0)*(NRPyLeakageET_gamma_0))*tmp_17/(tmp_18 + 1);
  const CCTK_REAL tmp_21 = -1.0/2.0*tmp_20 - 1;
  const CCTK_REAL tmp_23 = (1.0/((NRPyLeakageET_harm_hc3)*(NRPyLeakageET_harm_hc3)));
  const CCTK_REAL tmp_24 = pow(T, 8);
  const CCTK_REAL tmp_26 = (1.0/3.0)*((M_PI)*(M_PI)*(M_PI))*NRPyLeakageET_beta*pow(NRPyLeakageET_gamma_0, 6)*((tmp_17)*(tmp_17)*(tmp_17))*tmp_23*tmp_24*(tmp_18 + 1)*exp(-tmp_18)/NRPyLeakageET_harm_alpha_fs;
  const CCTK_REAL tmp_27 = (1.0/2.0)*T*(tmp_20 + 2);
  const CCTK_REAL tmp_28 = NRPyLeakageET_Fermi_Dirac_integrals(3, tmp_1);
  const CCTK_REAL tmp_29 = (1.0/(tmp_28));
  const CCTK_REAL tmp_30 = NRPyLeakageET_Fermi_Dirac_integrals(4, -tmp_1);
  const CCTK_REAL tmp_31 = NRPyLeakageET_Fermi_Dirac_integrals(3, -tmp_1);
  const CCTK_REAL tmp_32 = (1.0/(tmp_31));
  const CCTK_REAL tmp_33 = -1.0/2.0*tmp_2*tmp_29 - 1.0/2.0*tmp_30*tmp_32;
  const CCTK_REAL tmp_34 = tmp_15*tmp_23*tmp_28;
  const CCTK_REAL tmp_35 = (16.0/9.0)*NRPyLeakageET_beta*tmp_24*tmp_31*tmp_34;
  const CCTK_REAL tmp_36 = 32*pow(T, 9);
  const CCTK_REAL tmp_37 = (1.0/64.0)*((NRPyLeakageET_harm_hc3)*(NRPyLeakageET_harm_hc3))*tmp_29*tmp_32*(tmp_15*tmp_2*tmp_23*tmp_31*tmp_36 + tmp_30*tmp_34*tmp_36)/(tmp_15*tmp_24);
  const CCTK_REAL tmp_38 = tmp_11 + EnsureFinite(NRPyLeakageET_enable_pair_nue_anue*tmp_37*EnsureFinite(NRPyLeakageET_C1pC2_nue_anue*tmp_35/((exp(tmp_14 + tmp_33) + 1)*(exp(tmp_33 + tmp_6) + 1)))) + EnsureFinite(NRPyLeakageET_enable_plasmon_nue_anue*tmp_27*EnsureFinite(((NRPyLeakageET_C_V)*(NRPyLeakageET_C_V))*tmp_26/((exp(tmp_14 + tmp_21) + 1)*(exp(tmp_21 + tmp_6) + 1))));
  const CCTK_REAL tmp_39 = tmp_38 + EnsureFinite(NRPyLeakageET_enable_beta_nue*T*tmp_3*EnsureFinite(Y_pn*tmp_10*tmp_2/(exp(-tmp_3 + tmp_6) + 1)));
  const CCTK_REAL tmp_40 = NRPyLeakageET_Fermi_Dirac_integrals(5, tmp_6);
  const CCTK_REAL tmp_41 = NRPyLeakageET_Fermi_Dirac_integrals(3, tmp_6);
  const CCTK_REAL tmp_42 = NRPyLeakageET_harm_N_A*NRPyLeakageET_sigma_0*((T)*(T))*rho_b_cgs/((NRPyLeakageET_harm_m_e_c2)*(NRPyLeakageET_harm_m_e_c2));
  const CCTK_REAL tmp_43 = tmp_40*tmp_42/tmp_41;
  const CCTK_REAL tmp_45 = (1 - Y_e)*((5.0/24.0)*tmp_7 + 1.0/24.0)/((2.0/3.0)*MAX(mu_n*tmp_0, 0) + 1);
  const CCTK_REAL tmp_46 = ((NRPyLeakageET_C_V - 1)*(NRPyLeakageET_C_V - 1));
  const CCTK_REAL tmp_47 = Y_e*((1.0/6.0)*tmp_46 + (5.0/24.0)*tmp_7)/((2.0/3.0)*MAX(mu_p*tmp_0, 0) + 1);
  const CCTK_REAL tmp_48 = (3.0/4.0)*tmp_7 + 1.0/4.0;
  const CCTK_REAL tmp_49 = 4*((T)*(T)*(T)*(T))*tmp_8;
  const CCTK_REAL tmp_50 = 6/NRPyLeakageET_harm_units_geom_to_cgs_L;
  const CCTK_REAL tmp_52 = NRPyLeakageET_harm_units_cgs_to_geom_Q*W*((alpha)*(alpha))*sqrt(gammaxx*gammayy*gammazz - gammaxx*((gammayz)*(gammayz)) - ((gammaxy)*(gammaxy))*gammazz + 2*gammaxy*gammaxz*gammayz - ((gammaxz)*(gammaxz))*gammayy);
  const CCTK_REAL tmp_53 = NRPyLeakageET_Fermi_Dirac_integrals(5, -tmp_1)/tmp_30;
  const CCTK_REAL tmp_54 = tmp_38 + EnsureFinite(NRPyLeakageET_enable_beta_anue*T*tmp_53*EnsureFinite(Y_np*tmp_10*tmp_30/(exp(tmp_14 - tmp_53) + 1)));
  const CCTK_REAL tmp_55 = NRPyLeakageET_Fermi_Dirac_integrals(5, tmp_14);
  const CCTK_REAL tmp_56 = NRPyLeakageET_Fermi_Dirac_integrals(3, tmp_14);
  const CCTK_REAL tmp_57 = tmp_55/tmp_56;
  const CCTK_REAL tmp_60 = NRPyLeakageET_Fermi_Dirac_integrals(3, 0);
  const CCTK_REAL tmp_61 = NRPyLeakageET_Fermi_Dirac_integrals(5, 0)/tmp_60;
  *lum_nue = tmp_39*tmp_52/(((tau_nue[1])*(tau_nue[1]))*tmp_39*tmp_50/((EnsureFinite(tmp_43*tmp_45) + EnsureFinite(tmp_43*tmp_47) + EnsureFinite(Y_np*tmp_43*tmp_48/(exp(tmp_1 - tmp_40/NRPyLeakageET_Fermi_Dirac_integrals(4, tmp_6)) + 1)))*MAX(tmp_41*tmp_49, 1.0000000000000001e-15)) + 1);
  *lum_anue = tmp_52*tmp_54/(((tau_anue[1])*(tau_anue[1]))*tmp_50*tmp_54/((EnsureFinite(tmp_42*tmp_45*tmp_57) + EnsureFinite(tmp_42*tmp_47*tmp_57) + EnsureFinite(Y_pn*tmp_42*tmp_48*tmp_57/(exp(-tmp_1 - tmp_55/NRPyLeakageET_Fermi_Dirac_integrals(4, tmp_14)) + 1)))*MAX(tmp_49*tmp_56, 1.0000000000000001e-15)) + 1);
  *lum_nux = tmp_52*(tmp_11 + EnsureFinite(NRPyLeakageET_enable_pair_nux_anux*tmp_37*EnsureFinite(NRPyLeakageET_C1pC2_nux_anux*tmp_35/((exp(tmp_33) + 1)*(exp(tmp_33) + 1)))) + EnsureFinite(NRPyLeakageET_enable_plasmon_nux_anux*tmp_27*EnsureFinite(tmp_26*tmp_46/((exp(tmp_21) + 1)*(exp(tmp_21) + 1)))))/(((tau_nux[1])*(tau_nux[1]))*tmp_39*tmp_50/((EnsureFinite(tmp_42*tmp_45*tmp_61) + EnsureFinite(tmp_42*tmp_47*tmp_61))*MAX(tmp_49*tmp_60, 1.0000000000000001e-15)) + 1);
}
