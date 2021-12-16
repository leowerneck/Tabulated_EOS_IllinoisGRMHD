#include "Basic_defines.hh"
#include "Leakage.hh"
#include "WVU_EOS_Tabulated_headers.hh"

/*
 * (c) Leo Werneck
 * Compute free neutrino emission and cooling rates following
 * Ruffert et al. (1996)
 * https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
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
                REAL *restrict Q_free_total_nux) {


  // Step 1: Get chemical potentials and mass
  //         fractions using the EOS
  REAL mu_e, mu_n, mu_p, muhat, X_n, X_p;
  WVU_EOS_mue_mup_mun_muhat_Xn_and_Xp_from_rho_Ye_T(rho_b, Y_e, T,
                                                    &mu_e, &mu_p, &mu_n, &muhat, &X_p, &X_n);
                                                    
  // Step 2: Compute rho_b in cgs units
  const REAL rho_b_cgs = rho_b * Leakage::units_geom_to_cgs_density;

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

  // Step 4: The code below is generated by NRPy+
  const REAL tmp_0 = (1.0/(T));
  const REAL tmp_1 = mu_e*tmp_0;
  const REAL tmp_2 = Leakage_Fermi_Dirac_integrals(4, tmp_1);
  const REAL tmp_3 = 8*Leakage::beta*M_PI*rho_b_cgs*((3.0/8.0)*((Leakage::alpha)*(Leakage::alpha)) + 1.0/8.0)/Leakage::hc3;
  const REAL tmp_4 = ((T)*(T)*(T)*(T)*(T))*tmp_3;
  const REAL tmp_5 = Leakage_Fermi_Dirac_integrals(5, tmp_1);
  const REAL tmp_6 = exp(-tau_nue);
  const REAL tmp_10 = Leakage::eta_nue_0*tmp_6 + (1 - tmp_6)*(-Leakage::Q_npmass*tmp_0 - mu_n*tmp_0 + mu_p*tmp_0 + tmp_1);
  const REAL tmp_11 = Y_pn/(exp(tmp_10 - tmp_5/tmp_2) + 1);
  const REAL tmp_12 = Leakage::Brems_zeta*rho_b_cgs*(((X_n)*(X_n)) + (28.0/3.0)*X_n*X_p + ((X_p)*(X_p)));
  const REAL tmp_13 = Leakage_Fermi_Dirac_integrals(3, tmp_1);
  const REAL tmp_15 = Leakage_Fermi_Dirac_integrals(4, -tmp_1);
  const REAL tmp_16 = Leakage_Fermi_Dirac_integrals(3, -tmp_1);
  const REAL tmp_17 = -1.0/2.0*tmp_15/tmp_16 - 1.0/2.0*tmp_2/tmp_13;
  const REAL tmp_18 = exp(-tau_anue);
  const REAL tmp_19 = Leakage::eta_anue_0*tmp_18 + (1 - tmp_18)*(Leakage::Q_npmass*tmp_0 + mu_n*tmp_0 - mu_p*tmp_0 - tmp_1);
  const REAL tmp_20 = Leakage::C1pC2_nue_anue/((exp(tmp_10 + tmp_17) + 1)*(exp(tmp_17 + tmp_19) + 1));
  const REAL tmp_21 = ((M_PI)*(M_PI));
  const REAL tmp_23 = (1.0/((Leakage::hc3)*(Leakage::hc3)));
  const REAL tmp_24 = Leakage::beta*pow(T, 8)*tmp_23;
  const REAL tmp_25 = tmp_13*tmp_16*tmp_21*tmp_24;
  const REAL tmp_26 = (1.0/3.0)*tmp_21 + ((mu_e)*(mu_e))/((T)*(T));
  const REAL tmp_27 = Leakage::gamma_0*sqrt(tmp_26);
  const REAL tmp_29 = ((Leakage::gamma_0)*(Leakage::gamma_0))*tmp_26/(tmp_27 + 1);
  const REAL tmp_30 = -1.0/2.0*tmp_29 - 1;
  const REAL tmp_31 = pow(Leakage::gamma_0, 6)*((M_PI)*(M_PI)*(M_PI))*((tmp_26)*(tmp_26)*(tmp_26))*(tmp_27 + 1)*exp(-tmp_27)/Leakage::alpha_fs;
  const REAL tmp_32 = ((Leakage::C_V)*(Leakage::C_V))*tmp_31/((exp(tmp_10 + tmp_30) + 1)*(exp(tmp_19 + tmp_30) + 1));
  const REAL tmp_33 = Leakage::Brems_C1*pow(T, 4.5)*tmp_12 + (16.0/9.0)*tmp_20*tmp_25 + (1.0/3.0)*tmp_24*tmp_32;
  const REAL tmp_34 = Leakage_Fermi_Dirac_integrals(5, -tmp_1);
  const REAL tmp_35 = Y_np/(exp(tmp_19 - tmp_34/tmp_15) + 1);
  const REAL tmp_36 = Leakage::C1pC2_nux_anux/((exp(tmp_17) + 1)*(exp(tmp_17) + 1));
  const REAL tmp_37 = Leakage::Brems_C2*pow(T, 5.5)*tmp_12;
  const REAL tmp_38 = pow(T, 9)*tmp_23;
  const REAL tmp_39 = (1.0/6.0)*Leakage::beta*tmp_38*(tmp_29 + 2);
  const REAL tmp_40 = tmp_31*tmp_39*((Leakage::C_V - 1)*(Leakage::C_V - 1))/((exp(tmp_30) + 1)*(exp(tmp_30) + 1)) + tmp_37;
  const REAL tmp_41 = pow(T, 6)*tmp_3;
  const REAL tmp_43 = Leakage::beta*(32*tmp_13*tmp_15*tmp_21*tmp_38 + 32*tmp_16*tmp_2*tmp_21*tmp_38);
  const REAL tmp_44 = (1.0/36.0)*tmp_20*tmp_43 + tmp_32*tmp_39 + tmp_37;
  *R_free_total_nue = tmp_11*tmp_2*tmp_4 + tmp_33;
  *R_free_total_anue = tmp_15*tmp_35*tmp_4 + tmp_33;
  *R_free_total_nux = (64.0/9.0)*tmp_25*tmp_36 + tmp_40;
  *Q_free_total_nue = tmp_11*tmp_41*tmp_5 + tmp_44;
  *Q_free_total_anue = tmp_34*tmp_35*tmp_41 + tmp_44;
  *Q_free_total_nux = (1.0/9.0)*tmp_36*tmp_43 + tmp_40;
}