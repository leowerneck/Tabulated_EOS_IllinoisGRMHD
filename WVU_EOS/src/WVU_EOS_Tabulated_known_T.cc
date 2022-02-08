// Thorn      : WVU_EOS
// File       : WVU_EOS_Tabulated_known_T.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This file provides wrapper functions to compute various
//              tables quantities neeeded by WVUThorns from (rho,Ye,T).

#include "cctk.h"
#include "WVU_EOS_Tabulated_headers.hh"

// -------------------------------------
// -----------   P(rho,Ye,T) -----------
// ----------- eps(rho,Ye,T) -----------
// -------------------------------------
extern "C"
void WVU_EOS_P_and_eps_from_rho_Ye_T_impl( const CCTK_REAL rho,
                                           const CCTK_REAL Ye,
                                           const CCTK_REAL T,
                                           CCTK_REAL *restrict P,
                                           CCTK_REAL *restrict eps ) {
  // Number of interpolated quantities: 2 (P and eps)
  const CCTK_INT n = 2;
  // Table variables keys
  const CCTK_INT keys[n] = {WVU_EOS::press_key,WVU_EOS::eps_key};
  // Declare error variable
  WVU_EOS::eos_error_report report;
  // Set output variable array
  CCTK_REAL outvars[n];

  // Get P and eps
  WVU_EOS_from_rho_Ye_T_interpolate_n_quantities( n,rho,Ye,T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    CCTK_VInfo(CCTK_THORNSTRING,"Inside WVU_EOS_P_and_eps_from_rho_Ye_T. Error message: %s (key = %d)",report.message.c_str(),report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // Then update P and eps
  *P   = outvars[0];
  *eps = outvars[1];
}

// -------------------------------------
// -----------   P(rho,Ye,T) -----------
// ----------- eps(rho,Ye,T) -----------
// -----------   S(rho,Ye,T) -----------
// -------------------------------------
extern "C"
void WVU_EOS_P_eps_and_S_from_rho_Ye_T_impl( const CCTK_REAL rho,
                                             const CCTK_REAL Ye,
                                             const CCTK_REAL T,
                                             CCTK_REAL *restrict P,
                                             CCTK_REAL *restrict eps,
                                             CCTK_REAL *restrict S ) {
  // Number of interpolated quantities: 3 (P, eps, and S)
  const CCTK_INT n = 3;
  // Table variables keys
  const CCTK_INT keys[n] = {WVU_EOS::press_key,WVU_EOS::eps_key,WVU_EOS::entropy_key};
  // Declare error variable
  WVU_EOS::eos_error_report report;
  // Set output variable array
  CCTK_REAL outvars[n];

  // Get P, eps, and S
  WVU_EOS_from_rho_Ye_T_interpolate_n_quantities( n,rho,Ye,T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    CCTK_VInfo(CCTK_THORNSTRING,"Inside WVU_EOS_P_eps_and_S_from_rho_Ye_T. Error message: %s (key = %d)",report.message.c_str(),report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // Then update P, eps, and S
  *P   = outvars[0];
  *eps = outvars[1];
  *S   = outvars[2];
}

// -------------------------------------
// -----------   P(rho,Ye,T) -----------
// ----------- eps(rho,Ye,T) -----------
// -----------   S(rho,Ye,T) -----------
// ----------- cs2(rho,Ye,T) -----------
// -------------------------------------
extern "C"
void WVU_EOS_P_eps_S_and_cs2_from_rho_Ye_T_impl( const CCTK_REAL rho,
                                                 const CCTK_REAL Ye,
                                                 const CCTK_REAL T,
                                                 CCTK_REAL *restrict P,
                                                 CCTK_REAL *restrict eps,
                                                 CCTK_REAL *restrict S,
                                                 CCTK_REAL *restrict cs2 ) {
  // Number of interpolated quantities: 3 (P, eps, S, and cs2)
  const CCTK_INT n = 4;
  // Table variables keys
  const CCTK_INT keys[n] = {WVU_EOS::press_key,WVU_EOS::eps_key,WVU_EOS::entropy_key,WVU_EOS::cs2_key};
  // Declare error variable
  WVU_EOS::eos_error_report report;
  // Set output variable array
  CCTK_REAL outvars[n];

  // Get P, eps, S, and cs2
  WVU_EOS_from_rho_Ye_T_interpolate_n_quantities( n,rho,Ye,T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    CCTK_VInfo(CCTK_THORNSTRING,"Inside WVU_EOS_P_eps_and_S_from_rho_Ye_T. Error message: %s (key = %d)",report.message.c_str(),report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // Then update P, eps, S, and cs2
  *P   = outvars[0];
  *eps = outvars[1];
  *S   = outvars[2];
  *cs2 = outvars[3];
  // Must update cs2
  *cs2 = rho * (*cs2) / (rho + rho*(*eps) + (*P));
}

// ----------------------------------------
// -----------      P(rho,Ye,T) -----------
// -----------    eps(rho,Ye,T) -----------
// ----------- depsdT(rho,Ye,T) -----------
// ----------------------------------------
extern "C"
void WVU_EOS_P_eps_and_depsdT_from_rho_Ye_T_impl( const CCTK_REAL rho,
                                                  const CCTK_REAL Ye,
                                                  const CCTK_REAL T,
                                                  CCTK_REAL *restrict P,
                                                  CCTK_REAL *restrict eps,
                                                  CCTK_REAL *restrict depsdT ) {
  // Number of interpolated quantities: 3 (P, eps, and deps/dT)
  const CCTK_INT n = 3;
  // Table variables keys
  const CCTK_INT keys[n] = {WVU_EOS::press_key,WVU_EOS::eps_key,WVU_EOS::depsdT_key};
  // Declare error variable
  WVU_EOS::eos_error_report report;
  // Set output variable array
  CCTK_REAL outvars[n];

  // Get P and eps
  WVU_EOS_from_rho_Ye_T_interpolate_n_quantities( n,rho,Ye,T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    CCTK_VInfo(CCTK_THORNSTRING,"Inside WVU_EOS_P_eps_and_S_from_rho_Ye_T. Error message: %s (key = %d)",report.message.c_str(),report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // Then update P, eps, and deps/dT
  *P      = outvars[0];
  *eps    = outvars[1];
  *depsdT = outvars[2];
}

// ----------------------------------------
// ----------        P(rho,Ye,T) ----------
// ----------      eps(rho,Ye,T) ----------
// ----------   dPdrho(rho,Ye,T) ----------
// ----------     dPdT(rho,Ye,T) ----------
// ---------- depsdrho(rho,Ye,T) ----------
// ----------   depsdT(rho,Ye,T) ----------
// ----------------------------------------
extern "C"
void WVU_EOS_P_eps_dPdrho_dPdT_depsdrho_and_depsdT_from_rho_Ye_T_impl( const CCTK_REAL rho,
                                                                       const CCTK_REAL Ye,
                                                                       const CCTK_REAL T,
                                                                       CCTK_REAL *restrict P,
                                                                       CCTK_REAL *restrict eps,
                                                                       CCTK_REAL *restrict dPdrho,
                                                                       CCTK_REAL *restrict dPdT,
                                                                       CCTK_REAL *restrict depsdrho,
                                                                       CCTK_REAL *restrict depsdT ) {
  // This function is a little different than the others. We need
  // dP/dT and deps/drho, but we can only get from the table the
  // following quantities:
  //
  // -> dP/drho
  // -> dP/deps
  // -> deps/dT
  //
  // Therefore we get all three of the quantities above and then compute:
  // .----------------------------.
  // | dP/dT = (dP/deps)(deps/dT) |
  // .----------------------------.-------------------------.
  // | deps/drho = (deps/dP)(dP/drho) = (dP/drho)/(dP/deps) |
  // .------------------------------------------------------.
  
  // Number of interpolated quantities: 5 (P, eps, dPdrho, depsdrho, and depsdT)
  const CCTK_INT n = 5;
  // Table variables keys (we use the table order here)
  const CCTK_INT keys[n] = {WVU_EOS::press_key,WVU_EOS::eps_key,WVU_EOS::depsdT_key,WVU_EOS::dPdrho_key,WVU_EOS::dPdeps_key};
  // Declare error variable
  WVU_EOS::eos_error_report report;
  // Set output variable array
  CCTK_REAL outvars[n];

  // Get P, eps, dP/drho, dP/deps, and deps/dT
  WVU_EOS_from_rho_Ye_T_interpolate_n_quantities( n,rho,Ye,T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    CCTK_VInfo(CCTK_THORNSTRING,"Inside WVU_EOS_P_eps_and_S_from_rho_Ye_T. Error message: %s (key = %d)",report.message.c_str(),report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // Auxiliary variables

  // Then update P, eps, deps/dT, and dP/drho
  *P        = outvars[0];
  *eps      = outvars[1];
  *depsdT   = outvars[2];
  *dPdrho   = outvars[3];

  // Finally compute dP/dT
  *dPdT     = outvars[4]*(*depsdT);
  // and deps/drho
  *depsdrho = (*dPdrho)/outvars[4];
}

// -------------------------------------
// ----------  mu_e(rho,Ye,T) ----------
// ----------  mu_p(rho,Ye,T) ----------
// ----------  mu_n(rho,Ye,T) ----------
// ---------- muhat(rho,Ye,T) ----------
// ----------   X_p(rho,Ye,T) ----------
// ----------   X_n(rho,Ye,T) ----------
// -------------------------------------
extern "C"
void WVU_EOS_mue_mup_mun_muhat_Xn_and_Xp_from_rho_Ye_T_impl( const CCTK_REAL rho,
                                                             const CCTK_REAL Ye,
                                                             const CCTK_REAL T,
                                                             CCTK_REAL *restrict mu_e,
                                                             CCTK_REAL *restrict mu_p,
                                                             CCTK_REAL *restrict mu_n,
                                                             CCTK_REAL *restrict muhat,
                                                             CCTK_REAL *restrict X_n,
                                                             CCTK_REAL *restrict X_p) {
  // Number of interpolated quantities: 6 (mu_e, mu_p, mu_n, mu_hat, X_p, and X_n)
  const CCTK_INT n = 6;
  // Table variables keys
  const CCTK_INT keys[n] = {WVU_EOS::muhat_key, WVU_EOS::mu_e_key, WVU_EOS::mu_p_key, WVU_EOS::mu_n_key, WVU_EOS::Xn_key, WVU_EOS::Xp_key};
  // Declare error variable
  WVU_EOS::eos_error_report report;
  // Set output variable array
  CCTK_REAL outvars[n];

  // Get P and eps
  WVU_EOS_from_rho_Ye_T_interpolate_n_quantities( n,rho,Ye,T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    CCTK_VInfo(CCTK_THORNSTRING,"Inside WVU_EOS_mue_mup_mun_muhat_Xn_and_Xp_from_rho_Ye_T. Error message: %s (key = %d)",report.message.c_str(),report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // Then update mu_e, mu_p, mu_n, mu_hat, X_p, and X_n
  *muhat = outvars[0];
  *mu_e  = outvars[1];
  *mu_p  = outvars[2];
  *mu_n  = outvars[3];
  *X_n   = outvars[4];
  *X_p   = outvars[5];
}
