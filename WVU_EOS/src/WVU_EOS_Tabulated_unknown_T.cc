// Thorn      : WVU_EOS
// File       : WVU_EOS_Tabulated_unknown_T.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This file provides wrapper functions to compute various
//              tables quantities neeeded by WVUThorns from (rho,Ye,aux),
//              where aux is any table quantity other than (rho,Ye,T).

#include "cctk.h"
#include "WVU_EOS_Tabulated_headers.hh"

// -----------------------------------
// ---------- T(rho,Ye,eps) ----------
// ---------- P(rho,Ye,T)   ----------
// -----------------------------------
extern "C" void WVU_EOS_P_and_T_from_rho_Ye_eps_impl( const CCTK_REAL rho,
                                                      const CCTK_REAL Ye,
                                                      const CCTK_REAL eps,
                                                      CCTK_REAL *restrict P,
                                                      CCTK_REAL *restrict T ) {
  // Number of interpolated quantities: 1 (P)
  const CCTK_INT n = 1;
  // Set the key to the auxiliary variable (eps)
  const CCTK_INT auxkey = WVU_EOS::eps_key;
  const CCTK_REAL aux   = eps;
  // Table variables keys
  const CCTK_INT keys[n] = {WVU_EOS::press_key};
  // Declare error variable
  WVU_EOS::eos_error_report report;
  // Set output variable array
  CCTK_REAL outvars[n];
  // Set root-finding precision
  CCTK_REAL root_finding_precision = 1e-10;

  // Get T, then P
  WVU_EOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( n,root_finding_precision,
                                                               rho,Ye,aux,auxkey, keys,outvars,T, &report );

  // Error handling
  if( report.error ) {
    CCTK_VInfo(CCTK_THORNSTRING,"Inside WVU_EOS_P_and_T_from_rho_Ye_eps. Error message: %s (key = %d)",report.message.c_str(),report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // T has already been updated, so update P
  *P = outvars[0];
}

// -----------------------------------
// ---------- T(rho,Ye,eps) ----------
// ---------- P(rho,Ye,T)   ----------
// ---------- S(rho,Ye,T)   ----------
// -----------------------------------
extern "C" void WVU_EOS_P_S_and_T_from_rho_Ye_eps_impl( const CCTK_REAL rho,
                                                        const CCTK_REAL Ye,
                                                        const CCTK_REAL eps,
                                                        CCTK_REAL *restrict P,
                                                        CCTK_REAL *restrict S,
                                                        CCTK_REAL *restrict T ) {
  // Number of interpolated quantities: 2 (P and S)
  const CCTK_INT n = 2;
  // Set the key to the auxiliary variable (eps)
  const CCTK_INT auxkey = WVU_EOS::eps_key;
  const CCTK_REAL aux   = eps;
  // Table variables keys
  const CCTK_INT keys[n] = {WVU_EOS::press_key,WVU_EOS::entropy_key};
  // Declare error variable
  WVU_EOS::eos_error_report report;
  // Set output variable array
  CCTK_REAL outvars[n];
  // Set root-finding precision
  CCTK_REAL root_finding_precision = 1e-10;

  // Get T, then P and S
  WVU_EOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( n,root_finding_precision,
                                                               rho,Ye,aux,auxkey, keys,outvars,T, &report );

  // Error handling
  if( report.error ) {
    CCTK_VInfo(CCTK_THORNSTRING,"Inside WVU_EOS_P_S_and_T_from_rho_Ye_eps. Error message: %s (key = %d)",report.message.c_str(),report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // T has already been updated, so update P and S
  *P = outvars[0];
  *S = outvars[1];
}

// -----------------------------------
// ---------- T(rho,Ye,P)   ----------
// ---------- eps(rho,Ye,T) ----------
// ---------- S(rho,Ye,T)   ----------
// -----------------------------------
extern "C" void WVU_EOS_eps_S_and_T_from_rho_Ye_P_impl( const CCTK_REAL rho,
                                                        const CCTK_REAL Ye,
                                                        const CCTK_REAL P,
                                                        CCTK_REAL *restrict eps,
                                                        CCTK_REAL *restrict S,
                                                        CCTK_REAL *restrict T ) {
  // Number of interpolated quantities: 2 (eps and S)
  const CCTK_INT n = 2;
  // Set the key to the auxiliary variable (P)
  const CCTK_INT auxkey = WVU_EOS::press_key;
  const CCTK_REAL aux   = P;
  // Table variables keys
  const CCTK_INT keys[n] = {WVU_EOS::eps_key,WVU_EOS::entropy_key};
  // Declare error variable
  WVU_EOS::eos_error_report report;
  // Set output variable array
  CCTK_REAL outvars[n];
  // Set root-finding precision
  CCTK_REAL root_finding_precision = 1e-10;

  // Get T, then eps and S
  WVU_EOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( n,root_finding_precision,
                                                               rho,Ye,aux,auxkey, keys,outvars,T, &report );

  // Error handling
  if( report.error ) {
    CCTK_VInfo(CCTK_THORNSTRING,"Inside WVU_EOS_eps_S_and_T_from_rho_Ye_P. Error message: %s (key = %d)",report.message.c_str(),report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // T has already been updated, so update eps and S
  *eps = outvars[0];
  *S   = outvars[1];
}

// -----------------------------------
// ---------- T(rho,Ye,S)   ----------
// ---------- P(rho,Ye,T)   ----------
// ---------- eps(rho,Ye,T) ----------
// -----------------------------------
extern "C" void WVU_EOS_P_eps_and_T_from_rho_Ye_S_impl( const CCTK_REAL rho,
                                                        const CCTK_REAL Ye,
                                                        const CCTK_REAL S,
                                                        CCTK_REAL *restrict P,
                                                        CCTK_REAL *restrict eps,
                                                        CCTK_REAL *restrict T ) {
  // Number of interpolated quantities: 2 (P and eps)
  const CCTK_INT n = 2;
  // Set the key to the auxiliary variable (S)
  const CCTK_INT auxkey = WVU_EOS::entropy_key;
  const CCTK_REAL aux   = S;
  // Table variables keys
  const CCTK_INT keys[n] = {WVU_EOS::press_key,WVU_EOS::eps_key};
  // Declare error variable
  WVU_EOS::eos_error_report report;
  // Set output variable array
  CCTK_REAL outvars[n];
  // Set root-finding precision
  CCTK_REAL root_finding_precision = 1e-10;

  // Get T, then P and eps
  WVU_EOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( n,root_finding_precision,
                                                               rho,Ye,aux,auxkey, keys,outvars,T, &report );

  // Error handling
  if( report.error ) {
    CCTK_VInfo(CCTK_THORNSTRING,"Inside WVU_EOS_eps_S_and_T_from_rho_Ye_P. Error message: %s (key = %d)",report.message.c_str(),report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // T has already been updated, so update P and eps
  *P   = outvars[0];
  *eps = outvars[1];
}
