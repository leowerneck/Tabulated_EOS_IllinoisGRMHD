// Thorn      : NRPyEOS
// File       : NRPyEOS_Tabulated_unknown_T.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This file provides wrapper functions to compute various
//              tables quantities neeeded by WVUThorns from (rho,Ye,aux),
//              where aux is any table quantity other than (rho,Ye,T).

#include "NRPyEOS.h"

// -----------------------------------
// ---------- T(rho,Ye,eps) ----------
// ---------- P(rho,Ye,T)   ----------
// -----------------------------------
void NRPyEOS_P_and_T_from_rho_Ye_eps( const NRPyEOS_params *restrict eos_params,
                                      const REAL rho,
                                      const REAL Ye,
                                      const REAL eps,
                                      REAL *restrict P,
                                      REAL *restrict T ) {
  // Number of interpolated quantities: 1 (P)
  const int n = 1;
  // Set the key to the auxiliary variable (eps)
  const int auxkey = NRPyEOS_eps_key;
  const REAL aux   = eps;
  // Table variables keys
  const int keys[1] = {NRPyEOS_press_key};
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];
  // Set root-finding precision
  REAL root_finding_precision = 1e-10;

  // Get T, then P
  NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( eos_params, n,root_finding_precision,
                                                               rho,Ye,aux,auxkey, keys,outvars,T, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_P_and_T_from_rho_Ye_eps. Error message: %s (key = %d)",report.message,report.error_key);
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
void NRPyEOS_P_S_and_T_from_rho_Ye_eps( const NRPyEOS_params *restrict eos_params,
                                        const REAL rho,
                                        const REAL Ye,
                                        const REAL eps,
                                        REAL *restrict P,
                                        REAL *restrict S,
                                        REAL *restrict T ) {
  // Number of interpolated quantities: 2 (P and S)
  const int n = 2;
  // Set the key to the auxiliary variable (eps)
  const int auxkey = NRPyEOS_eps_key;
  const REAL aux   = eps;
  // Table variables keys
  const int keys[2] = {NRPyEOS_press_key,NRPyEOS_entropy_key};
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];
  // Set root-finding precision
  REAL root_finding_precision = 1e-10;

  // Get T, then P and S
  NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( eos_params, n,root_finding_precision,
                                                               rho,Ye,aux,auxkey, keys,outvars,T, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_P_S_and_T_from_rho_Ye_eps. Error message: %s (key = %d)",report.message,report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // T has already been updated, so update P and S
  *P = outvars[0];
  *S = outvars[1];
}

// --------------------------------------
// ----------    T(rho,Ye,eps) ----------
// ----------      P(rho,Ye,T) ----------
// ----------      S(rho,Ye,T) ----------
// ---------- depsdT(rho,Ye,T) ----------
// --------------------------------------
void NRPyEOS_P_S_T_and_depsdT_from_rho_Ye_eps( const NRPyEOS_params *restrict eos_params,
                                               const REAL rho,
                                               const REAL Ye,
                                               const REAL eps,
                                               REAL *restrict P,
                                               REAL *restrict S,
                                               REAL *restrict depsdT,
                                               REAL *restrict T ) {
  // Number of interpolated quantities: 3 (P, S, and depsdT)
  const int n = 3;
  // Set the key to the auxiliary variable (eps)
  const int auxkey = NRPyEOS_eps_key;
  const REAL aux   = eps;
  // Table variables keys
  const int keys[3] = {NRPyEOS_press_key,NRPyEOS_entropy_key,NRPyEOS_depsdT_key};
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];
  // Set root-finding precision
  REAL root_finding_precision = 1e-10;

  // Get T, then P and S
  NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( eos_params, n,root_finding_precision,
                                                               rho,Ye,aux,auxkey, keys,outvars,T, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_P_S_T_and_depsdT_from_rho_Ye_eps. Error message: %s (key = %d)",report.message,report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // T has already been updated, so update P and S
  *P      = outvars[0];
  *S      = outvars[1];
  *depsdT = outvars[2];
}

// -----------------------------------
// ---------- T(rho,Ye,P)   ----------
// ---------- eps(rho,Ye,T) ----------
// ---------- S(rho,Ye,T)   ----------
// -----------------------------------
void NRPyEOS_eps_S_and_T_from_rho_Ye_P( const NRPyEOS_params *restrict eos_params,
                                        const REAL rho,
                                        const REAL Ye,
                                        const REAL P,
                                        REAL *restrict eps,
                                        REAL *restrict S,
                                        REAL *restrict T ) {
  // Number of interpolated quantities: 2 (eps and S)
  const int n = 2;
  // Set the key to the auxiliary variable (P)
  const int auxkey = NRPyEOS_press_key;
  const REAL aux   = P;
  // Table variables keys
  const int keys[2] = {NRPyEOS_eps_key,NRPyEOS_entropy_key};
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];
  // Set root-finding precision
  REAL root_finding_precision = 1e-10;

  // Get T, then eps and S
  NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( eos_params, n,root_finding_precision,
                                                               rho,Ye,aux,auxkey, keys,outvars,T, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_eps_S_and_T_from_rho_Ye_P. Error message: %s (key = %d)",report.message,report.error_key);
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
void NRPyEOS_P_eps_and_T_from_rho_Ye_S( const NRPyEOS_params *restrict eos_params,
                                        const REAL rho,
                                        const REAL Ye,
                                        const REAL S,
                                        REAL *restrict P,
                                        REAL *restrict eps,
                                        REAL *restrict T ) {
  // Number of interpolated quantities: 2 (P and eps)
  const int n = 2;
  // Set the key to the auxiliary variable (S)
  const int auxkey = NRPyEOS_entropy_key;
  const REAL aux   = S;
  // Table variables keys
  const int keys[2] = {NRPyEOS_press_key,NRPyEOS_eps_key};
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];
  // Set root-finding precision
  REAL root_finding_precision = 1e-10;

  // Get T, then P and eps
  NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( eos_params, n,root_finding_precision,
                                                               rho,Ye,aux,auxkey, keys,outvars,T, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_eps_S_and_T_from_rho_Ye_P. Error message: %s (key = %d)",report.message,report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // T has already been updated, so update P and eps
  *P   = outvars[0];
  *eps = outvars[1];
}
