// Thorn      : NRPyEOS
// File       : NRPyEOS_Tabulated_known_T.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This file provides wrapper functions to compute various
//              tables quantities neeeded by WVUThorns from (rho,Ye,T).

#include "NRPyEOS_Tabulated_headers.h"

// -------------------------------------
// -----------  P(rho,Ye,T)  -----------
// -------------------------------------
void NRPyEOS_P_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                              const REAL rho,
                              const REAL Ye,
                              const REAL T,
                              REAL *restrict P ) {
  // Number of interpolated quantities: 1 (P)
  const int n = 1;
  // Table variables keys
  const int keys[1] = {NRPyEOS_press_key};
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];

  // Get P and eps
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, n,rho,Ye,T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_P_from_rho_Ye_T. Error message: %s (key = %d)\n",report.message,report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // Then update P
  *P = outvars[0];
}

// -------------------------------------
// -----------  S(rho,Ye,T)  -----------
// -------------------------------------
void NRPyEOS_S_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                              const REAL rho,
                              const REAL Ye,
                              const REAL T,
                              REAL *restrict S ) {
  // Number of interpolated quantities: 1 (P)
  const int n = 1;
  // Table variables keys
  const int keys[1] = {NRPyEOS_entropy_key};
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];

  // Get P and eps
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, n,rho,Ye,T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_P_from_rho_Ye_T. Error message: %s (key = %d)\n",report.message,report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // Then update P
  *S = outvars[0];
}

// -------------------------------------
// ----------- eps(rho,Ye,T) -----------
// -------------------------------------
void NRPyEOS_eps_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                const REAL rho,
                                const REAL Ye,
                                const REAL T,
                                REAL *restrict eps ) {
  // Number of interpolated quantities: 1 (eps)
  const int n = 1;
  // Table variables keys
  const int keys[1] = {NRPyEOS_eps_key};
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];

  // Get P and eps
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, n,rho,Ye,T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_eps_from_rho_Ye_T. Error message: %s (key = %d)\n",report.message,report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // Then update eps
  *eps = outvars[0];
}

// -------------------------------------
// -----------   P(rho,Ye,T) -----------
// ----------- eps(rho,Ye,T) -----------
// -------------------------------------
void NRPyEOS_P_and_eps_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                      const REAL rho,
                                      const REAL Ye,
                                      const REAL T,
                                      REAL *restrict P,
                                      REAL *restrict eps ) {
  // Number of interpolated quantities: 2 (P and eps)
  const int n = 2;
  // Table variables keys
  const int keys[2] = {NRPyEOS_press_key,NRPyEOS_eps_key};
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];

  // Get P and eps
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, n,rho,Ye,T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_P_and_eps_from_rho_Ye_T. Error message: %s (key = %d)\n",report.message,report.error_key);
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
void NRPyEOS_P_eps_and_S_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                        const REAL rho,
                                        const REAL Ye,
                                        const REAL T,
                                        REAL *restrict P,
                                        REAL *restrict eps,
                                        REAL *restrict S ) {
  // Number of interpolated quantities: 3 (P, eps, and S)
  const int n = 3;
  // Table variables keys
  const int keys[3] = {NRPyEOS_press_key,NRPyEOS_eps_key,NRPyEOS_entropy_key};
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];

  // Get P, eps, and S
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, n,rho,Ye,T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_P_eps_and_S_from_rho_Ye_T. Error message: %s (key = %d)\n",report.message,report.error_key);
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
void NRPyEOS_P_eps_S_and_cs2_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                            const REAL rho,
                                            const REAL Ye,
                                            const REAL T,
                                            REAL *restrict P,
                                            REAL *restrict eps,
                                            REAL *restrict S,
                                            REAL *restrict cs2 ) {
  // Number of interpolated quantities: 3 (P, eps, S, and cs2)
  const int n = 4;
  // Table variables keys
  const int keys[4] = {NRPyEOS_press_key,NRPyEOS_eps_key,NRPyEOS_entropy_key,NRPyEOS_cs2_key};
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];

  // Get P, eps, S, and cs2
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, n,rho,Ye,T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_P_eps_and_S_from_rho_Ye_T. Error message: %s (key = %d)\n",report.message,report.error_key);
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
void NRPyEOS_P_eps_and_depsdT_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                             const REAL rho,
                                             const REAL Ye,
                                             const REAL T,
                                             REAL *restrict P,
                                             REAL *restrict eps,
                                             REAL *restrict depsdT ) {
  // Number of interpolated quantities: 3 (P, eps, and deps/dT)
  const int n = 3;
  // Table variables keys
  const int keys[3] = {NRPyEOS_press_key,NRPyEOS_eps_key,NRPyEOS_depsdT_key};
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];

  // Get P and eps
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, n,rho,Ye,T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_P_eps_and_S_from_rho_Ye_T. Error message: %s (key = %d)\n",report.message,report.error_key);
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
void NRPyEOS_P_eps_dPdrho_dPdT_depsdrho_and_depsdT_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                                                  const REAL rho,
                                                                  const REAL Ye,
                                                                  const REAL T,
                                                                  REAL *restrict P,
                                                                  REAL *restrict eps,
                                                                  REAL *restrict dPdrho,
                                                                  REAL *restrict dPdT,
                                                                  REAL *restrict depsdrho,
                                                                  REAL *restrict depsdT ) {
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
  const int n = 5;
  // Table variables keys (we use the table order here)
  const int keys[5] = {NRPyEOS_press_key,NRPyEOS_eps_key,NRPyEOS_depsdT_key,NRPyEOS_dPdrho_key,NRPyEOS_dPdeps_key};
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];

  // Get P, eps, dP/drho, dP/deps, and deps/dT
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, n,rho,Ye,T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_P_eps_and_S_from_rho_Ye_T. Error message: %s (key = %d)\n",report.message,report.error_key);
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
// -------------------------------------
void NRPyEOS_mue_mup_mun_and_muhat_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                                  const REAL rho,
                                                  const REAL Ye,
                                                  const REAL T,
                                                  REAL *restrict mu_e,
                                                  REAL *restrict mu_p,
                                                  REAL *restrict mu_n,
                                                  REAL *restrict muhat ) {
  // Number of interpolated quantities: 4 (mu_e, mu_p, mu_n, and mu_hat)
  const int n = 4;
  // Table variables keys
  const int keys[4] = {NRPyEOS_muhat_key, NRPyEOS_mu_e_key, NRPyEOS_mu_p_key, NRPyEOS_mu_n_key};
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];

  // Get P and eps
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, n,rho,Ye,T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_mue_mup_mun_and_muhat_from_rho_Ye_T. Error message: %s (key = %d)\n",report.message,report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // Then update mu_e, mu_p, mu_n, mu_hat, X_p, and X_n
  *muhat = outvars[0];
  *mu_e  = outvars[1];
  *mu_p  = outvars[2];
  *mu_n  = outvars[3];
}

// -------------------------------------
// ----------  mu_e(rho,Ye,T) ----------
// ----------  mu_p(rho,Ye,T) ----------
// ----------  mu_n(rho,Ye,T) ----------
// ---------- muhat(rho,Ye,T) ----------
// ----------   X_p(rho,Ye,T) ----------
// ----------   X_n(rho,Ye,T) ----------
// -------------------------------------
void NRPyEOS_mue_mup_mun_muhat_Xn_and_Xp_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                                        const REAL rho,
                                                        const REAL Ye,
                                                        const REAL T,
                                                        REAL *restrict mu_e,
                                                        REAL *restrict mu_p,
                                                        REAL *restrict mu_n,
                                                        REAL *restrict muhat,
                                                        REAL *restrict X_p,
                                                        REAL *restrict X_n ) {
  // Number of interpolated quantities: 6 (mu_e, mu_p, mu_n, mu_hat, X_p, and X_n)
  const int n = 6;
  // Table variables keys
  const int keys[6] = {NRPyEOS_muhat_key, NRPyEOS_mu_e_key, NRPyEOS_mu_p_key, NRPyEOS_mu_n_key, NRPyEOS_Xn_key, NRPyEOS_Xp_key};
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];

  // Get P and eps
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, n,rho,Ye,T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_mue_mup_mun_muhat_Xn_and_Xp_from_rho_Ye_T. Error message: %s (key = %d)\n",report.message,report.error_key);
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

// -------------------------------------
// ----------  mu_e(rho,Ye,T) ----------
// ----------  mu_p(rho,Ye,T) ----------
// ----------  mu_n(rho,Ye,T) ----------
// ---------- muhat(rho,Ye,T) ----------
// ----------   X_p(rho,Ye,T) ----------
// ----------   X_n(rho,Ye,T) ----------
// -------------------------------------
void NRPyEOS_P_eps_mue_mup_mun_and_muhat_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                                        const REAL rho,
                                                        const REAL Ye,
                                                        const REAL T,
                                                        REAL *restrict P,
                                                        REAL *restrict eps,
                                                        REAL *restrict mu_e,
                                                        REAL *restrict mu_p,
                                                        REAL *restrict mu_n,
                                                        REAL *restrict muhat ) {
  // Number of interpolated quantities: 6 (P, eps, mu_e, mu_p, mu_n, and mu_hat)
  const int n = 6;
  // Table variables keys
  const int keys[6] = {NRPyEOS_press_key, NRPyEOS_eps_key, NRPyEOS_muhat_key, NRPyEOS_mu_e_key, NRPyEOS_mu_p_key, NRPyEOS_mu_n_key};
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];

  // Get P and eps
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, n,rho,Ye,T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_P_eps_mue_mup_mun_and_muhat_from_rho_Ye_T. Error message: %s (key = %d)\n",report.message,report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // Then update P, eps, mu_e, mu_p, mu_n, and mu_hat
  *P     = outvars[0];
  *eps   = outvars[1];
  *muhat = outvars[2];
  *mu_e  = outvars[3];
  *mu_p  = outvars[4];
  *mu_n  = outvars[5];
}

// -------------------------------------
// ----------  mu_e(rho,Ye,T) ----------
// ----------  mu_p(rho,Ye,T) ----------
// ----------  mu_n(rho,Ye,T) ----------
// ---------- muhat(rho,Ye,T) ----------
// ----------   X_p(rho,Ye,T) ----------
// ----------   X_n(rho,Ye,T) ----------
// -------------------------------------
void NRPyEOS_mue_mup_mun_muhat_Xn_Xp_and_Xh_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                                           const REAL rho,
                                                           const REAL Ye,
                                                           const REAL T,
                                                           REAL *restrict mu_e,
                                                           REAL *restrict mu_p,
                                                           REAL *restrict mu_n,
                                                           REAL *restrict muhat,
                                                           REAL *restrict X_p,
                                                           REAL *restrict X_n,
                                                           REAL *restrict X_h ) {
  // Number of interpolated quantities: 6 (mu_e, mu_p, mu_n, mu_hat, X_p, and X_n)
  const int n = 7;
  // Table variables keys
  const int keys[7] = {NRPyEOS_muhat_key, NRPyEOS_mu_e_key, NRPyEOS_mu_p_key, NRPyEOS_mu_n_key, NRPyEOS_Xn_key, NRPyEOS_Xn_key, NRPyEOS_Xp_key};
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];

  // Get P and eps
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, n,rho,Ye,T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_mue_mup_mun_muhat_Xn_and_Xp_from_rho_Ye_T. Error message: %s (key = %d)\n",report.message,report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // Then update mu_e, mu_p, mu_n, mu_hat, X_p, and X_n
  *muhat = outvars[0];
  *mu_e  = outvars[1];
  *mu_p  = outvars[2];
  *mu_n  = outvars[3];
  *X_h   = outvars[4];
  *X_n   = outvars[5];
  *X_p   = outvars[6];
}

// -------------------------------------
// -----------   P(rho,Ye,T) -----------
// ----------- eps(rho,Ye,T) -----------
// -------------------------------------
void nrpyeos_p_and_eps_from_rho_ye_t_fortran( const NRPyEOS_params *restrict eos_params,
                                              const REAL *restrict rho,
                                              const REAL *restrict Ye,
                                              const REAL *restrict T,
                                              REAL *restrict P,
                                              REAL *restrict eps ) {
  // Number of interpolated quantities: 2 (P and eps)
  const int n = 2;
  // Table variables keys
  const int keys[2] = {NRPyEOS_press_key,NRPyEOS_eps_key};
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];

  // Get P and eps
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, n,*rho,*Ye,*T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_P_and_eps_from_rho_Ye_T. Error message: %s (key = %d)\n",report.message,report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // Then update P and eps
  *P   = outvars[0];
  *eps = outvars[1];
}

void nrpyeos_full( const NRPyEOS_params *restrict eos_params,
                   const REAL *restrict rho,
                   const REAL *restrict Ye,
                   const REAL *restrict T,
                   REAL *restrict P,
                   REAL *restrict eps,
                   REAL *restrict S,
                   REAL *restrict cs2,
                   REAL *restrict depsdT,
                   REAL *restrict dPdeps,
                   REAL *restrict dPdrho,
                   REAL *restrict X_a,
                   REAL *restrict X_h,
                   REAL *restrict X_n,
                   REAL *restrict X_p,
                   REAL *restrict Abar,
                   REAL *restrict Zbar,
                   REAL *restrict mu_e,
                   REAL *restrict mu_n,
                   REAL *restrict mu_p,
                   REAL *restrict muhat ) {
  // Number of interpolated quantities: 17
  const int n = 17;
  // Table variables keys
  const int keys[17] = {
    NRPyEOS_press_key,NRPyEOS_eps_key,NRPyEOS_entropy_key,NRPyEOS_cs2_key,
    NRPyEOS_depsdT_key,NRPyEOS_dPdrho_key,NRPyEOS_dPdeps_key,
    NRPyEOS_muhat_key,NRPyEOS_mu_e_key,NRPyEOS_mu_p_key,NRPyEOS_mu_n_key,
    NRPyEOS_Xa_key,NRPyEOS_Xh_key,NRPyEOS_Xn_key,NRPyEOS_Xp_key,
    NRPyEOS_Abar_key,NRPyEOS_Zbar_key
  };
  // Declare error variable
  NRPyEOS_error_report report;
  // Set output variable array
  REAL outvars[n];

  // Get P and eps
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, n,*rho,*Ye,*T, keys,outvars, &report );

  // Error handling
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside nrpyeos_full. Error message: %s (key = %d)\n",report.message,report.error_key);
    // May want to terminate depending on the error. We'll just warn for now.
  }

  // Then update output
  *P      = outvars[ 0];
  *eps    = outvars[ 1];
  *S      = outvars[ 2];
  *cs2    = outvars[ 3];
  *depsdT = outvars[ 4];
  *dPdrho = outvars[ 5];
  *dPdeps = outvars[ 6];
  *muhat  = outvars[ 7];
  *mu_e   = outvars[ 8];
  *mu_p   = outvars[ 9];
  *mu_n   = outvars[10];
  *X_a    = outvars[11];
  *X_h    = outvars[12];
  *X_n    = outvars[13];
  *X_p    = outvars[14];
  *Abar   = outvars[15];
  *Zbar   = outvars[16];
}
