// Thorn      : WVU_EOS
// File       : WVU_EOS_Tabulated_headers.hh
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This file contains all the variables, structs, and
//              function prototypes we need in the rest of the thorn.
//              It also provides an interface with EOS_Omni variables.

#ifndef WVU_EOS_HEADERS_H
#define WVU_EOS_HEADERS_H

// ------------------------------------------------------
// ------ Functions where the temperature is known ------
// ------------------------------------------------------
void WVU_EOS_P_and_eps_from_rho_Ye_T_impl( const CCTK_REAL rho,
                                           const CCTK_REAL Ye,
                                           const CCTK_REAL T,
                                           CCTK_REAL *restrict P,
                                           CCTK_REAL *restrict eps );

void WVU_EOS_P_eps_and_S_from_rho_Ye_T_impl( const CCTK_REAL rho,
                                             const CCTK_REAL Ye,
                                             const CCTK_REAL T,
                                             CCTK_REAL *restrict P,
                                             CCTK_REAL *restrict eps,
                                             CCTK_REAL *restrict S );

void WVU_EOS_P_eps_S_and_cs2_from_rho_Ye_T_impl( const CCTK_REAL rho,
                                                 const CCTK_REAL Ye,
                                                 const CCTK_REAL T,
                                                 CCTK_REAL *restrict P,
                                                 CCTK_REAL *restrict eps,
                                                 CCTK_REAL *restrict S,
                                                 CCTK_REAL *restrict cs2 );

void WVU_EOS_P_eps_and_depsdT_from_rho_Ye_T_impl( const CCTK_REAL rho,
                                                  const CCTK_REAL Ye,
                                                  const CCTK_REAL T,
                                                  CCTK_REAL *restrict P,
                                                  CCTK_REAL *restrict eps,
                                                  CCTK_REAL *restrict depsdT );

void WVU_EOS_P_eps_dPdrho_dPdT_depsdrho_and_depsdT_from_rho_Ye_T_impl( const CCTK_REAL rho,
                                                                       const CCTK_REAL Ye,
                                                                       const CCTK_REAL T,
                                                                       CCTK_REAL *restrict P,
                                                                       CCTK_REAL *restrict eps,
                                                                       CCTK_REAL *restrict dPdrho,
                                                                       CCTK_REAL *restrict dPdT,
                                                                       CCTK_REAL *restrict depsdrho,
                                                                       CCTK_REAL *restrict depsdT );

void WVU_EOS_mue_mup_mun_muhat_Xn_and_Xp_from_rho_Ye_T_impl( const CCTK_REAL rho,
                                                             const CCTK_REAL Ye,
                                                             const CCTK_REAL T,
                                                             CCTK_REAL *restrict mu_e,
                                                             CCTK_REAL *restrict mu_p,
                                                             CCTK_REAL *restrict mu_n,
                                                             CCTK_REAL *restrict muhat,
                                                             CCTK_REAL *restrict X_n,
                                                             CCTK_REAL *restrict X_p );

// ------------------------------------------------------
// ---- Functions where the temperature is not known ----
// ------------------------------------------------------
void WVU_EOS_P_and_T_from_rho_Ye_eps_impl( const CCTK_REAL rho,
                                           const CCTK_REAL Ye,
                                           const CCTK_REAL eps,
                                           CCTK_REAL *restrict P,
                                           CCTK_REAL *restrict T );

void WVU_EOS_P_S_and_T_from_rho_Ye_eps_impl( const CCTK_REAL rho,
                                             const CCTK_REAL Ye,
                                             const CCTK_REAL eps,
                                             CCTK_REAL *restrict P,
                                             CCTK_REAL *restrict S,
                                             CCTK_REAL *restrict T );

void WVU_EOS_eps_S_and_T_from_rho_Ye_P_impl( const CCTK_REAL rho,
                                             const CCTK_REAL Ye,
                                             const CCTK_REAL P,
                                             CCTK_REAL *restrict eps,
                                             CCTK_REAL *restrict S,
                                             CCTK_REAL *restrict T );

void WVU_EOS_P_eps_and_T_from_rho_Ye_S_impl( const CCTK_REAL rho,
                                             const CCTK_REAL Ye,
                                             const CCTK_REAL S,
                                             CCTK_REAL *restrict P,
                                             CCTK_REAL *restrict eps,
                                             CCTK_REAL *restrict T );

#endif // WVU_EOS_HEADERS_H
