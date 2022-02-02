// Thorn      : NRPyEOS
// File       : NRPyEOS_Tabulated_headers.hh
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This file contains all the variables, structs, and
//              function prototypes we need in the rest of the thorn.
//              It also provides an interface with EOS_Omni variables.

#ifndef NRPyEOS_TABULATED_HEADERS_H
#define NRPyEOS_TABULATED_HEADERS_H

#include "Basic_defines.h"
#include "NRPyEOS.h"

// Table reader
void NRPyEOS_readtable_set_EOS_params(const char *nuceos_table_name, NRPyEOS_params *restrict eos_params);

// Free all memory allocated for the table
void NRPyEOS_free_memory(NRPyEOS_params *restrict eos_params);

// ------------------------------------------------------
// ------------- New general interpolators --------------
// ------------------------------------------------------
void NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( const NRPyEOS_params *restrict eos_params,
                                                     const int  n,
                                                     const REAL rho,
                                                     const REAL Ye,
                                                     const REAL T,
                                                     const int *restrict tablevars_keys,
                                                     REAL *restrict tablevars,
                                                     NRPyEOS_error_report *restrict report );

void NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( const NRPyEOS_params *restrict eos_params,
                                                                  const int  n,
                                                                  const REAL prec,
                                                                  const REAL rho,
                                                                  const REAL Ye,
                                                                  const REAL tablevar_in,
                                                                  const int  tablevar_in_key,
                                                                  const int *restrict tablevars_keys,
                                                                  REAL *restrict tablevars,
                                                                  REAL *restrict T,
                                                                  NRPyEOS_error_report *restrict report );

// ------------------------------------------------------
// ------ Functions where the temperature is known ------
// ------------------------------------------------------
void NRPyEOS_P_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                              const REAL rho,
                              const REAL Ye,
                              const REAL T,
                              REAL *restrict P );

void NRPyEOS_S_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                              const REAL rho,
                              const REAL Ye,
                              const REAL T,
                              REAL *restrict S );

void NRPyEOS_eps_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                const REAL rho,
                                const REAL Ye,
                                const REAL T,
                                REAL *restrict eps );

void NRPyEOS_P_and_eps_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                      const REAL rho,
                                      const REAL Ye,
                                      const REAL T,
                                      REAL *restrict P,
                                      REAL *restrict eps );

void NRPyEOS_P_eps_and_S_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                        const REAL rho,
                                        const REAL Ye,
                                        const REAL T,
                                        REAL *restrict P,
                                        REAL *restrict eps,
                                        REAL *restrict S );

void NRPyEOS_P_eps_S_and_cs2_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                            const REAL rho,
                                            const REAL Ye,
                                            const REAL T,
                                            REAL *restrict P,
                                            REAL *restrict eps,
                                            REAL *restrict S,
                                            REAL *restrict cs2 );

void NRPyEOS_P_eps_and_depsdT_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                             const REAL rho,
                                             const REAL Ye,
                                             const REAL T,
                                             REAL *restrict P,
                                             REAL *restrict eps,
                                             REAL *restrict depsdT );

void NRPyEOS_P_eps_dPdrho_dPdT_depsdrho_and_depsdT_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                                                  const REAL rho,
                                                                  const REAL Ye,
                                                                  const REAL T,
                                                                  REAL *restrict P,
                                                                  REAL *restrict eps,
                                                                  REAL *restrict dPdrho,
                                                                  REAL *restrict dPdT,
                                                                  REAL *restrict depsdrho,
                                                                  REAL *restrict depsdT );

void NRPyEOS_mue_mup_mun_muhat_Xn_and_Xp_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                                        const REAL rho,
                                                        const REAL Ye,
                                                        const REAL T,
                                                        REAL *restrict mu_e,
                                                        REAL *restrict mu_p,
                                                        REAL *restrict mu_n,
                                                        REAL *restrict muhat,
                                                        REAL *restrict X_p,
                                                        REAL *restrict X_n );

void NRPyEOS_P_eps_mue_mup_mun_and_muhat_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                                        const REAL rho,
                                                        const REAL Ye,
                                                        const REAL T,
                                                        REAL *restrict P,
                                                        REAL *restrict eps,
                                                        REAL *restrict mu_e,
                                                        REAL *restrict mu_p,
                                                        REAL *restrict mu_n,
                                                        REAL *restrict muhat );

void NRPyEOS_mue_mup_mun_and_muhat_from_rho_Ye_T( const NRPyEOS_params *restrict eos_params,
                                                  const REAL rho,
                                                  const REAL Ye,
                                                  const REAL T,
                                                  REAL *restrict mu_e,
                                                  REAL *restrict mu_p,
                                                  REAL *restrict mu_n,
                                                  REAL *restrict muhat );

// ------------------------------------------------------
// ---- Functions where the temperature is not known ----
// ------------------------------------------------------
void NRPyEOS_P_and_T_from_rho_Ye_eps( const NRPyEOS_params *restrict eos_params,
                                      const REAL rho,
                                      const REAL Ye,
                                      const REAL eps,
                                      REAL *restrict P,
                                      REAL *restrict T );

void NRPyEOS_P_S_and_T_from_rho_Ye_eps( const NRPyEOS_params *restrict eos_params,
                                        const REAL rho,
                                        const REAL Ye,
                                        const REAL eps,
                                        REAL *restrict P,
                                        REAL *restrict S,
                                        REAL *restrict T );

void NRPyEOS_eps_S_and_T_from_rho_Ye_P( const NRPyEOS_params *restrict eos_params,
                                        const REAL rho,
                                        const REAL Ye,
                                        const REAL P,
                                        REAL *restrict eps,
                                        REAL *restrict S,
                                        REAL *restrict T );

void NRPyEOS_P_eps_and_T_from_rho_Ye_S( const NRPyEOS_params *restrict eos_params,
                                        const REAL rho,
                                        const REAL Ye,
                                        const REAL S,
                                        REAL *restrict P,
                                        REAL *restrict eps,
                                        REAL *restrict T );

void NRPyEOS_Xa_Xh_Xn_Xp_and_T_from_rho_Ye_eps( const NRPyEOS_params *restrict eos_params,
                                                const REAL rho,
                                                const REAL Ye,
                                                const REAL eps,
                                                REAL *restrict X_a,
                                                REAL *restrict X_h,
                                                REAL *restrict X_n,
                                                REAL *restrict X_p,
                                                REAL *restrict T );

// ----------------------------------------

#endif // NRPyEOS_TABULATED_HEADERS_H
