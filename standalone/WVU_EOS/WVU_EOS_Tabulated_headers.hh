// Thorn      : WVU_EOS
// File       : WVU_EOS_Tabulated_headers.hh
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This file contains all the variables, structs, and
//              function prototypes we need in the rest of the thorn.
//              It also provides an interface with EOS_Omni variables.

#ifndef WVU_EOS_HEADERS_H
#define WVU_EOS_HEADERS_H

#include "Basic_defines.hh"

namespace WVU_EOS {

  // Number of entries in the EOS table
  constexpr int ntables = 19;

  // Keys for table entries
  enum table_var_key { press_key,eps_key,entropy_key,munu_key,cs2_key,depsdT_key,
                       dPdrho_key, dPdeps_key, muhat_key, mu_e_key, mu_p_key, mu_n_key,
                       Xa_key, Xh_key, Xn_key, Xp_key, Abar_key, Zbar_key, Gamma_key };

  // Name of the variables. This is only used to print
  // information about the keys during startup
  static std::string table_var_names[ntables] { "logpress","logenergy","entropy","munu","cs2","dedt",
                                                "dpdrhoe", "dpderho", "muhat", "mu_e", "mu_p", "mu_n",
                                                "Xa","Xh","Xn","Xp","Abar","Zbar","Gamma"};

  // Error handling struct
  struct eos_error_report {
    bool error;
    int error_key;
    std::string message;
  };

}

extern "C" {

// Table reader
void WVU_EOS_ReadTable(const char *nuceos_table_name);

// Free all memory allocated for the table
void WVU_EOS_free_memory();

// ------------------------------------------------------
// ------------- New general interpolators --------------
// ------------------------------------------------------
void WVU_EOS_from_rho_Ye_T_interpolate_n_quantities( const int& n,
                                                     const REAL& rho,
                                                     const REAL& Ye,
                                                     const REAL& T,
                                                     const int *restrict tablevars_keys,
                                                     REAL *restrict tablevars,
                                                     WVU_EOS::eos_error_report *restrict report );

void WVU_EOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( const int& n,
                                                                  const REAL& prec,
                                                                  const REAL& rho,
                                                                  const REAL& Ye,
                                                                  const REAL& tablevar_in,
                                                                  const int& tablevar_in_key,
                                                                  const int *restrict tablevars_keys,
                                                                  REAL *restrict tablevars,
                                                                  REAL *restrict T,
                                                                  WVU_EOS::eos_error_report *restrict report );

// ------------------------------------------------------
// ------ Functions where the temperature is known ------
// ------------------------------------------------------
void WVU_EOS_P_and_eps_from_rho_Ye_T_impl( const REAL rho,
                                           const REAL Ye,
                                           const REAL T,
                                           REAL *restrict P,
                                           REAL *restrict eps );

void WVU_EOS_P_eps_and_S_from_rho_Ye_T_impl( const REAL rho,
                                             const REAL Ye,
                                             const REAL T,
                                             REAL *restrict P,
                                             REAL *restrict eps,
                                             REAL *restrict S );

void WVU_EOS_P_eps_S_and_cs2_from_rho_Ye_T_impl( const REAL rho,
                                                 const REAL Ye,
                                                 const REAL T,
                                                 REAL *restrict P,
                                                 REAL *restrict eps,
                                                 REAL *restrict S,
                                                 REAL *restrict cs2 );

void WVU_EOS_P_eps_and_depsdT_from_rho_Ye_T_impl( const REAL rho,
                                                  const REAL Ye,
                                                  const REAL T,
                                                  REAL *restrict P,
                                                  REAL *restrict eps,
                                                  REAL *restrict depsdT );

void WVU_EOS_P_eps_dPdrho_dPdT_depsdrho_and_depsdT_from_rho_Ye_T_impl( const REAL rho,
                                                                       const REAL Ye,
                                                                       const REAL T,
                                                                       REAL *restrict P,
                                                                       REAL *restrict eps,
                                                                       REAL *restrict dPdrho,
                                                                       REAL *restrict dPdT,
                                                                       REAL *restrict depsdrho,
                                                                       REAL *restrict depsdT );

// ------------------------------------------------------
// ---- Functions where the temperature is not known ----
// ------------------------------------------------------
void WVU_EOS_P_and_T_from_rho_Ye_eps_impl( const REAL rho,
                                           const REAL Ye,
                                           const REAL eps,
                                           REAL *restrict P,
                                           REAL *restrict T );

void WVU_EOS_P_S_and_T_from_rho_Ye_eps_impl( const REAL rho,
                                             const REAL Ye,
                                             const REAL eps,
                                             REAL *restrict P,
                                             REAL *restrict S,
                                             REAL *restrict T );

void WVU_EOS_eps_S_and_T_from_rho_Ye_P_impl( const REAL rho,
                                             const REAL Ye,
                                             const REAL P,
                                             REAL *restrict eps,
                                             REAL *restrict S,
                                             REAL *restrict T );

void WVU_EOS_P_eps_and_T_from_rho_Ye_S_impl( const REAL rho,
                                             const REAL Ye,
                                             const REAL S,
                                             REAL *restrict P,
                                             REAL *restrict eps,
                                             REAL *restrict T );
} // extern "C"

// ----------------------------------------
// ---------- EOS_Omni interface ----------
// ----------------------------------------
namespace nuc_eos {
  extern double temp0, temp1;
  extern double energy_shift;

  extern double eos_rhomax, eos_rhomin;
  extern double eos_tempmin, eos_tempmax;
  extern double eos_yemin, eos_yemax;
}

namespace nuc_eos_private {
  extern int nrho;
  extern int ntemp;
  extern int nye;

  extern double * restrict alltables;
  extern double * restrict epstable;
  extern double * restrict logrho;
  extern double * restrict logtemp;
  extern double dlintemp,dlintempi;
  extern double drholintempi;
  extern double dlintempyei;
  extern double drholintempyei;
  extern double * restrict yes;
  extern double dtemp, dtempi;
  extern double drho, drhoi;
  extern double dye, dyei;
  extern double drhotempi;
  extern double drhoyei;
  extern double dtempyei;
  extern double drhotempyei;
}
// ----------------------------------------

#endif // WVU_EOS_HEADERS_H
