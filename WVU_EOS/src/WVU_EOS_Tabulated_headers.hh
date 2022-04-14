// Thorn      : WVU_EOS
// File       : WVU_EOS_Tabulated_headers.hh
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This file contains all the variables, structs, and
//              function prototypes we need in the rest of the thorn.
//              It also provides an interface with EOS_Omni variables.

#ifndef WVU_EOS_HEADERS_HH
#define WVU_EOS_HEADERS_HH

#include <string>

#ifdef MAX
#undef MAX
#endif

#ifdef MIN
#undef MIN
#endif

// MAX and MIN macros
#define MAX(a,b) ( (a) > (b) ? (a) : (b) )
#define MIN(a,b) ( (a) < (b) ? (a) : (b) )

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
    CCTK_INT error_key;
    std::string message;
  };

}

typedef struct WVU_EOS_slice_const_Ye_and_S {
  CCTK_REAL Y_e,S,*lr_arr,*lt_arr,*lh_arr;
} WVU_EOS_slice_const_Ye_and_S;

extern "C" {

// ------------------------------------------------------
// ------------- New general interpolators --------------
// ------------------------------------------------------
void WVU_EOS_from_rho_Ye_T_interpolate_n_quantities( const CCTK_INT& n,
                                                     const CCTK_REAL& rho,
                                                     const CCTK_REAL& Ye,
                                                     const CCTK_REAL& T,
                                                     const CCTK_INT *restrict tablevars_keys,
                                                     CCTK_REAL *restrict tablevars,
                                                     WVU_EOS::eos_error_report *restrict report );

void WVU_EOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( const CCTK_INT& n,
                                                                  const CCTK_REAL& prec,
                                                                  const CCTK_REAL& rho,
                                                                  const CCTK_REAL& Ye,
                                                                  const CCTK_REAL& tablevar_in,
                                                                  const CCTK_INT& tablevar_in_key,
                                                                  const CCTK_INT *restrict tablevars_keys,
                                                                  CCTK_REAL *restrict tablevars,
                                                                  CCTK_REAL *restrict T,
                                                                  WVU_EOS::eos_error_report *restrict report );

// ------------------------------------------------------
// ------ Functions where the temperature is known ------
// ------------------------------------------------------
void WVU_EOS_P_from_rho_Ye_T_impl( const CCTK_REAL rho,
                                   const CCTK_REAL Ye,
                                   const CCTK_REAL T,
                                   CCTK_REAL *restrict P );

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
                                                             CCTK_REAL *restrict X_p);

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

// Table slicer
WVU_EOS_slice_const_Ye_and_S WVU_EOS_slice_eos_table_constant_Ye_and_S( const double Y_e, const double S );
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

#endif // WVU_EOS_HEADERS_HH
