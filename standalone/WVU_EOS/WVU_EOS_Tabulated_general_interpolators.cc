// Thorn      : WVU_EOS
// File       : WVU_EOS_Tabulated_general_interpolators.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: In this file we provide a general function to retrieve
//              table quantities from (rho,Ye,T). We also provide a
//              general function which is able to get the temperature
//              from (rho,Ye,aux), where aux is any table quantity other
//              than (rho,Ye,T).
#include "Basic_defines.hh"
#include "WVU_EOS_Tabulated_headers.hh"
#include "WVU_EOS_Tabulated_helpers.hh"

// General interpolator to get table quantities from (rho,Ye,T)
extern "C" void WVU_EOS_from_rho_Ye_T_interpolate_n_quantities( const int& n,
                                                                const REAL& rho,
                                                                const REAL& Ye,
                                                                const REAL& T,
                                                                const int *restrict tablevars_keys,
                                                                REAL *restrict tablevars,
                                                                WVU_EOS::eos_error_report *restrict report ) {
  // This function will interpolate n table quantities from
  // (rho,Ye,T). It replaces EOS_Omni calls with keytemp = 1
  if( n > WVU_EOS::ntables ) {
    fprintf(stderr,"(WVU_EOS) WVU_EOS::from_rho_Ye_T_interpolate_n_quantities: number of quantities exceed maximum allowed: %d > %d. ABORTING.",
            n,WVU_EOS::ntables);
  }

  // Start by assuming no errors
  report->error = false;

  // Check table bounds for input variables
  report->error_key = WVU_EOS::checkbounds(rho,T,Ye);
  if( report->error_key != 0 ) {
    // This should never happen, because we enforce
    // limits before calling this function
    report->message = "WVU_EOS::from_rho_Ye_T_interpolate_n_quantities: problem with checkbounds";
    report->error   = true;
    return;
  }

  // Get interpolation spots
  int idx[8];
  REAL delx,dely,delz;
  const REAL lr = log(rho);
  const REAL lt = log(T);
  WVU_EOS::get_interp_spots(lr,lt,Ye,&delx,&dely,&delz,idx);

  for(int i=0;i<n;i++) {
    // Now perform the interpolations
    int key = tablevars_keys[i];
    REAL tablevar_out;
    WVU_EOS::nuc_eos_C_linterp_one(idx,delx,dely,delz,&tablevar_out,key);

    // We have the result, but we must convert appropriately.
    // The only edge cases are P and eps, for which we obtain
    // log(P) and log(eps+eps0). We must check for them here
    if( key == WVU_EOS::press_key ) {
      tablevar_out = exp(tablevar_out);
    }
    else if( key == WVU_EOS::eps_key ) {
      tablevar_out = exp(tablevar_out) - nuc_eos::energy_shift;
    }

    // Then update tablevars
    tablevars[i] = tablevar_out;
  }
}

// General interpolator to recover the temperature
// and then get table quantities from (rho,Ye,T)
extern "C" void WVU_EOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( const int& n,
                                                                             const REAL& prec,
                                                                             const REAL& rho,
                                                                             const REAL& Ye,
                                                                             const REAL& tablevar_in,
                                                                             const int& tablevar_in_key,
                                                                             const int *restrict tablevars_keys,
                                                                             REAL *restrict tablevars,
                                                                             REAL *restrict T,
                                                                             WVU_EOS::eos_error_report *restrict report ) {
  // This function will interpolate n table quantities from
  // (rho,Ye,aux). It replaces EOS_Omni calls with keytemp != 1
  if( n > WVU_EOS::ntables ) {
    fprintf(stderr,"(WVU_EOS) WVU_EOS::from_rho_Ye_aux_find_T_and_interpolate_n_quantities: number of quantities exceed maximum allowed: %d > %d. ABORTING.",
            n,WVU_EOS::ntables);
  }

  // Check table bounds for input variables
  report->error_key = WVU_EOS::checkbounds_kt0_noTcheck(rho,Ye);
  if( report->error_key != 0 ) {
    // This should never happen, because we enforce
    // limits before calling this function
    report->message = "WVU_EOS::from_rho_Ye_aux_find_T_and_interpolate_n_quantities: problem with checkbounds_kt0_noTcheck";
    report->error   = true;
    return;
  }

  // First step is to recover the temperature. The variable
  // tablevar_in is the one used in the temperature recovery.
  // For example, if tablevar_in = eps, then we recover T
  // using (rho,Ye,eps).
  REAL aux = tablevar_in;

  if( tablevar_in_key == WVU_EOS::press_key ) {
    // If aux = P, then we need log(P).
    aux = log(aux);
  }
  else if( tablevar_in_key == WVU_EOS::eps_key ) {
    // If aux = eps, then we need log(eps+eps0).
    // Compute eps+eps0
    aux += nuc_eos::energy_shift;
    // At this point, aux *must* be positive. If not, error out.
    if( aux < 0.0 ) {
      fprintf(stderr,"(WVU_EOS) WVU_EOS::from_rho_Ye_aux_find_T_and_interpolate_n_quantities: found eps+energy_shift < 0.0 (%e). ABORTING.",
              aux);
    }
    // Compute log(eps+eps0)
    aux = log(aux);
  }

  // Now compute the temperature
  const REAL lr  = log(rho);
  const REAL lt0 = log(*T);
  REAL lt        = 0.0;
  int keyerr=0;
  WVU_EOS::findtemp_from_any(tablevar_in_key,lr,lt0,Ye,aux,prec,&lt,&keyerr);

  // Now set the temperature
  *T = exp(lt);

  // Then interpolate the quantities we want from (rho,Ye,T)
  int anyerr=0;
  WVU_EOS_from_rho_Ye_T_interpolate_n_quantities(n,rho,Ye,*T,tablevars_keys,tablevars,report);
  report->error_key = keyerr;
  report->error     = anyerr;
}
