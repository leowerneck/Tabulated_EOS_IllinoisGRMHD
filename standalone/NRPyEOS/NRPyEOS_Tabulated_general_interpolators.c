// Thorn      : NRPyEOS
// File       : NRPyEOS_Tabulated_general_interpolators.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: In this file we provide a general function to retrieve
//              table quantities from (rho,Ye,T). We also provide a
//              general function which is able to get the temperature
//              from (rho,Ye,aux), where aux is any table quantity other
//              than (rho,Ye,T).
#include "Basic_defines.h"
#include "NRPyEOS_Tabulated_headers.h"
#include "NRPyEOS_Tabulated_helpers.h"

// General interpolator to get table quantities from (rho,Ye,T)
void NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( const NRPyEOS_params *restrict eos_params,
                                                     const int n,
                                                     const REAL rho,
                                                     const REAL Ye,
                                                     const REAL T,
                                                     const int *restrict tablevars_keys,
                                                     REAL *restrict tablevars,
                                                     NRPyEOS_error_report *restrict report ) {
  // This function will interpolate n table quantities from
  // (rho,Ye,T). It replaces EOS_Omni calls with keytemp = 1
  if( n > NRPyEOS_ntablekeys ) {
    fprintf(stderr,"(NRPyEOS) from_rho_Ye_T_interpolate_n_quantities: number of quantities exceed maximum allowed: %d > %d. ABORTING.",
            n,NRPyEOS_ntablekeys);
  }

  // Start by assuming no errors
  report->error = false;

  // Check table bounds for input variables
  report->error_key = NRPyEOS_checkbounds(eos_params,rho,T,Ye);
  if( report->error_key != 0 ) {
    // This should never happen, because we enforce
    // limits before calling this function
    sprintf(report->message,"from_rho_Ye_T_interpolate_n_quantities: problem with checkbounds");
    report->error = true;
    return;
  }

  // Get interpolation spots
  int idx[8];
  REAL delx,dely,delz;
  const REAL lr = log(rho);
  const REAL lt = log(T);
  NRPyEOS_get_interp_spots(eos_params,lr,lt,Ye,&delx,&dely,&delz,idx);

  for(int i=0;i<n;i++) {
    // Now perform the interpolations
    int key = tablevars_keys[i];
    REAL tablevar_out;
    NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&tablevar_out,key);

    // We have the result, but we must convert appropriately.
    // The only edge cases are P and eps, for which we obtain
    // log(P) and log(eps+eps0). We must check for them here
    if( key == NRPyEOS_press_key ) {
      tablevar_out = exp(tablevar_out);
    }
    else if( key == NRPyEOS_eps_key ) {
      tablevar_out = exp(tablevar_out) - eos_params->energy_shift;
    }

    // Then update tablevars
    tablevars[i] = tablevar_out;
  }
}

// General interpolator to recover the temperature
// and then get table quantities from (rho,Ye,T)
void NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( const NRPyEOS_params *restrict eos_params,
                                                                  const int n,
                                                                  const REAL prec,
                                                                  const REAL rho,
                                                                  const REAL Ye,
                                                                  const REAL tablevar_in,
                                                                  const int tablevar_in_key,
                                                                  const int *restrict tablevars_keys,
                                                                  REAL *restrict tablevars,
                                                                  REAL *restrict T,
                                                                  NRPyEOS_error_report *restrict report ) {
  // This function will interpolate n table quantities from
  // (rho,Ye,aux). It replaces EOS_Omni calls with keytemp != 1
  if( n > NRPyEOS_ntablekeys ) {
    fprintf(stderr,"(NRPyEOS) NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities: number of quantities exceed maximum allowed: %d > %d. ABORTING.",
            n,NRPyEOS_ntablekeys);
  }

  // Check table bounds for input variables
  report->error_key = NRPyEOS_checkbounds_kt0_noTcheck(eos_params,rho,Ye);
  if( report->error_key != 0 ) {
    // This should never happen, because we enforce
    // limits before calling this function
    sprintf(report->message,"NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities: problem with checkbounds_kt0_noTcheck");
    report->error = true;
    return;
  }

  // First step is to recover the temperature. The variable
  // tablevar_in is the one used in the temperature recovery.
  // For example, if tablevar_in = eps, then we recover T
  // using (rho,Ye,eps).
  REAL aux = tablevar_in;

  if( tablevar_in_key == NRPyEOS_press_key ) {
    // If aux = P, then we need log(P).
    aux = log(aux);
  }
  else if( tablevar_in_key == NRPyEOS_eps_key ) {
    // If aux = eps, then we need log(eps+eps0).
    // Compute eps+eps0
    aux += eos_params->energy_shift;
    // At this point, aux *must* be positive. If not, error out.
    if( aux < 0.0 ) {
      fprintf(stderr,"(NRPyEOS) NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities: found eps+energy_shift < 0.0 (%e). ABORTING.",
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
  NRPyEOS_findtemp_from_any(eos_params,tablevar_in_key,lr,lt0,Ye,aux,prec,&lt,&keyerr);

  // Now set the temperature
  *T = exp(lt);

  // Then interpolate the quantities we want from (rho,Ye,T)
  int anyerr=0;
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(eos_params,n,rho,Ye,*T,tablevars_keys,tablevars,report);
  report->error_key = keyerr;
  report->error     = anyerr;
}
