// Thorn      : IllinoisGRMHD
// File       : EOS_floors_ceilings.hh
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This file contains a function which applies floors
//              and ceilings to the appropriate hydro quantities
#include "cctk.h"
#include "cctk_Parameters.h"

#include "IllinoisGRMHD_headers.h"
#include "EOS_headers.hh"

void apply_floors_and_ceilings_to_prims__recompute_prims( const igm_eos_parameters eos,
                                                          const CCTK_REAL *restrict METRIC_LAP_PSI4,
                                                          CCTK_REAL *restrict PRIMS ) {

  DECLARE_CCTK_PARAMETERS;

  // The density floor and ceiling is always applied
  PRIMS[RHOB] = MIN(MAX(PRIMS[RHOB],eos.rho_atm),eos.rho_max);

  // Hybrid EOS specific floors and ceilings
  if( eos.is_Hybrid ) {
    // Pressure and epsilon must be recomputed
    // Compute P and eps
    CCTK_REAL prs_cold = 0.0;
    CCTK_REAL eps_cold = 0.0;
    compute_P_cold__eps_cold(eos,PRIMS[RHOB],prs_cold,eps_cold);
    // Set P_min and P_max
    CCTK_REAL P_min = 0.9*prs_cold;
    CCTK_REAL P_max = 100.0*prs_cold;
    // Set Psi6
    const CCTK_REAL Psi6  = METRIC_LAP_PSI4[PSI6];
    // Adjust P_max based on Psi6
    if(Psi6 > Psi6threshold) P_max = 1e5*prs_cold; // <-- better than 10.
    // Now apply floors and ceilings to P
    if(PRIMS[PRESSURE]<P_min) PRIMS[PRESSURE] = P_min;
    // Finally, perform the last check
    if((PRIMS[RHOB] < 100.0*eos.rho_atm || Psi6 > Psi6threshold) && PRIMS[PRESSURE]>P_max) {
      PRIMS[PRESSURE] = P_max;
    }
    // Now compute eps
    PRIMS[EPSILON ] = eps_cold + (PRIMS[PRESSURE]-prs_cold)/(eos.Gamma_th-1.0)/PRIMS[RHOB];
    // If needed, recompute the entropy function
    compute_entropy_function(eos,PRIMS[RHOB],PRIMS[PRESSURE],&PRIMS[ENTROPY]);
  }

  // Tabulated EOS specific floors and ceilings
  else if( eos.is_Tabulated ) {
    // Apply floors and ceilings to Y_e and T
    const CCTK_REAL xye   = MIN(MAX(PRIMS[YEPRIM     ],eos.Ye_min),eos.Ye_atm);
    const CCTK_REAL xtemp = MIN(MAX(PRIMS[TEMPERATURE],eos.T_atm ),eos.T_max );

    // Additional variables used for the EOS call
    const CCTK_REAL xrho  = PRIMS[RHOB];
    CCTK_REAL xprs        = 0.0;
    CCTK_REAL xeps        = 0.0;
    CCTK_REAL xent        = 0.0;
    get_P_eps_and_S_from_rho_Ye_and_T(eos,xrho,xye,xtemp, &xprs,&xeps,&xent);
    // Now update the primitives (rho has already been set)
    PRIMS[YEPRIM     ] = xye;
    PRIMS[TEMPERATURE] = xtemp;
    PRIMS[PRESSURE   ] = xprs;
    PRIMS[EPSILON    ] = xeps;
    PRIMS[ENTROPY    ] = xent;
  }  
  
}
