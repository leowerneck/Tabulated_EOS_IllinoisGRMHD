// Thorn      : IllinoisGRMHD
// File       : con2prim_set_cons_and_prim_from_CONSERVS_and_PRIMS.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This provides functions which 1. convert IllinoisGRMHD's set
//              of conservative variables into the appropriate variables
//              required by the C2P routines and 2. set appropriate primitive
//              guesses.

#include "cctk.h"

#include "IllinoisGRMHD_headers.h"
#include "con2prim_headers.h"

void set_cons_from_PRIMS_and_CONSERVS( const igm_eos_parameters eos,
                                       const CCTK_INT c2p_key,
                                       const CCTK_REAL *restrict METRIC,
                                       const CCTK_REAL *restrict METRIC_LAP_PSI4,
                                       const CCTK_REAL *restrict PRIMS,
                                       const CCTK_REAL *restrict CONSERVS,
                                       CCTK_REAL *restrict cons ) {

  const CCTK_REAL psim6 = 1.0/METRIC_LAP_PSI4[PSI6];

  cons[DD    ] = CONSERVS[RHOSTAR] * psim6;
  cons[S1_cov] = CONSERVS[STILDEX] * psim6;
  cons[S2_cov] = CONSERVS[STILDEY] * psim6;
  cons[S3_cov] = CONSERVS[STILDEZ] * psim6;
  cons[B1_con] = PRIMS[BX_CENTER ] * ONE_OVER_SQRT_4PI;
  cons[B2_con] = PRIMS[BY_CENTER ] * ONE_OVER_SQRT_4PI;
  cons[B3_con] = PRIMS[BZ_CENTER ] * ONE_OVER_SQRT_4PI;

  CCTK_REAL uu = -CONSERVS[TAUENERGY]*METRIC_LAP_PSI4[LAPSE] - (METRIC_LAP_PSI4[LAPSE]-1.0)*CONSERVS[RHOSTAR] +
    METRIC[SHIFTX]*CONSERVS[STILDEX] + METRIC[SHIFTY]*CONSERVS[STILDEY]  + METRIC[SHIFTZ]*CONSERVS[STILDEZ]; // note the minus sign on tau
  cons[UU] = (uu - CONSERVS[RHOSTAR]) * psim6;
  cons[TAU] = MAX(CONSERVS[TAUENERGY] * psim6,eos.tau_atm);
  cons[YE ] = CONSERVS[YESTAR   ] * psim6;
  cons[WS] = CONSERVS[ENTSTAR] * psim6;
}

void set_prim_from_PRIMS_and_CONSERVS( const igm_eos_parameters eos,
                                       const CCTK_INT c2p_key,
                                       const CCTK_INT which_guess,
                                       const CCTK_REAL *restrict METRIC,
                                       const CCTK_REAL *restrict METRIC_LAP_PSI4,
                                       const CCTK_REAL *restrict PRIMS,
                                       const CCTK_REAL *restrict CONSERVS,
                                       const CCTK_REAL *restrict cons,
                                       CCTK_REAL *restrict prim ) {

  if( which_guess == 1 ) {
    prim[TEMP  ] = eos.T_atm;
  }
  else {
    prim[TEMP  ] = eos.T_max;
  }
  prim[YE      ] = CONSERVS[YESTAR]/CONSERVS[RHOSTAR];
}
