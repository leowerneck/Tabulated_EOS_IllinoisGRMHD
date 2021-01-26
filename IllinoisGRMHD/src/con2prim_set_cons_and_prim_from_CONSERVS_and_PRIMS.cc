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
                                       CCTK_REAL *restrict cons) {

  // IllinoisGRMHD's variable is the "densitised" version of
  // the standard conservative variables (D,tau,S_{i}). In
  // other words, we have the relationships:
  //
  // rho_star   = alpha*sqrt(gamma) *  D
  // tilde(tau) = alpha*sqrt(gamma) * tau
  // tilde(S)_i = alpha*sqrt(gamma) * S_i
  //
  // Therefore the conversion between the two is straightfoward.
  //
  // One important note, however, is that the Noble2D con2prim
  // routine does not use the variable tau, but instead the
  // energy variable u which is related to IllinoisGRMHD's
  // conservatives via the relation:
  //
  // u = -alpha*tau - (alpha-1)*rho_star + beta^{i}tilde(S)_{i}
  //
  // The magnetic fields in IllinoisGRMHD also need to be
  // rescaled by a factor of sqrt(4pi). In the case of
  // the Noble2D routine, an extra factor of the lapse
  // is also required.
  //
  // Now let us begin the conversion.

  // First compute alpha*sqrt(gamma) = alpha*psi^(6)
  
  // Finally compute the remaining quantities (which are routine specific)
  if( (c2p_key == Noble2D         ) ||
      (c2p_key == Noble1D         ) ||
      (c2p_key == Noble1D_entropy ) ||
      (c2p_key == Noble1D_entropy2) ) {

    // Auxiliary variable (tau to u conversion)
    CCTK_REAL uu = -CONSERVS[TAUENERGY]*METRIC_LAP_PSI4[LAPSE] - (METRIC_LAP_PSI4[LAPSE]-1.0)*CONSERVS[RHOSTAR] +
      METRIC[SHIFTX]*CONSERVS[STILDEX] + METRIC[SHIFTY]*CONSERVS[STILDEY]  + METRIC[SHIFTZ]*CONSERVS[STILDEZ]; // note the minus sign on tau

    // Note that in the c2p routine we used to
    // multiply the cons vector by alpha/detg.
    // But this is equivalent to:
    //
    // alpha/detg = alpha/( alpha*sqrt(gamma) )
    //            = alpha/( alpha*psi^6 )
    //            = psi^{-6}
    const CCTK_REAL psim6 = 1.0/ METRIC_LAP_PSI4[PSI6];

    // Now convert from one system to the next
    cons[RHO   ] = CONSERVS[RHOSTAR];
    cons[UU    ] = (uu - cons[RHO])  * psim6;
    cons[RHO   ] = cons[RHO]         * psim6;
    cons[UTCON1] = CONSERVS[STILDEX] * psim6;
    cons[UTCON2] = CONSERVS[STILDEY] * psim6;
    cons[UTCON3] = CONSERVS[STILDEZ] * psim6;
    cons[BCON1 ] = PRIMS[BX_CENTER ] * ONE_OVER_SQRT_4PI;
    cons[BCON2 ] = PRIMS[BY_CENTER ] * ONE_OVER_SQRT_4PI;
    cons[BCON3 ] = PRIMS[BZ_CENTER ] * ONE_OVER_SQRT_4PI;
    if( (c2p_key == Noble1D_entropy ) || (c2p_key == Noble1D_entropy2) ) {
      // Note that:
      //
      // S_star = alpha * sqrt(gamma) * S * u^{0}
      //        = ( alpha * u^{0} ) * psi^{6} * S
      //        = W_Lorentz * psi^{6} * S
      //    .--------------------------------.
      // => | S_star/psi^{6} = W_Lorentz * S | ,
      //    .--------------------------------.
      // which is the variable expeected by the Noble1D entropy routines.
      cons[WS ]  = CONSERVS[ENTSTAR] * psim6;
    }

  }
  else if( (c2p_key == Palenzuela1D) || (c2p_key == Palenzuela1D_entropy) ) {

    // All that is required in this case is to "undensitize"
    // the conservative variables, i.e. divide by sqrt(-g).

    // Auxiliary variables 1/detg = 1/( alpha*sqrt(gamma) )
    const CCTK_REAL detg     = METRIC_LAP_PSI4[LAPSE]*METRIC_LAP_PSI4[PSI6];
    const CCTK_REAL inv_detg = 1.0/detg;

    // Now convert from one system to the next
    cons[DD    ] = CONSERVS[RHOSTAR  ] * inv_detg;
    cons[S1_cov] = CONSERVS[STILDEX  ] * inv_detg;
    cons[S2_cov] = CONSERVS[STILDEY  ] * inv_detg;
    cons[S3_cov] = CONSERVS[STILDEZ  ] * inv_detg;
    cons[TAU   ] = CONSERVS[TAUENERGY] * inv_detg;
    cons[B1_con] = PRIMS[BX_CENTER   ] * ONE_OVER_SQRT_4PI;
    cons[B2_con] = PRIMS[BY_CENTER   ] * ONE_OVER_SQRT_4PI;
    cons[B3_con] = PRIMS[BZ_CENTER   ] * ONE_OVER_SQRT_4PI;
    cons[YE    ] = CONSERVS[YESTAR   ] * inv_detg;
    if( c2p_key == Palenzuela1D_entropy ) {
      cons[DS  ] = CONSERVS[ENTSTAR  ] * inv_detg;
    }

  }
  else {
    CCTK_VError(VERR_DEF_PARAMS,"Unknown c2p_routine! (Key: %d)",c2p_key);
  }
      
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

  // Initialization is done to silence compiler warnings
  CCTK_INT  polytropic_index = 0;
  CCTK_REAL K_ppoly_tab      = 0.0;
  CCTK_REAL Gamma_ppoly_tab  = 0.0;
  CCTK_REAL u0L              = 1.0;
  CCTK_REAL rho_b_oldL       = PRIMS[RHOB    ];
  CCTK_REAL P_oldL           = PRIMS[PRESSURE];
  CCTK_REAL vxL              = PRIMS[VX      ];
  CCTK_REAL vyL              = PRIMS[VY      ];
  CCTK_REAL vzL              = PRIMS[VZ      ];
  // First the Noble et al. con2prim guesses, which are the more complicated ones
  if( (c2p_key == Noble2D         ) ||
      (c2p_key == Noble1D         ) ||
      (c2p_key == Noble1D_entropy ) ||
      (c2p_key == Noble1D_entropy2) ) {

    if(which_guess==1) {
      //Use a different initial guess:
      rho_b_oldL = CONSERVS[RHOSTAR]/METRIC_LAP_PSI4[PSI6];
      
      /**********************************
       * Piecewise Polytropic EOS Patch *
       *  Finding Gamma_ppoly_tab and K_ppoly_tab *
       **********************************/
      /* Here we use our newly implemented
       * find_polytropic_K_and_Gamma() function
       * to determine the relevant polytropic
       * Gamma and K parameters to be used
       * within this function.
       */
      polytropic_index = find_polytropic_K_and_Gamma_index(eos,rho_b_oldL);
      K_ppoly_tab     = eos.K_ppoly_tab[polytropic_index];
      Gamma_ppoly_tab = eos.Gamma_ppoly_tab[polytropic_index];

      // After that, we compute P_cold
      P_oldL = K_ppoly_tab*pow(rho_b_oldL,Gamma_ppoly_tab);

      u0L = METRIC_LAP_PSI4[LAPSEINV];
      vxL = -METRIC[SHIFTX];
      vyL = -METRIC[SHIFTY];
      vzL = -METRIC[SHIFTZ];
    }

    if(which_guess==2) {
      //Use atmosphere as initial guess:
      rho_b_oldL = 100.0*eos.rho_atm;

      /**********************************
       * Piecewise Polytropic EOS Patch *
       *  Finding Gamma_ppoly_tab and K_ppoly_tab *
       **********************************/
      /* Here we use our newly implemented
       * find_polytropic_K_and_Gamma() function
       * to determine the relevant polytropic
       * Gamma and K parameters to be used
       * within this function.
       */
      polytropic_index = find_polytropic_K_and_Gamma_index(eos,rho_b_oldL);
      K_ppoly_tab     = eos.K_ppoly_tab[polytropic_index];
      Gamma_ppoly_tab = eos.Gamma_ppoly_tab[polytropic_index];

      // After that, we compute P_cold
      P_oldL = K_ppoly_tab*pow(rho_b_oldL,Gamma_ppoly_tab);

      u0L = METRIC_LAP_PSI4[LAPSEINV];
      vxL = -METRIC[SHIFTX];
      vyL = -METRIC[SHIFTY];
      vzL = -METRIC[SHIFTZ];
    }

    CCTK_REAL uL   = P_oldL/(Gamma_ppoly_tab - 1.0);
    CCTK_REAL utxL = u0L*(vxL + METRIC[SHIFTX]);
    CCTK_REAL utyL = u0L*(vyL + METRIC[SHIFTY]);
    CCTK_REAL utzL = u0L*(vzL + METRIC[SHIFTZ]);

    prim[RHO   ]   = rho_b_oldL;
    prim[UU    ]   = uL;
    prim[UTCON1]   = utxL;
    prim[UTCON2]   = utyL;
    prim[UTCON3]   = utzL;
    prim[BCON1 ]   = cons[BCON1];
    prim[BCON2 ]   = cons[BCON2];
    prim[BCON3 ]   = cons[BCON3];

  }
  else if( c2p_key == Palenzuela1D || c2p_key == Palenzuela1D_entropy ) {
    // This one is very simple! The only guess required is the temperature
    prim[TEMP  ] = pow(10.0,(which_guess-1)) * eos.T_atm;
  }
  else {
    CCTK_VError(VERR_DEF_PARAMS,"Unknown c2p_routine! (Key: %d)",c2p_key);
  }
  
}
