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

  // IllinoisGRMHD's variable is the "densitised" version of
  // the standard conservative variables (D,tau,S_{i}). In
  // other words, we have the relationships:
  //
  // rho_star   = sqrt(gamma) *  D
  // tilde(tau) = sqrt(gamma) * tau
  // tilde(S)_i = sqrt(gamma) * S_i
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

  // Note that in the c2p routine we used to
  // multiply the cons vector by alpha/detg.
  // But this is equivalent to:
  //
  // alpha/detg = alpha/( alpha*sqrt(gamma) )
  //            = alpha/( alpha*psi^6 )
  //            = psi^{-6}
  const CCTK_REAL psim6 = 1.0/METRIC_LAP_PSI4[PSI6];

  //----------------------------------------
  //----------- D, S_{i}, B^{i} ------------
  //----------------------------------------
  cons[DD    ] = CONSERVS[RHOSTAR] * psim6;
  cons[S1_cov] = CONSERVS[STILDEX] * psim6;
  cons[S2_cov] = CONSERVS[STILDEY] * psim6;
  cons[S3_cov] = CONSERVS[STILDEZ] * psim6;
  cons[B1_con] = PRIMS[BX_CENTER ] * ONE_OVER_SQRT_4PI;
  cons[B2_con] = PRIMS[BY_CENTER ] * ONE_OVER_SQRT_4PI;
  cons[B3_con] = PRIMS[BZ_CENTER ] * ONE_OVER_SQRT_4PI;

  // First compute the auxiliary variable uu:
  CCTK_REAL uu = -CONSERVS[TAUENERGY]*METRIC_LAP_PSI4[LAPSE] - (METRIC_LAP_PSI4[LAPSE]-1.0)*CONSERVS[RHOSTAR] +
    METRIC[SHIFTX]*CONSERVS[STILDEX] + METRIC[SHIFTY]*CONSERVS[STILDEY]  + METRIC[SHIFTZ]*CONSERVS[STILDEZ]; // note the minus sign on tau
  // Compute the conserv needed by the Noble routines
  cons[UU] = (uu - CONSERVS[RHOSTAR]) * psim6;
  // Then set tau
  cons[TAU] = CONSERVS[TAUENERGY] * psim6;

  //----------------------------------------

  //----------------------------------------
  //----------- Electron fraction ----------
  //----------------------------------------
  if( eos.is_Tabulated ) {
    cons[YE ] = CONSERVS[YESTAR   ] * psim6;
  }
  //----------------------------------------

  //----------------------------------------
  //--------------- Entropy ----------------
  //----------------------------------------
  if( (c2p_key == Noble1D_entropy ) ||
      (c2p_key == Noble1D_entropy2) ||
      (c2p_key == Palenzuela1D    ) ) {

    // The entropy variable is given by
    //
    // S_star / psi^{6} = alpha * psi^{6} * S * u^{0} / psi^{6}
    //                  = ( alpha * u^{0} ) * S
    //                  = W * S
    cons[WS] = CONSERVS[ENTSTAR] * psim6;
  }
  //---------------------------------------

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
      (c2p_key == Noble1D_entropy2) ||
      (c2p_key == CerdaDuran2D    ) ) {

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

    prim[RHO     ] = rho_b_oldL;
    prim[UU      ] = uL;
    prim[UTCON1  ] = utxL;
    prim[UTCON2  ] = utyL;
    prim[UTCON3  ] = utzL;
    prim[BCON1   ] = cons[BCON1];
    prim[BCON2   ] = cons[BCON2];
    prim[BCON3   ] = cons[BCON3];
    prim[WLORENTZ] = METRIC_LAP_PSI4[LAPSE] * u0L;

  }
  
  if( eos.is_Tabulated ) {

    for(int i=0;i<numprims;i++) prim[i] = 1.0/0.0;

    // This one is very simple! The only guess required is the temperature
    if( which_guess == 1 ) {
      prim[TEMP  ] = eos.T_max;
    }
    else {
      prim[TEMP  ] = eos.T_atm;
    }
    prim[YE      ] = CONSERVS[YESTAR]/CONSERVS[RHOSTAR];
  }
  
}
