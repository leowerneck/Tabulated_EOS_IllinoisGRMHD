// Thorn      : IllinoisGRMHD
// File       : EOS_Hybrid.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: In this file we provide piecewise polytropic-based
//              hybrid EOS functions.

#include "cctk.h"
#include "cctk_Parameters.h"

#include "IllinoisGRMHD_headers.h"

/* Function    : initialize_igm_eos_parameters_from_input()
 * Authors     : Leo Werneck
 * Description : Initialize the eos struct from user
 *               input
 * Dependencies: setup_K_ppoly_tab__and__eps_integ_consts()
 *             : cctk_parameters.h (FIXME)
 *
 * Inputs      : eos             - a struct containing the following
 *                                 relevant quantities:
 *             : neos            - number of polytropic EOSs used
 *             : rho_ppoly_tab   - array of rho values that determine
 *                                 the polytropic EOS to be used.
 *             : Gamma_ppoly_tab - array of Gamma_cold values to be
 *                                 used in each polytropic EOS.
 *             : K_ppoly_tab     - array of K_ppoly_tab values to be used
 *                                 in each polytropic EOS.
 *             : eps_integ_const - array of C_{j} values, which are the
 *                                 integration constants that arrise when
 *                                 determining eps_{cold} for a piecewise
 *                                 polytropic EOS.
 *
 * Outputs     : eos             - fully initialized EOS struct
 */
void initialize_Hybrid_EOS_parameters_from_input(igm_eos_parameters &eos) {

  // Initialize: neos, {rho_{j}}, {K_{0}}, and {Gamma_{j}}
#ifndef ENABLE_STANDALONE_IGM_C2P_SOLVER
  DECLARE_CCTK_PARAMETERS;
#endif

  // In the Hybrid EOS case we allow reconstruction
  // of either the pressure or the entropy.
  if( CCTK_EQUALS(igm_PPM_reconstructed_variable,"pressure") ) {
    eos.PPM_reconstructed_var = PRESSURE;
  }
  else if( CCTK_EQUALS(igm_PPM_reconstructed_variable,"entropy") ) {
    eos.PPM_reconstructed_var = ENTROPY;
  }
  else {
    CCTK_VError(VERR_DEF_PARAMS,
                "PPM reconstruction of variable \"%s\" not supported with Hybrid EOS. "
                "Can only reconstruct: pressure, entropy. ABORTING!",
                igm_PPM_reconstructed_variable);
  }

  eos.neos     = neos;
  eos.Gamma_th = Gamma_th;
  eos.K_ppoly_tab[0] = K_ppoly_tab0;
  for(int j=0; j<=neos-2; j++) eos.rho_ppoly_tab[j]   = rho_ppoly_tab_in[j];
  for(int j=0; j<=neos-1; j++) eos.Gamma_ppoly_tab[j] = Gamma_ppoly_tab_in[j];

  // Initialize {K_{j}}, j>=1, and {eps_integ_const_{j}}
  setup_K_ppoly_tab__and__eps_integ_consts(eos);

  // --------- Atmospheric values ---------
  // Compute P and eps atm
  CCTK_REAL P_atm,eps_atm;
  compute_P_cold__eps_cold(eos,rho_b_atm,P_atm,eps_atm);
  // Set atmospheric values
  eos.rho_atm = rho_b_atm;
  eos.P_atm   = P_atm;
  eos.eps_atm = eps_atm;
  eos.tau_atm = tau_atm;
  // --------------------------------------

  // -------------- Ceilings --------------
  // Compute maximum P and eps
  CCTK_REAL P_max,eps_max;
  compute_P_cold__eps_cold(eos,rho_b_max,P_max,eps_max);
  // Set maximum values
  eos.rho_max = rho_b_max;
  eos.P_max   = P_max;
  eos.eps_max = eps_max;
  eos.W_max   = GAMMA_SPEED_LIMIT;
  // --------------------------------------

  // --------------- Floors ---------------
  // We'll choose these as the atmospheric values
  eos.rho_min = eos.rho_atm;
  eos.P_min   = eos.P_atm;
  eos.eps_min = eos.eps_atm;
  // --------------------------------------
  
}

/* Function    : setup_K_ppoly_tab__and__eps_integ_consts()
 * Authors     : Leo Werneck
 * Description : For a given set of EOS inputs, determine
 *               values of K_ppoly_tab that will result in a
 *               everywhere continuous P_cold function.
 * Dependencies: none
 *
 * Inputs      : eos             - a struct containing the following
 *                                 relevant quantities:
 *             : neos            - number of polytropic EOSs used
 *             : rho_ppoly_tab   - array of rho values that determine
 *                                 the polytropic EOS to be used.
 *             : Gamma_ppoly_tab - array of Gamma_cold values to be
 *                                 used in each polytropic EOS.
 *             : K_ppoly_tab     - array of K_ppoly_tab values to be used
 *                                 in each polytropic EOS. Only K_ppoly_tab[0]
 *                                 is known prior to the function call
 *             : eps_integ_const - array of C_{j} values, which are the
 *                                 integration constants that arrise when
 *                                 determining eps_{cold} for a piecewise
 *                                 polytropic EOS. This array should be
 *                                 uninitialized or contain absurd values
 *                                 prior to this function call.
 *
 * Outputs     : K_ppoly_tab     - fully populated array of K_ppoly_tab
 *                                 to be used in each polytropic EOS.
 *             : eps_integ_const - fully populated array of C_{j}'s,
 *                                 used to compute eps_cold for
 *                                 a piecewise polytropic EOS.
 */
void setup_K_ppoly_tab__and__eps_integ_consts(igm_eos_parameters &eos) {

  /* When neos = 1, we will only need the value K_ppoly_tab[0] and eps_integ_const[0].
   * Since our only polytropic EOS is given by
   *  -----------------------------------
   * | P_{0} = K_{0} * rho ^ (Gamma_{0}) | ,
   *  -----------------------------------
   * and, therefore,
   *  ---------------------------------------------------------------
   * | eps_{0} = K_{0} * rho ^ (Gamma_{0}-1) / (Gamma_{0}-1) + C_{0} | ,
   *  ---------------------------------------------------------------
   * we only need to set up K_{0} := K_ppoly_tab[0] and C_{0} := eps_integ_const[0].
   * K_{0} is a user input, so we need to do nothing. C_{0}, on the other hand,
   * is fixed by demanding that eps(rho) -> 0 as rho -> 0. Thus, C_{0} = 0.
   */
  eos.eps_integ_const[0] = 0.0;
  if(eos.neos==1) return;

  /********************
   * Setting up K_{j} *
   ********************/
  /* When neos > 1, we have the following structure
   *
   *           /      K_{0} * rho^(Gamma_{0})      ,                 rho <  rho_{0}
   *           |      K_{1} * rho^(Gamma_{1})      ,      rho_{0} <= rho <  rho_{1}
   *           |              ...
   * P(rho) = <       K_{j} * rho^(Gamma_{j})      ,    rho_{j-1} <= rho <  rho_{j}
   *           |              ...
   *           | K_{neos-2} * rho^(Gamma_{neos-2}) , rho_{neos-3} <= rho <  rho_{neos-2}
   *           \ K_{neos-1} * rho^(Gamma_{neos-1}) ,                 rho >= rho_{neos-2}
   *
   * Imposing that P(rho) be everywhere continuous, we have
   *  -------------------------------------------------
   * | K_{j} = K_{j-1} * rho^(Gamma_{j-1} - Gamma_{j}) |
   *  -------------------------------------------------
   */
  for(int j=1; j<eos.neos; j++) {
    // Set a useful auxiliary variable to keep things more compact:
    // First, (Gamma_{j-1} - Gamma_{j}):
    CCTK_REAL Gamma_diff = eos.Gamma_ppoly_tab[j-1] - eos.Gamma_ppoly_tab[j];

    // Implement the boxed equation above, using our auxiliary variable:
    eos.K_ppoly_tab[j] = eos.K_ppoly_tab[j-1] * pow(eos.rho_ppoly_tab[j-1],Gamma_diff);
  }


  /********************
   * Setting up C_{j} *
   ********************/
  /* When neos > 1, we have the following structure (let neos->N):
   *
   *             /      K_{0}*rho^(Gamma_{0}-1)/(Gamma_{0}-1)  + C_{0},                rho <  rho_{0}
   *             |      K_{1}*rho^(Gamma_{1}-1)/(Gamma_{1}-1)  + C_{1},     rho_{0} <= rho <  rho_{1}
   *             |                     ...
   * eps(rho) = <       K_{j}*rho^(Gamma_{j}-1)/(Gamma_{j}-1)  + C_{j},   rho_{j-1} <= rho <  rho_{j}
   *             |                     ...
   *             | K_{N-2}*rho^(Gamma_{N-2}-1)/(Gamma_{N-2}-1) + C_{N-2}, rho_{N-3} <= rho <  rho_{N-2}
   *             \ K_{N-1}*rho^(Gamma_{N-1}-1)/(Gamma_{N-1}-1) + C_{N-1},              rho >= rho_{N-2}
   *
   * Imposing that eps_{cold}(rho) be everywhere continuous, we have
   *  ---------------------------------------------------------------
   * | C_{j} = C_{j-1}                                               |
   * |       + ( K_{j-1}*rho_{j-1}^(Gamma_{j-1}-1) )/(Gamma_{j-1}-1) |
   * |       - ( K_{j+0}*rho_{j-1}^(Gamma_{j+0}-1) )/(Gamma_{j+0}-1) |
   *  ---------------------------------------------------------------
   */
  for(int j=1; j<eos.neos; j++) {
    // Set a few useful auxiliary variables to keep things more compact:
    // First, (Gamma_{j-1}-1):
    CCTK_REAL Gammajm1m1 = eos.Gamma_ppoly_tab[j-1] - 1.0;

    // Then, (Gamma_{j+0}-1):
    CCTK_REAL Gammajp0m1 = eos.Gamma_ppoly_tab[j+0] - 1.0;

    // Next, ( K_{j-1}*rho_{j-1}^(Gamma_{j-1}-1) )/(Gamma_{j-1}-1):
    CCTK_REAL aux_epsm1  = eos.K_ppoly_tab[j-1]*pow(eos.rho_ppoly_tab[j-1],Gammajm1m1)/Gammajm1m1;

    // Finally, ( K_{j+0}*rho_{j+0}^(Gamma_{j+0}-1) )/(Gamma_{j+0}-1):
    CCTK_REAL aux_epsp0  = eos.K_ppoly_tab[j+0]*pow(eos.rho_ppoly_tab[j-1],Gammajp0m1)/Gammajp0m1;

    // Implement the boxed equation above, using our auxiliary variables:
    eos.eps_integ_const[j] = eos.eps_integ_const[j-1] + aux_epsm1 - aux_epsp0;
  }
}

/* Function    : find_polytropic_K_and_Gamma_index()
 * Authors     : Leo Werneck & Zach Etienne
 * Description : For a given value of rho, find the
 *               appropriate values of Gamma_ppoly_tab
 *               and K_ppoly_tab by determining the appropriate
 *               index
 * Dependencies: initialize_igm_eos_parameters_from_input()
 *             : cctk_parameters.h (FIXME)
 *
 * Inputs      : eos           - a struct containing the following
 *                               relevant quantities:
 *             : neos          - number of polytropic EOSs used
 *             : rho_ppoly_tab - array of rho values that determine
 *
 * Outputs     : index         - the appropriate index for the K_ppoly_tab
 *                               and Gamma_ppoly_tab array
 */
int find_polytropic_K_and_Gamma_index(igm_eos_parameters eos, CCTK_REAL rho_in) {

  /* We want to find the appropriate polytropic EOS for the
   * input value rho_in. Remember that:
   *
   * if rho < rho_{0}:                P_{0} , index: 0
   * if rho >= rho_{0} but < rho_{1}: P_{1} , index: 1
   * if rho >= rho_{1} but < rho_{2}: P_{2} , index: 2
   *                      ...
   * if rho >= rho_{j-1} but < rho_{j}: P_{j} , index: j
   *
   * Then, a simple way of determining the index is through
   * the formula:
   *  ---------------------------------------------------------------------------
   * | index = (rho >= rho_{0}) + (rho >= rho_{1}) + ... + (rho >= rho_{neos-2}) |
   *  ---------------------------------------------------------------------------
   */
  if(eos.neos == 1) return 0;

  int polytropic_index = 0;
  for(int j=0; j<=eos.neos-2; j++) polytropic_index += (rho_in >= eos.rho_ppoly_tab[j]);

  return polytropic_index;
}

/* Function    : find_polytropic_K_and_Gamma_index_from_P()
 * Authors     : Leo Werneck
 * Description : For a given value of , find the
 *               appropriate values of Gamma_ppoly_tab
 *               and K_ppoly_tab by determining the appropriate
 *               index
 * Dependencies: initialize_igm_eos_parameters_from_input()
 *             : cctk_parameters.h (FIXME)
 *
 * Inputs      : eos               - a struct containing the following
 *                                   relevant quantities:
 *                 : neos          - number of polytropic EOSs used
 *                 : rho_ppoly_tab - array of rho values that determine
 *             : P_in              - Input pressure
 *
 * Outputs     : index             - the appropriate index for the K_ppoly_tab
 *                                   and Gamma_ppoly_tab array
 */
int find_polytropic_K_and_Gamma_index_from_P(const igm_eos_parameters eos, const CCTK_REAL P_in) {

  if(eos.neos == 1) return 0;

  int polytropic_index = 0;
  for(int j=0; j<=eos.neos-2; j++) {
    // This function is just slightly more involved than the previous function.
    // Instead of comparing rho_in against the values of rho that separate
    // one polytropic index from the next, we compute P(rho) at the separation
    // values and compare P_in against them. This is guaranteed to work because
    // P(rho) increases monotonically with rho.
    const CCTK_REAL P_local = eos.K_ppoly_tab[j] * pow( eos.rho_ppoly_tab[j], eos.Gamma_ppoly_tab[j] );
    polytropic_index += (P_in >= P_local);
  }

  return polytropic_index;
}

/* Function    : compute_P_cold__eps_cold()
 * Authors     : Leo Werneck
 * Description : Computes P_cold and eps_cold.
 * Dependencies: initialize_igm_eos_parameters_from_input()
 *             : find_polytropic_K_and_Gamma_index()
 *
 * Inputs      : P_cold           - cold pressure
 *             : eps_cold         - cold specific internal energy
 *             : eos              - a struct containing the following
 *                                  relevant quantities:
 *             : neos             - number of polytropic EOSs used
 *             : rho_ppoly_tab         - array of rho values that determine
 *                                  the polytropic EOS to be used.
 *             : Gamma_ppoly_tab       - array of Gamma_cold values to be
 *                                  used in each polytropic EOS.
 *             : K_ppoly_tab           - array of K_ppoly_tab values to be used
 *                                  in each polytropic EOS.
 *             : eps_integ_const  - array of C_{j} values, which are the
 *                                  integration constants that arrise when
 *                                  determining eps_{cold} for a piecewise
 *                                  polytropic EOS.
 *
 * Outputs     : P_cold           - cold pressure (supports SPEOS and PPEOS)
 *             : eps_cold         - cold specific internal energy (supports SPEOS and PPEOS)
 *             : polytropic_index - polytropic index used for P_cold and eps_cold
 *
 *             SPEOS: Single-Polytrope Equation of State
 *             PPEOS: Piecewise Polytrope Equation of State
 */
void compute_P_cold__eps_cold(igm_eos_parameters eos, CCTK_REAL rho_in,
                              CCTK_REAL &P_cold,CCTK_REAL &eps_cold) {
  // This code handles equations of state of the form defined
  // in Eqs 13-16 in http://arxiv.org/pdf/0802.0200.pdf
  if(rho_in==0) {
    P_cold   = 0.0;
    eps_cold = 0.0;
    return;
  }

  /*  --------------------------------------------------
   * | Single and Piecewise Polytropic EOS modification |
   *  --------------------------------------------------
   *
   * We now begin our modifications to this function so that
   * it supports both single and piecewise polytropic equations
   * of state.
   *
   * The modifications below currently assume that the user
   * has called the recently added function
   *
   * - initialize_igm_eos_parameters_from_input()
   *
   * *before* this function is called. We can add some feature
   * to check this automatically as well, but we'll keep that as
   * a TODO/FIXME for now.
   */

  /* First, we compute the pressure, which in the case of a
   * piecewise polytropic EOS is given by
   *
   *           /   P_{1}      /    K_{1} * rho^(Gamma_{1})       ,      rho_{0} <= rho < rho_{1}
   *           |    ...       |            ...
   * P(rho) = <    P_{j}   = <     K_{j} * rho^(Gamma_{j})       ,    rho_{j-1} <= rho < rho_{j}
   *           |    ...       |            ...
   *           \ P_{neos-2}   \ K_{neos-2} * rho^(Gamma_{neos-2}), rho_{neos-3} <= rho < rho_{neos-2}
   *
   * The index j is determined by the find_polytropic_K_and_Gamma_index() function.
   */
  // Set up useful auxiliary variables
  int polytropic_index      = find_polytropic_K_and_Gamma_index(eos,rho_in);
  CCTK_REAL K_ppoly_tab     = eos.K_ppoly_tab[polytropic_index];
  CCTK_REAL Gamma_ppoly_tab = eos.Gamma_ppoly_tab[polytropic_index];
  CCTK_REAL eps_integ_const = eos.eps_integ_const[polytropic_index];

  // Then compute P_{cold}
  P_cold = K_ppoly_tab*pow(rho_in,Gamma_ppoly_tab);

  /* Then we compute the cold component of the specific internal energy,
   * which in the case of a piecewise polytropic EOS is given by (neos -> N)
   *
   *             /   P_{1}/(rho*(Gamma_{1}-1))   + C_{1}  ,   rho_{0} <= rho < rho_{1}
   *             |                     ...
   * eps(rho) = <    P_{j}/(rho*(Gamma_{j}-1))   + C_{j}  , rho_{j-1} <= rho < rho_{j}
   *             |                     ...
   *             \ P_{N-2}/(rho*(Gamma_{N-2}-1)) + C_{N-2}, rho_{N-3} <= rho < rho_{N-2}
   */
  eps_cold = P_cold/(rho_in*(Gamma_ppoly_tab-1.0)) + eps_integ_const;

}


/* Function    : compute_P_cold__eps_cold__dPcold_drho__eps_th__h__Gamma_cold()
 * Authors     : Leo Werneck & Zach Etienne
 * Description : Compute basic quantities related to
 *             : the EOS, namely: P_cold, eps_cold,
 *             : dPcold/drho, eps_th, h, and Gamma_cold
 * Dependencies: initialize_igm_eos_parameters_from_input()
 *
 * Inputs      : U               - array containing primitives {rho,P,v^{i},B^i}
 *             : P_cold          - cold pressure
 *             : eps_cold        - cold specific internal energy
 *             : dPcold_drho     - derivative of P_cold
 *             : eps_th          - thermal specific internal energy
 *             : h               - enthalpy
 *             : Gamma_cold      - cold polytropic Gamma
 *             : eos             - a struct containing the following
 *                                 relevant quantities:
 *             : neos            - number of polytropic EOSs used
 *             : rho_ppoly_tab   - array of rho values that determine
 *                                 the polytropic EOS to be used.
 *             : Gamma_ppoly_tab - array of Gamma_cold values to be
 *                                 used in each polytropic EOS.
 *             : K_ppoly_tab     - array of K_ppoly_tab values to be used
 *                                 in each polytropic EOS.
 *             : eps_integ_const - array of C_{j} values, which are the
 *                                 integration constants that arrise when
 *                                 determining eps_{cold} for a piecewise
 *                                 polytropic EOS.
 *
 * Outputs     : P_cold          - cold pressure (supports SPEOS and PPEOS)
 *             : eps_cold        - cold specific internal energy (supports SPEOS and PPEOS)
 *             : dPcold_drho     - derivative of P_cold (supports SPEOS and PPEOS)
 *             : eps_th          - thermal specific internal energy (supports SPEOS and PPEOS)
 *             : h               - enthalpy (supports SPEOS and PPEOS)
 *             : Gamma_cold      - cold polytropic Gamma (supports SPEOS and PPEOS)
 *
 *             SPEOS: Single-Polytrope Equation of State
 *             PPEOS: Piecewise Polytrope Equation of State
 */
void compute_P_cold__eps_cold__dPcold_drho__eps_th__h__Gamma_cold(CCTK_REAL *U, const igm_eos_parameters eos,
                                                                  CCTK_REAL &P_cold,CCTK_REAL &eps_cold,CCTK_REAL &dPcold_drho,CCTK_REAL &eps_th,CCTK_REAL &h,
                                                                  CCTK_REAL &Gamma_cold) {
  // This code handles equations of state of the form defined
  // in Eqs 13-16 in http://arxiv.org/pdf/0802.0200.pdf

  if(U[RHOB]==0) {
    P_cold      = 0.0;
    eps_cold    = 0.0;
    dPcold_drho = 0.0;
    eps_th      = 0.0;
    h           = 0.0;
    Gamma_cold  = eos.Gamma_ppoly_tab[0];
    return;
  }

  int polytropic_index = find_polytropic_K_and_Gamma_index(eos, U[RHOB]);
  compute_P_cold__eps_cold(eos,U[RHOB], P_cold,eps_cold);
  CCTK_REAL Gamma_ppoly_tab = eos.Gamma_ppoly_tab[polytropic_index];

  // Set auxiliary variable rho_b^{-1}
  CCTK_REAL U_RHOB_inv = 1.0/U[RHOB];

  // Next compute dP/drho = Gamma * P / rho
  dPcold_drho = Gamma_ppoly_tab*P_cold*U_RHOB_inv;

  // Then we compute eps_th, h, and set Gamma_cold = Gamma_ppoly_tab[j].
  eps_th = (U[PRESSURE] - P_cold)/(eos.Gamma_th-1.0)*U_RHOB_inv;
  h = 1.0 + eps_cold + eps_th + U[PRESSURE]*U_RHOB_inv;
  Gamma_cold = Gamma_ppoly_tab;

}


/* Function    : print_EOS_table()
 * Authors     : Leo Werneck
 * Description : Prints out the EOS table, for diagnostic purposes
 *
 * Dependencies: initialize_igm_eos_parameters_from_input()
 *
 * Inputs      : eos             - a struct containing the following
 *                                 relevant quantities:
 *             : neos            - number of polytropic EOSs used
 *             : rho_ppoly_tab   - array of rho values that determine
 *                                 the polytropic EOS to be used.
 *             : Gamma_ppoly_tab - array of Gamma_cold values to be
 *                                 used in each polytropic EOS.
 *
 * Outputs     : CCTK_VInfo string with the EOS table used by IllinoisGRMHD
 */
void print_EOS_Hybrid( igm_eos_parameters eos ) {

  /* Start by printint a header t the table */
#ifndef ENABLE_STANDALONE_IGM_C2P_SOLVER
  CCTK_VInfo(CCTK_THORNSTRING,
             "\n.--------------------------------------------.\n");
#else
  printf(".--------------------------------------------.\n");
#endif
         
  printf("|             Hybrid EOS Details             |\n"
         ".--------------------------------------------.\n"
         "|              rho_ppoly_tab[j]              |\n"
         ".--------------------------------------------.\n");

  /* Adjust the maximum index of rhob to
   * allow for single polytropes as well
   */
  int max_rho_index;
  if( eos.neos==1 ) {
    max_rho_index = 0;
  }
  else {
    max_rho_index = eos.neos-2;
  }

  /* Print out rho_pppoly_tab */
  for(int jj=0; jj<=max_rho_index; jj++) {
    printf("|  rho_ppoly_tab[%d] = %.15e  |\n",jj,eos.rho_ppoly_tab[jj]);
    if(jj == eos.neos-2) {
      printf(".--------------------------------------------.\n");
    }
  }

  /* Print out Gamma_ppoly_tab */
  printf("|             Gamma_ppoly_tab[j]             |\n"
         ".--------------------------------------------.\n");
  for(int jj=0; jj<=eos.neos-1; jj++) {
    printf("| Gamma_ppoly_tab[%d] = %.15e |\n",jj,eos.Gamma_ppoly_tab[jj]);
    if(jj == eos.neos-1) {
      printf(".--------------------------------------------.\n");
    }
  }

  /* Print out K_ppoly_tab */
  printf("|               K_ppoly_tab[j]               |\n"
         ".--------------------------------------------.\n");
  for(int jj=0; jj<=eos.neos-1; jj++) {
    printf("|   K_ppoly_tab[%d] = %.15e   |\n",jj,eos.K_ppoly_tab[jj]);
    if(jj == eos.neos-1) {
      printf(".--------------------------------------------.\n\n");
    }
  }
}

void compute_entropy_function( const igm_eos_parameters eos,
                               const CCTK_REAL rho,
                               const CCTK_REAL P,
                               CCTK_REAL *restrict S ) {
  // This function sets the entropy funtion:
  //
  // S = P / rho^(Gamma-1)
  //
  // See eq. (20) in https://arxiv.org/pdf/0808.3140.pdf
  const CCTK_INT  index = find_polytropic_K_and_Gamma_index(eos,rho);
  const CCTK_REAL Gamma = eos.Gamma_ppoly_tab[index];
  // Now compute S
  *S = P / pow(rho,Gamma-1.0);
}
