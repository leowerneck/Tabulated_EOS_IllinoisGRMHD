// Thorn      : IllinoisGRMHD
// File       : EOS_Tabulated.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: In this file we provide tabulated EOS functions.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "IllinoisGRMHD_headers.h"

namespace nuc_eos_private {
  extern int nrho,nye,ntemp;
  extern double * restrict alltables;
}

//--------------------------------------

// EOS_Omni does not provide functions to obtain the
// maximum table value of e.g. the specific internal
// energy. This function does that.
CCTK_REAL get_EOS_table_max( const int which_var ) {

  // It's simply too annoying to keep appending the
  // name to all the variables. This namespace is
  // also very small, so we should be okay.
  using namespace nuc_eos_private;

  // Loop over the table, searching for the maximum value
  CCTK_INT  totalsize     = nrho * nye * ntemp;
  CCTK_REAL var_max_value = alltables[which_var];

  for(int i=0;i<totalsize;i++) {
    CCTK_REAL var_aux = alltables[which_var + NTABLES*i];
    if( var_aux > var_max_value ) var_max_value = var_aux;
  }
  return var_max_value;
  
}

// EOS_Omni does not provide functions to obtain the
// maximum table value of e.g. the specific internal
// energy. This function does that.
CCTK_REAL get_EOS_table_min( const int which_var ) {

  // It's simply too annoying to keep appending the
  // name to all the variables. This namespace is
  // also very small, so we should be okay.
  using namespace nuc_eos_private;

  // Loop over the table, searching for the maximum value
  CCTK_INT  totalsize     = nrho * nye * ntemp;
  CCTK_REAL var_min_value = alltables[which_var];

  for(int i=0;i<totalsize;i++) {
    CCTK_REAL var_aux = alltables[which_var + NTABLES*i];
    if( var_aux < var_min_value ) var_min_value = var_aux;
  }
  return var_min_value;
  
}

void initialize_Tabulated_EOS_parameters_from_input( const CCTK_REAL cctk_time,igm_eos_parameters& eos ) {

  DECLARE_CCTK_PARAMETERS;

  // It's simply too annoying to keep appending the
  // name to all the variables. This namespace is
  // also very small, so we should be okay.
  using namespace nuc_eos;

  // Which variable do we want to reconstruct during PPM?
  // The tabulated EOS case actually allows for all the
  // possibilities.
  if( CCTK_EQUALS(igm_PPM_reconstructed_variable,"pressure") ) {
    eos.PPM_reconstructed_var = PRESSURE;
  }
  else if( CCTK_EQUALS(igm_PPM_reconstructed_variable,"epsilon") ) {
    eos.PPM_reconstructed_var = EPSILON;
  }
  else if( CCTK_EQUALS(igm_PPM_reconstructed_variable,"entropy") ) {
    eos.PPM_reconstructed_var = ENTROPY;
  }
  else {
    CCTK_VError(VERR_DEF_PARAMS,
                "PPM reconstruction of variable \"%s\" not supported with Tabulated EOS. "
                "Can only reconstruct: pressure, epsilon, entropy. ABORTING!",
                igm_PPM_reconstructed_variable);
  }

  // Initialize tabulated EOS parameters

  // Check if it's time to begin temperature evolution
  if( igm_evolve_temperature && (cctk_time >= igm_freeze_T_evolution_until_cctk_time) ) {
    eos.evolve_T = true;
  }
  else {
    eos.evolve_T = false;
  }

  // Root-finding precision (for table inversions)
  eos.root_finding_precision = igm_eos_root_finding_precision;

  // --------- Atmospheric values ---------
  // Atmospheric rho
  eos.rho_atm = rho_b_atm;
  // Atmospheric electron fraction
  eos.Ye_atm  = igm_Ye_atm;
  // Atmospheric temperature fraction
  eos.T_atm   = igm_T_atm;
  // Compute P, eps, and S in the atmosphere
  if( eos.evolve_entropy ) {
    WVU_EOS_P_eps_and_S_from_rho_Ye_T( eos.rho_atm,eos.Ye_atm,eos.T_atm,
                                       &eos.P_atm,&eos.eps_atm,&eos.S_atm );
  }
  else {
    WVU_EOS_P_and_eps_from_rho_Ye_T( eos.rho_atm,eos.Ye_atm,eos.T_atm,
                                     &eos.P_atm,&eos.eps_atm );
  }
  // Atmospheric tau
  eos.tau_atm = 1.1*tau_atm;
  // --------------------------------------

  // -------------- Ceilings --------------
  // Get the maximum pressure. Remember that the alltables
  // array actually constains ln(press), so we must adjust
  // appropriately.
  const CCTK_REAL eos_prsmax = exp(get_EOS_table_max( table_key_pressure ));
  // Then get the maximum value of eps. Remember that the 
  // alltables array actually contains ln(eps + energy_shift),
  // so we must adjust appropriately.
  const CCTK_REAL eos_epsmax = exp(get_EOS_table_max( table_key_epsilon )) - energy_shift;
  // Finally, get the maximum entropy
  const CCTK_REAL eos_entmax = get_EOS_table_max( table_key_entropy );
  // Now set the EOS struct variables
  eos.rho_max = eos_rhomax  * igm_eos_table_ceiling_safety_factor;
  eos.Ye_max  = eos_yemax   * igm_eos_table_ceiling_safety_factor;
  eos.T_max   = eos_tempmax * igm_eos_table_ceiling_safety_factor;
  eos.P_max   = eos_prsmax  * igm_eos_table_ceiling_safety_factor;
  eos.eps_max = eos_epsmax  * igm_eos_table_ceiling_safety_factor;
  eos.S_max   = eos_entmax  * igm_eos_table_ceiling_safety_factor;
  // --------------------------------------

  // --------------- Floors ---------------
  // Get the miniimum pressure. Remember that the alltables
  // array actually constains ln(press), so we must adjust
  // appropriately.
  const CCTK_REAL eos_prsmin = exp(get_EOS_table_min( table_key_pressure ));
  // Then get the minimum value of eps. Remember that the 
  // alltables array actually contains ln(eps + energy_shift),
  // so we must adjust appropriately.
  const CCTK_REAL eos_epsmin = exp(get_EOS_table_min( table_key_epsilon )) - energy_shift;
  // Finally, get the minimum entropy
  const CCTK_REAL eos_entmin = get_EOS_table_min( table_key_entropy );
  // Now set the EOS struct variables
  eos.rho_min = eos_rhomin  * igm_eos_table_floor_safety_factor;
  eos.Ye_min  = eos_yemin   * igm_eos_table_floor_safety_factor;
  eos.T_min   = eos_tempmin * igm_eos_table_floor_safety_factor;
  eos.P_min   = eos_prsmin  * igm_eos_table_floor_safety_factor;
  eos.eps_min = eos_epsmin  * igm_eos_table_floor_safety_factor;
  eos.S_min   = eos_entmin  * igm_eos_table_floor_safety_factor;
  // --------------------------------------

  // ----- con2prim threshold values ------
  eos.depsdT_threshold = palenzuela_depsdT_threshold;
  // --------------------------------------

  // All done!

}

void compute_remaining_prims_on_right_and_left_face( const igm_eos_parameters eos,
                                                     const cGH *restrict cctkGH,
                                                     const CCTK_INT *restrict cctk_lsh,
                                                     const gf_and_gz_struct *restrict in_prims,
                                                     gf_and_gz_struct *restrict out_prims_r,
                                                     gf_and_gz_struct *restrict out_prims_l ) {

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) {
    for(int j=0;j<cctk_lsh[1];j++) {
      for(int i=0;i<cctk_lsh[0];i++) {        
        CCTK_INT index  = CCTK_GFINDEX3D(cctkGH,i,j,k);

        //---------- Left face ----------
        CCTK_REAL xrhoR  = out_prims_r[RHOB    ].gf[index];
        CCTK_REAL xyeR   = out_prims_r[YEPRIM  ].gf[index];
        CCTK_REAL xtempR = in_prims[TEMPERATURE].gf[index];
        CCTK_REAL xentR  = 0.0;
        if( eos.PPM_reconstructed_var == PRESSURE ) {
          CCTK_REAL xprsR  = out_prims_r[PRESSURE].gf[index];
          CCTK_REAL xepsR  = 0.0;
          WVU_EOS_eps_S_and_T_from_rho_Ye_P( xrhoR,xyeR,xprsR, &xepsR,&xentR,&xtempR );
          out_prims_r[EPSILON  ].gf[index] = xepsR;
        }
        else if( eos.PPM_reconstructed_var == EPSILON ) {
          CCTK_REAL xprsR  = 0.0;
          CCTK_REAL xepsR  = out_prims_r[EPSILON].gf[index];
          WVU_EOS_P_S_and_T_from_rho_Ye_and_eps( xrhoR,xyeR,xepsR, &xprsR,&xentR,&xtempR );
          out_prims_r[PRESSURE ].gf[index] = xprsR;
        }
        out_prims_r[TEMPERATURE].gf[index] = xtempR;
        out_prims_r[ENTROPY    ].gf[index] = xentR;
        //-------------------------------

        //---------- Right face ---------
        CCTK_REAL xrhoL  = out_prims_l[RHOB    ].gf[index];
        CCTK_REAL xyeL   = out_prims_l[YEPRIM  ].gf[index];
        CCTK_REAL xtempL = in_prims[TEMPERATURE].gf[index];
        CCTK_REAL xentL  = 0.0;
        if( eos.PPM_reconstructed_var == PRESSURE ) {
          CCTK_REAL xprsL  = out_prims_l[PRESSURE].gf[index];
          CCTK_REAL xepsL  = 0.0;
          WVU_EOS_eps_S_and_T_from_rho_Ye_P( xrhoL,xyeL,xprsL, &xepsL,&xentL,&xtempL );
          out_prims_l[EPSILON  ].gf[index] = xepsL;
        }
        else if( eos.PPM_reconstructed_var == EPSILON ) {
          CCTK_REAL xprsL  = 0.0;
          CCTK_REAL xepsL  = out_prims_l[EPSILON].gf[index];
          WVU_EOS_P_S_and_T_from_rho_Ye_eps( xrhoL,xyeL,xepsL, &xprsL,&xentL,&xtempL );
          out_prims_l[PRESSURE ].gf[index] = xprsL;
        }
        out_prims_l[TEMPERATURE].gf[index] = xtempL;
        out_prims_l[ENTROPY    ].gf[index] = xentL;
        //-------------------------------
      }
    }
  }

}
