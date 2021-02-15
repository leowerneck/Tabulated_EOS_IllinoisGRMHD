// Thorn      : IllinoisGRMHD
// File       : IllinoisGRMHD_parameter_check.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This file provides a function which checks
//              the IllinoisGRMHD input parameters.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "IllinoisGRMHD_headers.h"
#include "EOS_headers.hh"
#include "con2prim_headers.h"

extern "C"
void IllinoisGRMHD_parameter_check(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Initialize the IGM EOS parameter struct
  igm_eos_parameters eos;
  initialize_igm_eos_parameters_from_input(igm_eos_key,0,eos);

  // ------------------------------------
  // ---------- EOS parameters ----------
  // ------------------------------------
  if( eos.is_Hybrid ) {
    // Check if Gamma_th has been set
    if( Gamma_th == -1 ) CCTK_VError(VERR_DEF_PARAMS,"You must set Gamma_th in the parameter file. ABORTING.");

    // Check if K_ppoly_tab0 has been set
    if( K_ppoly_tab0 == -1 ) CCTK_VError(VERR_DEF_PARAMS,"You must set K_ppoly_tab0 in the parameter file. ABORTING.");

    // Check if rho_ppoly_tab have been set
    for(int i=0;i<neos-1;i++) {
      if( rho_ppoly_tab_in[i] == -1 ) {
        CCTK_VError(VERR_DEF_PARAMS,"neos was set to %d, but rho_ppoly_tab_in[%d] was not set. You must specify %d values of rho_ppoly_tab_in. ABORTING.",neos,i,neos-1);
      }
    }

    // Check if Gamma_ppoly_tab have been set
    for(int i=0;i<neos;i++) {
      if( Gamma_ppoly_tab_in[i] == -1 ) {
        CCTK_VError(VERR_DEF_PARAMS,"neos was set to %d, but Gamma_ppoly_tab_in[%d] was not set. You must specify %d values of Gamma_ppoly_tab_in. ABORTING.",neos,i,neos);
      }
    }
    
    // Print Hybrid EOS information
    CCTK_VInfo(CCTK_THORNSTRING,"EOS type: Hybrid");
    print_EOS_Hybrid(eos);
  }
  else if( eos.is_Tabulated ) {

    // Temperature evolution info
    if( igm_evolve_temperature == true ) {
      CCTK_VInfo(CCTK_THORNSTRING,"Temperature evolution is ENABLED!");
      CCTK_VInfo(CCTK_THORNSTRING,"Temperature evolution will begin at cctk_time: %g",igm_freeze_T_evolution_until_cctk_time);
    }
    else
      CCTK_VInfo(CCTK_THORNSTRING,"Temperature evolution is DISABLED!");

    CCTK_VInfo(CCTK_THORNSTRING,"EOS type: Tabulated");
    
  }

  // Entropy evolution info
  if( igm_evolve_entropy == true )
    CCTK_VInfo(CCTK_THORNSTRING,"Entropy evolution is ENABLED!");
  else
    CCTK_VInfo(CCTK_THORNSTRING,"Entropy evolution is DISABLED!");

  // ------------------------------------
  // ------- con2prim parameters --------
  // ------------------------------------
  if( eos.is_Hybrid ) {

    // Check if selected a con2prim routine which is supported by Hybrid EOS
    if( CCTK_EQUALS(igm_con2prim_routine,"CerdaDuran2D") ||
        CCTK_EQUALS(igm_con2prim_routine,"Palenzuela1D") ) {
      CCTK_VError(VERR_DEF_PARAMS,
                  "Hybrid EOS only supports the following con2prim routines: "
                  "Noble2D, Noble1D, Noble1D_entropy, and Noble1D_entropy2. ABORTING.");
    }
    else if( (igm_evolve_entropy == false) &&
             ( CCTK_EQUALS(igm_con2prim_routine,"Noble1D_entropy" ) ||
               CCTK_EQUALS(igm_con2prim_routine,"Noble1D_entropy2") ) ) {
      CCTK_VError(VERR_DEF_PARAMS,
                  "Routines Noble1D_entropy and Noble1D_entropy2 require enabling entropy evolution. "
                  "Please set igm_evolve_entropy=\"yes\" in the parameter file. ABORTING.");
    }

  }
  else if( eos.is_Tabulated ) {

    // Check if selected a con2prim routine which is supported by tabulated EOS
    if( CCTK_EQUALS(igm_con2prim_routine,"Noble1D_entropy" ) ||
        CCTK_EQUALS(igm_con2prim_routine,"Noble1D_entropy2") ) {
      CCTK_VError(VERR_DEF_PARAMS,
                  "Tabulated EOS only supports the following con2prim routines:"
                  "Palenzuela1D, Noble2D, Noble1D, and CerdaDuran2D. ABORTING.");
    }
    else if( (igm_evolve_entropy == false) && CCTK_EQUALS(igm_con2prim_routine,"Palenzuela1D" ) ) {
      CCTK_VInfo(CCTK_THORNSTRING,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING");
      CCTK_VInfo(CCTK_THORNSTRING,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING");
      CCTK_VInfo(CCTK_THORNSTRING,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING");
      CCTK_WARN(CCTK_WARN_ALERT, "*********** It is highly recommended to enable entropy evolution when using the Palenzuela1D routine **********");
      CCTK_VInfo(CCTK_THORNSTRING,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING");
      CCTK_VInfo(CCTK_THORNSTRING,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING");
      CCTK_VInfo(CCTK_THORNSTRING,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING");
    }
  }

  CCTK_VInfo(CCTK_THORNSTRING,"Primary conservative-to-primitive routine: %s",igm_con2prim_routine);
  if( CCTK_EQUALS(igm_con2prim_backup_routine[0],"None") ) {
    CCTK_VInfo(CCTK_THORNSTRING,"No backup conservative-to-primitive routine selected.");
  }
  else {
    CCTK_VInfo(CCTK_THORNSTRING,"Backup conservative-to-primitive routine #1: %s",igm_con2prim_backup_routine[0]);
    if( !CCTK_EQUALS(igm_con2prim_backup_routine[1],"None") ) {
      CCTK_VInfo(CCTK_THORNSTRING,"Backup conservative-to-primitive routine #2: %s",igm_con2prim_backup_routine[1]);
      if( !CCTK_EQUALS(igm_con2prim_backup_routine[2],"None") ) {
        CCTK_VInfo(CCTK_THORNSTRING,"Backup conservative-to-primitive routine #3: %s",igm_con2prim_backup_routine[2]);
      }
    }
  }

  // ------------------------------------
  // --------- Other parameters ---------
  // ------------------------------------
  if( GAMMA_SPEED_LIMIT > 10.0 ) {
    CCTK_VInfo(CCTK_THORNSTRING,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING");
    CCTK_VInfo(CCTK_THORNSTRING,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING");
    CCTK_VInfo(CCTK_THORNSTRING,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING");
    CCTK_WARN(CCTK_WARN_ALERT, "******************** It is not recommended to set GAMMA_SPEED_LIMIT to values larger than 10 *******************");
    CCTK_VInfo(CCTK_THORNSTRING,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING");
    CCTK_VInfo(CCTK_THORNSTRING,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING");
    CCTK_VInfo(CCTK_THORNSTRING,"WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING");
  }
}
