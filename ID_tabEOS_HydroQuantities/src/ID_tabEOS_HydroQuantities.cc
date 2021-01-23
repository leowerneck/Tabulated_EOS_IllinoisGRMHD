#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "ID_tabEOS_HydroQuantities.h"

extern "C" void ID_tabEOS_HydroQuantities__initial_Y_e( const CCTK_INT  npoints,
                                                        const CCTK_REAL *restrict rho,
                                                        CCTK_REAL *restrict Y_e ) {

  DECLARE_CCTK_PARAMETERS;

  // Open the Y_e file, which should countain Y_e(rho) for the EOS table slice
  FILE *Y_e_file = fopen(Y_e_filename,"r");
  
  // Check if everything is OK with the file
  if( (Y_e_file = fopen(Y_e_filename,"r")) == NULL ) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "File \"%s\" does not exist. ABORTING", Y_e_filename);
  }
  else {
    // Set nrho
    const CCTK_INT nrho = nuc_eos_private::nrho;

    // Set rho
    CCTK_REAL rho_arr[nrho];
    for(int i=0;i<nrho;i++) rho_arr[i] = exp(nuc_eos_private::logrho[i]);

    // Now read in the Y_e file
    CCTK_REAL Y_e_of_rho_arr[nrho];
    read_1dfile__set_array(Y_e_file,Y_e_of_rho_arr);
    // Close the file
    fclose(Y_e_file);

    // Set interpolation stencil size
    const CCTK_INT interp_stencil_size = 5;
    for(int i=0;i<npoints;i++) {
      if( rho[i] > id_rho_atm ) {
        // Interpolate Y_e(rho_i) at gridpoint i
        CCTK_REAL Y_eL;
        interpolate_1d_quantity_as_function_of_rho(interp_stencil_size,nrho,rho_arr,Y_e_of_rho_arr,rho[i],&Y_eL);
        // Finally, set the Y_e gridfunction
        Y_e[i] = MIN(MAX(Y_eL,nuc_eos::eos_yemin),nuc_eos::eos_yemax);
      }
      else {
        Y_e[i] = id_Y_e_atm;
      }
    }
  }
}

// Set initial temperature to be constant everywhere (TODO: add other options)
extern "C" void ID_tabEOS_HydroQuantities__initial_temperature( const CCTK_INT  npoints,
                                                                CCTK_REAL *restrict temperature ) {

  DECLARE_CCTK_PARAMETERS;

  // Loop over the grid, initializing the temperature
  for(int i=0;i<npoints;i++) temperature[i] = id_temperature;
}

// Now recompute all HydroQuantities, to ensure consistent initial data
extern "C" void ID_tabEOS_HydroQuantities__recompute_HydroQuantities( const CCTK_INT npoints,
                                                                      CCTK_REAL *restrict rho,
                                                                      CCTK_REAL *restrict Y_e,
                                                                      CCTK_REAL *restrict temperature,
                                                                      CCTK_REAL *restrict press,
                                                                      CCTK_REAL *restrict eps,
                                                                      CCTK_REAL *restrict entropy,
                                                                      CCTK_REAL *restrict vel ) {

  DECLARE_CCTK_PARAMETERS;

  // Prepare for EOS function calls
  const CCTK_INT  havetemp     = 1;
  const CCTK_INT  eoskey       = EOS_Omni_GetHandle("nuc_eos");
  const CCTK_REAL rf_precision = id_rf_precision;

  // Check whether or not we want to initialize the entropy in this thorn
  bool initialize_entropy;
  if( CCTK_EQUALS( initial_entropy,"ID_tabEOS_HydroQuantities" ) ) {
    CCTK_VInfo(CCTK_THORNSTRING,"Entropy initialization is ENABLED!");
    initialize_entropy = true;
  }
  else {
    CCTK_VInfo(CCTK_THORNSTRING,"Entropy initialization is DISABLED!");
    initialize_entropy = false;
  }

  CCTK_VInfo(CCTK_THORNSTRING,"Recomputing all HydroBase quantities ...");

  // Set atmospheric values
  CCTK_REAL id_temp_atm  = id_temperature_atm;
  CCTK_REAL id_press_atm = 0.0;
  CCTK_REAL id_eps_atm   = 0.0;
  CCTK_REAL id_ent_atm   = 0.0;
  CCTK_REAL dummy        = 0.0;
  CCTK_INT  keyerr       = 0;
  CCTK_INT  anyerr       = 0;
  if( initialize_entropy ) {
    // Only call EOS_Omni_short() if we need the entropy.
    EOS_Omni_short(eoskey,havetemp,rf_precision,1,
                   &id_rho_atm,&id_eps_atm,&id_temp_atm,&id_Y_e_atm,&id_press_atm,&id_ent_atm,
                   &dummy,&dummy,&dummy,&dummy,&dummy,
                   &keyerr,&anyerr);
  }
  else {
    // Otherwise use EOS_Omni_press(), which performs fewer
    // interpolations and therefore is more efficient.
    EOS_Omni_press(eoskey,havetemp,rf_precision,1,
                   &id_rho_atm,&id_eps_atm,&id_temp_atm,&id_Y_e_atm,&id_press_atm,
                   &keyerr,&anyerr);
  }
  // Loop over the grid, recomputing the HydroBase quantities
  for(int i=0;i<npoints;i++) {
    CCTK_REAL xrho    = rho[i];

    if( xrho > id_rho_atm ) {
      CCTK_REAL xye     = Y_e[i];
      CCTK_REAL xtemp   = temperature[i];
      CCTK_REAL xpress  = 0.0;
      CCTK_REAL xeps    = 0.0;
      CCTK_REAL xent    = 0.0;
      dummy             = 0.0;
      if( initialize_entropy ) {
        // Only call EOS_Omni_short() if we need the entropy.
        EOS_Omni_short(eoskey,havetemp,rf_precision,1,
                       &xrho,&xeps,&xtemp,&xye,&xpress,&xent,
                       &dummy,&dummy,&dummy,&dummy,&dummy,
                       &keyerr,&anyerr);
      }
      else {
        // Otherwise use EOS_Omni_press(), which performs fewer
        // interpolations and therefore is more efficient.
        EOS_Omni_press(eoskey,havetemp,rf_precision,1,
                       &xrho,&xeps,&xtemp,&xye,&xpress,
                     &keyerr,&anyerr);
      }
    // Now set press, eps, and entropy gridfunctions.
    press[    i] = xpress;
    eps[      i] = xeps;
    if( initialize_entropy )
      entropy[i] = xent;
    }
    else {
      // Reset to atmosphere
      rho[          i] = id_rho_atm;
      Y_e[          i] = id_Y_e_atm;
      temperature[  i] = id_temperature_atm;
      press[        i] = id_press_atm;
      eps[          i] = id_eps_atm;
      vel[          i] = 0.0;
      vel[npoints  +i] = 0.0;
      vel[2*npoints+i] = 0.0;
      if( initialize_entropy )
        entropy[    i] = id_ent_atm;
    }
  }
  
}

extern "C" void ID_tabEOS_HydroQuantities(CCTK_ARGUMENTS) {
  
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Set total number of grid points
  const CCTK_INT npoints = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];

  // This function sets Y_e from a file containing Y_e(rho),
  // typically in neutrino-free beta-equilibrium.
  if( CCTK_EQUALS( initial_Y_e,"ID_tabEOS_HydroQuantities") ) {
    CCTK_VInfo(CCTK_THORNSTRING,"Y_e initialization is ENABLED!");
    // Initialize Y_e
    ID_tabEOS_HydroQuantities__initial_Y_e( npoints, rho, Y_e );
  }
  else {
    CCTK_VInfo(CCTK_THORNSTRING,"Y_e initialization is DISABLED!");
  }

  // This function sets the temperature everywhere in the grid to
  // a constant value. We can generalize it later to allow for
  // other types of temperature initial data.
  if( CCTK_EQUALS( initial_temperature,"ID_tabEOS_HydroQuantities") ) {
    CCTK_VInfo(CCTK_THORNSTRING,"temperature initialization is ENABLED!");
    // Initialize the temperature
    ID_tabEOS_HydroQuantities__initial_temperature( npoints, temperature );
  }
  else {
    CCTK_VInfo(CCTK_THORNSTRING,"temperature initialization is DISABLED!");
  }

  // Now recompute the HydroBase quantities, to ensure consistency
  if( CCTK_EQUALS( initial_Y_e        ,"ID_tabEOS_HydroQuantities") ||
      CCTK_EQUALS( initial_temperature,"ID_tabEOS_HydroQuantities") ||
      CCTK_EQUALS( initial_entropy    ,"ID_tabEOS_HydroQuantities") ) {
    ID_tabEOS_HydroQuantities__recompute_HydroQuantities( npoints,
                                                          rho,Y_e,temperature,
                                                          press,eps,entropy,
                                                          vel );
  }

  CCTK_VInfo(CCTK_THORNSTRING,"All done!");
  
}
