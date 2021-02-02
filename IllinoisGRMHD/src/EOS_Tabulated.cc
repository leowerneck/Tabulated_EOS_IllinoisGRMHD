// Thorn      : IllinoisGRMHD
// File       : EOS_Tabulated.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: In this file we provide tabulated EOS functions.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "IllinoisGRMHD_headers.h"
#include "HARM_TabEOS_helpers.hh"

// These are useful for us, though not defined
// explicitly by EOS_Omni
const CCTK_INT table_key_pressure = 0;
const CCTK_INT table_key_epsilon  = 1;
const CCTK_INT table_key_entropy  = 2;

// These are also very useful
const CCTK_INT have_eps  = 0;
const CCTK_INT have_temp = 1;
const CCTK_INT have_ent  = 2;
const CCTK_INT have_prs  = 3;
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
    get_P_eps_and_S_from_rho_Ye_and_T( eos,
                                       eos.rho_atm,eos.Ye_atm,eos.T_atm,
                                       &eos.P_atm,&eos.eps_atm,&eos.S_atm );
  }
  else {
    get_P_and_eps_from_rho_Ye_and_T( eos,
                                     eos.rho_atm,eos.Ye_atm,eos.T_atm,
                                     &eos.P_atm,&eos.eps_atm );
  }
  // Atmospheric tau
  eos.tau_atm = tau_atm;
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

  // All done!

}

void get_P_and_eps_from_rho_Ye_and_T( const igm_eos_parameters eos,
                                      const CCTK_REAL rho,
                                      const CCTK_REAL Y_e,
                                      const CCTK_REAL T,
                                      CCTK_REAL *restrict P,
                                      CCTK_REAL *restrict eps ) {
  // Set up for table call
  CCTK_INT  npoints = 1;
  CCTK_INT  keyerr  = 0;
  CCTK_INT  anyerr  = 0;
  // The EOS_Omni function does not like some of its arguments
  // being constant (even when they are not supposed to change).
  // We declare auxiliary variables here to avoid errors.
  CCTK_REAL rho_in  = rho;
  CCTK_REAL Y_e_in  = Y_e;
  CCTK_REAL T_in    = T;
  CCTK_REAL prs_out = 0.0;
  CCTK_REAL eps_out = 0.0;

  // Perform the table interpolations
  EOS_Omni_press( eos.key,have_temp,eos.root_finding_precision,npoints,
                  &rho_in,&eps_out,&T_in,&Y_e_in,&prs_out,
                  &keyerr,&anyerr );

  // Now update P and eps
  *P   = prs_out;
  *eps = eps_out;

  // FIXME: add error handling!  
}

void get_P_and_T_from_rho_Ye_and_eps( const igm_eos_parameters eos,
                                      const CCTK_REAL rho,
                                      const CCTK_REAL Y_e,
                                      const CCTK_REAL eps,
                                      CCTK_REAL *restrict P,
                                      CCTK_REAL *restrict T ) {
  // Set up for table call
  CCTK_INT  npoints = 1;
  CCTK_INT  keyerr  = 0;
  CCTK_INT  anyerr  = 0;
  // The EOS_Omni function does not like some of its arguments
  // being constant (even when they are not supposed to change).
  // We declare auxiliary variables here to avoid errors.
  CCTK_REAL rho_in  = rho;
  CCTK_REAL Y_e_in  = Y_e;
  CCTK_REAL eps_in  = eps;
  CCTK_REAL prs_out = 0.0;
  CCTK_REAL T_out   = *T;

  // Perform the table interpolations
  EOS_Omni_press( eos.key,have_eps,eos.root_finding_precision,npoints,
                  &rho_in,&eps_in,&T_out,&Y_e_in,&prs_out,
                  &keyerr,&anyerr );

  // Now update P and eps
  *P = prs_out;
  *T = T_out;

  // FIXME: add error handling!  
}

void get_P_S_and_T_from_rho_Ye_and_eps( const igm_eos_parameters eos,
                                        const CCTK_REAL rho,
                                        const CCTK_REAL Y_e,
                                        const CCTK_REAL eps,
                                        CCTK_REAL *restrict P,
                                        CCTK_REAL *restrict S,
                                        CCTK_REAL *restrict T ) {
  // Set up for table call
  CCTK_INT  npoints = 1;
  CCTK_INT  keyerr  = 0;
  CCTK_INT  anyerr  = 0;
  // The EOS_Omni function does not like some of its arguments
  // being constant (even when they are not supposed to change).
  // We declare auxiliary variables here to avoid errors.
  CCTK_REAL rho_in  = rho;
  CCTK_REAL Y_e_in  = Y_e;
  CCTK_REAL eps_in  = eps;
  CCTK_REAL prs_out = 0.0;
  CCTK_REAL S_out   = 0.0;
  CCTK_REAL T_out   = *T;
  CCTK_REAL dummy   = 0.0;

  // Perform the table interpolations
  EOS_Omni_short( eos.key,have_eps,eos.root_finding_precision,npoints,
                  &rho_in,&eps_in,&T_out,&Y_e_in,&prs_out,&S_out,
                  &dummy,&dummy,&dummy,&dummy,&dummy,
                  &keyerr,&anyerr );

  // Now update P and eps
  *P = prs_out;
  *S = S_out;
  *T = T_out;

  // FIXME: add error handling!  
}


void get_P_eps_and_S_from_rho_Ye_and_T( const igm_eos_parameters eos,
                                        const CCTK_REAL rho,
                                        const CCTK_REAL Y_e,
                                        const CCTK_REAL T,
                                        CCTK_REAL *restrict P,
                                        CCTK_REAL *restrict eps,
                                        CCTK_REAL *restrict S ) {
  // Set up for table call
  CCTK_INT  npoints = 1;
  CCTK_INT  keyerr  = 0;
  CCTK_INT  anyerr  = 0;
  // The EOS_Omni function does not like some of its arguments
  // being constant (even when they are not supposed to change).
  // We declare auxiliary variables here to avoid errors.
  CCTK_REAL rho_in  = rho;
  CCTK_REAL Y_e_in  = Y_e;
  CCTK_REAL T_in    = T;
  CCTK_REAL prs_out = 0.0;
  CCTK_REAL eps_out = 0.0;
  CCTK_REAL S_out   = 0.0;
  CCTK_REAL dummy   = 0.0;

  // Perform the table interpolations
  EOS_Omni_short( eos.key,have_temp,eos.root_finding_precision,npoints,
                  &rho_in,&eps_out,&T_in,&Y_e_in,&prs_out,&S_out,
                  &dummy,&dummy,&dummy,&dummy,&dummy,
                  &keyerr,&anyerr );

  // Now update P, eps, and S
  *P   = prs_out;
  *eps = eps_out;
  *S   = S_out;

  // FIXME: add error handling!
}

void get_P_eps_and_T_from_rho_Ye_and_S( const igm_eos_parameters eos,
                                        const CCTK_REAL rho,
                                        const CCTK_REAL Y_e,
                                        const CCTK_REAL S,
                                        CCTK_REAL *restrict P,
                                        CCTK_REAL *restrict eps,
                                        CCTK_REAL *restrict T ) {
  // Set up for table call
  CCTK_INT  npoints = 1;
  CCTK_INT  keyerr  = 0;
  CCTK_INT  anyerr  = 0;
  // The EOS_Omni function does not like some of its arguments
  // being constant (even when they are not supposed to change).
  // We declare auxiliary variables here to avoid errors.
  CCTK_REAL rho_in  = rho;
  CCTK_REAL Y_e_in  = Y_e;
  CCTK_REAL S_in    = S;
  CCTK_REAL prs_out = 0.0;
  CCTK_REAL eps_out = 0.0;
  CCTK_REAL T_out   = *T;
  CCTK_REAL dummy   = 0.0;

  // Perform the table interpolations
  EOS_Omni_short( eos.key,have_ent,eos.root_finding_precision,npoints,
                  &rho_in,&eps_out,&T_out,&Y_e_in,&prs_out,&S_in,
                  &dummy,&dummy,&dummy,&dummy,&dummy,
                  &keyerr,&anyerr );

  // Now update P, eps, and S
  *P   = prs_out;
  *eps = eps_out;
  *T   = T_out;

  // FIXME: add error handling!
}

void get_P_eps_S_and_cs2_from_rho_Ye_and_T( const igm_eos_parameters eos,
                                            const CCTK_REAL rho,
                                            const CCTK_REAL Y_e,
                                            const CCTK_REAL T,
                                            CCTK_REAL *restrict P,
                                            CCTK_REAL *restrict eps,
                                            CCTK_REAL *restrict S,
                                            CCTK_REAL *restrict cs2 ) {
  // Set up for table call
  CCTK_INT  npoints = 1;
  CCTK_INT  keyerr  = 0;
  CCTK_INT  anyerr  = 0;
  // The EOS_Omni function does not like some of its arguments
  // being constant (even when they are not supposed to change).
  // We declare auxiliary variables here to avoid errors.
  CCTK_REAL rho_in  = rho;
  CCTK_REAL Y_e_in  = Y_e;
  CCTK_REAL T_in    = T;
  CCTK_REAL prs_out = 0.0;
  CCTK_REAL eps_out = 0.0;
  CCTK_REAL ent_out = 0.0;
  CCTK_REAL cs2_out = 0.0;
  CCTK_REAL dummy   = 0.0;

  // Perform the table interpolations
  EOS_Omni_short( eos.key,have_temp,eos.root_finding_precision,npoints,
                  &rho_in,&eps_out,&T_in,&Y_e_in,&prs_out,&ent_out,
                  &cs2_out,&dummy,&dummy,&dummy,&dummy,
                  &keyerr,&anyerr );

  // Now update P, cs2, and T
  *P   = prs_out;
  *eps = eps_out;
  *S   = ent_out;
  *cs2 = cs2_out;

  // FIXME: add error handling!
}

void get_eps_S_and_T_from_rho_Ye_and_P( const igm_eos_parameters eos,
                                        const CCTK_REAL rho,
                                        const CCTK_REAL Y_e,
                                        const CCTK_REAL P,
                                        CCTK_REAL *restrict eps,
                                        CCTK_REAL *restrict S,
                                        CCTK_REAL *restrict T ) {
  // Set up for table call
  CCTK_INT  npoints = 1;
  CCTK_INT  keyerr  = 0;
  CCTK_INT  anyerr  = 0;
  // The EOS_Omni function does not like some of its arguments
  // being constant (even when they are not supposed to change).
  // We declare auxiliary variables here to avoid errors.
  CCTK_REAL rho_in  = rho;
  CCTK_REAL Y_e_in  = Y_e;
  CCTK_REAL P_in    = P;
  CCTK_REAL eps_out = 0.0;
  CCTK_REAL ent_out = 0.0;
  CCTK_REAL T_out   = *T;
  CCTK_REAL dummy   = 0.0;

  // Perform the table interpolations
  EOS_Omni_short( eos.key,have_prs,eos.root_finding_precision,npoints,
                  &rho_in,&eps_out,&T_out,&Y_e_in,&P_in,&ent_out,
                  &dummy,&dummy,&dummy,&dummy,&dummy,
                  &keyerr,&anyerr );

  // Now update eps, S, and T
  *eps = eps_out;
  *S   = ent_out;
  *T   = T_out;

  // FIXME: add error handling!
}

void check_temperature_reconstruction( const igm_eos_parameters eos,
                                       const cGH *restrict cctkGH,
                                       const int *restrict cctk_lsh,
                                       gf_and_gz_struct *restrict prims_center,
                                       gf_and_gz_struct *restrict prims_right,
                                       gf_and_gz_struct *restrict prims_left ) {

  // Since the default PPM reconstruction (rho,P,vx,vy,vz) is
  // incredibly robust and has been extensively validated, we
  // want to preserve it as much as we can. However, the nature
  // of the EOS table is that the pressure and the specific
  // internal energy at high density regions have a very weak
  // temperature dependence. This can be incredibly frustrating,
  // because a small numerical perturbation on the pressure, for
  // example, can be greatly amplified during a table inversion
  // and result in a temperature that is massively different than
  // what we expect. While the difference in pressure may be of
  // order ~0.1%,, the resulting temperature amplification
  // may be of order ~10 or more!
  //
  // Our goal here is to make sure that the temperature doesn't go
  // berserk after it has been recovered from the pressure using
  // the EOS table interpolator. To this end, we will compare
  // the results between the PPM reconstruction of the temperature
  // and the one obtained from the interpolator. Choosing which
  // one is best can be tricky, though.

  const int imin = 0, imax = cctk_lsh[0];
  const int jmin = 0, jmax = cctk_lsh[1];
  const int kmin = 0, kmax = cctk_lsh[2];
#pragma omp parallel for
  for(int k=kmin;k<kmax;k++) {
    for(int j=jmin;j<jmax;j++) {
      for(int i=imin;i<imax;i++) {

        // Current index
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // Read in the density
        const CCTK_REAL rho_c   = prims_center[RHOB].gf[index];
        const CCTK_REAL rho_r   = prims_right[ RHOB].gf[index];
        const CCTK_REAL rho_l   = prims_left[  RHOB].gf[index];

        // For now, let us apply this "fix" only if we are in
        // a high density region. No need to worry about the
        // atmosphere.
        if( (rho_c > 1e-6*eos.rho_max) ||
            (rho_r > 1e-6*eos.rho_max) ||
            (rho_l > 1e-6*eos.rho_max) ) {

          // Read in the temperature at the cell centered, which
          // was computed during conservative-to-primitive.
          const CCTK_REAL T_c = prims_center[TEMPERATURE].gf[index];

          CCTK_REAL PRIMS[MAXNUMVARS];

          //--------------------------------------------------------
          //---------------------- Right face ----------------------
          //--------------------------------------------------------
          PRIMS[RHOB       ] = prims_right[RHOB       ].gf[index];
          PRIMS[YEPRIM     ] = prims_right[YEPRIM     ].gf[index];
          PRIMS[TEMPERATURE] = prims_right[TEMPERATURE].gf[index];
          PRIMS[PRESSURE   ] = prims_right[PRESSURE   ].gf[index];
          PRIMS[EPSILON    ] = prims_right[EPSILON    ].gf[index];
          PRIMS[ENTROPY    ] = prims_right[ENTROPY    ].gf[index];

          // Look for the best temperature
          find_lowest_gradient_temperature_recompute_prims(eos,T_c,PRIMS);

          // Update the gridfunctions
          prims_right[TEMPERATURE].gf[index] = PRIMS[TEMPERATURE];
          prims_right[PRESSURE   ].gf[index] = PRIMS[PRESSURE   ];
          prims_right[EPSILON    ].gf[index] = PRIMS[EPSILON    ];
          prims_right[ENTROPY    ].gf[index] = PRIMS[ENTROPY    ];

          //--------------------------------------------------------
          //---------------------- Left face -----------------------
          //--------------------------------------------------------
          PRIMS[RHOB       ] = prims_left[RHOB       ].gf[index];
          PRIMS[YEPRIM     ] = prims_left[YEPRIM     ].gf[index];
          PRIMS[TEMPERATURE] = prims_left[TEMPERATURE].gf[index];
          PRIMS[PRESSURE   ] = prims_left[PRESSURE   ].gf[index];
          PRIMS[EPSILON    ] = prims_left[EPSILON    ].gf[index];
          PRIMS[ENTROPY    ] = prims_left[ENTROPY    ].gf[index];

          // Look for the best temperature
          find_lowest_gradient_temperature_recompute_prims(eos,T_c,PRIMS);

          // Update the gridfunctions
          prims_left[TEMPERATURE].gf[index] = PRIMS[TEMPERATURE];
          prims_left[PRESSURE   ].gf[index] = PRIMS[PRESSURE   ];
          prims_left[EPSILON    ].gf[index] = PRIMS[EPSILON    ];
          prims_left[ENTROPY    ].gf[index] = PRIMS[ENTROPY    ];

        }
        else {
          // We won't fix low density regions
          continue;
        }
      } // for(int i=imin;i<imax;i++)
    } // for(int j=jmin;j<jmax;j++)
  } // for(int k=kmin;k<kmax;k++)
} // check_temperature_reconstruction()

void find_lowest_gradient_temperature_recompute_prims( const igm_eos_parameters eos,
                                                       const CCTK_REAL T_center,
                                                       CCTK_REAL *restrict PRIMS ) {

  // Start by computing the primitives using the table.
  // We will recover the temperature at the right and
  // left faces using (rho,Ye,P), using the value of
  // the temperature at the cell center as an initial
  // guess for the root finding method.
  CCTK_REAL xrho  = PRIMS[RHOB       ];
  CCTK_REAL xye   = PRIMS[YEPRIM     ];
  CCTK_REAL xprs  = PRIMS[PRESSURE   ];
  CCTK_REAL xtemp = T_center;
  CCTK_REAL xeps  = 0.0;
  CCTK_REAL xent  = 0.0;

  // Now perform the table inversion
  get_eps_S_and_T_from_rho_Ye_and_P(eos,xrho,xye,xprs, &xeps,&xent,&xtemp);

  // Then check the gradient
  const CCTK_REAL dT_PPM = fabs(T_center - PRIMS[TEMPERATURE]);
  const CCTK_REAL dT_EOS = fabs(T_center - xtemp);

  if( dT_PPM < dT_EOS ) {
    // In this case we use the value of T that was
    // reconstructed. We must recompute P, eps, and S.
    get_P_eps_and_S_from_rho_Ye_and_T(eos,
                                      PRIMS[RHOB],PRIMS[YEPRIM],PRIMS[TEMPERATURE],
                                      &PRIMS[PRESSURE],&PRIMS[EPSILON],&PRIMS[ENTROPY]);
  }
  else {
    // In this case we keep the values obtained from the EOS
    // interpolator, so we must update T, eps, and S
    PRIMS[TEMPERATURE] = xtemp;
    PRIMS[EPSILON    ] = xeps;
    PRIMS[ENTROPY    ] = xent;
  }

}


// This is from Dan Siegel's repo, modified to suit our needs
void get_P_eps_T_dPdrho_and_dPdeps_from_rho_Ye_and_h( const igm_eos_parameters eos,
                                                      const CCTK_REAL rho_in,
                                                      const CCTK_REAL Ye_in,
                                                      const CCTK_REAL h_in,
                                                      CCTK_REAL *restrict P,
                                                      CCTK_REAL *restrict T,
                                                      CCTK_REAL *restrict eps,
                                                      CCTK_REAL *restrict dPdrho, 
                                                      CCTK_REAL *restrict dPdeps ) {
  int keyerr    = 0;
  int anyerr    = 0;
  CCTK_REAL xtemp    = *T;  
  CCTK_REAL xeps     = 0.0;
  CCTK_REAL xprs     = 0.0;
  CCTK_REAL xdpdeps  = 0.0;
  CCTK_REAL xdpdrho  = 0.0;

  HARM_TabEOS_helpers::HARM_EOS_dpdrho_dpdeps_hinv(eos.root_finding_precision,1,
                                                   &rho_in,&h_in,&xtemp,&Ye_in,&xeps,&xprs,&xdpdrho,&xdpdeps,
                                                   &keyerr,&anyerr);
  
  *P      = xprs;
  *T      = xtemp;
  *eps    = xeps;
  *dPdeps = xdpdeps;
  *dPdrho = xdpdrho;
  
}
