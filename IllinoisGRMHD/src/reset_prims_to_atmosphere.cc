#include "cctk.h"

#include "IllinoisGRMHD_headers.h"

void reset_prims_to_atmosphere( const igm_eos_parameters eos,
                                CCTK_REAL *restrict PRIMS ) {

  // Just a simple reset to atmospheric values.
  // Velocities are set to zero. Keeping it
  // inside a single function ensures that
  // resets are consistent throughout the code.
  PRIMS[RHOB         ] = eos.rho_atm;
  PRIMS[PRESSURE     ] = eos.P_atm;
  PRIMS[EPSILON      ] = eos.eps_atm;
  PRIMS[ENTROPY      ] = eos.S_atm;
  if( eos.is_Tabulated ) {
    PRIMS[YEPRIM     ] = eos.Ye_atm;
    PRIMS[TEMPERATURE] = eos.T_atm;
  }
  PRIMS[VX           ] = 0;
  PRIMS[VY           ] = 0;
  PRIMS[VZ           ] = 0;  
  
}
