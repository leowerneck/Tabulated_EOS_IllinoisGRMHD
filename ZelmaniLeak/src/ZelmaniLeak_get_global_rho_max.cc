//----------------------------------------------------------
//----------------------------------------------------------
// Thorn      : ZelmaniLeak
// File       : ZelmaniLeak_get_global_rho_max.cc
// Author(s)  : Leonardo R. Werneck (wernecklr@gmail.com)
// Description: This function gets the maximum value of the
//              density on the grid in all processors. It
//              is a slight generalization of the function
//              GRHydro_Rho_Minima_Setup_Final_PUGH() from
//              the GRHydro thorn.
//----------------------------------------------------------
//----------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C"
void ZelmaniLeak_get_global_rho_max(CCTK_ARGUMENTS) {

  // Gives us access to the pointer ZL_global_rho_max
  DECLARE_CCTK_ARGUMENTS;
  // Gives us access to the ZL_rho_gf_VarString parameter
  DECLARE_CCTK_PARAMETERS;

  // Get handle for the "maximum" MPI reduction operation
  const CCTK_INT Reduction_Handle = CCTK_ReductionHandle("maximum");
  // Set density gridfunction
  const CCTK_INT input_array_variable_indices={CCTK_VarIndex(ZL_rho_gf_VarString)};
  // Variable which will hold the maximum density
  CCTK_REAL max_rho_in_grid;
  // Perform the reduction
  const CCTK_INT ierr = CCTK_Reduce(cctkGH,
                                    -1, // target processors; -1 -> all
                                    Reduction_Handle,
                                    1,  // number of output variables   
                                    CCTK_VARIABLE_REAL,
                                    &max_rho_in_grid,
                                    1,  // number of variables to be reduced
                                    input_array_variable_indices);
  // Check if everything was okay
  if(ierr != 0) CCTK_VWarn(0,__LINE__, __FILE__, CCTK_THORNSTRING,"Failed to compute the global maximum of rho from gridfunction %s",ZL_rho_gf_VarString);

  // Update the ZelmaniLeak variable global_rho_max
  *ZL_global_rho_max = max_rho_in_grid;
}
