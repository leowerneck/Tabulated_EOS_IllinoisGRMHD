#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <assert.h>

void ZelmaniLeak_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(do_pnu) {
    // check that GRHydro RHS are available when adding pressure terms
    int varindex, rhsindex, numtl;

    varindex = CCTK_VarIndex(ZL_tau_gf_VarString);
    rhsindex = MoLQueryEvolvedRHS (varindex);
    numtl = CCTK_ActiveTimeLevelsVI(cctkGH, rhsindex);
    if(numtl <= 0) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "No storage for %s. varindex = %d, rhsindex = %d, numtl = %d",
                  ZL_tau_gf_VarString, varindex, rhsindex, numtl);
    }

    varindex = CCTK_VarIndex(ZL_SD0_gf_VarString);
    rhsindex = MoLQueryEvolvedRHS (varindex);
    numtl = CCTK_ActiveTimeLevelsVI(cctkGH, rhsindex);
    if(numtl <= 0) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "No storage for %s. varindex = %d, rhsindex = %d, numtl = %d",
                  ZL_SD0_gf_VarString, varindex, rhsindex, numtl);
    }

    varindex = CCTK_VarIndex(ZL_SD1_gf_VarString);
    rhsindex = MoLQueryEvolvedRHS (varindex);
    numtl = CCTK_ActiveTimeLevelsVI(cctkGH, rhsindex);
    if(numtl <= 0) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "No storage for %s. varindex = %d, rhsindex = %d, numtl = %d",
                  ZL_SD1_gf_VarString, varindex, rhsindex, numtl);
    }

    varindex = CCTK_VarIndex(ZL_SD2_gf_VarString);
    rhsindex = MoLQueryEvolvedRHS (varindex);
    numtl = CCTK_ActiveTimeLevelsVI(cctkGH, rhsindex);
    if(numtl <= 0) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "No storage for %s. varindex = %d, rhsindex = %d, numtl = %d",
                  ZL_SD2_gf_VarString, varindex, rhsindex, numtl);
    }
  }

  if(do_tau) {
    if(!CCTK_IsImplementationActive("IllinoisGRMHD") && !CCTK_IsImplementationActive("GRHydro")) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Please enable a GRMHD evolution thorn when using do_tau. Current supported options: IllinoisGRMHD & GRHydro");
    }
  }

  // TODO: instead consider adding STORAGE statements to schedule.ccl?
  if(CCTK_ActiveTimeLevels(cctkGH, ZL_entropy_gf_VarString) == 0 ||
     CCTK_ActiveTimeLevels(cctkGH, ZL_temperature_gf_VarString) == 0) {
    CCTK_ERROR("Require storage for entropy and temperature. If using IllinoisGRMHD, one needs igm_entropy and igm_temperature. If using GRHydro, one needs the HydroBase counterparts.");
  }
}
