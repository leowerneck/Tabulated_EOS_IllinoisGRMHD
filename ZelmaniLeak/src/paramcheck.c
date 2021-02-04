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

    varindex = CCTK_VarIndex("GRHydro::tau");
    rhsindex = MoLQueryEvolvedRHS (varindex);
    numtl = CCTK_ActiveTimeLevelsVI(cctkGH, rhsindex);
    if(numtl <= 0) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "No storage for %s. varindex = %d, rhsindex = %d, numtl = %d",
                  "taurhs", varindex, rhsindex, numtl);
    }

    varindex = CCTK_VarIndex("GRHydro::scon[0]");
    rhsindex = MoLQueryEvolvedRHS (varindex);
    numtl = CCTK_ActiveTimeLevelsVI(cctkGH, rhsindex);
    if(numtl <= 0) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "No storage for %s. varindex = %d, rhsindex = %d, numtl = %d",
                  "sxrhs", varindex, rhsindex, numtl);
    }

    varindex = CCTK_VarIndex("GRHydro::scon[1]");
    rhsindex = MoLQueryEvolvedRHS (varindex);
    numtl = CCTK_ActiveTimeLevelsVI(cctkGH, rhsindex);
    if(numtl <= 0) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "No storage for %s. varindex = %d, rhsindex = %d, numtl = %d",
                  "syrhs", varindex, rhsindex, numtl);
    }

    varindex = CCTK_VarIndex("GRHydro::scon[2]");
    rhsindex = MoLQueryEvolvedRHS (varindex);
    numtl = CCTK_ActiveTimeLevelsVI(cctkGH, rhsindex);
    if(numtl <= 0) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "No storage for %s. varindex = %d, rhsindex = %d, numtl = %d",
                  "szrhs", varindex, rhsindex, numtl);
    }
  }

  if(do_tau) {
    if(!CCTK_IsImplementationActive("GRHydro")) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Please enable a thorn implementing GRHydro when using do_tau. While this *might* work without GRHydro if another prim2con is provided, I am not yet aware of sucha thorn");
    }
  }

  // TODO: instead consider adding STORAGE statements to schedule.ccl?
  if(CCTK_ActiveTimeLevels(cctkGH, "HydroBase::entropy") == 0 ||
     CCTK_ActiveTimeLevels(cctkGH, "HydroBase::temperature") == 0) {
    CCTK_ERROR("Require storage for entropy and temperature. Please make sure to set HydroBase::initial_entropy and HydroBase::initial_temperature to the correct values");
  }

}
