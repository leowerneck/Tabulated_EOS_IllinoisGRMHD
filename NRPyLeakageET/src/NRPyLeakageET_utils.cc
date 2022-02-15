#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "carpet.hh"
#include "NRPyLeakageET.h"

#ifdef CCTK_MPI
#  include <mpi.h>
#endif

using namespace Carpet;

extern "C"
void NRPyLeakageET_Prolongate(CCTK_ARGUMENTS, const char *varname) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(verbose) CCTK_VInfo(CCTK_THORNSTRING, "Prolongating %s...", varname);

  int ml=Carpet::mglevel;
  int ll=Carpet::reflevel;
  for (comm_state state; not state.done(); state.step()) {
    const int g = CCTK_GroupIndex(varname); 
    const int active_tl = CCTK_ActiveTimeLevelsGI (cctkGH, g);
    assert (active_tl>=0);
    const int tl = active_tl > 1 ? timelevel : 0;
    for (int m=0; m<(int)arrdata.AT(g).size(); ++m) {
      for (int v = 0; v < (int)arrdata.AT(g).AT(m).data.size(); ++v) {
        ggf *const gv = arrdata.AT(g).AT(m).data.AT(v);
        gv->ref_prolongate_all (state, tl, ll, ml, cctk_time);
        gv->ref_bnd_prolongate_all (state, tl, ll, ml, cctk_time);
      }
    }
  } // for state
}

extern "C"
void NRPyLeakageET_Restrict(CCTK_ARGUMENTS, const char *varname) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(verbose) CCTK_VInfo(CCTK_THORNSTRING, "Restricting %s...", varname);

  int ml=Carpet::mglevel;
  int ll=Carpet::reflevel;
  for (comm_state state; not state.done(); state.step()) {
    const int g = CCTK_GroupIndex(varname); 
    const int active_tl = CCTK_ActiveTimeLevelsGI (cctkGH, g);
    assert (active_tl>=0);
    const int tl = active_tl > 1 ? timelevel : 0;
    for (int m=0; m<(int)arrdata.AT(g).size(); ++m) {
      for (int v = 0; v < (int)arrdata.AT(g).AT(m).data.size(); ++v) {
        ggf *const gv = arrdata.AT(g).AT(m).data.AT(v);
        gv->ref_restrict_all (state, tl, ll, ml);
      }
    }
  } // for state
}

extern "C"
void NRPyLeakageET_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Step 1: Initialize all optical depth gridfunctions to zero
  if(verbose) CCTK_INFO("Initializing optical depths gridfunctions to zero...");
  BEGIN_REFLEVEL_LOOP(cctkGH) {
    BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
      BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        NRPyLeakageET_optical_depths_initialize_to_zero(CCTK_PASS_CTOC);
      } END_COMPONENT_LOOP;
    } END_MAP_LOOP;
  } END_REFLEVEL_LOOP;
  if(verbose) CCTK_INFO("Initialized all optical depths gridfunctions to zero");

  exit(1);
}
