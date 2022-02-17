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
int NRPyLeakageET_ProcessOwnsData() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int pc = vhh.AT(Carpet::mglevel)->processor(Carpet::reflevel,Carpet::component);
  return (rank == pc);
}

extern "C"
void NRPyLeakageET_Prolongate(CCTK_ARGUMENTS, const char *varname) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(verbose) CCTK_VInfo(CCTK_THORNSTRING, "Prolongating %s from ref. lvl. %d to ref. lvl. %d...", varname, GetRefinementLevel(cctkGH)-1, GetRefinementLevel(cctkGH));

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

  if(verbose) CCTK_VInfo(CCTK_THORNSTRING, "Restricting %s from ref. lvl. %d to ref. lvl. %d...", varname, GetRefinementLevel(cctkGH)+1, GetRefinementLevel(cctkGH));

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
  int RemainingIterations = 0;
  if(verbose) CCTK_INFO("Initializing optical depths gridfunctions to zero...");
  BEGIN_REFLEVEL_LOOP(cctkGH) {
    BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
      BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        NRPyLeakageET_optical_depths_initialize_to_zero(CCTK_PASS_CTOC);
        NRPyLeakageET_GetMaxSize(CCTK_PASS_CTOC,&RemainingIterations);
      } END_COMPONENT_LOOP;
    } END_MAP_LOOP;
  } END_REFLEVEL_LOOP;
  if(verbose) CCTK_INFO("Initialized all optical depths gridfunctions to zero");

  if( CCTK_EQUALS(optical_depth_evolution_type,"PathOfLeastResistance") ) {
    // Step 2: Now perform iterations of the path of least resistance algorithm
    RemainingIterations *= IterationFactor;
    int counter = 0;
    while( RemainingIterations ) {
      // Step 2.a: First perform the iteration on refinement level 0
      if( verbose ) CCTK_VInfo(CCTK_THORNSTRING,"Starting iteration %d...",counter+1);
      ENTER_LEVEL_MODE(cctkGH,0) {
        BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
          BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
            NRPyLeakageET_compute_opacities(CCTK_PASS_CTOC);
            NRPyLeakageET_optical_depths_PathOfLeastResistance(CCTK_PASS_CTOC);
          } END_COMPONENT_LOOP;
        } END_MAP_LOOP;
        CCTK_SyncGroup(cctkGH,"NRPyLeakageET::NRPyLeakageET_optical_depths");
      } LEAVE_LEVEL_MODE;

      // Step 2.b: Now loop over remaining refinement levels, performing the
      //           path of least resistance (POLR) algorithm
      for(int rl=1;rl<Carpet::reflevels;rl++) {
        ENTER_LEVEL_MODE(cctkGH,rl) {
          // Step 2.b.i: Prolongate results from the previous refinement level
          NRPyLeakageET_Prolongate(CCTK_PASS_CTOC,"NRPyLeakageET::NRPyLeakageET_optical_depths");
          BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
            BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
              // Step 2.b.i: Compute opacities
              NRPyLeakageET_compute_opacities(CCTK_PASS_CTOC);
              // Step 2.b.ii: Perform POLR algorithm on every grid component at this level
              NRPyLeakageET_optical_depths_PathOfLeastResistance(CCTK_PASS_CTOC);
            } END_COMPONENT_LOOP;
          } END_MAP_LOOP;
          CCTK_SyncGroup(cctkGH,"NRPyLeakageET::NRPyLeakageET_optical_depths");
        } LEAVE_LEVEL_MODE;
      }

      // Step 2.c: Now loop over the refinement levels backwards, restricting the result
      for(int rl=Carpet::reflevels-2;rl>=0;rl--) {
        ENTER_LEVEL_MODE(cctkGH,rl) {
          NRPyLeakageET_Restrict(CCTK_PASS_CTOC,"NRPyLeakageET::NRPyLeakageET_optical_depths");
        } LEAVE_LEVEL_MODE;
      }
      RemainingIterations--;
      counter++;
      if( verbose ) CCTK_VInfo(CCTK_THORNSTRING,"Completed iteration %d. Remaining iterations: %4d",counter,RemainingIterations);
    }
    if(verbose) CCTK_INFO("Completed path of least resistance algorithm!");
  }

  // Step 3: Now copy the optical depths and opacities to all time levels
  if(verbose) CCTK_INFO("Copying initial data to all time levels...");
  BEGIN_REFLEVEL_LOOP(cctkGH) {
    BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
      BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        NRPyLeakageET_copy_opacities_and_optical_depths_to_previous_time_levels(CCTK_PASS_CTOC);
      } END_COMPONENT_LOOP;
    } END_MAP_LOOP;
  } END_REFLEVEL_LOOP;
  if(verbose) CCTK_INFO("Finished copying initial data to all time levels");
}