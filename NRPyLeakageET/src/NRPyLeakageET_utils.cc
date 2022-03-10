#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "carpet.hh"
#include "loopcontrol.h"
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

  if(verbosity_level>1) CCTK_VINFO("Prolongating %s from ref. lvl. %d to ref. lvl. %d...", varname, GetRefinementLevel(cctkGH)-1, GetRefinementLevel(cctkGH));

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

  if(verbosity_level>1) CCTK_VINFO("Restricting %s from ref. lvl. %d to ref. lvl. %d...", varname, GetRefinementLevel(cctkGH)+1, GetRefinementLevel(cctkGH));

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
void NRPyLeakageET_SyncAuxiliaryOpticalDepths(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if( verbosity_level > 1 ) CCTK_VINFO("Synchronizing auxiliary optical depths at ref. lev. %d",GetRefinementLevel(cctkGH));

  const int status = CCTK_SyncGroup(cctkGH,"NRPyLeakageET::NRPyLeakageET_auxiliary_optical_depths");
  if( status < 0 ) CCTK_VError(__LINE__,__FILE__,CCTK_THORNSTRING,"Could not synchronize NRPyLeakageET::NRPyLeakageET_auxiliary_optical_depths. ");

  if( verbosity_level > 1 ) CCTK_VINFO("Finished synchronizing auxiliary optical depths at ref. lev. %d",GetRefinementLevel(cctkGH));
}

extern "C"
void NRPyLeakageET_getGlobalL2norm(CCTK_ARGUMENTS, const int startRefLev, const int endRefLev, CCTK_REAL *restrict l2norm_global) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Step 1: Initialize changes to zero everywhere
  BEGIN_REFLEVEL_LOOP(cctkGH) {
    BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
      BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        NRPyLeakageET_initialize_optical_depths_changes_to_zero(CCTK_PASS_CTOC);
      } END_COMPONENT_LOOP;
    } END_MAP_LOOP;
  } END_REFLEVEL_LOOP;

  // Step 2: Loop over refinement levels, maps, and components, computing the
  // differences in the optical depths between two consecutive iterations
  ENTER_LEVEL_MODE(cctkGH,startRefLev) {
    BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
      BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        NRPyLeakageET_compute_optical_depth_change(CCTK_PASS_CTOC);
      } END_COMPONENT_LOOP;
    } END_MAP_LOOP;
  } LEAVE_LEVEL_MODE;

  for(int rl=startRefLev+1;rl<=endRefLev;rl++) {
    ENTER_LEVEL_MODE(cctkGH,rl) {
      BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
        BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
          NRPyLeakageET_compute_optical_depth_change(CCTK_PASS_CTOC);
        } END_COMPONENT_LOOP;
      } END_MAP_LOOP;
    } LEAVE_LEVEL_MODE;
  }

  // Step 3: Get sum operation handle
  const int op_sum = CCTK_ReductionHandle("sum");

  // Step 4: Set variables we want to reduce
  const char varnames[6][30] = {"NRPyLeakageET::tau_0_nue_aux ",
                                "NRPyLeakageET::tau_1_nue_aux ",
                                "NRPyLeakageET::tau_0_anue_aux",
                                "NRPyLeakageET::tau_1_anue_aux",
                                "NRPyLeakageET::tau_0_nux_aux ",
                                "NRPyLeakageET::tau_1_nux_aux "};


  // Step 5: Get indices for the variables we want to reduce
  const int varindices[6] = {CCTK_VarIndex("NRPyLeakageET::tau_0_nue_aux"),
                             CCTK_VarIndex("NRPyLeakageET::tau_1_nue_aux"),
                             CCTK_VarIndex("NRPyLeakageET::tau_0_anue_aux"),
                             CCTK_VarIndex("NRPyLeakageET::tau_1_anue_aux"),
                             CCTK_VarIndex("NRPyLeakageET::tau_0_nux_aux"),
                             CCTK_VarIndex("NRPyLeakageET::tau_1_nux_aux")};

  // Step 6: Compute the global L2-norm:
  // l2norm = sqrt(dTau_0_nue + dTau_0_anue + dTau_0_nux
  //             + dTau_1_nue + dTau_1_anue + dTau_1_nux)
  CCTK_REAL l2norm = 0.0;
  for(int i=0;i<6;i++) {
    int varindex = varindices[i];
    assert(varindex>=0);
    CCTK_REAL l2normL = 0.0;
    const int ierr = CCTK_Reduce(cctkGH,
                                 -1, // target processors; -1 -> all
                                 op_sum,
                                 1,  // Number of outputs
                                 CCTK_VARIABLE_REAL,
                                 &l2normL,
                                 1,  // Number of inputs
                                 varindex);
    if( ierr ) CCTK_VError(__LINE__,__FILE__,CCTK_THORNSTRING,"Error in reduction of L2-norm of %s",varnames[i]);

    l2norm += l2normL;

    if( verbosity_level > 1 ) CCTK_VINFO("L2-norm of %s: %e",varnames[i],sqrt(l2normL));
  }
  *l2norm_global = sqrt(l2norm);
}

extern "C"
void NRPyLeakageET_CopyOpticalDepthsToAux(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(!NRPyLeakageET_ProcessOwnsData()) return;

  // Step 1: Copy optical depths to auxiliary gridfunctions
#pragma omp parallel for
  //for(int k=cctk_nghostzones[2];k<cctk_lsh[2]-cctk_nghostzones[2];k++) {
  //  for(int j=cctk_nghostzones[1];j<cctk_lsh[1]-cctk_nghostzones[1];j++) {
  //    for(int i=cctk_nghostzones[0];i<cctk_lsh[0]-cctk_nghostzones[0];i++) {
  for(int k=0;k<cctk_lsh[2];k++) {
    for(int j=0;j<cctk_lsh[1];j++) {
      for(int i=0;i<cctk_lsh[0];i++) {
        // Step 1.a: Set local gridfunction index
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // Step 1.b: Copy optical depths to auxiliary gridfunctions
        tau_0_nue_aux [index] = tau_0_nue [index];
        tau_1_nue_aux [index] = tau_1_nue [index];
        tau_0_anue_aux[index] = tau_0_anue[index];
        tau_1_anue_aux[index] = tau_1_anue[index];
        tau_0_nux_aux [index] = tau_0_nux [index];
        tau_1_nux_aux [index] = tau_1_nux [index];
      }
    }
  }
}

extern "C"
void NRPyLeakageET_CopyOpticalDepthsFromAux(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(!NRPyLeakageET_ProcessOwnsData()) return;

  // Step 1: Copy optical depths from auxiliary gridfunctions
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) {
    for(int j=0;j<cctk_lsh[1];j++) {
      for(int i=0;i<cctk_lsh[0];i++) {
        // Step 1.a: Set local gridfunction index
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // Step 1.b: Copy optical depths from auxiliary gridfunctions
        tau_0_nue [index] = tau_0_nue_aux [index];
        tau_1_nue [index] = tau_1_nue_aux [index];
        tau_0_anue[index] = tau_0_anue_aux[index];
        tau_1_anue[index] = tau_1_anue_aux[index];
        tau_0_nux [index] = tau_0_nux_aux [index];
        tau_1_nux [index] = tau_1_nux_aux [index];
      }
    }
  }
}

extern "C"
void NRPyLeakageET_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Step 1: Initialize all optical depth gridfunctions to zero
  if(verbosity_level>0) CCTK_INFO("Initializing optical depths gridfunctions to zero...");
  BEGIN_REFLEVEL_LOOP(cctkGH) {
    BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
      BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        NRPyLeakageET_optical_depths_initialize_to_zero(CCTK_PASS_CTOC);
      } END_COMPONENT_LOOP;
    } END_MAP_LOOP;
  } END_REFLEVEL_LOOP;
  if(verbosity_level>0) CCTK_INFO("Initialized all optical depths gridfunctions to zero");

  if( CCTK_EQUALS(optical_depth_evolution_type,"PathOfLeastResistance") ) {

    const int startRefLev = MIN(MAX(minInitRefLevel,0),Carpet::reflevels-1);
    const int endRefLev   = maxInitRefLevel == 0 ? Carpet::reflevels-1 : MIN(maxInitRefLevel,Carpet::reflevels-1);

    if(verbosity_level>0) {
      CCTK_VINFO("Optical depths will be initialized on levels %d through %d with the path of least resistance algorithm",startRefLev,endRefLev);
      CCTK_VINFO("Number of iterations to be performed on each refinement level: %d",numberOfIterations);
    }

    // Step 2: Now perform iterations of the path of least resistance algorithm
    int counter = 0;
    int RemainingIterations = numberOfIterations;

    for(int i=1;i<=numberOfIterations;i++) {
      // Step 2.a: First perform the iteration on refinement level startRefLev
      if( verbosity_level>1 ) CCTK_VINFO("Beginning iteration %d...",counter+1);
      ENTER_LEVEL_MODE(cctkGH,startRefLev) {
        BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
          BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
            NRPyLeakageET_compute_opacities(CCTK_PASS_CTOC);
            NRPyLeakageET_optical_depths_PathOfLeastResistance(CCTK_PASS_CTOC);
            NRPyLeakageET_CopyOpticalDepthsToAux(CCTK_PASS_CTOC);
          } END_COMPONENT_LOOP;
        } END_MAP_LOOP;
        // Step 2.a.i: Synchronize *auxiliary* gridfunctions, for which
        // prolongation is disabled. This means we only *copy* data from one
        // processor to the next, which is what we want.
        NRPyLeakageET_SyncAuxiliaryOpticalDepths(CCTK_PASS_CTOC);
        // Step 2.a.ii: Copy the synchronized data back to the main gridfunctions
        BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
          BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
            NRPyLeakageET_CopyOpticalDepthsFromAux(CCTK_PASS_CTOC);
          } END_COMPONENT_LOOP;
        } END_MAP_LOOP;
      } LEAVE_LEVEL_MODE;

      // Step 2.b: Now loop over remaining refinement levels, performing the
      //           path of least resistance (POLR) algorithm
      for(int rl=startRefLev+1;rl<=endRefLev;rl++) {
        ENTER_LEVEL_MODE(cctkGH,rl) {
          // Step 2.b.i: Prolongate results from the previous refinement level
          NRPyLeakageET_Prolongate(CCTK_PASS_CTOC,"NRPyLeakageET::NRPyLeakageET_optical_depths");
          BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
            BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
              // Step 2.b.i: Compute opacities
              NRPyLeakageET_compute_opacities(CCTK_PASS_CTOC);
              // Step 2.b.ii: Perform POLR algorithm on every grid component at this level
              NRPyLeakageET_optical_depths_PathOfLeastResistance(CCTK_PASS_CTOC);
              NRPyLeakageET_CopyOpticalDepthsToAux(CCTK_PASS_CTOC);
            } END_COMPONENT_LOOP;
          } END_MAP_LOOP;
          // Step 2.b.iii: Synchronize *auxiliary* gridfunctions, for which
          // prolongation is disabled. This means we only *copy* data from one
          // processor to the next, which is what we want.
          NRPyLeakageET_SyncAuxiliaryOpticalDepths(CCTK_PASS_CTOC);
          // Step 2.b.iv: Copy the synchronized data back to the main gridfunctions
          BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
            BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
              NRPyLeakageET_CopyOpticalDepthsFromAux(CCTK_PASS_CTOC);
            } END_COMPONENT_LOOP;
          } END_MAP_LOOP;
        } LEAVE_LEVEL_MODE;
      }

      // Step 2.c: Perform a reduction of the L2-norm among all processors
      CCTK_REAL l2norm_global;
      NRPyLeakageET_getGlobalL2norm(CCTK_PASS_CTOC, startRefLev, endRefLev, &l2norm_global);

      // Step 2.d: Add up L2-norm on all refinement levels that we have initialized
      if( verbosity_level>1 ) CCTK_VINFO("Global L2-norm: %e",l2norm_global);

      // Step 2.e: Increment counter; check convergence
      counter++;
      if( l2norm_global < tauChangeThreshold ) {
        // Step 2.e.i: Converged!
        RemainingIterations = 0;
        i = numberOfIterations*10;
        if( verbosity_level>0 ) CCTK_VINFO("Algorithm converged! l2norm_tau_change = %e (threshold = %e)",l2norm_global,tauChangeThreshold);
      }
      else {
        // Step 2.e.ii: Have not converged. Restrict current solution from the
        // finest to the coarsest initialized refinement levels; Decrement
        // number of iterations remaining; continue.
        if( RemainingIterations > 1 ) {
          for(int rl=endRefLev-1;rl>=startRefLev;rl--) {
            ENTER_LEVEL_MODE(cctkGH,rl) {
              NRPyLeakageET_Restrict(CCTK_PASS_CTOC,"NRPyLeakageET::NRPyLeakageET_optical_depths");
            } LEAVE_LEVEL_MODE;
          }
          RemainingIterations--;
          if( verbosity_level>1 ) CCTK_VINFO("Completed iteration %d. Remaining iterations: %4d",counter,RemainingIterations);
        }
        else {
          CCTK_WARN(CCTK_WARN_ALERT,"Maximum iterations reached, but convergence threshold not met.");
        }
      }
    }

    // Step 3: All done! Now loop over the refinement levels backwards, restricting the result
    for(int rl=endRefLev-1;rl>=0;rl--) {
      ENTER_LEVEL_MODE(cctkGH,rl) {
        NRPyLeakageET_Restrict(CCTK_PASS_CTOC,"NRPyLeakageET::NRPyLeakageET_optical_depths");
      } LEAVE_LEVEL_MODE;
    }
    if(verbosity_level>0) CCTK_INFO("Completed path of least resistance algorithm!");
  }

  // Step 4: Now copy the optical depths and opacities to all time levels
  if(verbosity_level>0) CCTK_INFO("Copying initial data to all time levels...");
  BEGIN_REFLEVEL_LOOP(cctkGH) {
    BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
      BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        NRPyLeakageET_copy_opacities_and_optical_depths_to_previous_time_levels(CCTK_PASS_CTOC);
      } END_COMPONENT_LOOP;
    } END_MAP_LOOP;
  } END_REFLEVEL_LOOP;
  if(verbosity_level>0) CCTK_INFO("Finished copying initial data to all time levels");
}
