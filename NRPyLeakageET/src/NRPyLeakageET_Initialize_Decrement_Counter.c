#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "NRPyLeakageET.h"

void NRPyLeakageET_InitializeIterationCounter(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Initialize iteration counter
  *NRPyLeakageET_IterationCounter = NRPyLeakageET_IterationFactor*MAX(cctk_lsh[0],MAX(cctk_lsh[1],cctk_lsh[2]));
}

void NRPyLeakageET_DecrementIterationCounter(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Decrement iteration counter
  (*NRPyLeakageET_IterationCounter)--;
}
