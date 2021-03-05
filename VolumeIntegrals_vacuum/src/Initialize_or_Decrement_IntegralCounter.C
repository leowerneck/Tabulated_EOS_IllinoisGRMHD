#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "util_Table.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>

void InitializeIntegralCounterToZero(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  *IntegralCounter = 0;

  if(verbose==2) printf("VolumeIntegrals: Just set IntegralCounter to %d\n",*IntegralCounter);
}

void InitializeIntegralCounter(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(cctk_iteration%VolIntegral_out_every==0) {
    *IntegralCounter = NumIntegrals;
    if(verbose==2) printf("VolumeIntegrals: Just set IntegralCounter to %d == NumIntegrals\n",*IntegralCounter);
  }
}

void DecrementIntegralCounter(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

    (*IntegralCounter)--;
    if(verbose==2) printf("VolumeIntegrals: Just decremented IntegralCounter to %d\n",*IntegralCounter);
}
