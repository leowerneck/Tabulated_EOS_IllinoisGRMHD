#include "Basic_defines.h"
#include "NRPyEOS.h"

void NRPyEOS_free_memory(NRPyEOS_params *restrict eos_params) {

  // Free memory allocated for the table
  free(eos_params->logrho);
  free(eos_params->logtemp);
  free(eos_params->yes);
  free(eos_params->alltables);
  free(eos_params->epstable);

}
