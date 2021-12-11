#include "Basic_defines.hh"
#include "WVU_EOS_Tabulated_headers.hh"

extern "C"
void WVU_EOS_free_memory() {

  // Free memory allocated for the table
  free(nuc_eos_private::logrho);
  free(nuc_eos_private::logtemp);
  free(nuc_eos_private::yes);
  free(nuc_eos_private::alltables);
  free(nuc_eos_private::epstable);

}
