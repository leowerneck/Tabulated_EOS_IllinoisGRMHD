#include "Basic_defines.hh"
#include "WVU_EOS_Tabulated_headers.hh"

int main(int argc, char **argv) {

  WVU_EOS_ReadTable(argv[1]);
  WVU_EOS_free_memory();

  return 0;
}
