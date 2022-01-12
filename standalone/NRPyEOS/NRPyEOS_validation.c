#include "NRPyEOS.h"

int main(int argc, char **argv) {

  if( argc != 2 ) {
    fprintf(stderr,"(NRPyEOS_validation) Correct usage is ./NRPyEOS_validation eos_file_path\n");
    exit(1);
  }
  
  NRPyEOS_params eos_params;
  NRPyEOS_readtable_set_EOS_params(argv[1],&eos_params);
  NRPyEOS_free_memory(&eos_params);

  return 0;
}
