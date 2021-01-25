#include "cctk.h"

#include "IllinoisGRMHD_headers.h"
#include "con2prim_headers.h"

CCTK_INT con2prim_get_key( const char* routine_name ) {

  // Get con2prim key
  if( CCTK_EQUALS(routine_name,"None") ) {
    return None;
  }
  else if( CCTK_EQUALS(routine_name,"Noble2D") ) {
    return Noble2D;
  }
  else if( CCTK_EQUALS(routine_name,"Noble1D") ) {
    return Noble1D;
  }
  else if( CCTK_EQUALS(routine_name,"Noble1D_entropy") ) {
    return Noble1D_entropy;
  }
  else if( CCTK_EQUALS(routine_name,"Noble1D_entropy2") ) {
    return Noble1D_entropy2;
  }
  else if( CCTK_EQUALS(routine_name,"Palenzuela1D") ) {
    return Palenzuela1D;
  }
  else if( CCTK_EQUALS(routine_name,"Palenzuela1D_entropy") ) {
    return Palenzuela1D_entropy;
  }
  else {
    CCTK_VError(VERR_DEF_PARAMS,"Unknown con2prim routine: %s. Please check your parameter file. ABORTING!",routine_name);
  }

}
