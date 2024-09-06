#include "cctk.h"

#include "IllinoisGRMHD_headers.h"
#include "con2prim_headers.h"

CCTK_INT con2prim_get_key( const char* routine_name ) {

  // Get con2prim key
  if( CCTK_EQUALS(routine_name,"None") ) {
    return igm_None;
  }
  else if( CCTK_EQUALS(routine_name,"Noble2D") ) {
    return igm_Noble2D;
  }
  else if( CCTK_EQUALS(routine_name,"Noble1D") ) {
    return igm_Noble1D;
  }
  else if( CCTK_EQUALS(routine_name,"Noble1D_entropy") ) {
    return igm_Noble1D_entropy;
  }
  else if( CCTK_EQUALS(routine_name,"Noble1D_entropy2") ) {
    return igm_Noble1D_entropy2;
  }
  else if( CCTK_EQUALS(routine_name,"CerdaDuran2D") ) {
    return igm_CerdaDuran2D;
  }
  else if( CCTK_EQUALS(routine_name,"CerdaDuran3D") ) {
    return igm_CerdaDuran3D;
  }
  else if( CCTK_EQUALS(routine_name,"Palenzuela1D") ) {
    return igm_Palenzuela1D;
  }
  else if( CCTK_EQUALS(routine_name,"Newman1D") ) {
    return igm_Newman1D;
  }
  else {
    CCTK_VERROR("Unknown con2prim routine: %s. Please check your parameter file. ABORTING!",routine_name);
  }

}
