#include "cctk.h"
#include "cctk_Parameters.h"

#include "IllinoisGRMHD_headers.h"
#include "con2prim_headers.h"

void con2prim_select( int (*con2prim)( const igm_eos_parameters,
                                       const CCTK_REAL[4][4],const CCTK_REAL[4][4],
                                       CCTK_REAL *restrict,CCTK_REAL *restrict ) ) {

  DECLARE_CCTK_PARAMETERS;
  
  if( CCTK_EQUALS(igm_con2prim_routine,"Noble2D") ) {
    con2prim = &con2prim_Noble2D;
  }
  else if( CCTK_EQUALS(igm_con2prim_routine,"Palenzuela1D") ) {
    con2prim = &con2prim_Palenzuela1D;
  }
  else if( CCTK_EQUALS(igm_con2prim_routine,"Palenzuela1D_entropy") ) {
    con2prim = &con2prim_Palenzuela1D_entropy;
  }
  else {
    CCTK_VError(VERR_DEF_PARAMS,"Unknown con2prim method: %s. ABORTING!",igm_con2prim_routine);
  }

}
