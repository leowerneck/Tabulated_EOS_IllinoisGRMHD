#include "cctk.h"
#include "cctk_Parameters.h"

#include "IllinoisGRMHD_headers.h"
#include "con2prim_headers.h"
#include "EOS_headers.hh"

int con2prim_select( const igm_eos_parameters eos,
                     const CCTK_INT c2p_key,
                     const CCTK_REAL *restrict adm_quantities,
                     const CCTK_REAL g4dn[4][4],
                     const CCTK_REAL g4up[4][4],
                     const CCTK_REAL *restrict cons,
                     CCTK_REAL *restrict prim,
                     output_stats& stats ) {

  switch( c2p_key ) {

    // Palenzuela 1D routine (see https://arxiv.org/pdf/1712.07538.pdf)
    case old_Palenzuela1D:
      return( con2prim_Palenzuela1D(eos,adm_quantities,cons,prim,stats) );
      break;

    // Newman 1D routine (see https://escholarship.org/content/qt0s53f84b/qt0s53f84b.pdf)
    case old_Newman1D:
      return( con2prim_Newman1D(eos,adm_quantities,cons,prim,stats) );
      break;

    default:
      CCTK_VError(VERR_DEF_PARAMS,"Unknown c2p key in con2prim_select (%d). ABORTING!",c2p_key);
      break;

  }
}
