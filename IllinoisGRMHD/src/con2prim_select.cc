#include "cctk.h"
#include "cctk_Parameters.h"

#include "IllinoisGRMHD_headers.h"
#include "con2prim_headers.h"
#include "EOS_headers.hh"

int con2prim_select( igm_eos_parameters *restrict eos,
                     const CCTK_INT c2p_key,
                     const CCTK_REAL *restrict adm_quantities,
                     const CCTK_REAL g4dn[4][4],
                     const CCTK_REAL g4up[4][4],
                     const CCTK_REAL *restrict cons,
                     CCTK_REAL *restrict prim ) {

  switch( c2p_key ) {

    // Noble2D routine (see https://arxiv.org/pdf/astro-ph/0512420.pdf)
    case Noble2D:
      return( con2prim_Noble2D(*eos,g4dn,g4up,cons,prim) );
      break;

    // Noble1D routine (see https://arxiv.org/pdf/astro-ph/0512420.pdf)
    case Noble1D:
      return( con2prim_Noble1D(*eos,g4dn,g4up,cons,prim) );
      break;

    // Noble1D entropy routine (see https://arxiv.org/pdf/0808.3140.pdf)
    case Noble1D_entropy:
      return( con2prim_Noble1D_entropy(*eos,g4dn,g4up,cons,prim) );
      break;

    // Noble1D entropy 2 routine (see https://arxiv.org/pdf/0808.3140.pdf)
    case Noble1D_entropy2:
      return( con2prim_Noble1D_entropy2(*eos,g4dn,g4up,cons,prim) );
      break;

    // Cerda-Duran et al. 2D routine (see https://arxiv.org/pdf/0804.4572.pdf)
    case CerdaDuran2D:
      return( con2prim_CerdaDuran2D(*eos,adm_quantities,cons,prim) );

    // Palenzuela 1D routine (see https://arxiv.org/pdf/1712.07538.pdf)
    case Palenzuela1D:
      return( con2prim_Palenzuela1D(eos,adm_quantities,cons,prim) );
      break;

    default:
      CCTK_VError(VERR_DEF_PARAMS,"Unknown c2p key in con2prim_select (%d). ABORTING!",c2p_key);
      break;

  }
}
