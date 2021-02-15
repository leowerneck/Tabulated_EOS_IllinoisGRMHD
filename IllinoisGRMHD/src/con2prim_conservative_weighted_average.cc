#include "cctk.h"
#include "cctk_Parameters.h"

#include "IllinoisGRMHD_headers.h"

void con2prim_conservative_weighted_average( const CCTK_INT i,
                                             const CCTK_INT j,
                                             const CCTK_INT k,
                                             const CCTK_INT index,
                                             const cGH  *restrict cctkGH,
                                             CCTK_INT& num_of_conservative_averagings_needed,
                                             unsigned short *restrict con2prim_failed_flag,
                                             CCTK_REAL *restrict rho_star,
                                             CCTK_REAL *restrict mhd_st_x,
                                             CCTK_REAL *restrict mhd_st_y,
                                             CCTK_REAL *restrict mhd_st_z,
                                             CCTK_REAL *restrict tau,
                                             CCTK_REAL *restrict Ye_star,
                                             CCTK_REAL *restrict S_star ) {

}
