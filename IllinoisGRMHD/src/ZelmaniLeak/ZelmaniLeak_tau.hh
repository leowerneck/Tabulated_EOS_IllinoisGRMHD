#ifndef TAU_HH
#define TAU_HH

#define PRESS_GF      (1.7982953469278e-39)
#define RHO_GF        (1.61620075314614e-18)
#define EPS_GF        (1.11265006e-21)
#define LENGTH_GF     (6.77140812e-06)
#define INV_RHO_GF    (6.18735016707159e17)
#define INV_PRESS_GF  (5.5608218177753538)
#define INV_LENGTH_GF (147679.77092481e0)

#include <cctk.h>
#include <cctk_Parameters.h>

namespace ZelmaniLeak {
 
  inline void spher2cart( const CCTK_REAL rad,
                          const CCTK_REAL theta,
                          const CCTK_REAL phi,
                          const CCTK_REAL *restrict x0,
                          const CCTK_REAL *restrict y0,
                          const CCTK_REAL *restrict z0,
                          CCTK_REAL *restrict x,
                          CCTK_REAL *restrict y,
                          CCTK_REAL *restrict z ) {

    DECLARE_CCTK_PARAMETERS;

#ifdef SYMMETRIC_OPERATORS
    // phi == 0 cannot be decided!
    *x = *x0 + rad*std::cos(theta)*std::cos(phi);
    *y = *y0 + rad*std::cos(theta)*std::sin(phi);
    *z = *z0 + rad*std::sin(theta);
#else
    *x = *x0 + rad*std::sin(theta)*std::cos(phi);
    *y = *y0 + rad*std::sin(theta)*std::sin(phi);
    *z = *z0 + rad*std::cos(theta);
#endif
  }
  
} // namespace ZelmaniLeak

#endif // TAU_HH
