#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "ZelmaniLeak_tau.hh"

#define SQR(a) ((a)*(a))

extern "C" void ZelmaniLeak_calc_grr(CCTK_ARGUMENTS);

namespace ZelmaniLeak {

  extern "C" void ZelmaniLeak_calc_grr(CCTK_ARGUMENTS) {

    DECLARE_CCTK_PARAMETERS;
    DECLARE_CCTK_ARGUMENTS;

    // most code copied from ADMAnalysis
    /* loop over the grid */
    const CCTK_INT npoints = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];

#pragma omp parallel for
    for(int i=0;i<npoints;i++) {
      CCTK_REAL rvalue = r[i];
      CCTK_REAL sxy    = sqrt( SQR(x[i]) + SQR(y[i]));

      //  be careful with r=0 and xy plane
      CCTK_REAL sint, cost, sinp, cosp;
      if( rvalue <= 1.0e-14 ) {
        cost = 1.0;
        sint = 0.0;
        sinp = 0.0;
        cosp = 1.0;
      }
      else if( sxy == 0 ) {
        cost = 1.0;
        sint = 0.0;
        sinp = 0.0;
        cosp = 1.0;
      }
      else {
        cost = z[i]/rvalue;
        sint = sxy/rvalue;
        sinp = y[i]/sxy;
        cosp = x[i]/sxy;
      }

      ZLtau_grr[i]=
        gyy[i]*SQR(sinp)*SQR(sint)+
        2*cosp*gxy[i]*sinp*SQR(sint)+
        SQR(cosp)*gxx[i]*SQR(sint)+
        2*cost*gyz[i]*sinp*sint+
        2*cosp*cost*gxz[i]*sint+
        SQR(cost)*gzz[i];

    } // end loop

  }

} // namespace ZelmaniLeak
