#include <stdio.h>
#include <cassert>
#include <cmath>
#include <mpi.h>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>
#include <cctk_Functions.h>
#include <cctk_Faces.h>
#include <Symmetry.h>
#include <loopcontrol.h>
#include <time.h>
#include "tau.hh"

#define SQR(a) ((a)*(a))

extern "C" {
  void ZLtau_calc_grr(CCTK_ARGUMENTS);
}


void ZLtau_calc_grr(CCTK_ARGUMENTS) {

  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  //  CCTK_Info(CCTK_THORNSTRING,"Calculating grr");

  // most code copied from ADMAnalysis
  /* loop over the grid */
#pragma omp parallel for
  for(int i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
    {
      CCTK_REAL cost;
      CCTK_REAL sint;
      CCTK_REAL sinp;
      CCTK_REAL cosp;
      CCTK_REAL rvalue;
      CCTK_REAL sxy;

      rvalue  = r[i];
      sxy = sqrt( SQR(x[i]) + SQR(y[i]));
      
      //  be careful with r=0 and xy plane 
      if (rvalue <= 1.0e-14)
	{
	  cost = 1.0;
	  sint = 0.0;
	  sinp = 0.0;
	  cosp = 1.0;
	}
      else if (sxy==0)
	{
	  cost = 1.0;
	  sint = 0.0;
	  sinp = 0.0;
	  cosp = 1.0;
	}
      else
	{
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
