#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#define RHOGF 1.61620075314614e-18
#define INVRHOGF 6.18735016707159e17

void CCTK_FCALL
CCTK_FNAME(ZelmaniLeak_CalcPnu)( CCTK_REAL      * taurhs,
				 CCTK_REAL      * sxrhs,
				 CCTK_REAL      * syrhs,
				 CCTK_REAL      * szrhs,
                                 CCTK_REAL      * rho,
                                 CCTK_REAL      * eps,
                                 CCTK_REAL      * temperature,
			         CCTK_REAL      * entropy,
                                 CCTK_REAL      * Y_e,
                                 CCTK_REAL      * munu,
			         CCTK_REAL      * pnu,
                                 CCTK_REAL      * wlorentz,
                                 CCTK_REAL      * vel,
			         CCTK_REAL      * alp,
			         CCTK_REAL      * gxx,
                                 CCTK_REAL      * gxy,
                                 CCTK_REAL      * gxz,
                                 CCTK_REAL      * gyy,
                                 CCTK_REAL      * gyz,
                                 CCTK_REAL      * gzz,
                                 CCTK_REAL      * x,
                                 CCTK_REAL      * y,
                                 CCTK_REAL      * z,
			         CCTK_REAL      * dx,
			         CCTK_REAL      * dy,
			         CCTK_REAL      * dz,
                                 CCTK_INT const * const nx,
                                 CCTK_INT const * const ny,
                                 CCTK_INT const * const nz,
			         CCTK_REAL      * dtime );


// This is a c++ wrapper around our beautiful
// F90 code. We use it to get the pointers to
// the rhs variables
void ZelmaniLeak_NeutrinoPressureWrap(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT varindex,rhsindex;
  CCTK_REAL* taurhsptr;
  CCTK_REAL* sxrhsptr;
  CCTK_REAL* syrhsptr;
  CCTK_REAL* szrhsptr;
  CCTK_REAL dtime = CCTK_DELTA_TIME;
  CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  if(*global_rho_max*INVRHOGF < pnu_rho_start) return;

  /* Original:
  varindex = CCTK_VarIndex("GRHydro::tau"); */
  varindex = CCTK_VarIndex("IllinoisGRMHD::tau");
  rhsindex = MoLQueryEvolvedRHS (varindex);
  taurhsptr = (CCTK_REAL*) CCTK_VarDataPtrI (cctkGH, 0, rhsindex);
  assert (taurhsptr);

  /* Original:
   varindex = CCTK_VarIndex("GRHydro::scon[0]"); */
  varindex = CCTK_VarIndex("IllinoisGRMHD::mhd_st_x");
  rhsindex = MoLQueryEvolvedRHS (varindex);
  sxrhsptr = (CCTK_REAL*) CCTK_VarDataPtrI(cctkGH, 0, rhsindex);
  assert (sxrhsptr);

  /* Original:
   varindex = CCTK_VarIndex("GRHydro::scon[1]"); */
  varindex = CCTK_VarIndex("IllinoisGRMHD::mhd_st_y");
  rhsindex = MoLQueryEvolvedRHS (varindex);
  syrhsptr = (CCTK_REAL*) CCTK_VarDataPtrI(cctkGH, 0, rhsindex);
  assert (syrhsptr);

  /* Original:
   varindex = CCTK_VarIndex("GRHydro::scon[2]"); */
  varindex = CCTK_VarIndex("IllinoisGRMHD::mhd_st_z");
  rhsindex = MoLQueryEvolvedRHS (varindex);
  szrhsptr = (CCTK_REAL*) CCTK_VarDataPtrI(cctkGH, 0, rhsindex);
  assert (szrhsptr);

  // Leo says: We don't need to modify this function call.
  //           However, we must have computed both ADMBase
  //           and HydroBase quantities from Baikal/IllinoisGRMHD
  //           before calling this function.
  CCTK_FNAME(ZelmaniLeak_CalcPnu)
    (taurhsptr,sxrhsptr,syrhsptr,szrhsptr,
     rho,eps,temperature,entropy,Y_e,munu,pnu,
     w_lorentz,vel,alp,
     gxx,gxy,gxz,gyy,gyz,gzz,
     x,y,z,&dx,&dy,&dz,
     &cctk_lsh[0], &cctk_lsh[1], &cctk_lsh[2], &dtime);

}
