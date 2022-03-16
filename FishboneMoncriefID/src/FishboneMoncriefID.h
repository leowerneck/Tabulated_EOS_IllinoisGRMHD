#ifndef FISHBONEMONCRIEFID_H_
#define FISHBONEMONCRIEFID_H_

#include <stdio.h>
#include <stdlib.h> // Needed for rand()
#include <stdbool.h>
#include <math.h>

// Alias for "vel" vector gridfunction:
#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

#define POLYTROPE_EOS 0
#define TABULATED_EOS 1

#ifdef __cplusplus
extern "C" {
#endif

void FishboneMoncriefID_velocities(const cGH* restrict const cctkGH,const CCTK_INT *cctk_lsh,
                                   const CCTK_INT i0,const CCTK_INT i1,const CCTK_INT i2,
                                   const CCTK_REAL *xcoordGF,const CCTK_REAL *ycoordGF,const CCTK_REAL *zcoordGF,
                                   CCTK_REAL *Valencia3velocityU0GF, CCTK_REAL *Valencia3velocityU1GF, CCTK_REAL *Valencia3velocityU2GF);

void FishboneMoncriefID_KerrSchild(const cGH* restrict const cctkGH,const CCTK_INT *cctk_lsh,
                                   const CCTK_INT i0,const CCTK_INT i1,const CCTK_INT i2,
                                   const CCTK_REAL *xcoordGF,const CCTK_REAL *ycoordGF,const CCTK_REAL *zcoordGF,
                                   CCTK_REAL *alphaGF,CCTK_REAL *betaU0GF,CCTK_REAL *betaU1GF,CCTK_REAL *betaU2GF,
                                   CCTK_REAL *gammaDD00GF,CCTK_REAL *gammaDD01GF,CCTK_REAL *gammaDD02GF,CCTK_REAL *gammaDD11GF,CCTK_REAL *gammaDD12GF,CCTK_REAL *gammaDD22GF,
                                   CCTK_REAL     *KDD00GF,CCTK_REAL     *KDD01GF,CCTK_REAL     *KDD02GF,CCTK_REAL     *KDD11GF,CCTK_REAL     *KDD12GF,CCTK_REAL     *KDD22GF);

void FishboneMoncriefID_initial_polytrope_eos(CCTK_ARGUMENTS);
void FishboneMoncriefID_initial_tabulated_eos(CCTK_ARGUMENTS);

#ifdef __cplusplus
}
#endif

#endif // FISHBONEMONCRIEFID_H_
