#ifndef BASIC_DEFINES_HH
#define BASIC_DEFINES_HH

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include "hdf5.h"
#define H5_USE_16_API 1

#define REAL double

#ifdef __cplusplus
#define restrict __restrict__
#endif

// MAX and MIN macros
#ifdef MAX
#undef MAX
#endif

#ifdef MIN
#undef MIN
#endif

#define MAX(a,b) ( (a) > (b) ? (a) : (b) )
#define MIN(a,b) ( (a) < (b) ? (a) : (b) )

// nuc_eos stuff
#define HAVEGR 1
#define NTABLES 19
#define LENGTHGF 6.77269222552442e-06
#define TIMEGF 2.03040204956746e05
#define RHOGF 1.61887093132742e-18
#define PRESSGF 1.80123683248503e-39
#define EPSGF 1.11265005605362e-21
#define INVRHOGF 6.17714470405638e17
#define INVEPSGF 8.98755178736818e20
#define INVPRESSGF 5.55174079257738e38

#endif // BASIC_DEFINES_HH
