
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#ifdef INTEL_HEADER
#include <mathimf.h>
#endif

#ifdef USE_GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>
#endif

#include "options.h"
#include "macros.h"
#include "globals.h"

