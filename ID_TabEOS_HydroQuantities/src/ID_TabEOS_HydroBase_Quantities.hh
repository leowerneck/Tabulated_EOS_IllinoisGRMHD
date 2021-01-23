#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"

//------ EOS_Omni stuff ------
namespace nuc_eos {
  extern double eos_yemin, eos_yemax;
}

namespace nuc_eos_private {
  extern int nrho;
  extern double *restrict logrho;
}
//----------------------------

void read_1dfile__set_array(FILE *in1D,
                            CCTK_REAL *f_of_rho_arr,
                            const int num_header_lines=0)
{
  char* line = NULL;
  size_t len = 0;
  ssize_t read;
  int which_line = 0;
  // Skip header line
  while( (read = getline(&line, &len, in1D)) != -1 ) {
    if( which_line >= num_header_lines ) f_of_rho_arr[which_line] = strtod(line, NULL);
    which_line++;
  }
  free(line);
}

// Find interpolation index using Bisection root-finding algorithm:
static inline int bisection_idx_finder(const CCTK_REAL rrbar, const int numlines_in_file, const CCTK_REAL *restrict rbar_arr) {
  int x1 = 0;
  int x2 = numlines_in_file-1;
  CCTK_REAL y1 = rrbar-rbar_arr[x1];
  CCTK_REAL y2 = rrbar-rbar_arr[x2];
  if(y1*y2 >= 0) {
    fprintf(stderr,"INTERPOLATION BRACKETING ERROR %e | %e %e\n",rrbar,y1,y2);
    exit(1);
  }
  for(int i=0;i<numlines_in_file;i++) {
    int x_midpoint = (x1+x2)/2;
    CCTK_REAL y_midpoint = rrbar-rbar_arr[x_midpoint];
    if(y_midpoint*y1 < 0) {
      x2 = x_midpoint;
      y2 = y_midpoint;
    } else {
      x1 = x_midpoint;
      y1 = y_midpoint;
    }
    if( abs(x2-x1) == 1 ) {
      // If rbar_arr[x1] is closer to rrbar than rbar_arr[x2] then return x1:
      // if(fabs(rrbar-rbar_arr[x1]) < fabs(rrbar-rbar_arr[x2])) return x1;
      // Otherwiser return x2:
      // return x2;
      // Always return the left value
      return x1;
    }
  }
  fprintf(stderr,"INTERPOLATION BRACKETING ERROR: DID NOT CONVERGE.\n");
  exit(1);
}

// Reference: https://en.wikipedia.org/wiki/Linear_interpolation
inline void linear_interp_1D( const int nx,
                              const CCTK_REAL *restrict x_arr,
                              const CCTK_REAL *restrict y_arr,
                              const CCTK_REAL x_star,
                              CCTK_REAL *restrict y_star ) {
  // Find the index to the left
  int idx = bisection_idx_finder(x_star,nx,x_arr);
  // Compute (x0,y0) and (x1,y1)
  const CCTK_REAL x0 = x_arr[idx  ];
  const CCTK_REAL y0 = y_arr[idx  ];
  const CCTK_REAL x1 = x_arr[idx+1];
  const CCTK_REAL y1 = y_arr[idx+1];
  // Perform the interpolation. Note that:
  //
  // y_star = y0 + (x_star-x0)*(y1-y0)/(x1-x0)
  //        = y0 - y0*(x_star-x0)/(x1-x0) + y1*(x_star-x0)/(x1-x0)
  //    .--------------------------------.
  // => | y_star = y0*(1 - aux) + y1*aux | ,
  //    .--------------------------------.
  // where
  //    .---------------------------.
  //    | aux = (x_star-x0)/(x1-x0) | .
  //    .---------------------------.
  const CCTK_REAL aux = (x_star-x0)/(x1-x0);
  *y_star = y0*(1.0-aux) + y1*aux;
}
                              

void interpolate_1d_quantity_as_function_of_rho(const int interp_stencil_size,
                                                const int numlines_in_file,
                                                const CCTK_REAL *restrict rho_arr,
                                                const CCTK_REAL *restrict f_of_rho_arr,
                                                const CCTK_REAL rho,
                                                CCTK_REAL *restrict f_of_rho)
{

  // First find the central interpolation stencil index:
  int idx = bisection_idx_finder(rho,numlines_in_file,rho_arr);


#ifdef MAX
#undef MAX
#endif
#define MAX(A, B) ( ((A) > (B)) ? (A) : (B) )

#ifdef MIN
#undef MIN
#endif
#define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
  
  int idxmin = MAX(0,idx-interp_stencil_size/2-1);
  idxmin     = MIN(idxmin,numlines_in_file - interp_stencil_size + 1);

  // Now perform the Lagrange polynomial interpolation:

  // First set the interpolation coefficients:
  CCTK_REAL rho_sample[interp_stencil_size];
  for(int i=idxmin;i<idxmin+interp_stencil_size;i++) {
    rho_sample[i-idxmin] = rho_arr[i];
  }
  CCTK_REAL l_i_of_r[interp_stencil_size];
  for(int i=0;i<interp_stencil_size;i++) {
    CCTK_REAL numer = 1.0;
    CCTK_REAL denom = 1.0;
    for(int j=0;j<i;j++) {
      numer *= rho - rho_sample[j];
      denom *= rho_sample[i] - rho_sample[j];
    }
    for(int j=i+1;j<interp_stencil_size;j++) {
      numer *= rho - rho_sample[j];
      denom *= rho_sample[i] - rho_sample[j];
    }
    l_i_of_r[i] = numer/denom;
  }

  // Then perform the interpolation:
  *f_of_rho = 0.0;
  for(int i=idxmin;i<idxmin+interp_stencil_size;i++) {
    *f_of_rho += l_i_of_r[i-idxmin] * f_of_rho_arr[i];
  }
}
