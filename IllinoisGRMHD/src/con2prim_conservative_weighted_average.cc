// Thorn      : IllinoisGRMHD
// File       : EOS_Tabulated.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
//              Zachariah B. Etienne
// Description: In this file we implement a function to compute
//              the weighted average of the conservative variables
//              around a point (i,j,k). It is intended to be used
//              after a con2prim failure has ocurred at the point
//              (i,j,k), and the weighted average should bring the
//              conservative variables at that point back within a
//              range in which the con2prim routines would work.

#include "cctk.h"
#include "cctk_Parameters.h"

#include "IllinoisGRMHD_headers.h"

void con2prim_average_neighbor_conservatives( const CCTK_INT i,
                                              const CCTK_INT j,
                                              const CCTK_INT k,
                                              const CCTK_INT index,
                                              const CCTK_INT *restrict cctk_lsh,
                                              const cGH  *restrict cctkGH,
                                              unsigned short *restrict con2prim_failed_flag,
                                              CCTK_REAL *restrict rho_star,
                                              CCTK_REAL *restrict mhd_st_x,
                                              CCTK_REAL *restrict mhd_st_y,
                                              CCTK_REAL *restrict mhd_st_z,
                                              CCTK_REAL *restrict tau,
                                              CCTK_REAL *restrict Ye_star,
                                              CCTK_REAL *restrict S_star,
                                              CCTK_REAL *restrict CONSERVS_avg_neighbors ) {

  // Auxiliary variables
  const int Nx0 = cctk_lsh[0]-1;
  const int Nx1 = cctk_lsh[1]-1;
  const int Nx2 = cctk_lsh[2]-1;

  // Set the index limits for the averaging loop
  // Minimum values
  const int iimin = (i!=0)*(i-1); // i=0 => iimin=0. Otherwise iimin = i-1
  const int jjmin = (j!=0)*(j-1); // j=0 => jjmin=0. Otherwise jjmin = j-1
  const int kkmin = (k!=0)*(k-1); // k=0 => kkmin=0. Otherwise kkmin = k-1

  // Maximum values  
  const int iimax = i + 1*(i!=Nx0); // i=Nx0 => iimax = Nx0. Otherwise iimax = i+1
  const int jjmax = j + 1*(j!=Nx1); // j=Nx1 => jjmax = Nx1. Otherwise jjmax = i+1
  const int kkmax = k + 1*(k!=Nx2); // k=Nx2 => kkmax = Nx2. Otherwise kkmax = i+1
  
  // Initialize the averages to zero
  for(int which_con=0;which_con<NUM_CONSERVS;which_con++)
    CONSERVS_avg_neighbors[which_con] = 0.0;

  // Now add up each of the neighbor conservative variables
  int number_of_neighbors = 0;
  for(int kk=kkmin;kk<=kkmax;kk++) {
    for(int jj=jjmin;jj<=jjmax;jj++) {
      for(int ii=iimin;ii<=iimax;ii++) {

        // Only count if at a neighbor point
        if( (ii!=i) && (jj!=j) && (kk!=k) ) {

          // Index of the neighbor gridfunction
          const int idx = CCTK_GFINDEX3D(cctkGH,ii,jj,kk);

          // Only count if con2prim succeeded at neighbor point
          if( con2prim_failed_flag[idx] == 0 ) {

            CONSERVS_avg_neighbors[RHOSTAR  ] += rho_star[idx];
            CONSERVS_avg_neighbors[STILDEX  ] += mhd_st_x[idx];
            CONSERVS_avg_neighbors[STILDEY  ] += mhd_st_y[idx];
            CONSERVS_avg_neighbors[STILDEZ  ] += mhd_st_z[idx];
            CONSERVS_avg_neighbors[TAUENERGY] += tau     [idx];
            CONSERVS_avg_neighbors[YESTAR   ] += Ye_star [idx];
            CONSERVS_avg_neighbors[ENTSTAR  ] += S_star  [idx];
            number_of_neighbors++;

          } // if( con2prim_failed_flag[idx] == 0 )
        } // if( (ii!=i) && (jj!=j) && (kk!=k) )

      } // for(int ii=iimax;ii<=iimax;ii++)
    } // for(int jj=jjmin;jj<=jjmax;jj++)
  } // for(int kk=kkmin;kk<=kkmax;kk++)

  // If not enough neighbors can be found, this mean that we have
  // many con2prim failures clustered in a region of the grid.
  // Probably best to just terminate the run.
  if( number_of_neighbors < 4 ) CCTK_VError(VERR_DEF_PARAMS,"Could not find enough neighbors to perform the conservative average: %d. ABORTING!",number_of_neighbors);

  // Now just compute the average of the neighboring conservs
  const CCTK_REAL inv_num_neighbors = 1.0/((CCTK_REAL)number_of_neighbors);
  for(int which_con=0;which_con<NUM_CONSERVS;which_con++)
    CONSERVS_avg_neighbors[which_con] *= inv_num_neighbors;

  // All done!

}
