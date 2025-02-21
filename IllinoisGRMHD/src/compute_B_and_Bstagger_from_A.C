#include "cctk.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "IllinoisGRMHD_headers.h"

#define LOOP_DEFINE_SIMPLE                                                     \
  _Pragma("omp parallel for") for (int k = 0; k < cctk_lsh[2];                 \
                                   k++) for (int j = 0; j < cctk_lsh[1];       \
                                             j++) for (int i = 0;              \
                                                       i < cctk_lsh[0]; i++)

extern "C" void IllinoisGRMHD_compute_B_and_Bstagger_from_A(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL dxi = 1.0 / CCTK_DELTA_SPACE(0);
  CCTK_REAL dyi = 1.0 / CCTK_DELTA_SPACE(1);
  CCTK_REAL dzi = 1.0 / CCTK_DELTA_SPACE(2);

  CCTK_REAL gridfunc_syms_A_x_tilde[3] = {-1, 1, Sym_Bz};
  IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, cctk_lsh, x, y, z, A_x_tilde,
                                           gridfunc_syms_A_x_tilde, 0, 1, 1);
  CCTK_REAL gridfunc_syms_A_y_tilde[3] = {1, -1, Sym_Bz};
  IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, cctk_lsh, x, y, z, A_y_tilde,
                                           gridfunc_syms_A_y_tilde, 1, 0, 1);
  CCTK_REAL gridfunc_syms_A_z_tilde[3] = {1, 1, -Sym_Bz};
  IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, cctk_lsh, x, y, z, A_z_tilde,
                                           gridfunc_syms_A_z_tilde, 1, 1, 0);
  CCTK_REAL gridfunc_syms_Phi_tilde[3] = {1, 1, 1};
  IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, cctk_lsh, x, y, z, Phi_tilde,
                                           gridfunc_syms_Phi_tilde, 1, 1, 1);

  LOOP_DEFINE_SIMPLE {
    int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
    psi_bssn[index] = exp(phi_bssn[index]);
  }

  LOOP_DEFINE_SIMPLE {
    // Look Mom, no if() statements!
    int shiftedim1 =
        (i - 1) * (i != 0); // This way, i=0 yields shiftedim1=0 and shiftedi=1,
                            // used below for our COPY boundary condition.
    int shiftedi = shiftedim1 + 1;

    int shiftedjm1 = (j - 1) * (j != 0);
    int shiftedj = shiftedjm1 + 1;

    int shiftedkm1 = (k - 1) * (k != 0);
    int shiftedk = shiftedkm1 + 1;

    int index, indexim1, indexjm1, indexkm1;

    int actual_index = CCTK_GFINDEX3D(cctkGH, i, j, k);

    CCTK_REAL Psi = psi_bssn[actual_index];
    CCTK_REAL Psim3 = 1.0 / (Psi * Psi * Psi);

    // For the lower boundaries, the following applies a "copy"
    //    boundary condition on Bi_stagger where needed.
    //    E.g., Bx_stagger(i,jmin,k) = Bx_stagger(i,jmin+1,k)
    //    We find the copy BC works better than extrapolation.
    // For the upper boundaries, we do the following copy:
    //    E.g., Psi(imax+1,j,k)=Psi(imax,j,k)
    /**************/
    /* Bx_stagger */
    /**************/

    index = CCTK_GFINDEX3D(cctkGH, i, shiftedj, shiftedk);
    indexjm1 = CCTK_GFINDEX3D(cctkGH, i, shiftedjm1, shiftedk);
    indexkm1 = CCTK_GFINDEX3D(cctkGH, i, shiftedj, shiftedkm1);
    // Set Bx_stagger = \partial_y A_z - partial_z A_y
    // "Grid" A_x_tilde(i,j,k) is actually A_x_tilde(i,j+1/2,k+1/2)
    // "Grid" A_y_tilde(i,j,k) is actually A_y_tilde(i+1/2,j,k+1/2)
    // "Grid" A_z_tilde(i,j,k) is actually A_z_tilde(i+1/2,j+1/2,k)
    // Therefore, the 2nd order derivative \partial_z A_y at (i+1/2,j,k) is:
    //          ["Grid" A_y_tilde(i,j,k) - "Grid" A_y_tilde(i,j,k-1)]/dZ
    Bx_stagger[actual_index] = (A_z_tilde[index] - A_z_tilde[indexjm1]) * dyi -
                               (A_y_tilde[index] - A_y_tilde[indexkm1]) * dzi;

    // Now multiply Bx and Bx_stagger by 1/sqrt(gamma(i+1/2,j,k)]) = 1/sqrt(1/2
    // [gamma + gamma_ip1]) = exp(-6 x 1/2 [phi + phi_ip1] )
    int imax_minus_i = (cctk_lsh[0] - 1) - i;
    int indexip1jk = CCTK_GFINDEX3D(
        cctkGH, i + ((imax_minus_i > 0) - (0 > imax_minus_i)), j, k);
    CCTK_REAL Psi_ip1 = psi_bssn[indexip1jk];
    Bx_stagger[actual_index] *= Psim3 / (Psi_ip1 * Psi_ip1 * Psi_ip1);

    /**************/
    /* By_stagger */
    /**************/

    index = CCTK_GFINDEX3D(cctkGH, shiftedi, j, shiftedk);
    indexim1 = CCTK_GFINDEX3D(cctkGH, shiftedim1, j, shiftedk);
    indexkm1 = CCTK_GFINDEX3D(cctkGH, shiftedi, j, shiftedkm1);
    // Set By_stagger = \partial_z A_x - \partial_x A_z
    By_stagger[actual_index] = (A_x_tilde[index] - A_x_tilde[indexkm1]) * dzi -
                               (A_z_tilde[index] - A_z_tilde[indexim1]) * dxi;

    // Now multiply By and By_stagger by 1/sqrt(gamma(i,j+1/2,k)]) = 1/sqrt(1/2
    // [gamma + gamma_jp1]) = exp(-6 x 1/2 [phi + phi_jp1] )
    int jmax_minus_j = (cctk_lsh[1] - 1) - j;
    int indexijp1k = CCTK_GFINDEX3D(
        cctkGH, i, j + ((jmax_minus_j > 0) - (0 > jmax_minus_j)), k);
    CCTK_REAL Psi_jp1 = psi_bssn[indexijp1k];
    By_stagger[actual_index] *= Psim3 / (Psi_jp1 * Psi_jp1 * Psi_jp1);

    /**************/
    /* Bz_stagger */
    /**************/

    index = CCTK_GFINDEX3D(cctkGH, shiftedi, shiftedj, k);
    indexim1 = CCTK_GFINDEX3D(cctkGH, shiftedim1, shiftedj, k);
    indexjm1 = CCTK_GFINDEX3D(cctkGH, shiftedi, shiftedjm1, k);
    // Set Bz_stagger = \partial_x A_y - \partial_y A_x
    Bz_stagger[actual_index] = (A_y_tilde[index] - A_y_tilde[indexim1]) * dxi -
                               (A_x_tilde[index] - A_x_tilde[indexjm1]) * dyi;

    // Now multiply Bz_stagger by 1/sqrt(gamma(i,j,k+1/2)]) = 1/sqrt(1/2 [gamma
    // + gamma_kp1]) = exp(-6 x 1/2 [phi + phi_kp1] )
    int kmax_minus_k = (cctk_lsh[2] - 1) - k;
    int indexijkp1 = CCTK_GFINDEX3D(
        cctkGH, i, j, k + ((kmax_minus_k > 0) - (0 > kmax_minus_k)));
    CCTK_REAL Psi_kp1 = psi_bssn[indexijkp1];
    Bz_stagger[actual_index] *= Psim3 / (Psi_kp1 * Psi_kp1 * Psi_kp1);
  }

  LOOP_DEFINE_SIMPLE {
    // Look Mom, no if() statements!
    int shiftedim1 =
        (i - 1) * (i != 0); // This way, i=0 yields shiftedim1=0 and shiftedi=1,
                            // used below for our COPY boundary condition.
    int shiftedi = shiftedim1 + 1;

    int shiftedjm1 = (j - 1) * (j != 0);
    int shiftedj = shiftedjm1 + 1;

    int shiftedkm1 = (k - 1) * (k != 0);
    int shiftedk = shiftedkm1 + 1;

    int index, indexim1, indexjm1, indexkm1;

    int actual_index = CCTK_GFINDEX3D(cctkGH, i, j, k);

    // For the lower boundaries, the following applies a "copy"
    //    boundary condition on Bi and Bi_stagger where needed.
    //    E.g., Bx(imin,j,k) = Bx(imin+1,j,k)
    //    We find the copy BC works better than extrapolation.
    /******/
    /* Bx */
    /******/
    index = CCTK_GFINDEX3D(cctkGH, shiftedi, j, k);
    indexim1 = CCTK_GFINDEX3D(cctkGH, shiftedim1, j, k);
    // Set Bx = 0.5 ( Bx_stagger + Bx_stagger_im1 )
    // "Grid" Bx_stagger(i,j,k) is actually Bx_stagger(i+1/2,j,k)
    Bx_center[actual_index] = 0.5 * (Bx_stagger[index] + Bx_stagger[indexim1]);

    /******/
    /* By */
    /******/
    index = CCTK_GFINDEX3D(cctkGH, i, shiftedj, k);
    indexjm1 = CCTK_GFINDEX3D(cctkGH, i, shiftedjm1, k);
    // Set By = 0.5 ( By_stagger + By_stagger_im1 )
    // "Grid" By_stagger(i,j,k) is actually By_stagger(i,j+1/2,k)
    By_center[actual_index] = 0.5 * (By_stagger[index] + By_stagger[indexjm1]);

    /******/
    /* Bz */
    /******/
    index = CCTK_GFINDEX3D(cctkGH, i, j, shiftedk);
    indexkm1 = CCTK_GFINDEX3D(cctkGH, i, j, shiftedkm1);
    // Set Bz = 0.5 ( Bz_stagger + Bz_stagger_im1 )
    // "Grid" Bz_stagger(i,j,k) is actually Bz_stagger(i,j+1/2,k)
    Bz_center[actual_index] = 0.5 * (Bz_stagger[index] + Bz_stagger[indexkm1]);
  }

  // Finish up by setting symmetry ghostzones on Bx, By, Bz, and their staggered
  // variants.
  CCTK_REAL gridfunc_syms_Bx[3] = {-1, 1, -Sym_Bz};
  IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, cctk_lsh, x, y, z, Bx_center,
                                           gridfunc_syms_Bx, 0, 0, 0);
  IllinoisGRMHD_set_symmetry_gzs_staggered(
      cctkGH, cctk_lsh, x, y, z, Bx_stagger, gridfunc_syms_Bx, 1, 0, 0);
  CCTK_REAL gridfunc_syms_By[3] = {1, -1, -Sym_Bz};
  IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, cctk_lsh, x, y, z, By_center,
                                           gridfunc_syms_By, 0, 0, 0);
  IllinoisGRMHD_set_symmetry_gzs_staggered(
      cctkGH, cctk_lsh, x, y, z, By_stagger, gridfunc_syms_By, 0, 1, 0);
  CCTK_REAL gridfunc_syms_Bz[3] = {1, 1, Sym_Bz};
  IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, cctk_lsh, x, y, z, Bz_center,
                                           gridfunc_syms_Bz, 0, 0, 0);
  IllinoisGRMHD_set_symmetry_gzs_staggered(
      cctkGH, cctk_lsh, x, y, z, Bz_stagger, gridfunc_syms_Bz, 0, 0, 1);
}
