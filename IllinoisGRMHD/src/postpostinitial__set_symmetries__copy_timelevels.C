//-------------------------------------------------
// Stuff to run right after initial data is set up
//-------------------------------------------------

#include "cctk.h"
#include <cstdio>
#include <cstdlib>
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"
#include "IllinoisGRMHD_headers.h"

extern "C" void
IllinoisGRMHD_PostPostInitial_Set_Symmetries__Copy_Timelevels(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /**********************************
   * Piecewise Polytropic EOS Patch *
   *     Printing the EOS table     *
   **********************************/
  /*
   * The short piece of code below takes care
   * of initializing the EOS parameters.
   * Please refer to the "inlined_functions.h"
   * source file for the documentation on the
   * function.
   */
  igm_eos_parameters eos;
  initialize_igm_eos_parameters_from_input(igm_eos_key, cctk_time, eos);

  // For emfields, we assume that you've set Bx, By, Bz (the UN-tilded B^i's)
  //  or A_x_tilde, A_y_tilde, A_z_tilde (if using constrained transport scheme
  //  of Del Zanna)

  if (CCTK_EQUALS(Symmetry, "equatorial")) {
    // SET SYMMETRY GHOSTZONES ON ALL CONSERVATIVE AND PRIMIIVE VARIABLES!
    int ierr;
    ierr = CartSymGN(cctkGH, "IllinoisGRMHD::grmhd_conservatives");
    if (ierr != 0)
      CCTK_VError(
          VERR_DEF_PARAMS,
          "Microsoft error code #1874109358120048. Grep it in the source code");
    ierr = CartSymGN(cctkGH, "IllinoisGRMHD::grmhd_primitives_allbutBi");
    if (ierr != 0)
      CCTK_VError(
          VERR_DEF_PARAMS,
          "Microsoft error code #1874109358120049. Grep it in the source code");

    // Finish up by setting symmetry ghostzones on Bx, By, Bz, and their
    // staggered variants.
    CCTK_REAL gridfunc_syms_Bx[3] = {-1, 1, -Sym_Bz};
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, cctk_lsh, x, y, z, Bx_center,
                                             gridfunc_syms_Bx, 0, 0, 0);
    IllinoisGRMHD_set_symmetry_gzs_staggered(
        cctkGH, cctk_lsh, x, y, z, Bx_stagger, gridfunc_syms_Bx, 1, 0, 0);
    CCTK_REAL gridfunc_syms_By[3] = {1, -1, -Sym_Bz};
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, cctk_lsh, x, y, z, By_center,
                                             gridfunc_syms_Bx, 0, 0, 0);
    IllinoisGRMHD_set_symmetry_gzs_staggered(
        cctkGH, cctk_lsh, x, y, z, By_stagger, gridfunc_syms_By, 0, 1, 0);
    CCTK_REAL gridfunc_syms_Bz[3] = {1, 1, Sym_Bz};
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, cctk_lsh, x, y, z, Bz_center,
                                             gridfunc_syms_Bz, 0, 0, 0);
    IllinoisGRMHD_set_symmetry_gzs_staggered(
        cctkGH, cctk_lsh, x, y, z, Bz_stagger, gridfunc_syms_Bz, 0, 0, 1);

    CCTK_REAL gridfunc_syms_Phi_tilde[3] = {1, 1, 1};
    IllinoisGRMHD_set_symmetry_gzs_staggered(
        cctkGH, cctk_lsh, x, y, z, Phi_tilde, gridfunc_syms_Phi_tilde, 1, 1, 1);
    CCTK_REAL gridfunc_syms_A_x_tilde[3] = {-1, 1, Sym_Bz};
    IllinoisGRMHD_set_symmetry_gzs_staggered(
        cctkGH, cctk_lsh, x, y, z, A_x_tilde, gridfunc_syms_A_x_tilde, 0, 1, 1);
    CCTK_REAL gridfunc_syms_A_y_tilde[3] = {1, -1, Sym_Bz};
    IllinoisGRMHD_set_symmetry_gzs_staggered(
        cctkGH, cctk_lsh, x, y, z, A_y_tilde, gridfunc_syms_A_y_tilde, 1, 0, 1);
    CCTK_REAL gridfunc_syms_A_z_tilde[3] = {1, 1, -Sym_Bz};
    IllinoisGRMHD_set_symmetry_gzs_staggered(
        cctkGH, cctk_lsh, x, y, z, A_z_tilde, gridfunc_syms_A_z_tilde, 1, 1, 0);
  }

  //------------------------------------------------------------------
  // FILL _p AND _p_p TIMELEVELS. Probably don't need to do this if
  // Carpet::init_fill_timelevels=yes  and
  // MoL::initial_data_is_crap = yes
  // NOTE: We don't fill metric data here.
  // FIXME: Do we really need this?

#pragma omp parallel for
  for (int k = 0; k < cctk_lsh[2]; k++)
    for (int j = 0; j < cctk_lsh[1]; j++)
      for (int i = 0; i < cctk_lsh[0]; i++) {
        int index = CCTK_GFINDEX3D(cctkGH, i, j, k);

        rho_tilde_p[index] = rho_tilde[index];
        tau_tilde_p[index] = tau_tilde[index];
        S_x_tilde_p[index] = S_x_tilde[index];
        S_y_tilde_p[index] = S_y_tilde[index];
        S_z_tilde_p[index] = S_z_tilde[index];

        Phi_tilde_p[index] = Phi_tilde[index];
        A_x_tilde_p[index] = A_x_tilde[index];
        A_y_tilde_p[index] = A_y_tilde[index];
        A_z_tilde_p[index] = A_z_tilde[index];

        rho_tilde_p_p[index] = rho_tilde[index];
        tau_tilde_p_p[index] = tau_tilde[index];
        S_x_tilde_p_p[index] = S_x_tilde[index];
        S_y_tilde_p_p[index] = S_y_tilde[index];
        S_z_tilde_p_p[index] = S_z_tilde[index];

        Phi_tilde_p_p[index] = Phi_tilde[index];
        A_x_tilde_p_p[index] = A_x_tilde[index];
        A_y_tilde_p_p[index] = A_y_tilde[index];
        A_z_tilde_p_p[index] = A_z_tilde[index];

        if (eos.is_Tabulated) {
          Y_e_tilde_p[index] = Y_e_tilde[index];
          Y_e_tilde_p_p[index] = Y_e_tilde[index];
        }

        if (eos.evolve_entropy) {
          ent_tilde_p[index] = ent_tilde[index];
          ent_tilde_p_p[index] = ent_tilde[index];
        }
      }
}
