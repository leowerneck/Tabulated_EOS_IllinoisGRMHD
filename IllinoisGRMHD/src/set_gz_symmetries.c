//-------------------------------------------------
// Stuff to run right after initial data is set up
//-------------------------------------------------

#include "IllinoisGRMHD.h"
#include "Symmetry.h"

void IllinoisGRMHD_set_gz_symmetries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  //For emfields, we assume that you've set Bx, By, Bz (the UN-tilded B^i's)
  // or Ax, Ay, Az (if using constrained transport scheme of Del Zanna)

  if(CCTK_EQUALS(Symmetry,"equatorial")) {
    // Set up symmetry ghostzones on Bx, By, Bz, and their staggered variants.
    const CCTK_REAL gridfunc_syms_phitilde[3] = { 1, 1, 1};
    const CCTK_REAL gridfunc_syms_Bx[3] = {-1,  1, -Sym_Bz};
    const CCTK_REAL gridfunc_syms_By[3] = { 1, -1, -Sym_Bz};
    const CCTK_REAL gridfunc_syms_Bz[3] = { 1,  1,  Sym_Bz};

    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z, phitilde, gridfunc_syms_phitilde, 1, 1, 1);

    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z, Bx_center,  gridfunc_syms_Bx, 0, 0, 0);
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z, Bx_stagger, gridfunc_syms_Bx, 1, 0, 0);
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z, Ax,         gridfunc_syms_Bx, 0, 1, 1);

    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z, By_center,  gridfunc_syms_By, 0, 0, 0);
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z, By_stagger, gridfunc_syms_By, 0, 1, 0);
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z, Ay,         gridfunc_syms_By, 1, 0, 1);

    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z, Bz_center,  gridfunc_syms_Bz, 0, 0, 0);
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z, Bz_stagger, gridfunc_syms_Bz, 0, 0, 1);
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH, x, y, z, Az,         gridfunc_syms_Bz, 1, 1, 0);
  }
}
