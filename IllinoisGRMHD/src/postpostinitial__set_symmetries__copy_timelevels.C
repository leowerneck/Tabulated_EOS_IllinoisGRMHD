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

extern "C" void IllinoisGRMHD_PostPostInitial_Set_Symmetries__Copy_Timelevels(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /**********************************
   * Piecewise Polytropic EOS Patch *
   *     Printing the EOS table     *
   **********************************/
  /*
   * The short piece of code below takes care
   * of initializing the EOS parameters.
   * Please refer to the "inlined_functions.C"
   * source file for the documentation on the
   * function.
   */
  igm_eos_parameters eos;
  initialize_igm_eos_parameters_from_input(igm_eos_key,cctk_time,eos);

  // Perform parameter checks

  // Hybrid EOS
  if( CCTK_EQUALS(igm_eos_type,"Hybrid") ) {
    if( Gamma_th<0 ) {
      CCTK_VError(VERR_DEF_PARAMS,"Default Gamma_th (=-1) detected. You must set Gamma_th to the appropriate value in your initial data thorn, or your .par file!");
    }
    if( !CCTK_EQUALS(igm_con2prim_routine,"Noble2D") &&
        !CCTK_EQUALS(igm_con2prim_routine,"Noble1D") &&
        !CCTK_EQUALS(igm_con2prim_routine,"Noble1D_entropy") &&
        !CCTK_EQUALS(igm_con2prim_routine,"Noble1D_entropy2")) {
      CCTK_VError(VERR_DEF_PARAMS,"IllinoisGRMHD only supports the Noble2D, Noble1D, Noble1D_entropy, and Noble1D_entropy2 con2prim routines with Hybrid EOS. ABORTING!");
    }
  }

  // Tabulated EOS
  if( CCTK_EQUALS(igm_eos_type,"Tabulated") || CCTK_EQUALS(igm_eos_type,"nuc_eos") ) {
    if( !CCTK_EQUALS(igm_con2prim_routine,"Palenzuela1D") ) {
      CCTK_VError(VERR_DEF_PARAMS,"IllinoisGRMHD only supports the Palenzuela1D con2prim routine with Tabulated EOS. ABORTING!");
    }
  }

  // Entropy
  if( igm_evolve_entropy == false ) {
    if( CCTK_EQUALS(igm_con2prim_routine,"Noble1D_entropy" ) ||
        CCTK_EQUALS(igm_con2prim_routine,"Noble1D_entropy2") ||
        CCTK_EQUALS(igm_con2prim_routine,"Palenzuela1D"    ) ) {
      CCTK_VError(VERR_DEF_PARAMS,"Cannot use the an entropy con2prim routine without evolving the entropy. Please set igm_evolve_entropy=\"yes\" in the parfile. ABORTING!");
    }
  }

  CCTK_VInfo(CCTK_THORNSTRING,"Primary conservative-to-primitive routine: %s",igm_con2prim_routine);
  if( CCTK_EQUALS(igm_con2prim_backup_routine[0],"None") ) {
    CCTK_VInfo(CCTK_THORNSTRING,"No backup conservative-to-primitive routine selected.");
  }
  else {
    CCTK_VInfo(CCTK_THORNSTRING,"Backup conservative-to-primitive routine #1: %s",igm_con2prim_backup_routine[0]);
    if( !CCTK_EQUALS(igm_con2prim_backup_routine[1],"None") ) {
      CCTK_VInfo(CCTK_THORNSTRING,"Backup conservative-to-primitive routine #2: %s",igm_con2prim_backup_routine[1]);
      if( !CCTK_EQUALS(igm_con2prim_backup_routine[2],"None") ) {
        CCTK_VInfo(CCTK_THORNSTRING,"Backup conservative-to-primitive routine #3: %s",igm_con2prim_backup_routine[2]);
      }
    }
  }

  //For emfields, we assume that you've set Bx, By, Bz (the UN-tilded B^i's)
  // or Ax, Ay, Az (if using constrained transport scheme of Del Zanna)

  if(CCTK_EQUALS(Symmetry,"equatorial")) {
    // SET SYMMETRY GHOSTZONES ON ALL CONSERVATIVE AND PRIMIIVE VARIABLES!
    int ierr;
    ierr=CartSymGN(cctkGH,"IllinoisGRMHD::grmhd_conservatives"); if(ierr!=0) CCTK_VError(VERR_DEF_PARAMS,"Microsoft error code #1874109358120048. Grep it in the source code");
    ierr=CartSymGN(cctkGH,"IllinoisGRMHD::grmhd_primitives_allbutBi"); if(ierr!=0) CCTK_VError(VERR_DEF_PARAMS,"Microsoft error code #1874109358120049. Grep it in the source code");

    // Finish up by setting symmetry ghostzones on Bx, By, Bz, and their staggered variants.
    CCTK_REAL gridfunc_syms_Bx[3] = {-1, 1,-Sym_Bz};
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,  Bx        , gridfunc_syms_Bx,0,0,0);
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,  Bx_stagger, gridfunc_syms_Bx,1,0,0);
    CCTK_REAL gridfunc_syms_By[3] = { 1,-1,-Sym_Bz};
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,  By        , gridfunc_syms_Bx,0,0,0);
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,  By_stagger, gridfunc_syms_By,0,1,0);
    CCTK_REAL gridfunc_syms_Bz[3] = { 1, 1, Sym_Bz};
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,  Bz        , gridfunc_syms_Bz,0,0,0);
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,  Bz_stagger, gridfunc_syms_Bz,0,0,1);

    CCTK_REAL gridfunc_syms_psi6phi[3] = { 1, 1,      1};
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,psi6phi     , gridfunc_syms_psi6phi,1,1,1);
    CCTK_REAL gridfunc_syms_Ax[3]      = {-1, 1, Sym_Bz};
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,  Ax        , gridfunc_syms_Ax,0,1,1);
    CCTK_REAL gridfunc_syms_Ay[3]      = { 1,-1, Sym_Bz};
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,  Ay        , gridfunc_syms_Ay,1,0,1);
    CCTK_REAL gridfunc_syms_Az[3]      = { 1, 1,-Sym_Bz};
    IllinoisGRMHD_set_symmetry_gzs_staggered(cctkGH,cctk_lsh,x,y,z,  Az        , gridfunc_syms_Az,1,1,0);
  }


  //------------------------------------------------------------------
  // FILL _p AND _p_p TIMELEVELS. Probably don't need to do this if
  // Carpet::init_fill_timelevels=yes  and
  // MoL::initial_data_is_crap = yes
  // NOTE: We don't fill metric data here.
  // FIXME: Do we really need this?

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        rho_star_p[index]    = rho_star[index];
        tau_p[index]         = tau[index];
        mhd_st_x_p[index]    = mhd_st_x[index];
        mhd_st_y_p[index]    = mhd_st_y[index];
        mhd_st_z_p[index]    = mhd_st_z[index];

        psi6phi_p[index]     = psi6phi[index];
        Ax_p[index]          = Ax[index];
        Ay_p[index]          = Ay[index];
        Az_p[index]          = Az[index];

        rho_star_p_p[index]  = rho_star[index];
        tau_p_p[index]       = tau[index];
        mhd_st_x_p_p[index]  = mhd_st_x[index];
        mhd_st_y_p_p[index]  = mhd_st_y[index];
        mhd_st_z_p_p[index]  = mhd_st_z[index];

        psi6phi_p_p[index]   = psi6phi[index];
        Ax_p_p[index]        = Ax[index];
        Ay_p_p[index]        = Ay[index];
        Az_p_p[index]        = Az[index];

        if( eos.is_Tabulated ) {
          Ye_star_p[index]   = Ye_star[index];
          Ye_star_p_p[index] = Ye_star[index];
        }

        if( eos.evolve_entropy ) {
          S_star_p[index]    = S_star[index];
          S_star_p_p[index]  = S_star[index];
        }
      }
}

