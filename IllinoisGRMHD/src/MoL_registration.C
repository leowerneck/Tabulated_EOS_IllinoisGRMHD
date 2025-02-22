//--------------------------------------------------------------------------
// Register with the time stepper
// (MoL thorn, found in arrangements/CactusBase/MoL)
// To understand this, read documentation in arrangements/CactusBase/MoL/doc
//--------------------------------------------------------------------------

#include "cctk.h"
#include <cstdio>
#include <cmath>
#include <cstddef>
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "Symmetry.h"

extern "C" void IllinoisGRMHD_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, group, rhs;

  // Register evolution & RHS gridfunction groups

  /* A_x_tilde and A_x_rhs */
  group = CCTK_GroupIndex("IllinoisGRMHD::A_x_tilde");
  rhs = CCTK_GroupIndex("IllinoisGRMHD::A_x_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* A_y_tilde and A_y_rhs */
  group = CCTK_GroupIndex("IllinoisGRMHD::A_y_tilde");
  rhs = CCTK_GroupIndex("IllinoisGRMHD::A_y_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* A_z_tilde and A_z_rhs */
  group = CCTK_GroupIndex("IllinoisGRMHD::A_z_tilde");
  rhs = CCTK_GroupIndex("IllinoisGRMHD::A_z_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Phi_tilde and Phi_rhs */
  group = CCTK_GroupIndex("IllinoisGRMHD::Phi_tilde");
  rhs = CCTK_GroupIndex("IllinoisGRMHD::Phi_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* ALL OTHER EVOLVED VARIABLES (rho_tilde,tau,S_x,S_y,S_z) */
  group = CCTK_GroupIndex("IllinoisGRMHD::grmhd_conservatives");
  rhs = CCTK_GroupIndex("IllinoisGRMHD::grmhd_rhss");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  if (ierr) CCTK_ERROR("Problems registering with MoL");
  //***********************************************

  //***********************************************
  // Next register ADMBase variables needed by
  //    IllinoisGRMHD as SaveAndRestore, so that
  //    they are not set to NaN at the start of
  //    each timestep (requiring that they be
  //    e.g., recomputed from BSSN variables
  //    in the BSSN solver, like Baikal or
  //    ML_BSSN)
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("admbase::lapse"));
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("admbase::shift"));
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("admbase::metric"));
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("admbase::curv"));
  if (ierr) CCTK_ERROR("Problems registering with MoLRegisterSaveAndRestoreGroup");
}
