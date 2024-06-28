//--------------------------------------------------------------------------
// Register with the time stepper
// (MoL thorn, found in arrangements/CactusBase/MoL)
// To understand this, read documentation in arrangements/CactusBase/MoL/doc
//--------------------------------------------------------------------------

#include "IllinoisGRMHD.h"
#include "Symmetry.h"

void IllinoisGRMHD_RegisterVars(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, group, rhs;

  //***********************************************
  // Register evolution & RHS gridfunction variables

  /* Ax and Ax_rhs */
  group = CCTK_VarIndex("IllinoisGRMHD::Ax");
  rhs = CCTK_VarIndex("IllinoisGRMHD::Ax_rhs");
  ierr += MoLRegisterEvolved(group, rhs);

  /* Ay and Ay_rhs */
  group = CCTK_VarIndex("IllinoisGRMHD::Ay");
  rhs = CCTK_VarIndex("IllinoisGRMHD::Ay_rhs");
  ierr += MoLRegisterEvolved(group, rhs);

  /* Az and Az_rhs */
  group = CCTK_VarIndex("IllinoisGRMHD::Az");
  rhs = CCTK_VarIndex("IllinoisGRMHD::Az_rhs");
  ierr += MoLRegisterEvolved(group, rhs);

  /* phitilde and phitilde_rhs */
  group = CCTK_VarIndex("IllinoisGRMHD::phitilde");
  rhs = CCTK_VarIndex("IllinoisGRMHD::phitilde_rhs");
  ierr += MoLRegisterEvolved(group, rhs);

  /* ALL OTHER EVOLVED VARIABLES (rho_star,tau,Stilde_x,Stilde_y,Stilde_z) */
  group = CCTK_GroupIndex("IllinoisGRMHD::grmhd_conservatives");
  rhs = CCTK_GroupIndex("IllinoisGRMHD::grmhd_conservatives_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("HydroBase::rho"));
  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("HydroBase::press"));
  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("HydroBase::eps"));
  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("HydroBase::Y_e"));
  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("HydroBase::temperature"));

  if(ghl_params->evolve_entropy) {
  //  group = CCTK_GroupIndex("IllinoisGRMHD::ent_star");
   // rhs = CCTK_GroupIndex("IllinoisGRMHD::ent_star_rhs");
   // ierr += MoLRegisterEvolvedGroup(group, rhs);
    ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("HydroBase::entropy"));
  }

  if(ghl_eos->eos_type == ghl_eos_tabulated) {
  //  group = CCTK_GroupIndex("IllinoisGRMHD::Ye_star");
    //rhs = CCTK_GroupIndex("IllinoisGRMHD::Ye_star_rhs");
  //  ierr += MoLRegisterEvolvedGroup(group, rhs);
    ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("HydroBase::Y_e"));
    ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("HydroBase::temperature"));
  }

  if(update_Tmunu) {
    MoLRegisterConstrainedGroup(CCTK_GroupIndex("TmunuBase::stress_energy_scalar"));
    MoLRegisterConstrainedGroup(CCTK_GroupIndex("TmunuBase::stress_energy_vector"));
    MoLRegisterConstrainedGroup(CCTK_GroupIndex("TmunuBase::stress_energy_tensor"));
  }

  if (ierr) CCTK_ERROR("Problems registering with MoL");

  //***********************************************
  // Next register ADMBase variables needed by
  //    IllinoisGRMHD as SaveAndRestore, so that
  //    they are not set to NaN at the start of
  //    each timestep (requiring that they be
  //    e.g., recomputed from BSSN variables
  //    in the BSSN solver, like Baikal or
  //    ML_BSSN)
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("ADMBase::lapse"));
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("ADMBase::shift"));
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("ADMBase::metric"));
  ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("ADMBase::curv"));
  if (ierr) CCTK_ERROR("Problems registering with MoLRegisterSaveAndRestoreGroup");
}
