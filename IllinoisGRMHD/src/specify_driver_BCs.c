#include "IllinoisGRMHD.h"
/*
 * Set all boundary conditions for PreSync.
 *
 * All variables used by this thorn have their boundary
 * conditions handled manually, so we register "none" to
 * remove warnings.
 */
void IllinoisGRMHD_specify_driver_BCs(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS_IllinoisGRMHD_specify_driver_BCs;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;

  // Primitive and conservative variables BCs are set by the
  // IllinoisGRMHD_EOS_hydro_outer_boundaries() function using
  // copy BCs with optional outflow BCs for the velocities.
  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "HydroBase::rho", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for HydroBase::rho!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "HydroBase::press", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for HydroBase::press!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "HydroBase::eps", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for HydroBase::eps!");

  ierr = Driver_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "IllinoisGRMHD::grmhd_velocities", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for IllinoisGRMHD::grmhd_velocities!");

  ierr = Driver_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "IllinoisGRMHD::grmhd_conservatives", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for IllinoisGRMHD::grmhd_conservatives!");



  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "HydroBase::entropy", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for HydroBase::entropy!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "IllinoisGRMHD::ent_star", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for IllinoisGRMHD::ent_star!");



    ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "HydroBase::Y_e", "none");
    if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for HydroBase::Y_e!");

    ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "HydroBase::temperature", "none");
    if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for HydroBase::temperature!");

    ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "IllinoisGRMHD::Ye_star", "none");
    if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for IllinoisGRMHD::Ye_star!");



  // A_i BCs are set by the IllinoisGRMHD_A_i_outer_boundaries
  // function using linear interpolation.
  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "IllinoisGRMHD::Ax", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for IllinoisGRMHD::Ax!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "IllinoisGRMHD::Ay", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for IllinoisGRMHD::Ay!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "IllinoisGRMHD::Az", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for IllinoisGRMHD::Az!");

  ierr = Driver_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "IllinoisGRMHD::phitilde", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for IllinoisGRMHD::phitilde!");



  // B boundaries are sety by the IllinoisGRMHD_compute_B_and_Bstagger_from_A
  // function using copy boundary conditions.
  ierr = Driver_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "IllinoisGRMHD::grmhd_B_center", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for IllinoisGRMHD::grmhd_B_center!");

  ierr = Driver_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "IllinoisGRMHD::grmhd_B_stagger", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC with Driver for IllinoisGRMHD::grmhd_B_stagger!");
printf("3\n");
}
