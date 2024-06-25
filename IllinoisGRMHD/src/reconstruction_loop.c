#include "IllinoisGRMHD.h"

// This reconstruction function is only used on B_stagger and reconstructed velocities, so it
// is independent of EOS
void IllinoisGRMHD_reconstruction_loop(const cGH *restrict cctkGH, const int flux_dir, const int num_vars,
                         const int *restrict var_indices,
                         const CCTK_REAL *pressure,
                         const CCTK_REAL *v_flux,
                         const CCTK_REAL **in_prims,
                         CCTK_REAL **out_prims_r,
                         CCTK_REAL **out_prims_l) {

  const int xdir = (flux_dir == 0);
  const int ydir = (flux_dir == 1);
  const int zdir = (flux_dir == 2);

  const int imin = cctkGH->cctk_nghostzones[0];
  const int imax = (cctkGH->cctk_lsh[0] - cctkGH->cctk_nghostzones[0]) + 1;
  const int jmin = cctkGH->cctk_nghostzones[1];
  const int jmax = (cctkGH->cctk_lsh[1] - cctkGH->cctk_nghostzones[1]) + 1;
  const int kmin = cctkGH->cctk_nghostzones[2];
  const int kmax = (cctkGH->cctk_lsh[2] - cctkGH->cctk_nghostzones[2]) + 1;

#pragma omp parallel for
  for(int k=kmin; k<kmax; k++) {
    for(int j=jmin; j<jmax; j++) {
      for(int i=imin; i<imax; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
        CCTK_REAL press_stencil[6], v_flux_stencil[6];
        CCTK_REAL var_data[num_vars][6], vars_r[num_vars], vars_l[num_vars];

        for(int ind=0; ind<6; ind++) {
          // Stencil from -3 to +2 reconstructs to e.g. i-1/2
          const int stencil = CCTK_GFINDEX3D(cctkGH, i+xdir*(ind-3), j+ydir*(ind-3), k+zdir*(ind-3));
          v_flux_stencil[ind] = v_flux[stencil]; // Could be smaller; doesn't use full stencil
          press_stencil[ind] = pressure[stencil];
          for(int var=0; var<num_vars; var++) {
            var_data[var][ind] = in_prims[var_indices[var]][stencil];
          }
        }

        CCTK_REAL ftilde[2];
        ghl_compute_ftilde(ghl_params, press_stencil, v_flux_stencil, ftilde);

        for(int var=0; var<num_vars; var++) {
          ghl_ppm_reconstruction(ftilde, var_data[var], &vars_r[var], &vars_l[var]);
          out_prims_r[var_indices[var]][index] = vars_r[var];
          out_prims_l[var_indices[var]][index] = vars_l[var];
        }
      }
    }
  }
}
