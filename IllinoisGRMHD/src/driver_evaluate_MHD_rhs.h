#ifndef DRIVER_EVALUATE_MHD_RHS_H_
#define DRIVER_EVALUATE_MHD_RHS_H_

/* PRIVATE FUNCTIONS, Called within driver_evaluate_MHD_rhs.C ONLY */
static void ftilde_gf_compute(const cGH *cctkGH,const int *cctk_lsh,const int flux_dirn,gf_and_gz_struct *input,CCTK_REAL *ftilde_gf);
static void reconstruct_set_of_prims_PPM(const cGH *cctkGH,const int *cctk_lsh,const int flux_dirn,const int num_prims_to_reconstruct,const int *which_prims_to_reconstruct,
                                         igm_eos_parameters &eosi,gf_and_gz_struct *in_prims,gf_and_gz_struct *out_prims_r,gf_and_gz_struct *out_prims_l,
                                         CCTK_REAL *ftilde_gf,CCTK_REAL *temporary);

static void compute_tau_rhs_extrinsic_curvature_terms_and_TUPmunu
(const cGH *cctkGH,const int *cctk_lsh,const int *cctk_nghostzones,CCTK_REAL *dX,CCTK_REAL **metric,gf_and_gz_struct *prims,
 CCTK_REAL **TUPmunu,igm_eos_parameters &eos,
 CCTK_REAL *gupxy,CCTK_REAL *gupxz,CCTK_REAL *gupyz,
 CCTK_REAL *kxx,CCTK_REAL *kxy,CCTK_REAL *kxz,CCTK_REAL *kyy,CCTK_REAL *kyz,CCTK_REAL *kzz,
 CCTK_REAL *tau_rhs);
static void A_i_rhs_no_gauge_terms(const int A_dirn, const cGH *cctkGH,const int *cctk_lsh,const int *cctk_nghostzones,gf_and_gz_struct *out_prims_r,gf_and_gz_struct *out_prims_l,
                                   CCTK_REAL *phi_interped,CCTK_REAL *cmax_1,CCTK_REAL *cmin_1,CCTK_REAL *cmax_2,CCTK_REAL *cmin_2, CCTK_REAL *A3_rhs);

static void Lorenz_psi6phi_rhs__add_gauge_terms_to_A_i_rhs(const cGH *cctkGH,const int *cctk_lsh,const int *cctk_nghostzones,CCTK_REAL *dX,CCTK_REAL **interp_vars,CCTK_REAL *psi6phi,
                                                           CCTK_REAL *shiftx_iphjphkph,CCTK_REAL *shifty_iphjphkph,CCTK_REAL *shiftz_iphjphkph,
                                                           CCTK_REAL *alpha_iphjphkph,CCTK_REAL *alpha_Phi_minus_betaj_A_j_iphjphkph,CCTK_REAL *alpha_sqrtg_Ax_interp,
                                                           CCTK_REAL *alpha_sqrtg_Ay_interp,CCTK_REAL *alpha_sqrtg_Az_interp,
                                                           CCTK_REAL *psi6phi_rhs,CCTK_REAL *Ax_rhs,CCTK_REAL *Ay_rhs,CCTK_REAL *Az_rhs);

static void add_fluxes_and_source_terms_to_hydro_rhss( const igm_eos_parameters eos,
                                                       const int flux_dirn,
                                                       const cGH *restrict cctkGH,
                                                       const int *restrict cctk_lsh,
                                                       const int *restrict cctk_nghostzones,
                                                       const CCTK_REAL *restrict dX,
                                                       CCTK_REAL **metric,
                                                       CCTK_REAL **TUPmunu,
                                                       const int numvars_reconstructed,
                                                       gf_and_gz_struct *restrict in_prims,
                                                       gf_and_gz_struct *restrict out_prims_r,
                                                       gf_and_gz_struct *restrict out_prims_l,
                                                       CCTK_REAL *restrict cmax,
                                                       CCTK_REAL *restrict cmin,
                                                       CCTK_REAL *restrict rho_star_flux,
                                                       CCTK_REAL *restrict tau_flux,
                                                       CCTK_REAL *restrict st_x_flux,
                                                       CCTK_REAL *restrict st_y_flux,
                                                       CCTK_REAL *restrict st_z_flux,
                                                       CCTK_REAL *restrict Ye_star_flux,
                                                       CCTK_REAL *restrict S_star_flux,
                                                       CCTK_REAL *restrict rho_star_rhs,
                                                       CCTK_REAL *restrict tau_rhs,
                                                       CCTK_REAL *restrict st_x_rhs,
                                                       CCTK_REAL *restrict st_y_rhs,
                                                       CCTK_REAL *restrict st_z_rhs,
                                                       CCTK_REAL *restrict Ye_star_rhs,
                                                       CCTK_REAL *restrict S_star_rhs );

#endif /* DRIVER_EVALUATE_MHD_RHS_H_ */

