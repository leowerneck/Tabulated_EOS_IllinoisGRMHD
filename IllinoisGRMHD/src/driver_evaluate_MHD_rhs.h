#ifndef DRIVER_EVALUATE_MHD_RHS_H_
#define DRIVER_EVALUATE_MHD_RHS_H_

/* PRIVATE FUNCTIONS, Called within driver_evaluate_MHD_rhs.C ONLY */
static void ftilde_gf_compute(
    const cGH *restrict cctkGH,
    const int *restrict cctk_lsh,
    const int flux_dirn,
    gf_and_gz_struct *restrict input,
    CCTK_REAL *restrict ftilde_gf
);

static void reconstruct_set_of_prims_PPM(
    const igm_eos_parameters eos,
    const cGH *restrict cctkGH,
    const int *restrict cctk_lsh,
    const int flux_dirn,
    const int num_prims_to_reconstruct,
    const int *restrict which_prims_to_reconstruct,
    gf_and_gz_struct *restrict in_prims,
    gf_and_gz_struct *restrict out_prims_r,
    gf_and_gz_struct *restrict out_prims_l,
    CCTK_REAL *restrict ftilde_gf,
    CCTK_REAL *restrict temporary
);

static void compute_tau_rhs_extrinsic_curvature_terms_and_TUPmunu(
    const igm_eos_parameters eos,
    const cGH *restrict cctkGH,
    const int *restrict cctk_lsh,
    const int *restrict cctk_nghostzones,
    CCTK_REAL *restrict dX,
    CCTK_REAL **metric,
    gf_and_gz_struct *restrict prims,
    CCTK_REAL **TUPmunu,
    CCTK_REAL *restrict gupxy,
    CCTK_REAL *restrict gupxz,
    CCTK_REAL *restrict gupyz,
    CCTK_REAL *restrict kxx,
    CCTK_REAL *restrict kxy,
    CCTK_REAL *restrict kxz,
    CCTK_REAL *restrict kyy,
    CCTK_REAL *restrict kyz,
    CCTK_REAL *restrict kzz,
    CCTK_REAL *restrict tau_rhs,
    CCTK_REAL *restrict s_tau
);

static void A_i_rhs_no_gauge_terms(
    const int A_dirn,
    const cGH *restrict cctkGH,
    const int *restrict cctk_lsh,
    const int *restrict cctk_nghostzones,
    gf_and_gz_struct *restrict out_prims_r,
    gf_and_gz_struct *restrict out_prims_l,
    CCTK_REAL *restrict phi_interped,
    CCTK_REAL *restrict cmax_1,
    CCTK_REAL *restrict cmin_1,
    CCTK_REAL *restrict cmax_2,
    CCTK_REAL *restrict cmin_2,
    CCTK_REAL *restrict A3_rhs
);

static void Lorenz_psi6phi_rhs__add_gauge_terms_to_A_i_rhs(
    const cGH *restrict cctkGH,
    const int *restrict cctk_lsh,
    const int *restrict cctk_nghostzones,
    CCTK_REAL *restrict dX,
    CCTK_REAL **in_vars,
    CCTK_REAL *restrict psi6phi,
    CCTK_REAL *restrict shiftx_iphjphkph,
    CCTK_REAL *restrict shifty_iphjphkph,
    CCTK_REAL *restrict shiftz_iphjphkph,
    CCTK_REAL *restrict alpha_iphjphkph,
    CCTK_REAL *restrict alpha_Phi_minus_betaj_A_j_iphjphkph,
    CCTK_REAL *restrict alpha_sqrtg_Ax_interp,
    CCTK_REAL *restrict alpha_sqrtg_Ay_interp,
    CCTK_REAL *restrict alpha_sqrtg_Az_interp,
    CCTK_REAL *restrict psi6phi_rhs,
    CCTK_REAL *restrict Ax_rhs,
    CCTK_REAL *restrict Ay_rhs,
    CCTK_REAL *restrict Az_rhs
);

static void add_fluxes_and_source_terms_to_hydro_rhss(
    const igm_eos_parameters eos,
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
    CCTK_REAL *restrict S_star_rhs,
    CCTK_REAL *restrict s_tau,
    CCTK_REAL *restrict s_sx,
    CCTK_REAL *restrict s_sy,
    CCTK_REAL *restrict s_sz,
    CCTK_REAL *restrict gxx_dx,
    CCTK_REAL *restrict gxy_dx,
    CCTK_REAL *restrict gxz_dx,
    CCTK_REAL *restrict gyy_dx,
    CCTK_REAL *restrict gyz_dx,
    CCTK_REAL *restrict gzz_dx,
    CCTK_REAL *restrict gxx_dy,
    CCTK_REAL *restrict gxy_dy,
    CCTK_REAL *restrict gxz_dy,
    CCTK_REAL *restrict gyy_dy,
    CCTK_REAL *restrict gyz_dy,
    CCTK_REAL *restrict gzz_dy,
    CCTK_REAL *restrict gxx_dz,
    CCTK_REAL *restrict gxy_dz,
    CCTK_REAL *restrict gxz_dz,
    CCTK_REAL *restrict gyy_dz,
    CCTK_REAL *restrict gyz_dz,
    CCTK_REAL *restrict gzz_dz
);

#endif /* DRIVER_EVALUATE_MHD_RHS_H_ */
