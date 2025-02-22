#include "GRHayLMHD.h"

#include <assert.h>

#define GRHAYLMHD_COMPUTE_ERROR_CON(name,con) \
    errors_numer_##name = fabs(cons.con - cons_orig.con); \
    errors_denom_##name = fabs(cons_orig.con);

#define GRHAYLMHD_COMPUTE_ERRORS \
    GRHAYLMHD_COMPUTE_ERROR_CON(rho, rho) \
    GRHAYLMHD_COMPUTE_ERROR_CON(tau, tau) \
    GRHAYLMHD_COMPUTE_ERROR_CON(S_x, SD[0]) \
    GRHAYLMHD_COMPUTE_ERROR_CON(S_y, SD[1]) \
    GRHAYLMHD_COMPUTE_ERROR_CON(S_z, SD[2]) \
    GRHAYLMHD_COMPUTE_ERROR_CON(Y_e, Y_e) \
    GRHAYLMHD_COMPUTE_ERROR_CON(entropy, entropy)


// TODO: add more diagnostic information
static void print_diagnostics(const int imin,
                              const int imax,
                              const int jmin,
                              const int jmax,
                              const int kmin,
                              const int kmax,
                              const int iteration,
                              const int reflevel,
                              const int atm_resets,
                              const ghl_conservative_quantities *errors)
{
    const int ntotal = (imax - imin) * (jmax - jmin) * (kmax - kmin);
    CCTK_VINFO("Con2Prim -- It. %d -- Ref. Lev. %d -- ATM Resets: %d / %d",
               iteration,
               reflevel,
               atm_resets,
               ntotal);
    CCTK_VINFO("  Errors -- rho %.2e -- tau %.2e -- S_x %.2e -- S_y %.2e -- S_z %.2e -- Y_e %.2e -- ent %.2e",
               errors->rho,
               errors->tau,
               errors->SD[0],
               errors->SD[1],
               errors->SD[2],
               errors->Y_e,
               errors->entropy);
}

void GRHayLMHD_Con2Prim(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // TODO: move this to initialization
    assert(ghl_params->calc_prim_guess == true);

    const int imin = cctk_nghostzones[0];
    const int jmin = cctk_nghostzones[1];
    const int kmin = cctk_nghostzones[2];

    const int imax = cctk_lsh[0] - cctk_nghostzones[0];
    const int jmax = cctk_lsh[1] - cctk_nghostzones[1];
    const int kmax = cctk_lsh[2] - cctk_nghostzones[2];

    int atm_resets = 0;

    CCTK_REAL errors_numer_rho = 0.0, errors_numer_tau = 0.0, errors_numer_S_x = 0.0, errors_numer_S_y = 0.0, errors_numer_S_z = 0.0, errors_numer_Y_e = 0.0, errors_numer_entropy = 0.0;
    CCTK_REAL errors_denom_rho = 0.0, errors_denom_tau = 0.0, errors_denom_S_x = 0.0, errors_denom_S_y = 0.0, errors_denom_S_z = 0.0, errors_denom_Y_e = 0.0, errors_denom_entropy = 0.0;


#define ERROR_ARGS \
    errors_numer_rho, errors_numer_tau, errors_numer_S_x, errors_numer_S_y, errors_numer_S_z, errors_numer_Y_e, errors_numer_entropy, \
    errors_denom_rho, errors_denom_tau, errors_denom_S_x, errors_denom_S_y, errors_denom_S_z, errors_denom_Y_e, errors_denom_entropy

#pragma omp parallel for reduction(+ : atm_resets, ERROR_ARGS)
    LOOP3D(imin, imax, jmin, jmax, kmin, kmax)
    {
        ghl_con2prim_diagnostics diagnostics = { 0 };

        ghl_metric_quantities adm_metric = { 0 };
        GRHAYLMHD_LOAD_METRIC_ENFORCE_DETGAMMAEQ1(adm_metric);

        ghl_ADM_aux_quantities aux_metric = { 0 };
        ghl_compute_ADM_auxiliaries(&adm_metric, &aux_metric);

        ghl_conservative_quantities cons_orig = { 0 };
        GRHAYLMHD_LOAD_CONS(cons_orig);
        ghl_conservative_quantities cons = cons_orig;

        ghl_conservative_quantities cons_undens = { 0 };
        ghl_undensitize_conservatives(adm_metric.sqrt_detgamma, &cons, &cons_undens);

        ghl_primitive_quantities prims = { 0 };
        prims.BU[0]                    = Bvec[ijkx];
        prims.BU[1]                    = Bvec[ijky];
        prims.BU[2]                    = Bvec[ijkz];

        ghl_error_codes_t error = ghl_con2prim_tabulated_multi_method(ghl_params,
                                                                      ghl_eos,
                                                                      &adm_metric,
                                                                      &aux_metric,
                                                                      &cons_undens,
                                                                      &prims,
                                                                      &diagnostics);

        // TODO: averaging algorithm
        if(error) {
            atm_resets++;
            ghl_set_prims_to_constant_atm(ghl_eos, &prims);
        }

        bool speed_limited = false;
        ghl_enforce_primitive_limits_and_compute_u0(ghl_params,
                                                    ghl_eos,
                                                    &adm_metric,
                                                    &prims,
                                                    &speed_limited);

        ghl_compute_conservs(&adm_metric, &aux_metric, &prims, &cons);

        GRHAYLMHD_WRITE_PRIMS(prims);
        GRHAYLMHD_WRITE_CONS(cons);
        GRHAYLMHD_COMPUTE_ERRORS;
    }
    ENDLOOP3D

    ghl_conservative_quantities cons_errors = { 0 };
    cons_errors.rho     = errors_denom_rho != 0.0 ? errors_numer_rho / errors_denom_rho : 0.0;
    cons_errors.tau     = errors_denom_tau != 0.0 ? errors_numer_tau / errors_denom_tau : 0.0;
    cons_errors.SD[0]   = errors_denom_S_x != 0.0 ? errors_numer_S_x / errors_denom_S_x : 0.0;
    cons_errors.SD[1]   = errors_denom_S_y != 0.0 ? errors_numer_S_y / errors_denom_S_y : 0.0;
    cons_errors.SD[2]   = errors_denom_S_z != 0.0 ? errors_numer_S_z / errors_denom_S_z : 0.0;
    cons_errors.Y_e     = errors_denom_Y_e != 0.0 ? errors_numer_Y_e / errors_denom_Y_e : 0.0;
    cons_errors.entropy = errors_denom_entropy != 0.0 ? errors_numer_entropy / errors_denom_entropy : 0.0;

    print_diagnostics(imin,
                      imax,
                      jmin,
                      jmax,
                      kmin,
                      kmax,
                      cctk_iteration,
                      GetRefinementLevel(cctkGH),
                      atm_resets,
                      &cons_errors);
}
