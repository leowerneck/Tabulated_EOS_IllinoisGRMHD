
static inline
void HydroQuantities(const cGH* restrict const cctkGH, const CCTK_INT i0,const CCTK_INT i1,const CCTK_INT i2,

                   const CCTK_REAL IDPressure, const CCTK_REAL IDrho_baryonic,
                   const CCTK_REAL IDrho__total_energy_density,

                   CCTK_REAL *PressureGF,CCTK_REAL *rho_baryonicGF,
                   CCTK_REAL *epsilonGF,
                   CCTK_REAL *Valencia3velocityU0GF,
                   CCTK_REAL *Valencia3velocityU1GF,
                   CCTK_REAL *Valencia3velocityU2GF) {
    DECLARE_CCTK_PARAMETERS;
    if(IDrho__total_energy_density <= 0 || IDrho_baryonic <= 0 || IDPressure <= 0) {
        rho_baryonicGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = rho_atmosphere;
        PressureGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)]     = K_atmosphere*pow(rho_atmosphere,Gamma_atmosphere);
        epsilonGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)]      = 0;
        Valencia3velocityU0GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 0;
        Valencia3velocityU1GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 0;
        Valencia3velocityU2GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 0;
    } else {
      /*
       * NRPy+ Finite Difference Code Generation, Step 1 of 1: Evaluate SymPy expressions and write to main memory:
       */
      PressureGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = IDPressure;
      rho_baryonicGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = IDrho_baryonic;
      epsilonGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = IDrho__total_energy_density/IDrho_baryonic - 1;
      Valencia3velocityU0GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 0;
      Valencia3velocityU1GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 0;
      Valencia3velocityU2GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 0;

        // Apply pressure depletion.
        PressureGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] *= (1.0 - Pressure_depletion_factor);
    }
}
    