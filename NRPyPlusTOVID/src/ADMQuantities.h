
static inline
void ADMQuantities(const cGH* restrict const cctkGH, const CCTK_INT i0,const CCTK_INT i1,const CCTK_INT i2,

                   const CCTK_REAL *restrict xx0GF,const CCTK_REAL *restrict xx1GF,const CCTK_REAL *restrict xx2GF,

                   const CCTK_REAL IDalpha,
                   const CCTK_REAL IDgammaDD00,const CCTK_REAL IDgammaDD01, const CCTK_REAL IDgammaDD02,
                   const CCTK_REAL IDgammaDD11,const CCTK_REAL IDgammaDD12, const CCTK_REAL IDgammaDD22,

                   CCTK_REAL *alphaGF,CCTK_REAL *betaU0GF,CCTK_REAL *betaU1GF,CCTK_REAL *betaU2GF,
                   CCTK_REAL *gammaDD00GF, CCTK_REAL *gammaDD01GF, CCTK_REAL *gammaDD02GF,
                   CCTK_REAL *gammaDD11GF, CCTK_REAL *gammaDD12GF, CCTK_REAL *gammaDD22GF,
                   CCTK_REAL *KDD00GF, CCTK_REAL *KDD01GF, CCTK_REAL *KDD02GF,
                   CCTK_REAL *KDD11GF, CCTK_REAL *KDD12GF, CCTK_REAL *KDD22GF) {
    const CCTK_REAL xx0 = xx0GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)];
    const CCTK_REAL xx1 = xx1GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)];
    const CCTK_REAL xx2 = xx2GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)];

   /*
    * NRPy+ Finite Difference Code Generation, Step 1 of 1: Evaluate SymPy expressions and write to main memory:
    */
   const double FDPart3_0 = ((xx1)*(xx1));
   const double FDPart3_1 = ((xx0)*(xx0));
   const double FDPart3_2 = FDPart3_0 + FDPart3_1;
   const double FDPart3_3 = IDgammaDD22/((FDPart3_2)*(FDPart3_2));
   const double FDPart3_4 = ((xx2)*(xx2));
   const double FDPart3_5 = FDPart3_2 + FDPart3_4;
   const double FDPart3_6 = (1.0/(FDPart3_5));
   const double FDPart3_7 = FDPart3_6*IDgammaDD00;
   const double FDPart3_8 = (1.0/(FDPart3_2));
   const double FDPart3_9 = (1.0/sqrt(FDPart3_5));
   const double FDPart3_11 = FDPart3_8*FDPart3_9*IDgammaDD02;
   const double FDPart3_12 = xx0*xx1;
   const double FDPart3_13 = 2*FDPart3_12;
   const double FDPart3_15 = (1.0/((FDPart3_5)*(FDPart3_5)));
   const double FDPart3_17 = -FDPart3_4*FDPart3_6 + 1;
   const double FDPart3_18 = (1.0/sqrt(FDPart3_17));
   const double FDPart3_19 = FDPart3_18*IDgammaDD01;
   const double FDPart3_20 = 2*FDPart3_19*xx2;
   const double FDPart3_22 = IDgammaDD11/FDPart3_17;
   const double FDPart3_23 = FDPart3_22*FDPart3_4/((FDPart3_5)*(FDPart3_5)*(FDPart3_5));
   const double FDPart3_25 = FDPart3_18*FDPart3_8*IDgammaDD12;
   const double FDPart3_26 = pow(FDPart3_5, -3.0/2.0);
   const double FDPart3_27 = FDPart3_26*xx2;
   const double FDPart3_29 = FDPart3_13*FDPart3_25*FDPart3_27;
   const double FDPart3_35 = -FDPart3_26*FDPart3_4 + FDPart3_9;
   const double FDPart3_36 = FDPart3_35*xx1;
   alphaGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = IDalpha;
   betaU0GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 0;
   betaU1GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 0;
   betaU2GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 0;
   gammaDD00GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_0*FDPart3_3 + FDPart3_1*FDPart3_15*FDPart3_20 + FDPart3_1*FDPart3_23 + FDPart3_1*FDPart3_7 - FDPart3_11*FDPart3_13 - FDPart3_29;
   gammaDD01GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = -FDPart3_0*FDPart3_11 - FDPart3_0*FDPart3_25*FDPart3_27 + FDPart3_1*FDPart3_18*FDPart3_27*FDPart3_8*IDgammaDD12 + FDPart3_1*FDPart3_8*FDPart3_9*IDgammaDD02 + FDPart3_12*FDPart3_23 - FDPart3_12*FDPart3_3 + FDPart3_12*FDPart3_7 + FDPart3_13*FDPart3_15*FDPart3_19*xx2;
   gammaDD02GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = -FDPart3_11*xx1*xx2 + FDPart3_15*FDPart3_19*FDPart3_4*xx0 - FDPart3_19*FDPart3_35*FDPart3_9*xx0 - FDPart3_22*FDPart3_27*FDPart3_35*xx0 + FDPart3_25*FDPart3_36 + FDPart3_7*xx0*xx2;
   gammaDD11GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_0*FDPart3_15*FDPart3_20 + FDPart3_0*FDPart3_23 + FDPart3_0*FDPart3_7 + FDPart3_1*FDPart3_3 + FDPart3_11*FDPart3_13 + FDPart3_29;
   gammaDD12GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_11*xx0*xx2 + FDPart3_15*FDPart3_19*FDPart3_4*xx1 - FDPart3_19*FDPart3_36*FDPart3_9 - FDPart3_22*FDPart3_27*FDPart3_36 - FDPart3_25*FDPart3_35*xx0 + FDPart3_7*xx1*xx2;
   gammaDD22GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = -FDPart3_20*FDPart3_35*FDPart3_9 + FDPart3_22*((FDPart3_35)*(FDPart3_35)) + FDPart3_4*FDPart3_6*IDgammaDD00;
   KDD00GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 0;
   KDD01GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 0;
   KDD02GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 0;
   KDD11GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 0;
   KDD12GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 0;
   KDD22GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 0;

}
    