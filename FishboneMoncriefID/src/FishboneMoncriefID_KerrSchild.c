#include "cctk.h"
#include "cctk_Parameters.h"

#include "FishboneMoncriefID.h"

void FishboneMoncriefID_KerrSchild(const cGH* restrict const cctkGH,const CCTK_INT *cctk_lsh,
                                   const CCTK_INT i0,const CCTK_INT i1,const CCTK_INT i2,
                                   const CCTK_REAL *xcoordGF,const CCTK_REAL *ycoordGF,const CCTK_REAL *zcoordGF,
                                   CCTK_REAL *alphaGF,CCTK_REAL *betaU0GF,CCTK_REAL *betaU1GF,CCTK_REAL *betaU2GF,
                                   CCTK_REAL *gammaDD00GF,CCTK_REAL *gammaDD01GF,CCTK_REAL *gammaDD02GF,CCTK_REAL *gammaDD11GF,CCTK_REAL *gammaDD12GF,CCTK_REAL *gammaDD22GF,
                                   CCTK_REAL     *KDD00GF,CCTK_REAL     *KDD01GF,CCTK_REAL     *KDD02GF,CCTK_REAL     *KDD11GF,CCTK_REAL     *KDD12GF,CCTK_REAL     *KDD22GF) {

  DECLARE_CCTK_PARAMETERS;

   /*
    * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
    */
   const double xcoord = xcoordGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
   const double ycoord = ycoordGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
   const double zcoord = zcoordGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
   /*
    * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
    */
   const double FDPart3_0 = ((a)*(a));
   const double FDPart3_1 = ((zcoord)*(zcoord));
   const double FDPart3_2 = ((xcoord)*(xcoord));
   const double FDPart3_3 = ((ycoord)*(ycoord));
   const double FDPart3_4 = FDPart3_2 + FDPart3_3;
   const double FDPart3_5 = FDPart3_1 + FDPart3_4;
   const double FDPart3_6 = (1.0/(FDPart3_5));
   const double FDPart3_7 = FDPart3_1*FDPart3_6;
   const double FDPart3_8 = FDPart3_0*FDPart3_7 + FDPart3_5;
   const double FDPart3_9 = (1.0/(FDPart3_8));
   const double FDPart3_10 = sqrt(FDPart3_5);
   const double FDPart3_11 = 2*FDPart3_10;
   const double FDPart3_12 = FDPart3_11*M;
   const double FDPart3_13 = FDPart3_12*FDPart3_9;
   const double FDPart3_14 = FDPart3_13 + 1;
   const double FDPart3_15 = (1.0/(FDPart3_14));
   const double FDPart3_16 = M*xcoord;
   const double FDPart3_19 = M*ycoord;
   const double FDPart3_20 = 2*FDPart3_15*FDPart3_9;
   const double FDPart3_22 = 1 - FDPart3_7;
   const double FDPart3_23 = -FDPart3_13 - 1;
   const double FDPart3_24 = FDPart3_0*((FDPart3_22)*(FDPart3_22))*((FDPart3_23)*(FDPart3_23));
   const double FDPart3_26 = FDPart3_22*(FDPart3_0*FDPart3_13*FDPart3_22 + FDPart3_0 + FDPart3_5);
   const double FDPart3_27 = FDPart3_14*FDPart3_26;
   const double FDPart3_28 = -FDPart3_24*FDPart3_8 + FDPart3_27*FDPart3_8;
   const double FDPart3_29 = -FDPart3_24 + FDPart3_27;
   const double FDPart3_30 = FDPart3_29*((FDPart3_8)*(FDPart3_8))/((FDPart3_28)*(FDPart3_28)*(FDPart3_28));
   const double FDPart3_31 = (1.0/(-FDPart3_24*FDPart3_30 + FDPart3_27*FDPart3_30));
   const double FDPart3_32 = (1.0/((FDPart3_28)*(FDPart3_28)));
   const double FDPart3_33 = FDPart3_29*FDPart3_31*FDPart3_32*FDPart3_8;
   const double FDPart3_34 = FDPart3_14*FDPart3_33;
   const double FDPart3_35 = xcoord*ycoord;
   const double FDPart3_36 = FDPart3_22/FDPart3_4;
   const double FDPart3_37 = (1.0/(FDPart3_10));
   const double FDPart3_39 = FDPart3_23*FDPart3_37*a;
   const double FDPart3_41 = 2*FDPart3_33*FDPart3_35*FDPart3_36*FDPart3_39;
   const double FDPart3_43 = (1.0/((FDPart3_4)*(FDPart3_4)));
   const double FDPart3_46 = (1.0/(FDPart3_22));
   const double FDPart3_47 = FDPart3_32*((FDPart3_8)*(FDPart3_8));
   const double FDPart3_48 = FDPart3_31*(-FDPart3_24*FDPart3_47 + FDPart3_27*FDPart3_47);
   const double FDPart3_50 = FDPart3_46*FDPart3_48/((FDPart3_5)*(FDPart3_5)*(FDPart3_5));
   const double FDPart3_51 = FDPart3_34*FDPart3_6;
   const double FDPart3_53 = FDPart3_26*FDPart3_33*FDPart3_43;
   const double FDPart3_56 = FDPart3_37*zcoord;
   const double FDPart3_58 = FDPart3_23*FDPart3_33*FDPart3_36*FDPart3_56*a;
   const double FDPart3_59 = pow(FDPart3_5, -3.0/2.0);
   const double FDPart3_60 = FDPart3_1*FDPart3_59;
   const double FDPart3_61 = FDPart3_37 - FDPart3_60;
   const double FDPart3_62 = FDPart3_46*FDPart3_48*FDPart3_59*FDPart3_61;
   const double FDPart3_64 = FDPart3_46*((FDPart3_61)*(FDPart3_61));
   const double FDPart3_65 = (1.0/((FDPart3_5)*(FDPart3_5)));
   const double FDPart3_66 = FDPart3_1*FDPart3_2*FDPart3_65;
   const double FDPart3_67 = sqrt(FDPart3_14);
   const double FDPart3_69 = 2*FDPart3_2;
   const double FDPart3_70 = 2*FDPart3_3;
   const double FDPart3_72 = acos(FDPart3_56);
   const double FDPart3_73 = FDPart3_0*cos(2*FDPart3_72);
   const double FDPart3_75 = FDPart3_0 + 2*FDPart3_1 + FDPart3_69 + FDPart3_70 + FDPart3_73;
   const double FDPart3_76 = FDPart3_67/(4*FDPart3_10*M + FDPart3_75);
   const double FDPart3_77 = FDPart3_76*M;
   const double FDPart3_79 = 4*FDPart3_46*FDPart3_77;
   const double FDPart3_80 = (1.0/(FDPart3_75));
   const double FDPart3_82 = 16*FDPart3_0*FDPart3_80;
   const double FDPart3_83 = FDPart3_77*FDPart3_82;
   const double FDPart3_84 = (1.0/((FDPart3_75)*(FDPart3_75)));
   const double FDPart3_85 = FDPart3_84*(FDPart3_0 - 2*FDPart3_1 - FDPart3_69 - FDPart3_70 + FDPart3_73);
   const double FDPart3_86 = FDPart3_36*FDPart3_67*FDPart3_85;
   const double FDPart3_87 = FDPart3_37*FDPart3_86*a;
   const double FDPart3_88 = 4*FDPart3_16*FDPart3_87*ycoord;
   const double FDPart3_89 = ((a)*(a)*(a));
   const double FDPart3_90 = FDPart3_16*FDPart3_76;
   const double FDPart3_91 = FDPart3_90*ycoord;
   const double FDPart3_92 = 16*FDPart3_36*FDPart3_60*FDPart3_80*FDPart3_89*FDPart3_91;
   const double FDPart3_93 = FDPart3_85*(FDPart3_12 + FDPart3_75);
   const double FDPart3_94 = 4*FDPart3_77*FDPart3_93;
   const double FDPart3_95 = ((a)*(a)*(a)*(a));
   const double FDPart3_98 = FDPart3_22*FDPart3_43*FDPart3_84*(4*FDPart3_0*FDPart3_5*(FDPart3_11 - M) + 4*FDPart3_10*FDPart3_73*(FDPart3_0 + FDPart3_10*(FDPart3_11 + M)) + 8*pow(FDPart3_5, 5.0/2.0) + FDPart3_95*(FDPart3_10 - M)*cos(4*FDPart3_72) + FDPart3_95*(3*FDPart3_10 + M));
   const double FDPart3_99 = FDPart3_10*FDPart3_77*FDPart3_98;
   const double FDPart3_102 = FDPart3_1*FDPart3_65*FDPart3_91;
   const double FDPart3_104 = 8*FDPart3_80*FDPart3_89;
   const double FDPart3_105 = FDPart3_104*FDPart3_60*FDPart3_77;
   const double FDPart3_106 = 4*FDPart3_6*FDPart3_93;
   const double FDPart3_107 = 8*FDPart3_0*FDPart3_80;
   const double FDPart3_109 = FDPart3_65*((zcoord)*(zcoord)*(zcoord));
   const double FDPart3_110 = FDPart3_56*FDPart3_61;
   const double FDPart3_111 = 4*FDPart3_110*FDPart3_46;
   const double FDPart3_112 = FDPart3_56*FDPart3_86*a;
   const double FDPart3_113 = FDPart3_19*FDPart3_76;
   const double FDPart3_114 = FDPart3_104*FDPart3_36*FDPart3_61*zcoord;
   const double FDPart3_116 = FDPart3_1*FDPart3_3*FDPart3_65;
   alphaGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = sqrt(FDPart3_15);
   betaU0GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 2*FDPart3_15*FDPart3_16*FDPart3_9;
   betaU1GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_19*FDPart3_20;
   betaU2GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_20*M*zcoord;
   gammaDD00GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_1*FDPart3_2*FDPart3_50 + FDPart3_2*FDPart3_34*FDPart3_6 + FDPart3_26*FDPart3_3*FDPart3_33*FDPart3_43 - FDPart3_41;
   gammaDD01GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_1*FDPart3_35*FDPart3_50 + FDPart3_2*FDPart3_33*FDPart3_36*FDPart3_39 - FDPart3_3*FDPart3_33*FDPart3_36*FDPart3_39 + FDPart3_35*FDPart3_51 - FDPart3_35*FDPart3_53;
   gammaDD02GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_51*xcoord*zcoord - FDPart3_58*ycoord - FDPart3_62*xcoord*zcoord;
   gammaDD11GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_1*FDPart3_3*FDPart3_50 + FDPart3_2*FDPart3_53 + FDPart3_3*FDPart3_51 + FDPart3_41;
   gammaDD12GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_51*ycoord*zcoord + FDPart3_58*xcoord - FDPart3_62*ycoord*zcoord;
   gammaDD22GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_34*FDPart3_7 + FDPart3_48*FDPart3_64;
   KDD00GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_2*FDPart3_6*FDPart3_94 + FDPart3_66*FDPart3_79 + FDPart3_66*FDPart3_83 + FDPart3_70*FDPart3_99 + FDPart3_88 + FDPart3_92;
   KDD01GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 4*FDPart3_102*FDPart3_46 + FDPart3_102*FDPart3_82 - FDPart3_105*FDPart3_2*FDPart3_36 + FDPart3_105*FDPart3_3*FDPart3_36 + FDPart3_106*FDPart3_91 - FDPart3_12*FDPart3_35*FDPart3_76*FDPart3_98 - FDPart3_69*FDPart3_87*M + FDPart3_70*FDPart3_87*M;
   KDD02GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_106*FDPart3_90*zcoord + FDPart3_107*FDPart3_109*FDPart3_90 - FDPart3_107*FDPart3_110*FDPart3_90 - FDPart3_111*FDPart3_90 + 2*FDPart3_112*FDPart3_19 - FDPart3_113*FDPart3_114;
   KDD11GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_116*FDPart3_79 + FDPart3_116*FDPart3_83 + FDPart3_3*FDPart3_6*FDPart3_94 + FDPart3_69*FDPart3_99 - FDPart3_88 - FDPart3_92;
   KDD12GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_106*FDPart3_113*zcoord + FDPart3_107*FDPart3_109*FDPart3_113 - FDPart3_107*FDPart3_110*FDPart3_113 - FDPart3_111*FDPart3_113 - 2*FDPart3_112*FDPart3_16 + FDPart3_114*FDPart3_90;
   KDD22GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = -FDPart3_1*FDPart3_37*FDPart3_61*FDPart3_83 + 4*FDPart3_5*FDPart3_64*FDPart3_77 + FDPart3_7*FDPart3_94;

}
