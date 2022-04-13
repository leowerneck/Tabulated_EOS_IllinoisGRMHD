#ifndef CON2PRIM_HELPERS_H_
#define CON2PRIM_HELPERS_H_

static inline void set_gammaDD_and_gammaUU_from_ADM_quantities( const CCTK_REAL *restrict adm_quantities,
                                                                CCTK_REAL gammaDD[3][3],
                                                                CCTK_REAL gammaUU[3][3] ) {

  // Set gamma_{ij}
  gammaDD[0][0] = adm_quantities[GXX];
  gammaDD[0][1] = gammaDD[1][0] = adm_quantities[GXY];
  gammaDD[0][2] = gammaDD[2][0] = adm_quantities[GXZ];
  gammaDD[1][1] = adm_quantities[GYY];
  gammaDD[1][2] = gammaDD[2][1] = adm_quantities[GYZ];
  gammaDD[2][2] = adm_quantities[GZZ];

  // Set gamma^{ij}
  gammaUU[0][0] = adm_quantities[GUPXX];
  gammaUU[0][1] = gammaUU[1][0] = adm_quantities[GUPXY];
  gammaUU[0][2] = gammaUU[2][0] = adm_quantities[GUPXZ];
  gammaUU[1][1] = adm_quantities[GUPYY];
  gammaUU[1][2] = gammaUU[2][1] = adm_quantities[GUPYZ];
  gammaUU[2][2] = adm_quantities[GUPZZ];

}

static inline void raise_or_lower_indices_3d( const CCTK_REAL *restrict vecD_or_U,
                                              const CCTK_REAL gammaUU_or_DD[3][3],
                                              CCTK_REAL *restrict vecU_or_D ) {
  for(int i=0;i<3;i++) {
    vecU_or_D[i] = 0;
    for(int j=0;j<3;j++) {
      vecU_or_D[i] += gammaUU_or_DD[i][j] * vecD_or_U[j];
    }
  }
}

inline CCTK_REAL simple_rel_err( const CCTK_REAL a, const CCTK_REAL b ) {
  if     ( a != 0.0 ) return( fabs(1.0 - b/a) );
  else if( b != 0.0 ) return( fabs(1.0 - a/b) );
  else                return(       0.0       );
}

#endif // CON2PRIM_HELPERS_H_
