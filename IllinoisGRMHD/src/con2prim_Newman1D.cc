#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "IllinoisGRMHD_headers.h"
#include "con2prim_headers.h"

void newman( const igm_eos_parameters eos,
             const CCTK_REAL tol_x,
             const CCTK_REAL sqrt_detgamma,
             const CCTK_REAL S_squared,
             const CCTK_REAL BdotS,
             const CCTK_REAL B_squared,
             const CCTK_REAL *restrict SU,
             const CCTK_REAL *restrict con,
             CCTK_REAL *restrict prim,
             bool& c2p_failed );

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

static inline CCTK_REAL compute_sqrt_detgamma( const CCTK_REAL gammaDD[3][3] ) {
  /*               / a b c \
  ** gamma_{ij} = |  b d e |
  **               \ c e f /
  */
  const CCTK_REAL a = gammaDD[0][0];
  const CCTK_REAL b = gammaDD[0][1];
  const CCTK_REAL c = gammaDD[0][2];
  const CCTK_REAL d = gammaDD[1][1];
  const CCTK_REAL e = gammaDD[1][2];
  const CCTK_REAL f = gammaDD[2][2];

  // Compute the deterrminant
  const CCTK_REAL detgamma = a*(d*f-e*e) + b*(c*e-b*f) + c*(b*e-d*c);

  // Return sqrt(detgamma)
  return( sqrt(detgamma) );
}

int con2prim_Newman1D( const igm_eos_parameters eos,
                       const CCTK_REAL *restrict adm_quantities,
                       const CCTK_REAL *restrict con,
                       CCTK_REAL *restrict prim,
                       output_stats& stats ) {

  // Set gamma_{ij} and gamma^{ij}
  CCTK_REAL gammaDD[3][3],gammaUU[3][3];
  set_gammaDD_and_gammaUU_from_ADM_quantities(adm_quantities,gammaDD,gammaUU);

  // Get sqrt(detgamma)
  CCTK_REAL sqrt_detgamma = compute_sqrt_detgamma(gammaDD);

  // Read in B^{i}
  CCTK_REAL BU[3];
  for(int i=0;i<3;i++) BU[i] = con[B1_con+i];

  // Read in S_{i}
  CCTK_REAL SD[3];
  for(int i=0;i<3;i++) SD[i] = con[S1_cov+i];

  // Compute B_{i} = gamma_{ij}B^{j}
  CCTK_REAL BD[3];
  raise_or_lower_indices_3d( BU,gammaDD, BD );

  // Compute S^{i} = gamma^{ij}S_{j}
  CCTK_REAL SU[3];
  raise_or_lower_indices_3d( SD,gammaUU, SU );

  // S^2 = S^i * S_i
  CCTK_REAL S_squared = 0.0;
  for(int i=0;i<3;i++) S_squared += SU[i] * SD[i];

  // Need to calculate for (21) and (22) in Cerda-Duran 2008
  // B * S = B^i * S_i
  CCTK_REAL BdotS = 0.0;
  for(int i=0;i<3;i++) BdotS += BU[i] * SD[i];
  
  // B^2 = B^i * B_i
  CCTK_REAL B_squared = 0.0;
  for(int i=0;i<3;i++) B_squared += BU[i] * BD[i];

  const CCTK_REAL tol_x = 5e-9;
  bool c2p_failed = false;
  newman(eos,tol_x,sqrt_detgamma,S_squared,BdotS,B_squared,SU,con,prim,c2p_failed);
  
  return( c2p_failed );
}

void newman( const igm_eos_parameters eos,
             const CCTK_REAL tol_x,
             const CCTK_REAL sqrt_detgamma,
             const CCTK_REAL S_squared,
             const CCTK_REAL BdotS,
             const CCTK_REAL B_squared,
             const CCTK_REAL *restrict SU,
             const CCTK_REAL *restrict con,
             CCTK_REAL *restrict prim,
             bool& c2p_failed ) {
 
  bool conacc = 0;
  // 2) Compute useful temporary variables
  const CCTK_REAL cE    = con[TAU] + con[DD];
  const CCTK_REAL cMsqr = S_squared/sqrt_detgamma;
  const CCTK_REAL cBsqr = B_squared;
  const CCTK_REAL cT    = BdotS/sqrt_detgamma;
  CCTK_REAL prec        = tol_x;

  // We always guess W = 1, since we don't keep
  // track of the velocities in between time
  // steps and therefore do not have a good
  // guess for them
  CCTK_REAL xrho  = con[DD];
  CCTK_REAL xye   = con[YE]/con[DD];
  CCTK_REAL xtemp = eos.T_atm;
  CCTK_REAL cP    = 0.0;
  CCTK_REAL xeps  = 0.0;
  get_P_and_eps_from_rho_Ye_and_T( eos,xrho,xye,xtemp, &cP,&xeps );

  //Now, begin iterative procedure to derive the primitive variables
  CCTK_REAL cPold = cP; // -> setup pressure initial guess
  int step=0;
  CCTK_REAL cL;

  const int maxsteps = 300;
  CCTK_REAL AtP[maxsteps]; // not length 3 in case of extrap. probs
  CCTK_REAL AtR;
  int AtStep=0;
  AtP[0]=cP;

  // Leo mod: compute auxiliary variables so we can find eps
  const CCTK_REAL q = con[TAU]/con[DD];
  const CCTK_REAL s = B_squared/con[DD];
  const CCTK_REAL t = BdotS/(pow(con[DD],1.5));
  bool use_entropy  = false;
  
  CCTK_REAL cD = 0.5*(cMsqr*cBsqr-cT*cT); 
  if(cD<0.) cD=0.;

  do {
    cPold = cP;
    step++;
    //This check is required to ensure that |cos(cPhi)|<1 below
    CCTK_REAL cA = cE + cP + cBsqr/2.;
    //Check that the values are physical...
    //What's the right thing to do if cD<0.?
    CCTK_REAL cPhi = acos(sqrt(27./4.*cD/cA)/cA);
    CCTK_REAL cEps = cA/3.*(1.-2.*cos(2./3.*cPhi+2./3.*M_PI));
    cL = cEps-cBsqr;

    CCTK_REAL cVsqr = (cMsqr*cL*cL+cT*cT*(cBsqr+2.*cL))/(cL*(cL+cBsqr)*cL*(cL+cBsqr));
    const CCTK_REAL cWsqr = 1./(1.-cVsqr);
    CCTK_REAL cH = cL/cWsqr;

    // We'll need a 1D root solve to find the temperature from density 
    // and enthalpy.
    CCTK_REAL W    = sqrt(cWsqr);
    CCTK_REAL xrho = con[RHO]/W;
    CCTK_REAL xye  = con[YE]/con[DD];
    CCTK_REAL xprs = 0.0;
    CCTK_REAL xent = 0.0;
    CCTK_REAL xeps = 0.0;

    prim[RHO     ] = con[DD]/W; //    rho[s] = tildeD[s]/(sqrtDetg[s]*W[s]);
    prim[WLORENTZ] = W;
    
    // define other dummy variables as needed for input arguments
    // EOS_press_from_rhoenthalpy(eoskey,keytemp,precEOS,prim[RHO],&xeps,&xtemp,xye,&xprs,&cH,&anyerr,&keyerr);
    if( use_entropy == true ) {
      xent = con[WS]/W;
      get_P_eps_and_T_from_rho_Ye_and_S( eos,xrho,xye,xent, &xprs,&xeps,&xtemp );
    }
    else {
      const CCTK_REAL x = cH * W / xrho;
      xeps = W - 1.0 + (1.0-W*W)*x/W + W*(q - s + t*t/(2*x*x) + s/(2*W*W)  );
      xeps=fmax(xeps, eos.eps_min);
      get_P_S_and_T_from_rho_Ye_and_eps( eos, xrho,xye,xeps, &xprs,&xent,&xtemp );
    }

    cP = xprs;
    prim[EPS] = xeps;

#if DEBUG
    printf("xprs= %e,cPhi = %e, cEps = %e, cH = %e, cVsqr = %e, cL = %e, cP = %e,prim[WLORENTZ]=%e, prim[RHO]=%e xeps = %e\n", xprs, cPhi, cEps, cH, cVsqr, cL, cP, prim[WLORENTZ], prim[RHO],xeps); 
#endif
 
#if DEBUG
    printf("At step %i, P = %0.12e, dP = %e, cP+cPold = %e, threshold = %e \n",step,cP,fabs(cP-cPold),fabs(cP+cPold),fabs(cP-cPold)/(cP+cPold));
#endif

    AtStep++;
    AtP[AtStep]=cP;

    if(AtStep>=2) {   //Aitken extrapolation

      AtR = (AtP[AtStep]-AtP[AtStep-1])/(AtP[AtStep-1]-AtP[AtStep-2]);
#if DEBUG    
      printf("At step %i, AtP[0] = %e, AtP[1] = %e, AtP[2]=%e, AtR=%e \n",step,AtP[0], AtP[1], AtP[2],AtR);
#endif
      if(AtR<1. && AtR>0.) {
        cP=AtP[AtStep-1]+(AtP[AtStep]-AtP[AtStep-1])/(1.-AtR);
#if DEBUG         
        printf("At step in acceleration  %i, P = %e, dP = %e \n",AtStep,cP,fabs(cP-cPold));
#endif        
        AtStep=0;
        conacc = 1;
        AtP[0]=cP;   //starting value for next Aitken extrapolation
      }
    }
  }
  while(fabs(cP-cPold)>prec*(cP+cPold) && step<maxsteps);

  if (conacc==1) {     //converged on an extrap. so recompute vars
    const CCTK_REAL cA = cE + cP + cBsqr/2.;
    const CCTK_REAL cPhi = acos(sqrt(27./4.*cD/cA)/cA);
    const CCTK_REAL cEps = cA/3.*(1.-2.*cos(2./3.*cPhi+2./3.*M_PI));
    cL = cEps-cBsqr;
    const CCTK_REAL cVsqr = 
      (cMsqr*cL*cL+cT*cT*(cBsqr+2.*cL))/(cL*(cL+cBsqr)*cL*(cL+cBsqr));
    if(cVsqr<0. || cVsqr>=1.) {
    }

    const CCTK_REAL cWsqr = 1./(1.-cVsqr);

    prim[WLORENTZ] = sqrt(cWsqr); //    W[s] = sqrt(cWsqr);
    prim[RHO     ] = con[DD]/prim[WLORENTZ]; //    rho[s] = tildeD[s]/(sqrtDetg[s]*W[s]);
    prim[PRESS   ] = cP;

    prim[B1_con  ] = con[B1_con];
    prim[B2_con  ] = con[B2_con];
    prim[B3_con  ] = con[B3_con];
  }

  //Compute v^i
  const CCTK_REAL cS = cT/cL;
   
  prim[B1_con] = con[B1_con];
  prim[B2_con] = con[B2_con];
  prim[B3_con] = con[B3_con];

  prim[UTCON1] = prim[WLORENTZ]*(cS * con[B1_con] + SU[0]) / (cL+cBsqr); // v[i][s] = cS*evolvedVariables.tildeB(i)[s];
  prim[UTCON2] = prim[WLORENTZ]*(cS * con[B2_con] + SU[1]) / (cL+cBsqr); // v[i][s] = cS*evolvedVariables.tildeB(i)[s];
  prim[UTCON3] = prim[WLORENTZ]*(cS * con[B3_con] + SU[2]) / (cL+cBsqr); // v[i][s] = cS*evolvedVariables.tildeB(i)[s];

  prim[YE  ] = con[YE]/con[DD];
  prim[RHO ] = con[DD]/prim[WLORENTZ];
  prim[ENT ] = con[WS]/prim[WLORENTZ];
  get_P_eps_and_T_from_rho_Ye_and_S( eos, prim[RHO],prim[YE],prim[ENT],
                                     &prim[PRESS],&prim[EPS],&prim[TEMP] );

  if (step >= maxsteps) c2p_failed = true;

}
