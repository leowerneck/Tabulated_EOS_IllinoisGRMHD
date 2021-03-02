#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "IllinoisGRMHD_headers.h"
#include "con2prim_headers.h"

static const CCTK_INT senergyvar = 0;
static const CCTK_INT entropyvar = 1;

void newman( const igm_eos_parameters eos,
             const CCTK_REAL tol_x,
             const CCTK_REAL S_squared,
             const CCTK_REAL BdotS,
             const CCTK_REAL B_squared,
             const CCTK_REAL *restrict SU,
             const CCTK_REAL *restrict con,
             CCTK_REAL *restrict prim,
             CCTK_INT& got_temp_from,
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

static inline CCTK_REAL simple_rel_err( const CCTK_REAL a, const CCTK_REAL b ) {
  if     ( a != 0.0 ) return( fabs(1.0 - b/a) );
  else if( b != 0.0 ) return( fabs(1.0 - a/b) );
  else                return(       0.0       );
}

int con2prim_Newman1D( const igm_eos_parameters eos,
                       const CCTK_REAL *restrict adm_quantities,
                       const CCTK_REAL *restrict con,
                       CCTK_REAL *restrict prim,
                       output_stats& stats ) {

  // Set gamma_{ij} and gamma^{ij}
  CCTK_REAL gammaDD[3][3],gammaUU[3][3];
  set_gammaDD_and_gammaUU_from_ADM_quantities(adm_quantities,gammaDD,gammaUU);

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

  // Enforce ceiling on S^{2} (A5 of Palenzuela et al. https://arxiv.org/pdf/1505.01607.pdf)
  CCTK_REAL S_squared_max = SQR(con[DD] + con[TAU]);
  if( S_squared > 0.9999 * S_squared_max ) {
    // Compute rescaling factor
    CCTK_REAL rescale_factor_must_be_less_than_one = sqrt(0.9999*S_squared_max/S_squared);
    // Rescale S_{i}
    for(int i=0;i<3;i++) SD[i] *= rescale_factor_must_be_less_than_one;
    // S_{i} has been rescaled. Recompute S^{i}.
    raise_or_lower_indices_3d( SD,gammaUU, SU );
    // Now recompute S^{2} := gamma^{ij}S^{i}S_{j}.
    S_squared = 0.0;
    for(int i=0;i<3;i++) S_squared += SU[i] * SD[i];
    // Check if the fix was successful
    if( simple_rel_err(S_squared,0.9999*S_squared_max) > 1e-12 ) CCTK_VError(VERR_DEF_PARAMS,"Incompatible values of S_squared after rescaling: %.15e %.15e\n",S_squared,0.9999*S_squared_max);
  }

  // Need to calculate for (21) and (22) in Cerda-Duran 2008
  // B * S = B^i * S_i
  CCTK_REAL BdotS = 0.0;
  for(int i=0;i<3;i++) BdotS += BU[i] * SD[i];
  
  // B^2 = B^i * B_i
  CCTK_REAL B_squared = 0.0;
  for(int i=0;i<3;i++) B_squared += BU[i] * BD[i];

  const CCTK_REAL tol_x = 1e-10;
  bool c2p_failed = false;
  CCTK_INT got_temp_from = None;
  newman(eos,tol_x,S_squared,BdotS,B_squared,SU,con,prim,got_temp_from,c2p_failed);

  if( c2p_failed && ( got_temp_from == entropyvar ) ) {
    // If the recovery failed when using the entropy, try
    // again using the specific internal energy.
    got_temp_from = senergyvar;
    newman(eos,tol_x,S_squared,BdotS,B_squared,SU,con,prim,got_temp_from,c2p_failed);
  }
  
  return( c2p_failed );
}

void newman( const igm_eos_parameters eos,
             const CCTK_REAL tol_x,
             const CCTK_REAL S_squared,
             const CCTK_REAL BdotS,
             const CCTK_REAL B_squared,
             const CCTK_REAL *restrict SU,
             const CCTK_REAL *restrict con,
             CCTK_REAL *restrict prim,
             CCTK_INT& got_temp_from,
             bool& c2p_failed ) {
 
  bool conacc    = false;
  CCTK_REAL prec = tol_x;

  // We always guess W = 1, since we don't keep
  // track of the velocities in between time
  // steps and therefore do not have a good
  // guess for them
  CCTK_REAL xrho   = con[DD];
  CCTK_REAL xye    = con[YE]/con[DD];
  CCTK_REAL xtemp  = eos.T_atm;
  CCTK_REAL xprs   = 0.0;
  CCTK_REAL xeps   = 0.0;
  CCTK_REAL xent   = 0.0;
  CCTK_REAL depsdT = 0.0;
  bool use_entropy = false;

  if( got_temp_from == None ) {

    // Only need to compute depsdT if we are not sure whether or not to use the entropy
    enforce_table_bounds_rho_Ye_T( eos,&xrho,&xye,&xtemp );
    WVU_EOS_P_eps_and_depsdT_from_rho_Ye_T( xrho,xye,xtemp, &xprs,&xeps,&depsdT );

    if( eos.evolve_entropy && (depsdT < eos.depsdT_threshold) ) {
      // We should use the entropy instead
      got_temp_from = entropyvar;
      use_entropy   = true;
    }
    else {
      got_temp_from = senergyvar;
    }

  }
  else {
    // If got_temp_from != None, then that means we tried to use
    // the entropy and the recovery failed. So it is safe to
    // skip the depsdT check and use a more efficient EOS call.
    enforce_table_bounds_rho_Ye_T( eos,&xrho,&xye,&xtemp );
    WVU_EOS_P_and_eps_from_rho_Ye_T( xrho,xye,xtemp, &xprs,&xeps );
    
  }

  //Now, begin iterative procedure to derive the primitive variables
  CCTK_REAL P_old = xprs; // -> setup pressure initial guess
  int step=0;

  const int maxsteps = 300;
  CCTK_REAL AtP[maxsteps]; // not length 3 in case of extrap. probs
  CCTK_REAL AtR;
  int AtStep=0;
  AtP[0]=xprs;

  // Leo mod: compute auxiliary variables so we can find eps
  const CCTK_REAL q = con[TAU]/con[DD];
  const CCTK_REAL s = B_squared/con[DD];
  const CCTK_REAL t = BdotS/(pow(con[DD],1.5));

  // d = 0.5( S^{2}*B^{2} - (B.S)^{2} ) (eq. 5.7 in Newman & Hamlin 2014)
  const CCTK_REAL d = fmax(0.5*(S_squared*B_squared-BdotS*BdotS),0.0);
  // e = tau + D
  const CCTK_REAL e = con[TAU] + con[DD];
  // z = rho*h*W^{2} = D*h*W.
  // Initialize to zero, for now.
  CCTK_REAL z = 0;

  // Useful auxiliary variable
  const CCTK_REAL BdotSsq = BdotS * BdotS;

  do {
    P_old = xprs;
    step++;
    //This check is required to ensure that |cos(cPhi)|<1 below

    // a = e + P + B^{2}/2 (eq. 5.7 in Newman & Hamlin 2014)
    CCTK_REAL a = e + xprs + 0.5*B_squared;

    // phi = acos( sqrt(27d/4a^{3}) ) (eq. 5.10 in Newman & Hamlin 2014)
    CCTK_REAL phi = acos( sqrt(27.0*d/(4.0*a))/a );
    // Eps = a/3 - 2a/3 cos( 2phi/3 + 2pi/3 )
    CCTK_REAL Eps = (a/3.0)*(1.0 - 2.0*cos((2.0/3.0)*(phi + M_PI)));

    // From the definition of Eps: Eps = z + B^{2} => z = Eps - B^{2}
    z = Eps-B_squared;

    // Now compute (eq. 5.2 in Newman & Hamlin 2014)
    //
    // v^{2} = ( z^{2}S^{2} + (2z + B^{2})(B.S)^{2}/(z^{2}(z+B^{2})^{2})
    const CCTK_REAL zBsq   = z + B_squared;
    const CCTK_REAL zBsqsq = zBsq*zBsq;
    const CCTK_REAL zsq    = z*z;
    const CCTK_REAL vsq    = (zsq * S_squared + (z+zBsq)*BdotSsq)/(zsq*zBsqsq);

    // Then compute W^{2} = 1/(1-v^{2})
    const CCTK_REAL Wsq    = 1.0/(1.0-vsq);

    // Now set W
    CCTK_REAL W = sqrt(Wsq);

    // Impose physical limits on W
    W = fmin(fmax(W,1.0),eos.W_max);
    
    // Then compute rho = D/W
    CCTK_REAL xrho = con[RHO]/W;

    // Initialize P, eps, and S to zero
    xprs = 0.0;
    xeps = 0.0;
    xent = 0.0;

    prim[RHO     ] = con[DD]/W; // rho[s] = tildeD[s]/(sqrtDetg[s]*W[s]);
    prim[WLORENTZ] = W;
    
    if( use_entropy == true ) {
      // If using the entropy, compute S = (WS)/W
      xent = con[WS]/W;
      // Then compute P, eps, and T using (rho,Ye,S)
      enforce_table_bounds_rho_Ye_S( eos,&xrho,&xye,&xent );
      WVU_EOS_P_eps_and_T_from_rho_Ye_S( xrho,xye,xent, &xprs,&xeps,&xtemp );
    }
    else {
      // If using the specific internal energy, remember that
      //
      // z = rho * h * W^{2} = D * h * W
      //
      // and therefore:
      //
      // x = h * W = z / D
      const CCTK_REAL x = z / con[DD];

      // Now use the Palenzuela formula:
      //
      // eps = - 1.0 + x(1-W^{2})/W + W( 1 + q - s + 0.5*( s/W^{2} + t^{2}/x^{2} ) )
      xeps = - 1.0 + (1.0-W*W)*x/W + W*( 1.0 + q - s + 0.5*( s/(W*W) + (t*t)/(x*x) ) );
      // Then compute P, S, and T using (rho,Ye,eps)
      enforce_table_bounds_rho_Ye_eps( eos,&xrho,&xye,&xeps );
      WVU_EOS_P_S_and_T_from_rho_Ye_eps( xrho,xye,xeps, &xprs,&xent,&xtemp );
    }

    prim[EPS] = xeps;

    AtStep++;
    AtP[AtStep]=xprs;

    if(AtStep>=2) {   //Aitken extrapolation

      AtR = (AtP[AtStep]-AtP[AtStep-1])/(AtP[AtStep-1]-AtP[AtStep-2]);
      if(AtR<1. && AtR>0.) {
        xprs=AtP[AtStep-1]+(AtP[AtStep]-AtP[AtStep-1])/(1.-AtR);
        AtStep=0;
        conacc = 1;
        AtP[0]=xprs;   //starting value for next Aitken extrapolation
      }
    }
  }
  while(fabs(xprs-P_old)>prec*(xprs+P_old) && step<maxsteps);

  if (conacc==1) {     //converged on an extrap. so recompute vars
    const CCTK_REAL a      = e + xprs + 0.5*B_squared;
    const CCTK_REAL phi    = acos(sqrt(27.0*d/(4.0*a))/a);
    const CCTK_REAL Eps    = a/3.0*( 1.0 - 2.0*cos( (2.0/3.0)*(phi + M_PI) ) );
    z = Eps-B_squared;
    const CCTK_REAL zBsq   = z + B_squared;
    const CCTK_REAL zBsqsq = zBsq*zBsq;
    const CCTK_REAL zsq    = z*z;
    const CCTK_REAL vsq    = (zsq * S_squared + (z+zBsq)*BdotSsq)/(zsq*zBsqsq);
    if(vsq<0. || vsq>=1.) {
    }

    const CCTK_REAL Wsq = 1.0/(1.0-vsq);
    const CCTK_REAL W   = fmin(fmax(sqrt(Wsq),1.0),eos.W_max);

    prim[WLORENTZ] = W;         //    W[s] = sqrt(cWsqr);
    prim[RHO     ] = con[DD]/W; //    rho[s] = tildeD[s]/(sqrtDetg[s]*W[s]);
    prim[PRESS   ] = xprs;

    prim[B1_con  ] = con[B1_con];
    prim[B2_con  ] = con[B2_con];
    prim[B3_con  ] = con[B3_con];
  }

  //Compute v^i
  const CCTK_REAL W = prim[WLORENTZ];
   
  prim[B1_con] = con[B1_con];
  prim[B2_con] = con[B2_con];
  prim[B3_con] = con[B3_con];

  prim[UTCON1] = W*(SU[0] + (BdotS)*con[B1_con]/z)/(z+B_squared);
  prim[UTCON2] = W*(SU[1] + (BdotS)*con[B2_con]/z)/(z+B_squared);
  prim[UTCON3] = W*(SU[2] + (BdotS)*con[B3_con]/z)/(z+B_squared);

  prim[RHO   ] = con[DD]/W;
  prim[YE    ] = xye;
  prim[ENT   ] = 0.0;

  if( use_entropy == true ) {
    // If using the entropy, compute S = (WS)/W
    prim[ENT] = con[WS]/W;
    // Then compute P, eps, and T using (rho,Ye,S)
    enforce_table_bounds_rho_Ye_S( eos,&prim[RHO],&prim[YE],&prim[ENT] );
    WVU_EOS_P_eps_and_T_from_rho_Ye_S( prim[RHO],prim[YE],prim[ENT],
                                       &prim[PRESS],&prim[EPS],&prim[TEMP] );
  }
  else {
    // If using the specific internal energy, remember that
    //
    // z = rho * h * W^{2} = D * h * W
    //
    // and therefore:
    //
    // x = h * W = z / D
    const CCTK_REAL x = z / con[DD];

    // Now use the Palenzuela formula:
    //
    // eps = - 1.0 + x(1-W^{2})/W + W( 1 + q - s + 0.5*( s/W^{2} + t^{2}/x^{2} ) )
    prim[EPS] = - 1.0 + (1.0-W*W)*x/W + W*( 1.0 + q - s + 0.5*( s/(W*W) + (t*t)/(x*x) ) );
    // Then compute P, S, and T using (rho,Ye,eps)
    prim[TEMP] = eos.T_atm;
    enforce_table_bounds_rho_Ye_eps( eos,&prim[RHO],&prim[YE],&prim[EPS] );
    WVU_EOS_P_S_and_T_from_rho_Ye_eps( prim[RHO],prim[YE],prim[EPS],
                                       &prim[PRESS],&prim[ENT],&prim[TEMP] );
  }

  if (step >= maxsteps) c2p_failed = true;

}
