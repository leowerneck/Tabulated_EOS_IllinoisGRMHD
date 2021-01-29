#include <stdio.h>
#include <stdlib.h>

#include "driver.h"
#include "nuc_eos.h"

int eoskey                 = 4;
int treat_floor_as_failure = 0;
int temp_guess             = 5.0;
double gam                 = 4.0/3.0;
double gam_m1_o_gam        = 0.25;
double Kpoly               = 100.0;

#define c2p_routine_3D   (0)
#define c2p_routine_2D   (1)
#define c2p_routine_SG   (2)
#define c2p_routine_BU   (3)
#define c2p_num_routines (4)

const char c2p_routine_names[c2p_num_routines][22] = {"3D routine",
                                                      "2D routine",
                                                      "2D safe guess routine",
                                                      "2D backup routine"};

void collect_metric_noble(struct of_geom *h, const double hxx, const double hxy, const double hxz,
        const double hyy, const double hyz, const double hzz){

  // get the metric and inverse metric components
  
  for (int i=0;i<4;i++) {
    h->gcov[0][i] = 0.0;
    h->gcov[i][0] = 0.0;
    h->gcon[0][i] = 0.0;
    h->gcon[i][0] = 0.0;
  }
  
  h->gcov[1][1] = hxx;
  h->gcov[1][2] = hxy;
  h->gcov[1][3] = hxz;
  h->gcov[2][2] = hyy;
  h->gcov[2][3] = hyz;
  h->gcov[3][3] = hzz;
  h->gcov[0][0] = -1.;
        
  h->gcov[2][1] = h->gcov[1][2];
  h->gcov[3][1] = h->gcov[1][3];
  h->gcov[3][2] = h->gcov[2][3];

  

  h->g = sqrt(hxx*hyy*hzz + 2.0*hxy*hyz*hxz - hxz*hxz*hyy - hxx*hyz*hyz - hxy*hxy*hzz);

  h-> g_inv = 1./h->g;


  double up_det = h-> g_inv;

  h->gcon[1][1] = (-hyz*hyz + hyy*hzz) * up_det;
  h->gcon[1][2] = (hxz*hyz - hxy*hzz) * up_det;
  h->gcon[2][2] = (-hxz*hxz + hxx*hzz) * up_det;
  h->gcon[1][3] = (-hxz*hyy + hxy*hyz) * up_det;
  h->gcon[2][3] = (hxy*hxz - hxx*hyz) * up_det;
  h->gcon[3][3] = (-hxy*hxy + hxx*hyy) *up_det;
  h->gcon[0][0] = -1.;
  
  h->gcon[2][1] = h->gcon[1][2];
  h->gcon[3][1] = h->gcon[1][3];
  h->gcon[3][2] = h->gcon[2][3];


  double lo_det = hxx*hyy*hzz + 2.0*hxy*hyz*hxz - hxz*hxz*hyy - hxx*hyz*hyz - hxy*hxy*hzz;
  
  h->alpha = 1.0 / lo_det;
  h->ncon[0] = 1./h->alpha;//sqrt(h->gcon[0][0]) ;
  //h->alpha = 1.0 / h->ncon[0];
  h->ncon[1] = h->alpha * h->gcon[0][1];
  h->ncon[2] = h->alpha * h->gcon[0][2];
  h->ncon[3] = h->alpha * h->gcon[0][3];
  h->beta[0] =  h->gcon[0][1]/h->gcon[0][0];
  h->beta[1] =  h->gcon[0][2]/h->gcon[0][0];
  h->beta[2] =  h->gcon[0][3]/h->gcon[0][0];
          
 //fprintf(stderr, "alpha%15.6E\n",h->alpha);
 //fprintf(stderr, "g_inv%15.6E\n",h-> g_inv);
}

inline double relative_error( const double a, const double b ) {

  if( a != 0.0 ) {
    return( fabs(1.0-b/a) );
  }
  else if( b != 0.0 ) {
    return( fabs(1.0-a/b) );
  }
  else {
    return 0.0;
  }
  
}

void prim2con_MHD(const double * prim, double * U, const double g_con[4][4],
                  const double g_cov[4][4], double alpha, double * bsq){
  

  
  const double lorentz = prim[WLORENTZ];
  
  U[RHO] = prim[RHO]*lorentz/alpha;

  U[BCON1] = prim[B1_con];
  U[BCON2] = prim[B2_con];
  U[BCON3] = prim[B3_con];
  
  const double uvec0_con=lorentz/alpha;
  const double uvec1_con=prim[u1_con]-alpha*lorentz*g_con[0][1];
  const double uvec2_con=prim[u2_con]-alpha*lorentz*g_con[0][2];
  const double uvec3_con=prim[u3_con]-alpha*lorentz*g_con[0][3];

  const double B1_cov = g_cov[1][1]*prim[B1_con]+g_cov[1][2]*prim[B2_con]+g_cov[1][3]*prim[B3_con];
  const double B2_cov = g_cov[2][1]*prim[B1_con]+g_cov[2][2]*prim[B2_con]+g_cov[2][3]*prim[B3_con];
  const double B3_cov = g_cov[3][1]*prim[B1_con]+g_cov[3][2]*prim[B2_con]+g_cov[3][3]*prim[B3_con];

  const double uvec0_cov=g_cov[0][0]*uvec0_con+g_cov[0][1]*uvec1_con+g_cov[0][2]*uvec2_con+g_cov[0][3]*uvec3_con;
  const double uvec1_cov=g_cov[1][0]*uvec0_con+g_cov[1][1]*uvec1_con+g_cov[1][2]*uvec2_con+g_cov[1][3]*uvec3_con;
  const double uvec2_cov=g_cov[2][0]*uvec0_con+g_cov[2][1]*uvec1_con+g_cov[2][2]*uvec2_con+g_cov[2][3]*uvec3_con;
  const double uvec3_cov=g_cov[3][0]*uvec0_con+g_cov[3][1]*uvec1_con+g_cov[3][2]*uvec2_con+g_cov[3][3]*uvec3_con;

  const double u0_cov=g_cov[0][1]*prim[u1_con]+g_cov[0][2]*prim[u2_con]+g_cov[0][3]*prim[u3_con];
  const double u1_cov=g_cov[1][1]*prim[u1_con]+g_cov[1][2]*prim[u2_con]+g_cov[1][3]*prim[u3_con];
  const double u2_cov=g_cov[2][1]*prim[u1_con]+g_cov[2][2]*prim[u2_con]+g_cov[2][3]*prim[u3_con];
  const double u3_cov=g_cov[3][1]*prim[u1_con]+g_cov[3][2]*prim[u2_con]+g_cov[3][3]*prim[u3_con];



  const double b0_con=1.0/alpha * (prim[B1_con]*u1_cov+prim[B2_con]*u2_cov+prim[B3_con]*u3_cov);

  //const double b0_cov=1./lorentz * alpha*u0_cov*b0_con;
  const double b1_cov=B1_cov/lorentz+1./lorentz * alpha*u1_cov*b0_con;
  const double b2_cov=B2_cov/lorentz+1./lorentz * alpha*u2_cov*b0_con;
  const double b3_cov=B3_cov/lorentz+1./lorentz * alpha*u3_cov*b0_con;

  const double b1_con=prim[B1_con]/lorentz+b0_con*(alpha*prim[u1_con]/lorentz-pow(alpha,2)*g_con[0][1]);
  const double b2_con=prim[B2_con]/lorentz+b0_con*(alpha*prim[u2_con]/lorentz-pow(alpha,2)*g_con[0][2]);
  const double b3_con=prim[B3_con]/lorentz+b0_con*(alpha*prim[u3_con]/lorentz-pow(alpha,2)*g_con[0][3]);

  const double Bsquared=prim[B1_con]*B1_cov + prim[B2_con]*B2_cov + prim[B3_con]*B3_cov;
  const double bsquared=(Bsquared+pow((alpha*b0_con),2))/pow(lorentz,2);
  
  const double w_enthalpy=prim[RHO]+prim[UU]+prim[PRESS];

  const double b0_cov=-b0_con;//for minkowski//1./b0_con * (b1_cov*b1_con+b2_cov*b2_con+b3_cov*b3_con);

  if (bsquared==0.){
  U[UU]=lorentz*(w_enthalpy)*uvec0_cov+alpha*(prim[PRESS])+prim[RHO]*lorentz/alpha;
  }

  else {
  U[UU]=lorentz*(w_enthalpy+bsquared)*uvec0_cov+alpha*(prim[PRESS]+bsquared/2.)+(-alpha*b0_con)*b0_cov+prim[RHO]*lorentz/alpha;
  }

  U[QCOV1]=lorentz*(w_enthalpy+bsquared)*uvec1_cov+(-alpha*b0_con)*b1_cov;
  U[QCOV2]=lorentz*(w_enthalpy+bsquared)*uvec2_cov+(-alpha*b0_con)*b2_cov;
  U[QCOV3]=lorentz*(w_enthalpy+bsquared)*uvec3_cov+(-alpha*b0_con)*b3_cov;

  U[YE] = prim[YE]*U[RHO];

  *bsq = bsquared;
}

int main(int argc,char** argv) {

  if( argc != 3 ) {
    fprintf(stderr,"(ERROR) Correct usage is ./standalone <eos_table> <which_routine_to_test>\n");
    fprintf(stderr,"(INFO) eos_table should be an HDF5 file (e.g. one from stellarcollapse.org)\n");
    for(int i=0;i<c2p_num_routines;i++ )
      fprintf(stderr,"(INFO) which_routine_to_test = %d -> %s\n",i,c2p_routine_names[i]);
    exit(1);
  }
  else {

    // int npoints  = pow(2,12);
    int numprims = 13; // rho, u, v^{x,y,z}, B^{x,y,z}, ye, T, eps, P, W
    int numcons  = 9;

    nuc_eos_C_ReadTable(argv[1]);
    printf("Successfully read EOS table from file %s\n",argv[1]);
    printf("rho_min = %e | rho_max = %e\n",eos_rhomin ,eos_rhomax );
    printf("Y_e_min = %e | Y_e_max = %e\n",eos_yemin  ,eos_yemax  );
    printf("  T_min = %e |   T_max = %e\n",eos_tempmin,eos_tempmax);

    const int c2p_routine = atoi(argv[2]);

    char prims_names[numprims][4];
    sprintf(prims_names[RHO     ],"RHO");
    sprintf(prims_names[YE      ],"Y_e");
    sprintf(prims_names[TEMP    ]," T ");
    sprintf(prims_names[PRESS   ]," P ");
    sprintf(prims_names[EPS     ],"eps");
    sprintf(prims_names[UU      ]," u ");
    sprintf(prims_names[u1_con  ],"v^x");
    sprintf(prims_names[u2_con  ],"v^y");
    sprintf(prims_names[u3_con  ],"v^z");
    sprintf(prims_names[B1_con  ],"B^x");
    sprintf(prims_names[B2_con  ],"B^y");
    sprintf(prims_names[B3_con  ],"B^z");
    sprintf(prims_names[WLORENTZ]," W ");

    // We'll keep W and Ye constant throughout
    const double W       = 2.0;
    const double ye      = 0.1;

    // We will try a constant guess for *all*
    // the temperatures we look at. This will
    // be useful in IGM because we do not keep
    // track of the temperature in between MHD
    // time steps.
    const double T_guess = 1e-2;

    // No need to compute this over and over during the test
    struct of_geom geom_noble;
    collect_metric_noble(&geom_noble, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0);

    //-----------------------------------------------------------
    //---------------------- con2prim test ----------------------
    //-----------------------------------------------------------

    const int npoints = pow(2,12);
    const int N       = (int)sqrt(npoints);

    // Now let us discretize rho and T evenly in log space
    const double lrmin = log(1.1*eos_rhomin);
    const double lrmax = log(0.5*eos_rhomax);
    const double dlr   = (lrmax-lrmin)/N;

    const double ltmin = log(1.1*eos_tempmin);
    const double ltmax = log(0.9*eos_tempmax);
    const double dlt   = (ltmax-ltmin)/N;

    int failures = 0;

    char filename[40];
    switch(c2p_routine) {
      case c2p_routine_3D:
        sprintf(filename,"out_c2p_test_3D.asc");
        break;
      case c2p_routine_2D:
        sprintf(filename,"out_c2p_test_2D.asc");
        break;
      case c2p_routine_SG:
        sprintf(filename,"out_c2p_test_SG.asc");
        break;
      case c2p_routine_BU:
        sprintf(filename,"out_c2p_test_BU.asc");
        break;
    }
    FILE* outfile = fopen(filename,"w");
    
    for(int i=0;i<N;i++) {
      for(int j=0;j<N;j++) {

        // (rho,Ye,T,P,eps)
        double xrho  = exp(lrmin + dlr*i);
        double xtemp = exp(ltmin + dlt*j);
        double xye   = ye;
        double xprs  = 0.0;
        double xeps  = 0.0;

        EOS_press_with_T(xrho,&xeps,&xtemp,xye,&xprs);
        
        // Velocity magnitude
        const double v = sqrt(1.0-1.0/(W*W));
        // Randomly choose vx, vy, vz
        const double vx = v*((double)rand())/((double)RAND_MAX);
        const double vy = sqrt(v*v - vx*vx)*((double)rand())/((double)RAND_MAX);
        const double vz = sqrt(v*v - vx*vx - vy*vy);

        // Now set the prim array
        double prims[numprims];
        prims[RHO     ] = xrho;
        prims[YE      ] = xye;
        prims[TEMP    ] = xtemp;
        prims[PRESS   ] = xprs;
        prims[EPS     ] = xeps;
        prims[UU      ] = xrho*xeps;
        prims[u1_con  ] = W*vx; // utilde^{x}
        prims[u2_con  ] = W*vy; // utilde^{y}
        prims[u3_con  ] = W*vz; // utilde^{z}
        prims[B1_con  ] = 0;
        prims[B2_con  ] = 0;
        prims[B3_con  ] = 0;   
        prims[WLORENTZ] = W;

        // Copy original prim values
        double orig_prims[numprims];
        for(int i=0;i<numprims;i++) orig_prims[i] = prims[i];

        // Now compute the conservs
        double cons[numcons], bsq;
        prim2con_MHD( prims, cons, geom_noble.gcon, geom_noble.gcov, geom_noble.alpha, &bsq );

        // These are IGM's guesses
        prims[RHO     ] = cons[RHO]/geom_noble.g; // here g = sqrt(gamma) = psi^{6}
        prims[YE      ] = cons[YE]/cons[RHO];
        prims[TEMP    ] = T_guess;

        // Now compute P and eps
        EOS_press_with_T(prims[RHO],&prims[EPS],&prims[TEMP],prims[YE],&prims[PRESS]);

        // Then compute u
        prims[UU      ] = prims[EPS]/prims[RHO];

        // Finally u^{i} and W
        prims[u1_con  ] = 0.0;
        prims[u2_con  ] = 0.0;
        prims[u3_con  ] = 0.0;
        prims[WLORENTZ] = 1.0;

        // con2prim
        int check;

        switch(c2p_routine) {
          case c2p_routine_3D:
            check = Utoprim_3d_eos(cons, &geom_noble, prims, &prims[WLORENTZ], &bsq);
            break;
          case c2p_routine_2D:
            check = Utoprim_2d_eos(cons, &geom_noble, prims, &prims[WLORENTZ], &bsq);
            break;
          case c2p_routine_SG:
            check = Utoprim_2d_eos_safe_guess(cons, &geom_noble, prims, &prims[WLORENTZ], &bsq);
            break;
          case c2p_routine_BU:
            check = Utoprim_2d_eos_backup(cons, &geom_noble, prims, &prims[WLORENTZ], &bsq);
            break;
        }

        // Check results
        double acc_err=0,acc_err_siegel=0;
        if( check != 0 ) {
          fprintf(stderr,"(ERROR) con2prim FAILURE! Error code: %d\n",check);
          failures++;
          acc_err=acc_err_siegel=1e300;
        }
        else {

          // The con2prim routine sets:
          //
          // rho, u, u^{x,y,z}, T
          //
          // So we must recompute P and eps
          EOS_press_with_T(prims[RHO],&prims[EPS],&prims[TEMP],prims[YE],&prims[PRESS]);

          // double gijuiuj = geom_noble.gcov[1][1]*prims[u1_con]*prims[u1_con] +
          //   geom_noble.gcov[2][2]*prims[u2_con]*prims[u2_con] +
          //   geom_noble.gcov[3][3]*prims[u3_con]*prims[u3_con];

          // prims[WLORENTZ] = gijuiuj/( 1.0 + sqrt(1.0+gijuiuj) ) + 1.0;

          printf("(INFO) con2prim SUCCESS!\n");
          double err[numprims];
          for(int i=0;i<numprims;i++) {
            err[i]   = relative_error(orig_prims[i],prims[i]);
            acc_err += err[i];
            printf("(INFO) con2prim error for %s   : %e\n",prims_names[i],err[i]);
          }
          acc_err_siegel = (err[RHO] + err[EPS] + err[u1_con] + err[u2_con] + err[u3_con])/5.0;
          acc_err        = MAX(1e-16,acc_err);
          acc_err_siegel = MAX(1e-16,acc_err_siegel);

          printf("(INFO) Siegel  accumulated error: %e\n",acc_err_siegel);
          printf("(INFO) Average accumulated error: %e\n",acc_err/numprims);
          printf("(INFO) Total   accumulated error: %e\n\n",acc_err);
        }
        fprintf(outfile,"%d    %e %e    %e %e %e\n",check,xrho,xtemp,acc_err_siegel,acc_err/numprims,acc_err);
      }
      fprintf(outfile,"\n");
    }

    fclose(outfile);

    printf("(INFO) Final report\n\n");
    printf("(INFO) Routine tested         : %s\n",c2p_routine_names[c2p_routine]);
    printf("(INFO) Number of points tested: %d\n",npoints);
    printf("(INFO) Number of failures     : %d\n",failures);
    printf("(INFO) Failure rate           : %.2lf%%\n\n",((double)failures/(double)npoints)*100.0);
      
    printf("(INFO) All done!\n");
    
    if( alltables  != NULL ) free(alltables );
    if( epstable   != NULL ) free(epstable  );
    if( presstable != NULL ) free(presstable);
    if( logrho     != NULL ) free(logrho    );
    if( logtemp    != NULL ) free(logtemp   );

    return 0;
  }

}
