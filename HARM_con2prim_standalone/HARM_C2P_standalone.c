/* .--------------------------------------------------------------.
 * | Tabulated EOS and neutrino leakage support for IllinoisGRMHD |
 * | Author: Leo Werneck (wernecklr@gmail.com)                    |
 * .--------------------------------------------------------------.
 *     This file is a modification of a HARM3D standalone code
 * written by Ariadna Murguia Berthier. It is intended to test
 * IllinoisGRMHD's conservative-to-primitive (C2P) routines and
 * demonstrate that it is compatible with HARM3D's routines.
 */

/*********************************************************
 * Defines and globals which are needed by either IGM's
 * version of this standalone of HARM3D's version of this
 * standalone
 *********************************************************/
#define cGH int
#define CCTK_REAL double
int update_Tmunu;
// For HARM3D
int numprims = 13;
int numcons  = 9;

#define GAMMAMAX (50.0)

#define MASS_BH (1.0)

#define SQR(x) ( (x)*(x) )

#define DEBUG_RUN (0) // 1 disables OpenMP

double temp_guess=5.;
double gam=4./3.;
double gam_m1_o_gam=0.25;//(gam-1.)/gam;!em
int treat_floor_as_failure = 0;
double ye=0.5;
double Kpoly=1.e2;
int eoskey=2;
int myid = 0;
int n_substep = 0;
int nstep = 0;

/*********************************************************/
/*********************************************************/

#include "driver.h"
// #include "nuc_eos.h"

/*********************************************************
 * IGM stuff
 *********************************************************/
// The order here MATTERS, as we assume that GUPXX+1=GUPYY, etc.
static const int PHI=0,PSI=1,GXX=2,GXY=3,GXZ=4,GYY=5,GYZ=6,GZZ=7,
  LAPM1=8,SHIFTX=9,SHIFTY=10,SHIFTZ=11,GUPXX=12,GUPYY=13,GUPZZ=14,
  NUMVARS_FOR_METRIC_FACEVALS=15; //<-- Be _sure_ to set this correctly, or you'll have memory access bugs!

// These are not used for facevals in the reconstruction step, but boy are they useful anyway.
static const int GUPXY=15,GUPXZ=16,GUPYZ=17,
  NUMVARS_FOR_METRIC=18; //<-- Be _sure_ to set this correctly, or you'll have memory access bugs!

// The order here MATTERS, and must be consistent with the order in the in_prims[] array in driver_evaluate_MHD_rhs.C.
static const int RHOB=0,PRESSURE=1,VX=2,VY=3,VZ=4,
  BX_CENTER=5,BY_CENTER=6,BZ_CENTER=7,BX_STAGGER=8,BY_STAGGER=9,BZ_STAGGER=10,
  VXR=11,VYR=12,VZR=13,VXL=14,VYL=15,VZL=16,YEPRIM=17,EPSILON=18,TEMPERATURE=19,MAXNUMVARS=20;  //<-- Be _sure_ to define MAXNUMVARS appropriately!

// Again, the order here MATTERS, since we assume in the code that, CONSERV[STILDEX+1] = \tilde{S}_y
static const int RHOSTAR=0,STILDEX=1,STILDEY=2,STILDEZ=3,TAUENERGY=4,YESTAR=5,NUM_CONSERVS=6;

// Auxiliary vars for the metric
static const int LAPSE=0,PSI2=1,PSI4=2,PSI6=3,PSIM4=4,LAPSEINV=5,NUMVARS_METRIC_AUX=6;
#define SET_LAPSE_PSI4(array_name,METRIC)   {                   \
      array_name[LAPSE] = METRIC[LAPM1]+1.0;                    \
      array_name[PSI2]  = exp(2.0*METRIC[PHI]);                 \
      array_name[PSI4]  = SQR(array_name[PSI2]);                \
      array_name[PSI6]  = array_name[PSI4]*array_name[PSI2];    \
      array_name[PSIM4]  = 1.0/array_name[PSI4];                \
      array_name[LAPSEINV]  = 1.0/array_name[LAPSE];            \
  }

// These functions are implemented
// right after the main function below
void set_ADM_3metric( const CCTK_REAL perturbation_strength, CCTK_REAL *ADM_3METRIC );
void set_velocities( const CCTK_REAL gamma, CCTK_REAL *PRIMS );
void set_magnetic_fields( const CCTK_REAL log_P_o_Pmag, CCTK_REAL *PRIMS );
void IllinoisGRMHD_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij( CCTK_REAL *ADM_METRIC, CCTK_REAL *BSSN_METRIC );
void compute_ADM_4metric_and_inverse( const CCTK_REAL *ADM_METRIC, const CCTK_REAL *METRIC_LAP_PSI4, CCTK_REAL g4dn[4][4], CCTK_REAL g4up[4][4] );
void set_harm_prims_from_IGM_prims( const CCTK_REAL *IGM_prims,
                                    const CCTK_REAL T_local,
                                    const CCTK_REAL Ye_local,
                                    const CCTK_REAL eps_local,
                                    const CCTK_REAL gamma_local,
                                    CCTK_REAL *harm_prims);
void set_harm_metric_from_IGM_metric( const CCTK_REAL *METRIC,
                                      const CCTK_REAL *METRIC_LAP_PSI4,
                                      const CCTK_REAL g4dn[4][4],
                                      const CCTK_REAL g4up[4][4],
                                      struct of_geom *geom );
void set_harm_cons_and_bsq_from_harm_prims( const struct of_geom geom,
                                            const CCTK_REAL *prims,
                                            CCTK_REAL *cons,
                                            CCTK_REAL *bsq );

void set_IGM_primitives_from_HARM3D_primitives( const CCTK_REAL *METRIC_PHYS,
                                                const CCTK_REAL *METRIC_LAP_PSI4,
                                                const CCTK_REAL *harm_cons,
                                                const CCTK_REAL *harm_prims,
                                                CCTK_REAL *igm_prims,
                                                CCTK_REAL *igm_eps );
/*********************************************************/
/*********************************************************/

// Extra defines
#define maybe_unused __attribute__((unused))
#define MeV_to_K (1.0/8.617330337217213e-11) // From Ari's script
#define random_number_between_m1_and_p1 ( ((CCTK_REAL)rand())/((CCTK_REAL)RAND_MAX) )

inline CCTK_REAL rel_err( const CCTK_REAL a, const CCTK_REAL b ) {

  // If both a and b are zero, return absolute error instead
  if( a != 0 )     return( fabs((a-b)/a) );
  else if( b!= 0 ) return( fabs((a-b)/b) );
  else             return 0.0;
  
}

// .----------------------------.
// | BEGINNING OF MAIN FUNCTION |
// .----------------------------.
int main( int argc, char *argv[] ) {

  // Start by checking for correct usage
  if( argc != 2 ) {
    fprintf(stderr,"(ERROR) Correct usage is ./HARM_C2P_standalone [EOS table path]\n");
    fflush(stderr);
    exit(1);
  }

  // Now read in the EOS table. Note that by doing this we
  // need to free the allocated memory in the end of this
  // function
  nuc_eos_C_ReadTable( argv[1] );

  fprintf(stdout,"(INFO) Basic info for table \"%s\" ...\n",argv[1]);
  fprintf(stdout,"(INFO) rho_min = %e | rho_max = %e\n",exp(logrho[0]) ,exp(logrho[nrho-1]) );
  fprintf(stdout,"(INFO) T_min   = %e | T_max   = %e\n",eos_tempmin,eos_tempmax);
  fprintf(stdout,"(INFO) Ye_min  = %e | Ye_max  = %e\n",eos_yemin  ,eos_yemax  );
  fflush(stdout);

  // We know before hand that for the particular stellarcollapse.opg table
  // "SLy4_NSE_rho391_temp163_ye66.h5"
  // we have
  // .-----------------.-------.
  // |    Quantity     | Value |
  // .-----------------.-------.
  // | log10(T_min[K]) |  ~ 7  |
  // | log10(T_max[K]) | ~12.5 |
  // | log10(rho_min ) |  ~ 3  |
  // | log10(rho_max ) |  ~16  |
  // .-----------------.-------.
  // and so we can sample rho and T such that we are always within table
  // bounds. Here we approximate what was done by Ari.
  const CCTK_REAL logrho_sampling_min =  5.0 + log10(RHOGF);
  const CCTK_REAL logrho_sampling_max = 14.0 + log10(RHOGF);
  const CCTK_REAL logT_sampling_min   =  8.0 - log10(MeV_to_K);
  const CCTK_REAL logT_sampling_max   = 10.5 - log10(MeV_to_K);

  // Define the sampling points
  const int nsamples = 100;
  const int ntotal   = nsamples+1;

  // Define the "step sizes"
  const CCTK_REAL dlogrho_sampling = (logrho_sampling_max-logrho_sampling_min)/nsamples;
  const CCTK_REAL dlogT_sampling   = (  logT_sampling_max-logT_sampling_min  )/nsamples;

  // Define and populate the sampling arrays for log(rho) and log(T)
  CCTK_REAL logrho_sampling[ntotal], logT_sampling[ntotal];
  for(int i=0;i<ntotal;i++) {
    logrho_sampling[i] = logrho_sampling_min + i*dlogrho_sampling;
    logT_sampling[i]   = logT_sampling_min   + i*dlogT_sampling;
  }

  // Now that we have the samples, define basic quantities that will be used for all points
  // const CCTK_REAL Ye_local     =  0.1;
  const CCTK_REAL gamma_local  =  2.0;
  const CCTK_REAL log_Pmag_o_P = -5.0;

  // Set the random seed and perturbation strength, which will be
  // used to perturb metric quantities about Minkowski spacetime
  const CCTK_REAL perturbation_strength = 0.0;
  const unsigned int random_seed        = 0;
  srand(random_seed);

  // Set output file
  FILE *outfile = fopen("c2p_errors.asc","w");
  fprintf(outfile,"# IlliinoisGRMHD - C2P test\n");
  fprintf(outfile,"# Columns:\n");
  fprintf(outfile,"# 1. Density in g/cm^3\n");
  fprintf(outfile,"# 2. Temperature in MeV\n");
  fprintf(outfile,"# 3. Sum of relative errors (P,vx,vy,vz,T) from Utoprim_3d_eos()\n");
  fprintf(outfile,"# 5. Sum of relative errors (P,vx,vy,vz,T) from Utoprim_2d_eos()\n");
  fprintf(outfile,"# 6. Sum of relative errors (P,vx,vy,vz,T) from Utoprim_2d_eos_safe_guess()\n");
  fprintf(outfile,"# 7. Sum of relative errors (P,vx,vy,vz,T) from Utoprim_2d_eos_backup()\n");

  // Now we are ready to start. We loop over the number of sampling points,
  // computing the pressure P from rho, Ye, and T, and setting values for
  // v^{i}, B^{i}, and metric quantities randomly.
  int fails_3d=0;
  int fails_2d=0;
  int fails_sg=0;
  int fails_bu=0;
  FILE *out3d = fopen("out3d.asc","w");
  FILE *out2d = fopen("out2d.asc","w");
  FILE *outsg = fopen("outsg.asc","w");
  FILE *outbu = fopen("outbu.asc","w");
  //#pragma omp parallel for schedule(static)
  for(int jj=0;jj<ntotal;jj++) {
    for(int ii=0;ii<ntotal;ii++) {
      // Get the local values of rho and T
      CCTK_REAL rho_local = pow(10.0,logrho_sampling[ii]);
      CCTK_REAL T_local   = pow(10.0,logT_sampling[jj]);
      CCTK_REAL Ye_local  = 1e-1;

      // Set the pressure. The specific internal energy will
      // be used below to compute u.
      CCTK_REAL P_local;
      CCTK_REAL eps_local;
      // This function call gets P and eps from (rho,Ye,T).
      EOS_press_with_T(rho_local, &eps_local, &T_local, Ye_local, &P_local);

      /********************************************/
      /********************************************/
      // .-----------------------.
      // | Set up IGM quantities |
      // .-----------------------.
      // Then initialize the PRIMS array following
      // IllinoisGRMHD's standards
      CCTK_REAL PRIMS[MAXNUMVARS];
      
      // Set P and rho
      PRIMS[RHOB]     = rho_local;
      PRIMS[PRESSURE] = P_local;

      // Set the velocities
      // set_velocities( gamma_local, PRIMS );
      PRIMS[VX] = PRIMS[VY] = PRIMS[VZ] = 0.0;

      // Set magnetic fields
      set_magnetic_fields( log_Pmag_o_P, PRIMS );
      // PRIMS[BX_CENTER] = PRIMS[BY_CENTER] = PRIMS[BZ_CENTER] = 0.0;

      // Now set the original prims
      CCTK_REAL prims_orig[numprims];
      prims_orig[UU      ] = rho_local * eps_local;
      prims_orig[YE      ] = Ye_local;
      prims_orig[RHO     ] = rho_local;
      prims_orig[TEMP    ] = T_local;
      prims_orig[PRESS   ] = P_local;
      prims_orig[BCON1   ] = PRIMS[BX_CENTER];
      prims_orig[BCON2   ] = PRIMS[BY_CENTER];
      prims_orig[BCON3   ] = PRIMS[BZ_CENTER];
      prims_orig[UTCON1  ] = PRIMS[VX];
      prims_orig[UTCON2  ] = PRIMS[VY];
      prims_orig[UTCON3  ] = PRIMS[VZ];
      prims_orig[WLORENTZ] = 1.0;//gamma_local;

      // Set the physical metric gamma_{ij}, the lapse, and the shift
      CCTK_REAL METRIC_PHYS[NUMVARS_FOR_METRIC];
      set_ADM_3metric( perturbation_strength, METRIC_PHYS );
      
      // Set BSSN quantities. Also recompute gamma_{ij} and compute gamma^{ij}.
      CCTK_REAL METRIC[NUMVARS_FOR_METRIC];
      IllinoisGRMHD_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij( METRIC_PHYS, METRIC );

      // Set auxiliary functions associated with the metric
      CCTK_REAL METRIC_LAP_PSI4[NUMVARS_METRIC_AUX];
      SET_LAPSE_PSI4(METRIC_LAP_PSI4,METRIC);

      // Compute ADM 4-metric, its inverse, and its determinant.
      CCTK_REAL g4dn[4][4],g4up[4][4];
      compute_ADM_4metric_and_inverse( METRIC_PHYS, METRIC_LAP_PSI4, g4dn, g4up );
      
      /********************************************/
      /********************************************/
      // .--------------------------.
      // | Set up HARM3D quantities |
      // .--------------------------.
      // Translate IGM's prims to HARM3D
      CCTK_REAL prims_3d[numprims], prims_2d[numprims], prims_sg[numprims], prims_bu[numprims];
      for(int i=0;i<numprims;i++) prims_3d[i] = prims_2d[i] = prims_sg[i] = prims_bu[i] = prims_orig[i];
      // set_harm_prims_from_IGM_prims( PRIMS,T_local,Ye_local,eps_local,gamma_local, prims_3d );

      /* // Now set HARM's struct of geom */
      struct of_geom geom;
      set_harm_metric_from_IGM_metric( METRIC,METRIC_LAP_PSI4,g4dn,g4up, &geom );
      
      // Set HARM3D conservs from HARM3D prims
      // This is based on Ari's script
      CCTK_REAL cons[numcons],bsq;
      set_harm_cons_and_bsq_from_harm_prims( geom, prims_3d, cons, &bsq );

      // Make copies of bsq
      CCTK_REAL bsq_3d, bsq_2d, bsq_sg, bsq_bu;
      bsq_3d = bsq_2d = bsq_sg = bsq_bu = bsq;

      // And make copies of gamma
      CCTK_REAL gamma_3d,gamma_2d,gamma_sg,gamma_bu;
      gamma_3d = gamma_2d = gamma_sg = gamma_bu = 1.0;//gamma_local;

      // printf("(INFO) Struct of geom information.\n(INFO) gcon:\n");
      // for(int j=0;j<4;j++) {
      //   for(int i=0;i<4;i++) {
      //     printf(" %23.16e",geom.gcon[i][j]);
      //   }
      //   printf("\n");
      // }
      // printf("(INFO) gcov:\n");
      // for(int j=0;j<4;j++) {
      //   for(int i=0;i<4;i++) {
      //     printf(" %23.16e",geom.gcov[i][j]);
      //   }
      //   printf("\n");
      // }

      // printf("(INFO   ) Primitive guesses (rho Ye T U vx vy vz Wlorentz):\n");
      // printf("(INFO 3D) %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n",
      //        prims_3d[RHO],prims_3d[YE],prims_3d[TEMP],prims_3d[UU],prims_3d[UTCON1],prims_3d[UTCON2],prims_3d[UTCON3],prims_3d[WLORENTZ]);
      // printf("(INFO 2D) %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n",
      //        prims_2d[RHO],prims_2d[YE],prims_2d[TEMP],prims_2d[UU],prims_2d[UTCON1],prims_2d[UTCON2],prims_2d[UTCON3],prims_2d[WLORENTZ]);
      // printf("(INFO SG) %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n",
      //        prims_sg[RHO],prims_sg[YE],prims_sg[TEMP],prims_sg[UU],prims_sg[UTCON1],prims_sg[UTCON2],prims_sg[UTCON3],prims_sg[WLORENTZ]);
      // printf("(INFO BU) %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n",
      //        prims_bu[RHO],prims_bu[YE],prims_bu[TEMP],prims_bu[UU],prims_bu[UTCON1],prims_bu[UTCON2],prims_bu[UTCON3],prims_bu[WLORENTZ]);
      // printf("(CONS   ) %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n",
      //        cons[RHO],cons[YE],0.0,cons[UU],cons[UTCON1],cons[UTCON2],cons[UTCON3]);

      // Call the different C2P routines
      int check_3d = Utoprim_3d_eos(            cons, &geom, prims_3d, &gamma_3d, &bsq_3d );
      int check_2d = Utoprim_2d_eos(            cons, &geom, prims_2d, &gamma_2d, &bsq_2d );
      int check_sg = Utoprim_2d_eos_safe_guess( cons, &geom, prims_sg, &gamma_sg, &bsq_sg );
      int check_bu = Utoprim_2d_eos_backup(     cons, &geom, prims_bu, &gamma_bu, &bsq_bu );

      // printf("(INFO) Function return values. 3d,2d,sg,bu: %d,%d,%d,%d\n",check_3d,check_2d,check_sg,check_bu);
      // printf("(INFO) 2d: %d\n",check_2d);

      // printf("(INFO) Printing results for those routines which didn't fail:\n");
      // printf("(INFO) Legend: OR=original, 3D=utoprim_3d_eos, 2D=utoprim_2d_eos, SG=utoprim_2d_eos_safe_guess, BU=utoprim_2d_eos_backup\n");
      // printf("(INFO) Original prims (rho Ye T U vx vy vz Wlorentz):\n");
      // printf("(INFO OR) %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n",
      //        prims_orig[RHO],prims_orig[YE],prims_orig[TEMP],prims_orig[UU],prims_orig[UTCON1],prims_orig[UTCON2],prims_orig[UTCON3],prims_orig[WLORENTZ]);
      // if( check_3d == 0 ) {
      //   printf("(INFO 3D) %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n",
      //        prims_3d[RHO],prims_3d[YE],prims_3d[TEMP],prims_3d[UU],prims_3d[UTCON1],prims_3d[UTCON2],prims_3d[UTCON3],prims_3d[WLORENTZ]);
      // }
      // if( check_2d == 0 ) {
      //   printf("(INFO 2D) %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n",
      //        prims_2d[RHO],prims_2d[YE],prims_2d[TEMP],prims_2d[UU],prims_2d[UTCON1],prims_2d[UTCON2],prims_2d[UTCON3],prims_2d[WLORENTZ]);
      // }
      // if( check_sg == 0 ) {
      //   printf("(INFO SG) %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n",
      //        prims_sg[RHO],prims_sg[YE],prims_sg[TEMP],prims_sg[UU],prims_sg[UTCON1],prims_sg[UTCON2],prims_sg[UTCON3],prims_sg[WLORENTZ]);
      // }
      // if( check_bu == 0 ) {
      //   printf("(INFO BU) %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n",
      //        prims_bu[RHO],prims_bu[YE],prims_bu[TEMP],prims_bu[UU],prims_bu[UTCON1],prims_bu[UTCON2],prims_bu[UTCON3],prims_bu[WLORENTZ]);
      // }

      // return 0;

      // Now compute the errors
      CCTK_REAL rho_rel_err,eps_rel_err,vx_rel_err,vy_rel_err,vz_rel_err;
      // Errors from the 3d routine
      if( check_3d==0 ) {
        CCTK_REAL eps_3d = prims_3d[UU]/prims_3d[RHO];
        rho_rel_err = rel_err( prims_3d[RHO   ],prims_orig[RHO   ]);
        eps_rel_err = rel_err( eps_3d          ,prims_orig[UU    ]/prims_orig[RHO]  );
        vx_rel_err  = 0.0;//rel_err( prims_3d[UTCON1],prims_orig[UTCON1]);
        vy_rel_err  = 0.0;//rel_err( prims_3d[UTCON2],prims_orig[UTCON2]);
        vz_rel_err  = 0.0;//rel_err( prims_3d[UTCON3],prims_orig[UTCON3]);
        /* printf("Orig  3d: rho=%e eps=%e vx=%e vy=%e vz=%e\n",rho_b_orig,eps_orig,vx_orig,vy_orig,vz_orig); */
        /* printf("Outp  3d: rho=%e eps=%e vx=%e vy=%e vz=%e\n",prims_3d[RHO],eps_3d,prims_3d[UTCON1],prims_3d[UTCON2],prims_3d[UTCON3]); */
        /* printf("Error 3d: rho=%e eps=%e vx=%e vy=%e vz=%e\n",rho_rel_err,eps_rel_err,vx_rel_err,vy_rel_err,vz_rel_err); */
      }
      else {
        rho_rel_err = 1e100;
        eps_rel_err = 1e100;
        vx_rel_err  = 1e100;
        vy_rel_err  = 1e100;
        vz_rel_err  = 1e100;
        fails_3d++;
      }
      // fprintf(out3d,"(Original) %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",rho_b_orig,eps_orig,vx_orig,vy_orig,vz_orig,T_local,Ye_local);
      // fprintf(out3d,"( Result ) %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n\n",prims_3d[RHO],prims_3d[UU]/prims_3d[RHO],prims_3d[UTCON1],prims_3d[UTCON2],prims_3d[UTCON3],prims_3d[TEMP],prims_3d[YE]);
      const CCTK_REAL sum_3d_errors = 0.2*(rho_rel_err + eps_rel_err + vx_rel_err + vy_rel_err + vz_rel_err);

      // Errors from the 2d routine
      if( check_2d==0 ) {
        CCTK_REAL eps_2d;
        set_IGM_primitives_from_HARM3D_primitives( METRIC_PHYS,METRIC_LAP_PSI4,cons,prims_2d,PRIMS,&eps_2d );
        rho_rel_err = rel_err( prims_2d[RHO   ],prims_orig[RHO   ]);
        eps_rel_err = rel_err( eps_2d          ,prims_orig[UU    ]/prims_orig[RHO]  );
        vx_rel_err  = 0.0;//rel_err( prims_2d[UTCON1],prims_orig[UTCON1]);
        vy_rel_err  = 0.0;//rel_err( prims_2d[UTCON2],prims_orig[UTCON2]);
        vz_rel_err  = 0.0;//rel_err( prims_2d[UTCON3],prims_orig[UTCON3]);
        /* printf("Orig  2d: rho=%e eps=%e vx=%e vy=%e vz=%e\n",rho_b_orig,eps_orig,vx_orig,vy_orig,vz_orig); */
        /* printf("Outp  2d: rho=%e eps=%e vx=%e vy=%e vz=%e\n",prims_2d[RHO],eps_2d,prims_2d[UTCON1],prims_2d[UTCON2],prims_2d[UTCON3]); */
        /* printf("Error 2d: rho=%e eps=%e vx=%e vy=%e vz=%e\n",rho_rel_err,eps_rel_err,vx_rel_err,vy_rel_err,vz_rel_err); */
      }
      else {
        rho_rel_err = 1e100;
        eps_rel_err = 1e100;
        vx_rel_err  = 1e100;
        vy_rel_err  = 1e100;
        vz_rel_err  = 1e100;
        fails_2d++;
      }
      // fprintf(out2d,"(Original) %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",rho_b_orig,eps_orig,vx_orig,vy_orig,vz_orig,T_local,Ye_local);
      // fprintf(out2d,"( Result ) %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n\n",prims_2d[RHO],prims_2d[UU]/prims_2d[RHO],prims_2d[UTCON1],prims_2d[UTCON2],prims_2d[UTCON3],prims_2d[TEMP],prims_2d[YE]);
      const CCTK_REAL sum_2d_errors = 0.2*(rho_rel_err + eps_rel_err + vx_rel_err + vy_rel_err + vz_rel_err);

      // Errors from the safe guess routine
      if( check_sg==0 ) {
        CCTK_REAL eps_sg;
        set_IGM_primitives_from_HARM3D_primitives( METRIC_PHYS,METRIC_LAP_PSI4,cons,prims_sg,PRIMS,&eps_sg );
        rho_rel_err = rel_err( prims_sg[RHO   ],prims_orig[RHO   ]);
        eps_rel_err = rel_err( eps_sg          ,prims_orig[UU    ]/prims_orig[RHO]  );
        vx_rel_err  = 0.0;//rel_err( prims_sg[UTCON1],prims_orig[UTCON1]);
        vy_rel_err  = 0.0;//rel_err( prims_sg[UTCON2],prims_orig[UTCON2]);
        vz_rel_err  = 0.0;//rel_err( prims_sg[UTCON3],prims_orig[UTCON3]);
        /* printf("Orig  sg: rho=%e eps=%e vx=%e vy=%e vz=%e\n",rho_b_orig,eps_orig,vx_orig,vy_orig,vz_orig); */
        /* printf("Outp  sg: rho=%e eps=%e vx=%e vy=%e vz=%e\n",prims_sg[RHO],eps_sg,prims_sg[UTCON1],prims_sg[UTCON2],prims_sg[UTCON3]); */
        /* printf("Error sg: rho=%e eps=%e vx=%e vy=%e vz=%e\n",rho_rel_err,eps_rel_err,vx_rel_err,vy_rel_err,vz_rel_err); */
      }
      else {
        rho_rel_err = 1e100;
        eps_rel_err = 1e100;
        vx_rel_err  = 1e100;
        vy_rel_err  = 1e100;
        vz_rel_err  = 1e100;
        fails_sg++;
      }
      // fprintf(outsg,"(Original) %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",rho_b_orig,eps_orig,vx_orig,vy_orig,vz_orig,T_local,Ye_local);
      // fprintf(outsg,"( Result ) %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n\n",prims_sg[RHO],prims_sg[UU]/prims_sg[RHO],prims_sg[UTCON1],prims_sg[UTCON2],prims_sg[UTCON3],prims_sg[TEMP],prims_sg[YE]);
      const CCTK_REAL sum_sg_errors = 0.2*(rho_rel_err + eps_rel_err + vx_rel_err + vy_rel_err + vz_rel_err);

      // Errors from the backup routine
      if( check_bu==0 ) {
        CCTK_REAL eps_bu;
        set_IGM_primitives_from_HARM3D_primitives( METRIC_PHYS,METRIC_LAP_PSI4,cons,prims_bu,PRIMS,&eps_bu );
        rho_rel_err = rel_err( prims_bu[RHO   ],prims_orig[RHO   ]);
        eps_rel_err = rel_err( eps_bu          ,prims_orig[UU    ]/prims_orig[RHO]  );
        vx_rel_err  = 0.0;//rel_err( prims_bu[UTCON1],prims_orig[UTCON1]);
        vy_rel_err  = 0.0;//rel_err( prims_bu[UTCON2],prims_orig[UTCON2]);
        vz_rel_err  = 0.0;//rel_err( prims_bu[UTCON3],prims_orig[UTCON3]);
        /* printf("Orig  bu: rho=%e eps=%e vx=%e vy=%e vz=%e\n",rho_b_orig,eps_orig,vx_orig,vy_orig,vz_orig); */
        /* printf("Outp  bu: rho=%e eps=%e vx=%e vy=%e vz=%e\n",prims_bu[RHO],eps_bu,prims_bu[UTCON1],prims_bu[UTCON2],prims_bu[UTCON3]); */
        /* printf("Error bu: rho=%e eps=%e vx=%e vy=%e vz=%e\n",rho_rel_err,eps_rel_err,vx_rel_err,vy_rel_err,vz_rel_err); */
      }
      else {
        rho_rel_err = 1e100;
        eps_rel_err = 1e100;
        vx_rel_err  = 1e100;
        vy_rel_err  = 1e100;
        vz_rel_err  = 1e100;
        fails_bu++;
      }
      // fprintf(outbu,"(Original) %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",rho_b_orig,eps_orig,vx_orig,vy_orig,vz_orig,T_local,Ye_local);
      // fprintf(outbu,"( Result ) %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n\n",prims_bu[RHO],prims_bu[UU]/prims_bu[RHO],prims_bu[UTCON1],prims_bu[UTCON2],prims_bu[UTCON3],prims_bu[TEMP],prims_bu[YE]);
      const CCTK_REAL sum_bu_errors = 0.2*(rho_rel_err + eps_rel_err + vx_rel_err + vy_rel_err + vz_rel_err);

      // If all routines fail, then the product of the return codes
      // we be different from zero. In that case, notify the user
      // that a major failure has occured. Such a point would have
      // to be handled by fixes or backups in an actual run.
      if( check_3d*check_2d*check_sg*check_bu != 0 ) {
        fprintf(stderr,"(ERROR) Major C2P failure! 3d,2d,sg,bu: %d,%d,%d,%d\n",check_3d,check_2d,check_sg,check_bu);
      }

      fprintf(outfile,"%.15e %.15e %.15e %.15e %.15e %.15e\n",rho_local*INVRHOGF,T_local,sum_3d_errors,sum_2d_errors,sum_sg_errors,sum_bu_errors);

    }
    fprintf(outfile,"\n");
  }

  // printf("Finished, closing files\n");

  fclose(outfile);
  // fclose(out3d);
  // fclose(out2d);
  // fclose(outsg);
  // fclose(outbu);

  const int totalpts = ntotal*ntotal;
  printf("Failures: 3d: %d/%d(%.0lf%%) | 2d: %d/%d(%.0lf%%) | sg: %d/%d(%.0lf%%) | bu: %d/%d(%.0lf%%)\n",
         fails_3d,totalpts,100.0*( ((CCTK_REAL)fails_3d)/((CCTK_REAL)totalpts) ),
         fails_2d,totalpts,100.0*( ((CCTK_REAL)fails_2d)/((CCTK_REAL)totalpts) ),
         fails_sg,totalpts,100.0*( ((CCTK_REAL)fails_sg)/((CCTK_REAL)totalpts) ),
         fails_bu,totalpts,100.0*( ((CCTK_REAL)fails_bu)/((CCTK_REAL)totalpts) ));

  return 0;
}
// .----------------------.
// | END OF MAIN FUNCTION |
// .----------------------.

// This function sets IGM primitives based on HARM3D primitives
void set_IGM_primitives_from_HARM3D_primitives( const CCTK_REAL *METRIC_PHYS,
                                                const CCTK_REAL *METRIC_LAP_PSI4,
                                                const CCTK_REAL *harm_cons,
                                                const CCTK_REAL *harm_prims,
                                                CCTK_REAL *igm_prims,
                                                CCTK_REAL *igm_eps ) {

  // Now that we have found some solution, we first limit velocity:
  // FIXME: Probably want to use exactly the same velocity limiter function here as in mhdflux.C
  CCTK_REAL utx_new = harm_prims[UTCON1];
  CCTK_REAL uty_new = harm_prims[UTCON2];
  CCTK_REAL utz_new = harm_prims[UTCON3];

  // Velocity limiter:
  CCTK_REAL gijuiuj = METRIC_PHYS[GXX]*SQR(utx_new) + 2.0*METRIC_PHYS[GXY]*utx_new*uty_new
                    + METRIC_PHYS[GYY]*SQR(uty_new) + 2.0*METRIC_PHYS[GXZ]*utx_new*utz_new
                    + METRIC_PHYS[GZZ]*SQR(utz_new) + 2.0*METRIC_PHYS[GYZ]*uty_new*utz_new;
  CCTK_REAL au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
  CCTK_REAL u0L = (au0m1+1.0)*METRIC_LAP_PSI4[LAPSEINV];

  CCTK_REAL rho = harm_prims[RHO];

  // Limit velocity
  if (au0m1 > 0.9999999*(GAMMAMAX-1.0)) {
    CCTK_REAL fac = sqrt((SQR(GAMMAMAX)-1.0)/(SQR(1.0+au0m1) - 1.0));
    utx_new *= fac;
    uty_new *= fac;
    utz_new *= fac;
    gijuiuj = gijuiuj * SQR(fac);
    au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
    // Reset rho_b and u0
    u0L = (au0m1+1.0)*METRIC_LAP_PSI4[LAPSEINV];
    rho =  harm_cons[RHO]/(METRIC_LAP_PSI4[LAPSE]*u0L*METRIC_LAP_PSI4[PSI6]);
  } //Finished limiting velocity

  // Now set prims. rho
  igm_prims[RHOB] = rho;
  // Pressure
  CCTK_REAL harm_Ye = harm_prims[YE];
  CCTK_REAL harm_T  = harm_prims[TEMP];
  CCTK_REAL eps_local;
  EOS_press_with_T(rho, &eps_local, &harm_T, harm_Ye, &igm_prims[PRESSURE]);
  // eps
  *igm_eps = eps_local;
  // Velocities
  igm_prims[VX] = utx_new/u0L - METRIC_PHYS[SHIFTX];
  igm_prims[VY] = uty_new/u0L - METRIC_PHYS[SHIFTY];
  igm_prims[VZ] = utz_new/u0L - METRIC_PHYS[SHIFTZ];
  
}

// Set HARM3D conservs and b^2
void set_harm_cons_and_bsq_from_harm_prims( const struct of_geom geom,
                                            const CCTK_REAL *prims,
                                            CCTK_REAL *cons,
                                            CCTK_REAL *bsq ) {

  // Lorentz factor
  CCTK_REAL lorentz  = prims[WLORENTZ];
  const double alpha = geom.alpha;

  // D = rho * W
  cons[RHO] = prims[RHO] * lorentz / alpha;

  // printf("cons[RHO] = %.15e ( %e %e %e )\n",cons[RHO],prims[RHO],lorentz,alpha);

  // B^i
  cons[BCON1] = prims[B1_con];
  cons[BCON2] = prims[B2_con];
  cons[BCON3] = prims[B3_con];

  // printf("cons[BCONI] = %.15e %.15e %.15e ( %d=%d , %d=%d , %d=%d )\n",cons[BCON1],cons[BCON2],cons[BCON3],
         // BCON1,B1_con, BCON2,B2_con, BCON3,B3_con);

  // Bunch of auxiliary vars
  const double uvec0_con=lorentz/alpha;
  const double uvec1_con=prims[u1_con]-alpha*lorentz*geom.gcon[0][1];
  const double uvec2_con=prims[u2_con]-alpha*lorentz*geom.gcon[0][2];
  const double uvec3_con=prims[u3_con]-alpha*lorentz*geom.gcon[0][3];

  const double B1_cov = geom.gcov[1][1]*prims[B1_con]+geom.gcov[1][2]*prims[B2_con]+geom.gcov[1][3]*prims[B3_con];
  const double B2_cov = geom.gcov[2][1]*prims[B1_con]+geom.gcov[2][2]*prims[B2_con]+geom.gcov[2][3]*prims[B3_con];
  const double B3_cov = geom.gcov[3][1]*prims[B1_con]+geom.gcov[3][2]*prims[B2_con]+geom.gcov[3][3]*prims[B3_con];

  const double uvec0_cov=geom.gcov[0][0]*uvec0_con+geom.gcov[0][1]*uvec1_con+geom.gcov[0][2]*uvec2_con+geom.gcov[0][3]*uvec3_con;
  const double uvec1_cov=geom.gcov[1][0]*uvec0_con+geom.gcov[1][1]*uvec1_con+geom.gcov[1][2]*uvec2_con+geom.gcov[1][3]*uvec3_con;
  const double uvec2_cov=geom.gcov[2][0]*uvec0_con+geom.gcov[2][1]*uvec1_con+geom.gcov[2][2]*uvec2_con+geom.gcov[2][3]*uvec3_con;
  const double uvec3_cov=geom.gcov[3][0]*uvec0_con+geom.gcov[3][1]*uvec1_con+geom.gcov[3][2]*uvec2_con+geom.gcov[3][3]*uvec3_con;

  const double u0_cov=geom.gcov[0][1]*prims[u1_con]+geom.gcov[0][2]*prims[u2_con]+geom.gcov[0][3]*prims[u3_con];
  const double u1_cov=geom.gcov[1][1]*prims[u1_con]+geom.gcov[1][2]*prims[u2_con]+geom.gcov[1][3]*prims[u3_con];
  const double u2_cov=geom.gcov[2][1]*prims[u1_con]+geom.gcov[2][2]*prims[u2_con]+geom.gcov[2][3]*prims[u3_con];
  const double u3_cov=geom.gcov[3][1]*prims[u1_con]+geom.gcov[3][2]*prims[u2_con]+geom.gcov[3][3]*prims[u3_con];

  const double b0_con=1.0/alpha * (prims[B1_con]*u1_cov+prims[B2_con]*u2_cov+prims[B3_con]*u3_cov);

  //const double b0_cov=1./lorentz * alpha*u0_cov*b0_con;
  const double b1_cov=B1_cov/lorentz+1./lorentz * alpha*u1_cov*b0_con;
  const double b2_cov=B2_cov/lorentz+1./lorentz * alpha*u2_cov*b0_con;
  const double b3_cov=B3_cov/lorentz+1./lorentz * alpha*u3_cov*b0_con;

  const double b1_con=prims[B1_con]/lorentz+b0_con*(alpha*prims[u1_con]/lorentz-pow(alpha,2)*geom.gcon[0][1]);
  const double b2_con=prims[B2_con]/lorentz+b0_con*(alpha*prims[u2_con]/lorentz-pow(alpha,2)*geom.gcon[0][2]);
  const double b3_con=prims[B3_con]/lorentz+b0_con*(alpha*prims[u3_con]/lorentz-pow(alpha,2)*geom.gcon[0][3]);

  const double Bsquared=prims[B1_con]*B1_cov + prims[B2_con]*B2_cov + prims[B3_con]*B3_cov;
  const double bsquared=(Bsquared+pow((alpha*b0_con),2))/pow(lorentz,2);
  
  const double w_enthalpy=prims[RHO]+prims[UU]+prims[PRESS];

  const double b0_cov=-b0_con;//for minkowski//1./b0_con * (b1_cov*b1_con+b2_cov*b2_con+b3_cov*b3_con);

  // b^2
  *bsq = bsquared;

  // u
  if (bsquared==0.){
    cons[UU]=lorentz*(w_enthalpy)*uvec0_cov+alpha*(prims[PRESS])+prims[RHO]*lorentz/alpha;
  }
  else {
    cons[UU]=lorentz*(w_enthalpy+bsquared)*uvec0_cov+alpha*(prims[PRESS]+bsquared/2.)+(-alpha*b0_con)*b0_cov+prims[RHO]*lorentz/alpha;
  }

  // printf("cons[UU] = %.15e\n",cons[UU]);

  // Q_{mu}
  cons[QCOV1]=lorentz*(w_enthalpy+bsquared)*uvec1_cov+(-alpha*b0_con)*b1_cov;
  cons[QCOV2]=lorentz*(w_enthalpy+bsquared)*uvec2_cov+(-alpha*b0_con)*b2_cov;
  cons[QCOV3]=lorentz*(w_enthalpy+bsquared)*uvec3_cov+(-alpha*b0_con)*b3_cov;

  // Ye
  cons[YE] = prims[YE]*cons[RHO];

  // printf("cons[YE] = %.15e\n",cons[YE]);

}

// Set HARM3D metric from IGM metric
void set_harm_metric_from_IGM_metric( const CCTK_REAL *METRIC,
                                      const CCTK_REAL *METRIC_LAP_PSI4,
                                      const CCTK_REAL g4dn[4][4],
                                      const CCTK_REAL g4up[4][4],
                                      struct of_geom *geom ) {
  // First g_{munu} and g^{munu}
  for(int mu=0;mu<4;mu++)
    for(int nu=0;nu<4;nu++) {
      geom->gcov[mu][nu] = g4dn[mu][nu];
      geom->gcon[mu][nu] = g4up[mu][nu];
    }
  // Then the determinant of g_{munu}
  geom->g     = METRIC_LAP_PSI4[LAPSE]*METRIC_LAP_PSI4[PSI6];
  // It's inverse
  geom->g_inv = 1.0/geom->g;
  // The lapse
  geom->alpha = METRIC_LAP_PSI4[LAPSE];
  // And the shift
  for(int i=0;i<3;i++)
    geom->beta[i] = METRIC[SHIFTX+i];
  // Then n^{mu} = (1/alpha,-beta^{i}/alpha)
  geom->ncon[0] = 1.0/geom->alpha;
  for(int i=0;i<3;i++)
    geom->ncon[i+1] = -geom->beta[i]/geom->alpha;
}

// Set HARM3D prims from IGM prims
void set_harm_prims_from_IGM_prims( const CCTK_REAL *IGM_prims,
                                    const CCTK_REAL T_local,
                                    const CCTK_REAL Ye_local,
                                    const CCTK_REAL eps_local,
                                    const CCTK_REAL gamma_local,
                                    CCTK_REAL *harm_prims) {

  // rho
  harm_prims[RHO     ] = IGM_prims[RHOB];
  // eps
  harm_prims[EPS     ] = eps_local;
  // u = rho*eps
  harm_prims[UU      ] = eps_local * IGM_prims[RHOB];
  // P
  harm_prims[PRESS   ] = IGM_prims[PRESSURE];
  // B^i
  harm_prims[B1_con  ] = IGM_prims[BX_CENTER];
  harm_prims[B2_con  ] = IGM_prims[BY_CENTER];
  harm_prims[B3_con  ] = IGM_prims[BZ_CENTER];
  // u^i. This is not the actual relation between
  // the velocities in IGM and HARM3D, but since
  // we have set up initial data for IGM's velocities
  // based on Ari's script, they should be the same.
  harm_prims[u1_con  ] = IGM_prims[VX];
  harm_prims[u2_con  ] = IGM_prims[VY];
  harm_prims[u3_con  ] = IGM_prims[VZ];
  // T
  harm_prims[TEMP    ] = T_local;
  // Ye
  harm_prims[YE      ] = Ye_local;
  // Lorentz factor
  harm_prims[WLORENTZ] = gamma_local;
}

// This function sets the ADM 3-metric as a perturbation around Minkoswski spacetime
void set_ADM_3metric( const CCTK_REAL perturbation_strength, CCTK_REAL *ADM_METRIC ) {

  ADM_METRIC[GXX  ]  = 1.0 + perturbation_strength*random_number_between_m1_and_p1;
  ADM_METRIC[GXY  ]  = 0.0 + perturbation_strength*random_number_between_m1_and_p1;
  ADM_METRIC[GXZ  ]  = 0.0 + perturbation_strength*random_number_between_m1_and_p1;
  ADM_METRIC[GYY  ]  = 1.0 + perturbation_strength*random_number_between_m1_and_p1;
  ADM_METRIC[GYZ  ]  = 0.0 + perturbation_strength*random_number_between_m1_and_p1;
  ADM_METRIC[GZZ  ]  = 1.0 + perturbation_strength*random_number_between_m1_and_p1;
  ADM_METRIC[LAPM1]  = 0.0 + perturbation_strength*random_number_between_m1_and_p1;
  ADM_METRIC[SHIFTX] = 0.0 + perturbation_strength*random_number_between_m1_and_p1;
  ADM_METRIC[SHIFTY] = 0.0 + perturbation_strength*random_number_between_m1_and_p1;
  ADM_METRIC[SHIFTZ] = 0.0 + perturbation_strength*random_number_between_m1_and_p1;

}

// This function is an adaptation of the IllinoisGRMHD function
// IllinoisGRMHD_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij()
void IllinoisGRMHD_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij( CCTK_REAL *ADM_METRIC, CCTK_REAL *BSSN_METRIC ) {

  CCTK_REAL gxx_physL = ADM_METRIC[GXX];
  CCTK_REAL gxy_physL = ADM_METRIC[GXY];
  CCTK_REAL gxz_physL = ADM_METRIC[GXZ];
  CCTK_REAL gyy_physL = ADM_METRIC[GYY];
  CCTK_REAL gyz_physL = ADM_METRIC[GYZ];
  CCTK_REAL gzz_physL = ADM_METRIC[GZZ];

  /**********************************************************************
   * Compute \tilde{\gamma_{ij}}, phi, and psi (BSSN) from g_{ij} (ADM) *
   **********************************************************************/
  CCTK_REAL gijdet = gxx_physL * gyy_physL * gzz_physL + gxy_physL * gyz_physL * gxz_physL + gxz_physL * gxy_physL * gyz_physL
    - gxz_physL * gyy_physL * gxz_physL - gxy_physL * gxy_physL * gzz_physL - gxx_physL * gyz_physL * gyz_physL;

  gijdet = fabs(gijdet);

   CCTK_REAL phiL = (1.0/12.0) * log(gijdet);
   CCTK_REAL psiL = exp(phiL);

   CCTK_REAL Psim4 = 1.0/(psiL*psiL*psiL*psiL);
   CCTK_REAL gtxxL = gxx_physL*Psim4;
   CCTK_REAL gtxyL = gxy_physL*Psim4;
   CCTK_REAL gtxzL = gxz_physL*Psim4;
   CCTK_REAL gtyyL = gyy_physL*Psim4;
   CCTK_REAL gtyzL = gyz_physL*Psim4;
   CCTK_REAL gtzzL = gzz_physL*Psim4;

   /*********************************
    * Apply det gtij = 1 constraint *
    *********************************/
   CCTK_REAL gtijdet = gtxxL * gtyyL * gtzzL + gtxyL * gtyzL * gtxzL + gtxzL * gtxyL * gtyzL -
     gtxzL * gtyyL * gtxzL - gtxyL * gtxyL * gtzzL - gtxxL * gtyzL * gtyzL;

   CCTK_REAL gtijdet_Fm1o3 = fabs(1.0/cbrt(gtijdet));

   gtxxL = gtxxL * gtijdet_Fm1o3;
   gtxyL = gtxyL * gtijdet_Fm1o3;
   gtxzL = gtxzL * gtijdet_Fm1o3;
   gtyyL = gtyyL * gtijdet_Fm1o3;
   gtyzL = gtyzL * gtijdet_Fm1o3;
   gtzzL = gtzzL * gtijdet_Fm1o3;

   CCTK_REAL Psi4 = psiL*psiL*psiL*psiL;
   /*****************************************
    * Set all the needed BSSN gridfunctions *
    *****************************************/
   // Set phi and psi
   BSSN_METRIC[PHI]    = phiL;
   BSSN_METRIC[PSI]    = psiL;

   // Set the lapse and shift
   BSSN_METRIC[LAPM1]  = ADM_METRIC[LAPM1];
   BSSN_METRIC[SHIFTX] = ADM_METRIC[SHIFTX];
   BSSN_METRIC[SHIFTY] = ADM_METRIC[SHIFTY];
   BSSN_METRIC[SHIFTZ] = ADM_METRIC[SHIFTZ];

   // \tilde{g}_{ij}
   BSSN_METRIC[GXX]    = gtxxL;
   BSSN_METRIC[GXY]    = gtxyL;
   BSSN_METRIC[GXZ]    = gtxzL;
   BSSN_METRIC[GYY]    = gtyyL;
   BSSN_METRIC[GYZ]    = gtyzL;
   BSSN_METRIC[GZZ]    = gtzzL;

   // \tilde{g}^{ij}
   BSSN_METRIC[GUPXX]  =   ( gtyyL * gtzzL - gtyzL * gtyzL );
   BSSN_METRIC[GUPXY]  = - ( gtxyL * gtzzL - gtyzL * gtxzL );
   BSSN_METRIC[GUPXZ]  =   ( gtxyL * gtyzL - gtyyL * gtxzL );
   BSSN_METRIC[GUPYY]  =   ( gtxxL * gtzzL - gtxzL * gtxzL );
   BSSN_METRIC[GUPYZ]  = - ( gtxxL * gtyzL - gtxyL * gtxzL );
   BSSN_METRIC[GUPZZ]  =   ( gtxxL * gtyyL - gtxyL * gtxyL );

      // Recompute ADM 3-metric
   ADM_METRIC[GXX]     = gtxxL*Psi4;
   ADM_METRIC[GXY]     = gtxyL*Psi4;
   ADM_METRIC[GXZ]     = gtxzL*Psi4;
   ADM_METRIC[GYY]     = gtyyL*Psi4;
   ADM_METRIC[GYZ]     = gtyzL*Psi4;
   ADM_METRIC[GZZ]     = gtzzL*Psi4;

   // Compute inverse ADM 3-metric
   ADM_METRIC[GUPXX]   = BSSN_METRIC[GUPXX]*Psim4;
   ADM_METRIC[GUPXY]   = BSSN_METRIC[GUPXY]*Psim4;
   ADM_METRIC[GUPXZ]   = BSSN_METRIC[GUPXZ]*Psim4;
   ADM_METRIC[GUPYY]   = BSSN_METRIC[GUPYY]*Psim4;
   ADM_METRIC[GUPYZ]   = BSSN_METRIC[GUPYZ]*Psim4;
   ADM_METRIC[GUPZZ]   = BSSN_METRIC[GUPZZ]*Psim4;

}

// This function sets the fluid 3-velocity
// and is an adaptation of Ari's code
void set_velocities( const CCTK_REAL gamma, CCTK_REAL *PRIMS ) {
  // Velocity magnitude
  const CCTK_REAL v = sqrt(1.0-1.0/(gamma*gamma));
  // Randomly choose vx, vy, vz based on v
  const CCTK_REAL vx = v*random_number_between_m1_and_p1;
  const CCTK_REAL vy = sqrt(v*v - vx*vx)*random_number_between_m1_and_p1;
  const CCTK_REAL vz = sqrt(v*v - vx*vx - vy*vy);
  // Now check if everything is okay
  const CCTK_REAL vsq = vx*vx + vy*vy + vz*vz;
  if( fabs(v*v - vsq)/vsq > 1e-10 ) {
    fprintf(stderr,"(WARNING) Possible problem with the velocities: v*v = %.15e | vsq = %.15e\n",v*v,vsq);
  }
  PRIMS[VX] = vx;
  PRIMS[VY] = vy;
  PRIMS[VZ] = vz;
}

// This function sets the magnetic fields
// and is an adaptation of Ari's code
void set_magnetic_fields( const CCTK_REAL log_Pmag_o_P, CCTK_REAL *PRIMS ) {

  // Set velocities
  const CCTK_REAL vx = PRIMS[VX];
  const CCTK_REAL vy = PRIMS[VY];
  const CCTK_REAL vz = PRIMS[VZ];

  // The magnitude of the magnetic field can be computed using
  //    log_Pmag_o_P = log(Pmag/P)
  // => 10^(log_Pmag_o_P) = Pmag/P
  // => P_mag = B^2/2 = 10^(log_P_o_Pmag) P
  //    .----------------------------------.
  // => | B = sqrt( 10^(log_Pmag_o_P) 2P ) |
  //    .----------------------------------.
  //
  // This is not entirely accurate since P_mag := b^2/2, but
  // it is a good way of getting the magnitude of B reasonably
  // close to what we want.
  const CCTK_REAL B = sqrt( 2.0*PRIMS[PRESSURE]*pow(10.0,log_Pmag_o_P) );
  // Randomly choose Bx, By, Bz based on B
  CCTK_REAL Bxhat = random_number_between_m1_and_p1;
  CCTK_REAL Byhat = sqrt(1.0 - Bxhat*Bxhat)*random_number_between_m1_and_p1;
  CCTK_REAL Bzhat = (-Bxhat*vx - Byhat*vy)/vz;
  CCTK_REAL Bnorm = sqrt(Bxhat*Bxhat + Byhat*Byhat + Bzhat*Bzhat);
  Bxhat /= Bnorm;
  Byhat /= Bnorm;
  Bzhat /= Bnorm;
  // Check orthogonality between B^i and v^i
  if( fabs(Bxhat*vx + Byhat*vy + Bzhat*vz) > 1e-10 ) {
    fprintf(stderr,"(WARNING) Orthogonality problem: |Bhat.v| = %.15e\n",fabs(Bxhat*vx + Byhat*vy + Bzhat*vz));
  }
  // Check if Bhat is unitary
  if(fabs(1.0 -(Bxhat*Bxhat+Byhat*Byhat+Bzhat*Bzhat))>1.0e-10) {
    fprintf(stderr,"(WARNING) Bhat is not unitary: Bhat.Bhat = %.15e\n",Bxhat*Bxhat+Byhat*Byhat+Bzhat*Bzhat);
  }
  PRIMS[BX_CENTER] = -B*Bxhat;
  PRIMS[BY_CENTER] = -B*Byhat;
  PRIMS[BZ_CENTER] = -B*Bzhat;
}

// Adaptation of the code in IllinoisGRMHD C2P driver
void compute_ADM_4metric_and_inverse( const CCTK_REAL *ADM_METRIC, const CCTK_REAL *METRIC_LAP_PSI4, CCTK_REAL g4dn[4][4], CCTK_REAL g4up[4][4] ) {

  CCTK_REAL shift_xL = ADM_METRIC[GXX]*ADM_METRIC[SHIFTX] + ADM_METRIC[GXY]*ADM_METRIC[SHIFTY] + ADM_METRIC[GXZ]*ADM_METRIC[SHIFTZ];
  CCTK_REAL shift_yL = ADM_METRIC[GXY]*ADM_METRIC[SHIFTX] + ADM_METRIC[GYY]*ADM_METRIC[SHIFTY] + ADM_METRIC[GYZ]*ADM_METRIC[SHIFTZ];
  CCTK_REAL shift_zL = ADM_METRIC[GXZ]*ADM_METRIC[SHIFTX] + ADM_METRIC[GYZ]*ADM_METRIC[SHIFTY] + ADM_METRIC[GZZ]*ADM_METRIC[SHIFTZ];
  CCTK_REAL beta2L   = shift_xL*ADM_METRIC[SHIFTX] + shift_yL*ADM_METRIC[SHIFTY] + shift_zL*ADM_METRIC[SHIFTZ];


  // Compute 4-metric, both g_{\mu \nu} and g^{\mu \nu}.
  // This is for computing T_{\mu \nu} and T^{\mu \nu}. Also the HARM con2prim lowlevel function requires them.
  g4dn[0][0] = -SQR(METRIC_LAP_PSI4[LAPSE]) + beta2L;
  g4dn[0][1] = g4dn[1][0] = shift_xL;
  g4dn[0][2] = g4dn[2][0] = shift_yL;
  g4dn[0][3] = g4dn[3][0] = shift_zL;
  g4dn[1][1]              = ADM_METRIC[GXX];
  g4dn[1][2] = g4dn[2][1] = ADM_METRIC[GXY];
  g4dn[1][3] = g4dn[3][1] = ADM_METRIC[GXZ];
  g4dn[2][2]              = ADM_METRIC[GYY];
  g4dn[2][3] = g4dn[3][2] = ADM_METRIC[GYZ];
  g4dn[3][3]              = ADM_METRIC[GZZ];

  CCTK_REAL alpha_inv_squared=SQR(METRIC_LAP_PSI4[LAPSEINV]);
  g4up[0][0] = -1.0*alpha_inv_squared;
  g4up[0][1] = g4up[1][0] = ADM_METRIC[SHIFTX]*alpha_inv_squared;
  g4up[0][2] = g4up[2][0] = ADM_METRIC[SHIFTY]*alpha_inv_squared;
  g4up[0][3] = g4up[3][0] = ADM_METRIC[SHIFTZ]*alpha_inv_squared;
  g4up[1][1]              = ADM_METRIC[GUPXX] - ADM_METRIC[SHIFTX]*ADM_METRIC[SHIFTX]*alpha_inv_squared;
  g4up[1][2] = g4up[2][1] = ADM_METRIC[GUPXY] - ADM_METRIC[SHIFTX]*ADM_METRIC[SHIFTY]*alpha_inv_squared;
  g4up[1][3] = g4up[3][1] = ADM_METRIC[GUPXZ] - ADM_METRIC[SHIFTX]*ADM_METRIC[SHIFTZ]*alpha_inv_squared;
  g4up[2][2]              = ADM_METRIC[GUPYY] - ADM_METRIC[SHIFTY]*ADM_METRIC[SHIFTY]*alpha_inv_squared;
  g4up[2][3] = g4up[3][2] = ADM_METRIC[GUPYZ] - ADM_METRIC[SHIFTY]*ADM_METRIC[SHIFTZ]*alpha_inv_squared;
  g4up[3][3]              = ADM_METRIC[GUPZZ] - ADM_METRIC[SHIFTZ]*ADM_METRIC[SHIFTZ]*alpha_inv_squared;

}

      /* printf("PRIMS before: rho = %e | Ye = %.1lf | T = %lf | P = %e | vi = %e,%e,%e | Bi = %e,%e,%e\n", */
      /*        PRIMS[RHOB],Ye_local,T_local,PRIMS[PRESSURE],PRIMS[VX],PRIMS[VY],PRIMS[VZ],PRIMS[BX_CENTER],PRIMS[BY_CENTER],PRIMS[BZ_CENTER]); */

      /* printf("PRIMS after : rho = %e | Ye = %.1lf | T = %lf | P = %e | vi = %e,%e,%e | Bi = %e,%e,%e\n", */
      /*        PRIMS[RHOB],Ye_local,T_local,PRIMS[PRESSURE],PRIMS[VX],PRIMS[VY],PRIMS[VZ],PRIMS[BX_CENTER],PRIMS[BY_CENTER],PRIMS[BZ_CENTER]); */

      /* printf("CONSERVS    : rho_star = %e | Ye_star = %e | tau = %e | Si = %e,%e,%e\n", */
      /*        CONSERVS[RHOSTAR],Ye_star_local,CONSERVS[TAUENERGY],CONSERVS[STILDEX],CONSERVS[STILDEY],CONSERVS[STILDEZ]); */

      /* printf("Lapse and shift: alpha = %e,%e,%e | betai = %e,%e %e,%e %e,%e\n", */
      /*        METRIC_LAP_PSI4[LAPSE],METRIC[LAPM1]+1.0,METRIC_PHYS[LAPM1]+1.0, */
      /*        METRIC_PHYS[SHIFTX],METRIC[SHIFTX],METRIC_PHYS[SHIFTY],METRIC[SHIFTY],METRIC_PHYS[SHIFTZ],METRIC[SHIFTZ]); */

      /* printf("phi and psi: phi = %e | psi = %e\n",METRIC[PHI],METRIC[PSI]); */
      
      /* printf("ADM  : g_xx = %e | g_xy = %e | g_xz = %e | g_yy = %e | g_yz = %e | g_zz = %e\n", */
      /*        METRIC_PHYS[GXX  ],METRIC_PHYS[GXY  ],METRIC_PHYS[GXZ  ],METRIC_PHYS[GYY  ],METRIC_PHYS[GYZ  ],METRIC_PHYS[GZZ  ]); */
      /* printf("ADM  : gxx  = %e | gxy  = %e | gxz  = %e | gyy  = %e | gyz  = %e | gzz  = %e\n", */
      /*        METRIC_PHYS[GUPXX],METRIC_PHYS[GUPXY],METRIC_PHYS[GUPXZ],METRIC_PHYS[GUPYY],METRIC_PHYS[GUPYZ],METRIC_PHYS[GUPZZ]); */

      /* printf("BSSN : g_xx = %e | g_xy = %e | g_xz = %e | g_yy = %e | g_yz = %e | g_zz = %e\n", */
      /*        METRIC[GXX  ],METRIC[GXY  ],METRIC[GXZ  ],METRIC[GYY  ],METRIC[GYZ  ],METRIC[GZZ  ]); */
      /* printf("BSSN : gxx  = %e | gxy  = %e | gxz  = %e | gyy  = %e | gyz  = %e | gzz  = %e\n", */
      /*        METRIC[GUPXX],METRIC[GUPXY],METRIC[GUPXZ],METRIC[GUPYY],METRIC[GUPYZ],METRIC[GUPZZ]); */
