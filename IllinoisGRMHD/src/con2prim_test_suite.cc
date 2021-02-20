// Thorn      : IllinoisGRMHD
// File       : con2prim_test_suite.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: In this file we provide an extensive test suite of
//              the con2prim routines available in the code.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "IllinoisGRMHD_headers.h"
#include "EOS_headers.hh"
#include "con2prim_headers.h"

#include <fstream>
using namespace std;

inline CCTK_REAL relative_error( const CCTK_REAL a, const CCTK_REAL b ) {
  if     ( a != 0 ) return( fabs(1.0-b/a) );
  else if( b != 0 ) return( fabs(1.0-a/b) );
  else              return( 0.0 );
}


extern "C"
void IllinoisGRMHD_con2prim_test_suit( CCTK_ARGUMENTS ) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Initialize the EOS parameters
  igm_eos_parameters eos;
  initialize_igm_eos_parameters_from_input(igm_eos_key,cctk_time,eos);

  // Print information
  if( CCTK_EQUALS(igm_eos_type,"Tabulated") || CCTK_EQUALS(igm_eos_type,"nuc_eos") )
    CCTK_VInfo(CCTK_THORNSTRING,"EOS type        : Tabulated");

  // Count number of routines tested
  CCTK_INT num_routines_tested = 1;
  CCTK_INT con2prim_test_keys[4];
  char con2prim_test_names[4][100];
  con2prim_test_keys[0] = eos.c2p_routine;
  sprintf(con2prim_test_names[0],"%s",igm_con2prim_routine);
  if( eos.c2p_backup[0] != None ) {
    num_routines_tested++;
    con2prim_test_keys[1] = eos.c2p_backup[0];
    sprintf(con2prim_test_names[1],"%s",igm_con2prim_backup_routine[0]);
    if( eos.c2p_backup[1] != None ) {
      num_routines_tested++;
      con2prim_test_keys[2] = eos.c2p_backup[1];
      sprintf(con2prim_test_names[2],"%s",igm_con2prim_backup_routine[1]);
      if( eos.c2p_backup[2] != None ) {
        num_routines_tested++;
        con2prim_test_keys[3] = eos.c2p_backup[2];
        sprintf(con2prim_test_names[3],"%s",igm_con2prim_backup_routine[2]);
      }
    }
  }

  CCTK_VInfo(CCTK_THORNSTRING,"Tested routine 1: %s",igm_con2prim_routine);
  for(int i=0;i<num_routines_tested-1;i++)
    CCTK_VInfo(CCTK_THORNSTRING,"Tested routine %d: %s",i+2,igm_con2prim_backup_routine[i]);

  // We will be performing the tabulated EOS test in the following way:
  //
  //      Y_e      = 0.1
  //       W       = 2
  // log10(Pmag/P) = -5
  //
  // rho will vary between rho_min and rho_max (uniformly in log space)
  //  T  will vary between  T_min  and  T_max  (uniformly in log space)

  // Number of points in the discretization of rho and T
  const CCTK_INT npoints       = igm_con2prim_standalone_npoints;

  const CCTK_REAL test_rho_min = igm_con2prim_standalone_rho_min;
  const CCTK_REAL test_rho_max = igm_con2prim_standalone_rho_max;

  const CCTK_REAL test_T_min   = igm_con2prim_standalone_T_min;
  const CCTK_REAL test_T_max   = igm_con2prim_standalone_T_max;

  // Compute the density step size
  const CCTK_REAL lrmin        = log(test_rho_min);
  const CCTK_REAL lrmax        = log(test_rho_max);
  const CCTK_REAL dlr          = (lrmax - lrmin)/npoints;

  // Compute the temperature step size
  const CCTK_REAL ltmin        = log(test_T_min);
  const CCTK_REAL ltmax        = log(test_T_max);
  const CCTK_REAL dlt          = (ltmax - ltmin)/npoints;
  
  // Fix Y_e
  const CCTK_REAL Ye_test = 0.1;
  // Fix W
  const CCTK_REAL W_test  = 2.0;
  // Fix log10(Pmag/P)
  const CCTK_REAL logPmoP = -5.0;

  // Useful for printing information
  char primnames[MAXNUMVARS][4];
  sprintf(primnames[RHOB       ],"rho");
  sprintf(primnames[YEPRIM     ],"Ye ");
  sprintf(primnames[TEMPERATURE],"T  ");
  sprintf(primnames[PRESSURE   ],"P  ");
  sprintf(primnames[EPSILON    ],"eps");
  sprintf(primnames[VX         ],"v^x");
  sprintf(primnames[VY         ],"v^y");
  sprintf(primnames[VZ         ],"v^z");
  sprintf(primnames[BX_CENTER  ],"B^x");
  sprintf(primnames[BY_CENTER  ],"B^y");
  sprintf(primnames[BZ_CENTER  ],"B^z");

  // Which variables we want to compute the errors
  int which_prims_in_error[MAXNUMVARS];
  int ww = 0;
  which_prims_in_error[ww++] = RHOB;
  which_prims_in_error[ww++] = YEPRIM;
  which_prims_in_error[ww++] = TEMPERATURE;
  which_prims_in_error[ww++] = PRESSURE;
  which_prims_in_error[ww++] = EPSILON;
  which_prims_in_error[ww++] = VX;
  which_prims_in_error[ww++] = VY;
  which_prims_in_error[ww++] = VZ;
  which_prims_in_error[ww++] = BX_CENTER;
  which_prims_in_error[ww++] = BY_CENTER;
  which_prims_in_error[ww++] = BZ_CENTER;
  int num_prims_in_error = ww;

  // tau is given by (see :
  //
  // tau := hW^{2} + B^{2} - P - 0.5*( (B.v)^{2} + (B/W)^{2} )
  // Absolutely minimum allowed tau

  // Now perform one test for each of the selected routines
  for(int which_routine=0;which_routine<num_routines_tested;which_routine++) {

    CCTK_VInfo(CCTK_THORNSTRING,"Beginning test for routine %s",con2prim_test_names[which_routine]);
    
    char filename[100];
    sprintf(filename,"con2prim_out_%s.asc",con2prim_test_names[which_routine]);
    FILE* outfile = fopen(filename,"w");

    int failures = 0;

    srand(0);
    
    for(int i=0;i<npoints;i++) {
      for(int j=0;j<npoints;j++) {

        // Start by setting the prims (rho,Ye,T,P,eps)
        CCTK_REAL xrho  = exp(lrmin + dlr*i);
        CCTK_REAL xtemp = exp(ltmin + dlt*j);
        CCTK_REAL xye   = Ye_test;
        CCTK_REAL xprs  = 0.0;
        CCTK_REAL xeps  = 0.0;
        WVU_EOS_P_and_eps_from_rho_Ye_T( xrho,xye,xtemp, &xprs,&xeps );

        // Now set the velocities
        // Velocity magnitude
        const CCTK_REAL v = sqrt(1.0-1.0/(W_test*W_test));
        const CCTK_REAL vx = v*((CCTK_REAL)rand())/((CCTK_REAL)RAND_MAX);
        const CCTK_REAL vy = sqrt(v*v - vx*vx)*((CCTK_REAL)rand())/((CCTK_REAL)RAND_MAX);
        const CCTK_REAL vz = sqrt(v*v - vx*vx - vy*vy);

        // Now the magnetic fields. We'll set them aligned
        // with the velocities, for simplicity.
        const CCTK_REAL Bhatx = vx/v;
        const CCTK_REAL Bhaty = vy/v;
        const CCTK_REAL Bhatz = vz/v;
        const CCTK_REAL B     = sqrt(2.0*pow(10.0,logPmoP)*xprs);
        const CCTK_REAL Bx    = -Bhatx * B;
        const CCTK_REAL By    = -Bhaty * B;
        const CCTK_REAL Bz    = -Bhatz * B;
        CCTK_REAL CONSERVS[NUM_CONSERVS];
        CCTK_REAL METRIC[NUMVARS_FOR_METRIC];
        CCTK_REAL METRIC_PHYS[NUMVARS_FOR_METRIC];
        CCTK_REAL METRIC_LAP_PSI4[NUMVARS_METRIC_AUX];
        
        // {
        //   ifstream myfile;
        //   char filename[100];
        //   sprintf(filename,"con2prim_debug_high_density_region.asc");
        //   myfile.open(filename, ios::in | ios::binary);
        //   if(myfile.fail()) {
        //     fprintf(stderr,"Error: file %s cannot be opened.\n",filename);
        //     exit(1);
        //   }
        //   myfile.read((char*)CONSERVS       ,(   NUM_CONSERVS   )*sizeof(CCTK_REAL));
        //   myfile.read((char*)METRIC         ,(NUMVARS_FOR_METRIC)*sizeof(CCTK_REAL));
        //   myfile.read((char*)METRIC_PHYS    ,(NUMVARS_FOR_METRIC)*sizeof(CCTK_REAL));
        //   myfile.read((char*)METRIC_LAP_PSI4,(NUMVARS_METRIC_AUX)*sizeof(CCTK_REAL));
        //   myfile.close();
        // }

        // Now set the primitive variables array, following IGM's standards
        CCTK_REAL PRIMS[MAXNUMVARS];
        PRIMS[RHOB       ] = xrho;
        PRIMS[YEPRIM     ] = xye;
        PRIMS[TEMPERATURE] = xtemp;
        PRIMS[PRESSURE   ] = xprs;
        PRIMS[EPSILON    ] = xeps;
        PRIMS[VX         ] = vx;
        PRIMS[VY         ] = vy;
        PRIMS[VZ         ] = vz;
        PRIMS[BX_CENTER  ] = Bx;
        PRIMS[BY_CENTER  ] = By;
        PRIMS[BZ_CENTER  ] = Bz;

        // Store original prims
        CCTK_REAL PRIMS_ORIG[MAXNUMVARS];
        for(int i=0;i<MAXNUMVARS;i++) PRIMS_ORIG[i] = PRIMS[i];

        // Set the metric to flat space
        METRIC[PHI   ] = 0.0;
        METRIC[LAPM1 ] = 0.0;
        METRIC[SHIFTX] = METRIC[SHIFTY] = METRIC[SHIFTZ] = 0.0;
        METRIC[GXX   ] = METRIC[GYY   ] = METRIC[GZZ   ] = 1.0;
        METRIC[GXY   ] = METRIC[GXZ   ] = METRIC[GYZ   ] = 0.0;
        METRIC[GUPXX ] = METRIC[GUPYY ] = METRIC[GUPZZ ] = 1.0;
        METRIC[GUPXY ] = METRIC[GUPXZ ] = METRIC[GUPYZ ] = 0.0;

        // We'll also need the auxilary metric variables array
        SET_LAPSE_PSI4(METRIC_LAP_PSI4,METRIC);

        // Now the physical metric
        METRIC_PHYS[GXX  ] = METRIC[GXX  ]*METRIC_LAP_PSI4[PSI4 ];
        METRIC_PHYS[GXY  ] = METRIC[GXY  ]*METRIC_LAP_PSI4[PSI4 ];
        METRIC_PHYS[GXZ  ] = METRIC[GXZ  ]*METRIC_LAP_PSI4[PSI4 ];
        METRIC_PHYS[GYY  ] = METRIC[GYY  ]*METRIC_LAP_PSI4[PSI4 ];
        METRIC_PHYS[GYZ  ] = METRIC[GYZ  ]*METRIC_LAP_PSI4[PSI4 ];
        METRIC_PHYS[GZZ  ] = METRIC[GZZ  ]*METRIC_LAP_PSI4[PSI4 ];
        METRIC_PHYS[GUPXX] = METRIC[GUPXX]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPXY] = METRIC[GUPXY]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPXZ] = METRIC[GUPXZ]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPYY] = METRIC[GUPYY]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPYZ] = METRIC[GUPYZ]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPZZ] = METRIC[GUPZZ]*METRIC_LAP_PSI4[PSIM4];

        // Then set the conservative variables array
        const int already_computed_physical_metric_and_inverse=0;
        CCTK_REAL TUPMUNU[10],TDNMUNU[10];
        CCTK_REAL g4dn[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
        CCTK_REAL g4up[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
        struct output_stats stats;
        IllinoisGRMHD_enforce_limits_on_primitives_and_recompute_conservs(already_computed_physical_metric_and_inverse,
                                                                          PRIMS,stats,eos,METRIC,g4dn,g4up,TUPMUNU,TDNMUNU,CONSERVS);

        // The con2prim routines require different conservative variables than
        // those that IllinoisGRMHD evolve. Therefore, we must first convert
        // them into a new set of variables, suitable for primitive recovery.
        CCTK_REAL cons[numcons];
        set_cons_from_PRIMS_and_CONSERVS(eos,con2prim_test_keys[which_routine],
                                         METRIC,METRIC_LAP_PSI4,
                                         PRIMS,CONSERVS, cons);

        // printf("Input conservs:\n");
        // for(int i=0;i<NUM_CONSERVS;i++) printf("%2d: %e\n",i,CONSERVS[i]);

        // printf("Computed cons:\n");
        // for(int i=0;i<numcons;i++) printf("%2d: %e\n",i,cons[i]);

        // printf("(Conformal) Metric quantities:\n");
        // for(int i=0;i<NUMVARS_FOR_METRIC;i++) printf("%2d: %e\n",i,METRIC[i]);

        // printf("(Physical)  Metric quantities:\n");
        // for(int i=0;i<NUMVARS_FOR_METRIC;i++) printf("%2d: %e\n",i,METRIC_PHYS[i]);

        // printf("(Auxiliary) Metric quantities:\n");
        // for(int i=0;i<NUMVARS_METRIC_AUX;i++) printf("%2d: %e\n",i,METRIC_LAP_PSI4[i]);

        // The con2prim routines require primitive guesses in order to perform
        // the recovery. In IllinoisGRMHD, we do not keep track of the primitives
        // in between time steps, and therefore our guesses are *not* the
        // values of the primitives in the previous time level. Instead, we provide
        // guesses based on the conservative variables.
        CCTK_INT  check = 0;
        CCTK_REAL prim[numprims];
        for(int which_guess=1;which_guess<=2;which_guess++) {
          set_prim_from_PRIMS_and_CONSERVS(eos,con2prim_test_keys[which_routine],which_guess,
                                           METRIC,METRIC_LAP_PSI4,
                                           PRIMS,CONSERVS,cons, prim);
          for(int i=0;i<numprims;i++) prim[i] = 1e300;
          // prim[TEMP    ] = eos.T_max;
          // prim[RHO     ] = PRIMS[RHOB       ];
          // prim[YE      ] = PRIMS[YEPRIM     ];
          prim[TEMP    ] = PRIMS[TEMPERATURE]*0.95;
          // prim[UTCON1  ] = PRIMS[VX         ];
          // prim[UTCON2  ] = PRIMS[VY         ];
          // prim[UTCON3  ] = PRIMS[VZ         ];
          prim[WLORENTZ] = W_test*0.95;
          // prim[B1_con  ] = PRIMS[BX_CENTER  ];
          // prim[B2_con  ] = PRIMS[BY_CENTER  ];
          // prim[B3_con  ] = PRIMS[BZ_CENTER  ];
          check = con2prim_select(eos,con2prim_test_keys[which_routine],METRIC_PHYS,g4dn,g4up,cons,prim,stats);
          if( check == 0 ) which_guess = 4;
        }

        CCTK_REAL ERRORS[MAXNUMVARS],accumulated_error = 0.0;
        if( check != 0 ) {
          failures++;
          CCTK_VInfo(CCTK_THORNSTRING,"Recovery FAILED!\n");
          accumulated_error = 1e300;
        }
        else {

          //Now that we have found some solution, we first limit velocity:
          //FIXME: Probably want to use exactly the same velocity limiter function here as in mhdflux.C
          CCTK_REAL utx_new = prim[UTCON1];
          CCTK_REAL uty_new = prim[UTCON2];
          CCTK_REAL utz_new = prim[UTCON3];
          
          // //Velocity limiter:
          CCTK_REAL gijuiuj = METRIC_PHYS[GXX]*SQR(utx_new ) +
            2.0*METRIC_PHYS[GXY]*utx_new*uty_new + 2.0*METRIC_PHYS[GXZ]*utx_new*utz_new +
            METRIC_PHYS[GYY]*SQR(uty_new) + 2.0*METRIC_PHYS[GYZ]*uty_new*utz_new +
            METRIC_PHYS[GZZ]*SQR(utz_new);
          CCTK_REAL au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
          CCTK_REAL u0L   = (au0m1+1.0)*METRIC_LAP_PSI4[LAPSEINV];
          PRIMS[YEPRIM     ] = prim[YE   ];
          PRIMS[TEMPERATURE] = prim[TEMP ];
          PRIMS[PRESSURE   ] = prim[PRESS];
          PRIMS[EPSILON    ] = prim[EPS  ];
          PRIMS[VX         ] = utx_new/u0L - METRIC[SHIFTX];
          PRIMS[VY         ] = uty_new/u0L - METRIC[SHIFTY];
          PRIMS[VZ         ] = utz_new/u0L - METRIC[SHIFTZ];

          CCTK_VInfo(CCTK_THORNSTRING,"Recovery SUCCEEDED!");
          for(int which_prim=0;which_prim<num_prims_in_error;which_prim++) {
            int primL = which_prims_in_error[which_prim];
            ERRORS[primL] = relative_error(PRIMS[primL],PRIMS_ORIG[primL]);
            CCTK_VInfo(CCTK_THORNSTRING,"Relative error for prim %s: %.3e (%e -> %e)",primnames[primL],ERRORS[primL],PRIMS_ORIG[primL],PRIMS[primL]);
            accumulated_error += ERRORS[primL];
          }
          CCTK_VInfo(CCTK_THORNSTRING,"Total accumulated error    : %e\n",accumulated_error);
        }
        fprintf(outfile,"%e %e %e\n",log10(PRIMS_ORIG[RHOB]),log10(PRIMS_ORIG[TEMPERATURE]),log10(MAX(accumulated_error,1e-16)));
      }
      fprintf(outfile,"\n");
    }

    fclose(outfile);

    int ntotal = npoints*npoints;
    
    CCTK_VInfo(CCTK_THORNSTRING,"Completed test for routine %s",con2prim_test_names[which_routine]);
    CCTK_VInfo(CCTK_THORNSTRING,"Final report:");
    CCTK_VInfo(CCTK_THORNSTRING,"    Number of recovery attempts: %d",ntotal);
    CCTK_VInfo(CCTK_THORNSTRING,"    Number of failed recoveries: %d",failures);
    CCTK_VInfo(CCTK_THORNSTRING,"    Recovery failure rate      : %.2lf%%",((CCTK_REAL)failures)/((CCTK_REAL)ntotal)*100.0);
    
  }

  CCTK_VInfo(CCTK_THORNSTRING,"All done!");
  // CCTK_VInfo(CCTK_THORNSTRING,"All done! Terminating the run.");
  // exit(1);
}
