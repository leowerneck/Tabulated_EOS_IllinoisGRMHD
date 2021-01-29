
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************

 macros.h : 
 ---------------

    Purpose: To contain all the constant macros.  These include
              constant parameters (e.g., enum-like sequences), in-line
              function macros, mnemonics, or defines that are
              determined by macros set in "options.h".
              
              Users should not need to modify this file unless they
              are adding functionality.  However, users may need to
              reference the values of these in order to understand
              output.

              It is meant to be included in decs.h and therefore in
              most every source file.

******************************************************************************/

/*************************************************************************************/
/*************************************************************************************/
/****   USERS SHOULD NOT NEED TO CHANGE ANY THING BELOW THIS LINE:   *****************/
/*************************************************************************************/
/*************************************************************************************/

/* Grid specifications : */
/* In the following, "ghost zones/cells" are those used to impose boundary conditions 
   and do not lie in the real, physical domain.  "Physical zones/cells" ARE
   those that lie in the real, physical domain.  
   x1,x2,x3 are typically thought of as 
*/
#define N0	  ( 2)	      /* number of time levels to store           */
#define N1	  (N1_glob)        /* number of zones in X1 direction          */
#define N2	  (N2_glob)	      /* number of zones in X2 direction          */
#define N3	  (N3_glob)	      /* number of zones in X3 direction          */
#define NG	  ( 3)	      /* number of ghost zones, same on each side */
#define NDIM      ( 4)        /* number of spacetime dimensions, 0=time   */
#define NPH       ( 5)        /* Number of hydro variables (i.e. rho,u,v^i) and assuming first NPH primitive variables are rho,uu,v^i */
#define NCELLS    (NCELLS_glob)  /* Number of physical cells  */

#define NP        ( 8 + N_PASSIVE_SCALARS )  /* number of EOM and # of primitive  variables */
#if( EVOLVE_ELECTRON_FRACTION )
#undef NP
#define NP ( 9 + N_PASSIVE_SCALARS )
#endif
#if(TABULATED_EOS_TABLES)
#undef NP
#define NP ( 10 + N_PASSIVE_SCALARS)
#endif

/* Specify whether we are using a binary black hole spacetime : */
#if( (METRIC_TYPE_CHOICE==METRIC_GENERAL_DYNAMIC||METRIC_TYPE_CHOICE==METRIC_GENERAL_PHI_AVG||METRIC_TYPE_CHOICE==METRIC_GENERAL_PHI_AVG2)&&((METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_35PN_NZ)||(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ)||(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_FAST)||(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ)||(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ)||(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_NEW)||METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_2ND||METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_2ND||(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_DROP)||(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NR_FZ)||(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_SIMPLE_NEWTONIAN)||(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_ASPIN_DROP)))
# define BBH_SPACETIME (1)   
#else 
# define BBH_SPACETIME (0)   
#endif

/* Specify whether we are using a binary black hole spacetime : */
#if( (METRIC_TYPE_CHOICE==METRIC_GENERAL_DYNAMIC||METRIC_TYPE_CHOICE==METRIC_GENERAL_PHI_AVG||METRIC_TYPE_CHOICE==METRIC_GENERAL_PHI_AVG2)&&((METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_35PN_NZ)||(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ)||(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_FAST)||(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ)||(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ)||(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_NEW)||METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_2ND||METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_2ND||(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_DROP)||(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NR_FZ)))
# define USE_ISOTROPIC_COORDS (1)   
#else 
# define USE_ISOTROPIC_COORDS (0)   
#endif


#if(METRIC_TYPE_CHOICE==METRIC_MINK_CARTESIAN)
#define METRIC_DIM   (0)
#elif(METRIC_TYPE_CHOICE==METRIC_KS_SPHERICAL||METRIC_TYPE_CHOICE==METRIC_BL_SPHERICAL||METRIC_TYPE_CHOICE==METRIC_GENERAL_PHI_AVG||METRIC_TYPE_CHOICE==METRIC_GENERAL_PHI_AVG2)
#define METRIC_DIM   (2)
#else
#define METRIC_DIM   (3)
#endif

/* This sets the topology type of the physical/numerical coordinates, which is useful to know : t*/
#define TOP_CARTESIAN   (0)
#define TOP_SPHERICAL  	(1)
#define TOP_CYLINDRICAL (2)

#if(METRIC_TYPE_CHOICE==METRIC_MINK_CARTESIAN||METRIC_TYPE_CHOICE==METRIC_KS_CARTESIAN||(BBH_SPACETIME&&COORD_TYPE_CHOICE==COORD_IDENTITY))
#define TOP_TYPE_CHOICE (TOP_CARTESIAN)
#elif(METRIC_TYPE_CHOICE==METRIC_MINK_SPHERICAL||METRIC_TYPE_CHOICE==METRIC_KS_SPHERICAL||METRIC_TYPE_CHOICE==METRIC_BL_SPHERICAL||METRIC_TYPE_CHOICE==METRIC_GENERAL_PHI_AVG||METRIC_TYPE_CHOICE==METRIC_GENERAL_PHI_AVG2)
#define TOP_TYPE_CHOICE (TOP_SPHERICAL)
#else
#define TOP_TYPE_CHOICE (TOP_SPHERICAL)
#endif

#if((METRIC_TYPE_CHOICE == METRIC_GENERAL_DYNAMIC || METRIC_TYPE_CHOICE == METRIC_GENERAL_STATIC) && ( METRIC_DYNAMIC_TYPE_CHOICE == METRIC_DYNAMIC_MINK_CARTESIAN) )
#undef TOP_TYPE_CHOICE
#define TOP_TYPE_CHOICE (TOP_CARTESIAN)
#endif
#if(METRIC_TYPE_CHOICE == METRIC_GENERAL_DYNAMIC || METRIC_TYPE_CHOICE == METRIC_GENERAL_STATIC || METRIC_TYPE_CHOICE == METRIC_GENERAL_PHI_AVG || METRIC_TYPE_CHOICE == METRIC_GENERAL_PHI_AVG2)
# if(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_KS_SPHERICAL||METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_MINK_SPHERICAL)
#   define METRIC_TOP_TYPE (TOP_SPHERICAL)
# elif(METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_35PN_NZ||METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ||METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_FAST||METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ||METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ||METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_NEW||METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_2ND||METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_2ND||METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NR_FZ||METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_DROP||METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_SIMPLE_NEWTONIAN||METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_MINK_CARTESIAN||METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_ASPIN_DROP)
#   define METRIC_TOP_TYPE (TOP_CARTESIAN)
# else 
#   define METRIC_TOP_TYPE (TOP_TYPE_CHOICE)
# endif
#else 
# define METRIC_TOP_TYPE (TOP_TYPE_CHOICE)
#endif


#if((COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD)||(COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL)||(COORD_TYPE_CHOICE==COORD_WARPED_CARTESIAN))
#  define DYNAMIC_COORDINATES (1)
#else 
#  define DYNAMIC_COORDINATES (0)
#endif

#if((METRIC_TYPE_CHOICE==METRIC_GENERAL_DYNAMIC)||(DYNAMIC_COORDINATES)||(COMOVING_METRIC))
#  define DYNAMIC_SPACETIME (1)
#else 
#  define DYNAMIC_SPACETIME (0)
#endif

 /* Number of time levels to store geometry at once :*/
#if(DYNAMIC_SPACETIME)
#  define  N0_GEOM  (3)   
#else
#  define  N0_GEOM  (1)   
#endif

 /* Number of time levels to store coordinate data at once :*/
#if( DYNAMIC_COORDINATES )
#define N0_COORD (3)
#else
#define N0_COORD (1)
#endif

/* Set the dimensionality of the coordinate radius */
#if( TOP_TYPE_CHOICE == TOP_SPHERICAL )
# if((COORD_TYPE_CHOICE==COORD_DIAGONAL         )||\
     (COORD_TYPE_CHOICE==COORD_MIXED            )||\
     (COORD_TYPE_CHOICE==COORD_DIAGONAL2        )||\
     (COORD_TYPE_CHOICE==COORD_DIAGONAL3        )||\
     (COORD_TYPE_CHOICE==COORD_DIAGONAL3_CUTOUT )||\
     (COORD_TYPE_CHOICE==COORD_DIAGONAL4        )||\
     (COORD_TYPE_CHOICE==COORD_DIAGONAL3_DYN_RAD)   )
#define COORD_RADIUS_DIM  (1)
# else
#define COORD_RADIUS_DIM  (3)
# endif
#else 
#define COORD_RADIUS_DIM  (3)
#endif

#if( RECON_TYPE_CHOICE >= RECON_PPM_MM && RECON_TYPE_CHOICE <= RECON_PPM_MC )
#define RECON_USE_PPM (1)
#else
#define RECON_USE_PPM (0)
#endif


/*************************************************************************************/
/******   MNEMONICS & CONSTANTS (NOT TO BE CHANGED)  *********************************/
/*************************************************************************************/

/* mnemonics for primitive vars */
#define RHO	(0)	  /* rest-mass density */
#define UU	(1)       /* internal energy density */ 
#define U1	(2)       /* spatial velocity component */
#define U2	(3)       /* spatial velocity component */
#define U3	(4)       /* spatial velocity component */
#define B1	(5)       /* spatial magnetic field component */
#define B2	(6)       /* spatial magnetic field component */
#define B3 	(7)       /* spatial magnetic field component */
#define YE      (8)       /* electron fraction component  */
#define TEMP    (9)   /* Temperature (in MeV)  */

/* Other grid parameters : */
#define N0TOT     (N0     )   /* number of numerical zone in X0 direction */
#define N1TOT     (N1TOT_glob)   /* number of numerical zone in X1 direction */
#define N2TOT     (N2TOT_glob)   /* number of numerical zone in X2 direction */
#define N3TOT     (N3TOT_glob)   /* number of numerical zone in X3 direction */
#define NTOT      (NTOT_glob)  /* total number of cells in numerical domain */

#define N0S       (0)         /* index of first physical cell in direction X0 */
#define N1S       (NG)        /* index of first physical cell in direction X1 */
#define N2S       (NG)        /* index of first physical cell in direction X2 */
#define N3S       (NG)        /* index of first physical cell in direction X3 */

#define N0E       (N0   -1)   /* index of last physical cell in X0 direction */
#define N1E       (N1E_glob)   /* index of last physical cell in X1 direction */
#define N2E       (N2E_glob)   /* index of last physical cell in X2 direction */
#define N3E       (N3E_glob)   /* index of last physical cell in X3 direction */

/* Determine the max. dimension of the data arrays : */
#define MAX_NTOT  (MAX_NTOT_glob)

/* Position types within a cell */
#define FACE1 (0)   /* Do not change values of FACE1-3, this order is used implicitly */
#define FACE2 (1)   /*  throughout the code. */
#define FACE3 (2)
#define CENT  (3)
#define CORN  (4)   
#define NPOS  (5)   /* Number of unique positions in a cell  */ 

/* mnemonics for dimensional indices */
#define TT      (0)     
#define RR      (1)
#define TH      (2)
#define PH      (3)

/* mnemonics for dimensional indices */
#define XX      (1)
#define YY      (2)
#define ZZ      (3)

/* mnemonics for (black hole) MASS_TYPE */
#define MTOT     (0)
#define MBH1     (1)
#define MBH2     (2)

/* Failure states (used by fail() ): */
#define FAIL_FLOOR           ( 1)   /* rho or u dropped below the floor values  */
#define FAIL_TMAX            ( 2)   /* u  >  rho * Tmax  */
#define FAIL_GAMMA_MAX       ( 3)   /* gamma became too large */ 
#define FAIL_GAMMA_CALC      ( 4)   /* A bad value of gamme was calculated  */
#define FAIL_USE_FULL_INV    ( 5)   /* Using the full inversion solution  */
#define FAIL_USE_EE          ( 6)   /* Using the entropy equation fix  */
#define FAIL_UTOPRIM         ( 7)   /* Regular Inversion failed to converge or solution was bad */
#define FAIL_UTOPRIM_EE      ( 8)   /* Entropy equation Inversion failed to converge or solution was bad */
#define FAIL_USE_INTERP_PRIM ( 9)   /* interpolated values from fixup_interp_prim() were used */
#define FAIL_USE_INTERP_V    (10)   /* interpolated values from fixup_interp_v() were used */
#define FAIL_INTERP_PRIM     (11)   /* No good stencils found to interpolate over prim. vars */
#define FAIL_INTERP_V        (12)   /* No good stencils found to interpolate over v^i */
#define FAIL_CUTOUT_PRESSURE (13)   /* Whether the pressure was change near the cutout */	    
#define FAIL_CUTOUT_GAMMA    (14)   /* Whether the Lorentz factor was changed along the cutout */
#define FAIL_CUTOUT_BC       (15)   /* Whether set_cutout_boundary*() changed anything */
#define FAIL_VCHAR_NEG       (16)   /* A Complex wave speed was calculated */
#define FAIL_VCHAR_SUPER     (17)   /* A super-luminal wave speed was calculated */
#define FAIL_VCHAR_DISCR     (18)   /* Discriminant is negative in wave speed calculation */
#define FAIL_BASIC           (19)   /* Basic failures, e.g. missing routine/file, I/O failure */
#define FAIL_METRIC          (20)   /* Failures with the metric routines */
#define FAIL_RESTART         (21)   /* Some problem with I/O of restart files */
#define FAIL_MPI_BASIC       (22)   /* Generic problem with a routine in harm_mpi.c  */
#define FAIL_HDF             (23)   /* Fundamental problem with a hdf routine   */
#define FAIL_EXCISE          (24)   /* Cell should not be evolved */
#define N_FAIL               (25)   /* Count null value that is not listed here "0"  */


/* nfail[]'s id numbers for different types of failures : */
#define NFAIL_FLOOR	       ( 0)   /* floor is reached in either rho or uu  */
#define NFAIL_TMAX             ( 1)   /* u surpassed  rho*Tmax                 */
#define NFAIL_GAMMA_MAX	       ( 2)   /* gamma reaches GAMMAMAX  */
#define NFAIL_GAMMA_CALC       ( 3)   /* gamma_calc() returns with a failure  */
#define NFAIL_USE_FULL_INV     ( 4)   /* using solution from full inversion routine */
#define NFAIL_USE_EE           ( 5)   /* using solution from fixup_entropy_eq() */
#define NFAIL_UTOPRIM          ( 6)   /* Utoprim() fails to converge to a solution */
#define NFAIL_UTOPRIM_EE       ( 7)   /* failed inversion for entropy equation */
#define NFAIL_USE_INTERP_PRIM  ( 8)   /* interpolated values from fixup_interp_prim() was used   */
#define NFAIL_USE_INTERP_V     ( 9)   /* interpolated values from fixup_interp_v() was used */
#define NFAIL_INTERP_PRIM      (10)   /* fixup_interp_prim() fails to provide adequate fix  */
#define NFAIL_INTERP_V         (11)   /* fixup_interp_v() fails to provide adequate fix */
#define NFAIL_CUTOUT_PRESSURE  (12)   /* Whether the pressure was change near the cutout */	    
#define NFAIL_CUTOUT_GAMMA     (13)   /* Whether the Lorentz factor was changed along the cutout */
#define NFAIL_CUTOUT_BC        (14)   /* Whether set_cutout_boundaryN() changed any values along the cutout boundary */
#define N_NFAIL                (15)   /* number of failure modes tallied by nfail[]  */


/* pflag codes : must be sequential starting with PFLAG_NONE=0, 
   the largest valued one should be the most conservative (or worst?) failure option */
#define PFLAG_NONE        (0)  /* No problems */
#define PFLAG_INTERP_V    (1)  /* Just interpolate v^i   */
#define PFLAG_INTERP_PRIM (2)  /* Interpolate rho,u,v^i  */
#define N_PFLAG_TYPES     (3)  

/* Mnemonics for the "evolution mask" which determines whether or a not certain 
   parts of the update are performed:  */
#define MASK_EXCISED  0    /* Absolutely nothing is calculated here */
#define MASK_BUFFER   1    /* Reconstruction is performed here      */
#define MASK_NORMAL   2    /* Regular evolution, everything is calculatesd */
#define MASK_CUTOUT   3    /* Masked cells for cutout */

/* Output types, RESTART needs to be on end before INIT&FINAL */
#define OUT_HISTORY ( 0)
#define OUT_ASCII   ( 1)
#define OUT_SDF     ( 2)
#define OUT_HDF5    ( 3)
#define OUT_IMAGE   ( 4)
#define OUT_STAT    ( 5)
#define OUT_STAT2   ( 6)
#define OUT_RADFLUX ( 7)
#define OUT_MIN_DT  ( 8)
#define OUT_SURFACE ( 9)
#define OUT_RESTART (10)
#define OUT_INIT    (11)
#define OUT_EVERY   (12)
#define OUT_FINAL   (13)
#define N_OUT_TYPES (OUT_RESTART)  /* number of output types excluding RESTART,INIT,EVERY,FINAL */


/* The number of gridfunctions written to hdf5 file. See "set_hdf5_gfuncs()" below */
/* Just MHD stuff plus EM current functions : */
#define N_HDF_FUNCS (NP+4+10+1+1+1)


/* Number of history and surface functions */
#define N_HISTORY_BASE (25)
#define N_SURFACE_BASE (21)

#if( CALC_CURRENT )
#define N_HISTORY_CURRENT (5)
#else 
#define N_HISTORY_CURRENT (0)
#endif
#define N_SURFACE_CURRENT (N_HISTORY_CURRENT)

#if( MAKE_RADFLUX )
#define N_HISTORY_RADFLUX (3)
#else
#define N_HISTORY_RADFLUX (0)
#endif
#define N_SURFACE_RADFLUX (N_HISTORY_RADFLUX)

#if( MAKE_EHT_DATA )
#define N_HISTORY_EHT_DATA (12)
#else
#define N_HISTORY_EHT_DATA (0)
#endif
#define N_SURFACE_EHT_DATA (N_HISTORY_EHT_DATA)

#if( METRIC_TYPE_CHOICE == METRIC_GENERAL_DYNAMIC )
#define N_HISTORY_DYNMET  (6)
#else
#define N_HISTORY_DYNMET  (0)
#endif
#define N_SURFACE_DYNMET  (N_HISTORY_DYNMET)

#define N_SURFACE  (N_SURFACE_BASE+N_SURFACE_CURRENT+N_SURFACE_RADFLUX+N_SURFACE_EHT_DATA+N_SURFACE_DYNMET)
#define N_HISTORY  (N_HISTORY_BASE+N_HISTORY_CURRENT+N_HISTORY_RADFLUX+N_HISTORY_EHT_DATA+N_HISTORY_DYNMET)

/* Number of groups or subdirectories in the root group (e.g. Bound, Unbound...): */
#define UNBOUND  (0)
#define BOUND    (1)
#define JET      (2)
#define DISK     (3)
#define CORONA   (4)
#define N_HIST_TYPES  (5)

#define N_HIST_POINTS (N_HIST_POINTS_glob)

#define N_SURF_TYPES  (N_HIST_TYPES)
#define N_SURF_POINTS (N_SURF_POINTS_glob)


/* To differential viscous stress, we need connection in ghost zones; 
     may need to increase this number if we want to increase the order of 
     accuracy of the finite differencing in phys.c  */
#if( USE_KINEMATIC_VISCOSITY ) 
# define NDC (1) 
#else 
# define NDC (0) 
#endif
#if( NDC > NG ) 
NDC cannot be greater than NG 
#endif


/* Determine the length of geom_arr, conn[] and arrays over the physical domain */
#define N_GEOM   (N_GEOM_glob)
#define N_CONN   (N_CONN_glob)

#define N_COORD  (N_COORD_glob)


/* for-loop aliases : */
#define LOOP       for(i=N1S;i<=N1E ;i++) for(j=N2S;j<=N2E;j++) for(k=N3S;k<=N3E;k++)
#define N1_LOOP    for(i=N1S;i<=N1E ;i++)
#define N2_LOOP    for(j=N2S;j<=N2E ;j++)
#define N3_LOOP    for(k=N3S;k<=N3E ;k++)
#define N1ALL_LOOP for(i=0  ;i<N1TOT;i++)
#define N2ALL_LOOP for(j=0  ;j<N2TOT;j++)
#define N3ALL_LOOP for(k=0  ;k<N3TOT;k++)
#define ALL_LOOP   for(i=0  ;i<N1TOT;i++) for(j=0  ;j<N2TOT;j++) for(k=0  ;k<N3TOT;k++)
#define PLOOP      for(l=0  ;l<NP   ;l++)
#define BLOOP      for(l=B1 ;l<=B3  ;l++)
#define GLOOP      for(g=0  ;g<NG   ;g++)
#define DLOOP1     for(i=0  ;i<NDIM ;i++)
#define DLOOP2     for(i=0  ;i<NDIM ;i++) for(j=0  ;j<NDIM ;j++)
#define DLOOP3     for(i=0  ;i<NDIM ;i++) for(j=0  ;j<NDIM ;j++) for(k=0  ;k<NDIM ;k++)
#define SDLOOP1    for(i=1  ;i<NDIM ;i++)
#define POSLOOP    for(pos=0;pos<NPOS;pos++) 
#define FACE_LOOP  for(d=0  ;d<(NDIM-1);d++)
#define NPH_LOOP   for(l=0;l<NPH;l++)

#if( N_PASSIVE_SCALARS )
#define PASSIVE_SCALAR_BEG  (NP-N_PASSIVE_SCALARS)
#define PASSIVE_SCALAR_END  (NP-1)
#define PASSIVE_SCALAR_LOOP        for(i_pscalar=0                   ; i_pscalar < N_PASSIVE_SCALARS    ; i_pscalar++)
#define PASSIVE_SCALAR_PRIM_LOOP   for(i_pscalar=(PASSIVE_SCALAR_BEG); i_pscalar <= (PASSIVE_SCALAR_END); i_pscalar++)
#endif


/* Prim var Fixup loop:  Averaging is not done over B^i since they have been updated fine */ 
#define FLOOP for(l=0;l<B1;l++)

/* Veloctiy Fixup loop */
#define VLOOP for(l=U1;l<=U3;l++)

#define OUT_LOOP   for(i=0; i<N_OUT_TYPES; i++) 
#define DUMP_LOOP   LOOP 
#define RDUMP_LOOP  DUMP_LOOP


#define CONN_LOOP  for(l=0; l<N_CONN; l++) 

#define CONN2_BEG1  (N1S-NDC)
#define CONN2_BEG2  (N2S-NDC)
#define CONN2_BEG3  (N3S-NDC)

#define CONN2_END1  (CONN2_END1_glob)
#define CONN2_END2  (CONN2_END2_glob)
#define CONN2_END3  (CONN2_END3_glob)

#define N1_LOOP_LOC   for(ii = ibeg; ii <= iend; ii++)
#define N2_LOOP_LOC   for(jj = jbeg; jj <= jend; jj++)
#define N3_LOOP_LOC   for(kk = kbeg; kk <= kend; kk++)

#if(METRIC_DIM==3) 
#define GEOM_LOOP   POSLOOP  N1ALL_LOOP N2ALL_LOOP N3ALL_LOOP 
#define GDUMP_LOOP     N1_LOOP    N2_LOOP    N3_LOOP 
#define CONN2_LOOP  for(i=CONN2_BEG1;i<=CONN2_END1 ;i++) for(j=CONN2_BEG2;j<=CONN2_END2 ;j++) for(k=CONN2_BEG3;k<=CONN2_END3 ;k++)
#elif(METRIC_DIM==2) 
#define GEOM_LOOP   for(k=N3S;k<=N3S;k++) POSLOOP N1ALL_LOOP N2ALL_LOOP
#define GDUMP_LOOP  for(k=N3S;k<=N3S;k++)    N1_LOOP    N2_LOOP 
#define CONN2_LOOP  for(k=N3S;k<=N3S;k++) for(i=CONN2_BEG1;i<=CONN2_END1 ;i++) for(j=CONN2_BEG2;j<=CONN2_END2 ;j++) 
#elif(METRIC_DIM==1) 
#define GEOM_LOOP   for(j=N2S;j<=N2S;j++) for(k=N3S;k<=N3S;k++) POSLOOP N1ALL_LOOP 
#define GDUMP_LOOP  for(j=N2S;j<=N2S;j++) for(k=N3S;k<=N3S;k++)    N1_LOOP 
#define CONN2_LOOP  for(j=N2S;j<=N2S;j++) for(k=N3S;k<=N3S;k++) for(i=CONN2_BEG1;i<=CONN2_END1 ;i++) 
#else
#define GEOM_LOOP   for(pos=0;pos<1;pos++) for(i=N1S;i<=N1S;i++) for(j=N2S;j<=N2S;j++) for(k=N3S;k<=N3S;k++)
#define GDUMP_LOOP  for(i=N1S;i<=N1S;i++) for(j=N2S;j<=N2S;j++) for(k=N3S;k<=N3S;k++) 
#define CONN2_LOOP  GDUMP_LOOP
#endif


#define SMALL	  (1.e-14)   /* numerical convenience */
#define SAFE      (1.3)      /* Max. factor to increase dt by */

/* This needs to be an impossible value of MPI "pid" */
#define BAD_PID    (-123456)
#define BC_PHYS    (BAD_PID)

/* Boundary condition mnemonics */
#define BCDN    (0)
#define BCUP    (1)


/* Number of instances we calculate U(P) in the timestep process for statistical reasons */
#define N_STAT (11) 

/* ASCII output formats : */ 
#define FMT_DBL_OUT "%28.18e "
#define FMT_INT_OUT "%14d "

 /* macros used for TRACE_CALLS : */


/* "Inline function" type  macros: */
#define   MAX(_a,_b) ((_a) > (_b) ? (_a) : (_b))
#define   MIN(_a,_b) ((_a) < (_b) ? (_a) : (_b))
#define  MINK(_a,_b) ((double)  ( (_a==_b)*( -1*(_a==0) + (_a!=0) ) ) ) 
#define   DOT(_a,_b) ((_a[0])*(_b[0]) + (_a[1])*(_b[1]) + (_a[2])*(_b[2]) + (_a[3])*(_b[3]))
#define  SIGN(_a)   ( ((_a) <0.) ? -1. : 1. )
#define DELTA(_i,_j) (((_i) == (_j)) ? 1. : 0.)
#define REL_DIFF_FUNC(_fval1,_fval2) ( ( ( (_fval1)==0. ) && ( (_fval2)==0. ) ) ? (0.) : ( fabs( 2.*((_fval1)-(_fval2))/((_fval1)+(_fval2)) ) ) )
#define IS_DIFFERENT_DOUBLE(_fval1,_fval2) ( (REL_DIFF_FUNC((_fval1),(_fval2)) > SMALL ) )
#define FREE(_p) { free(_p); _p = NULL; } 
#define FCLOSE(_f) { fclose(_f); _f = NULL; }
#define NINT( _f )  ((int) ( (_f) + 0.5 ) )

/* Should be precise to 1e-17 level : */
#define MY_SQRT_ONE_PLUS_EPS(_a) ( ((_a) < 1.e-1) ? (1.0+(0.5+(-0.125+(0.625E-1+(-0.390625E-1+(0.2734375E-1+(-0.205078125E-1+(0.1611328125E-1+(-0.13092041015625E-1+(0.109100341796875E-1+(-0.9273529052734375E-2+(0.8008956909179688E-2+(-0.7007837295532227E-2+0.6199240684509277E-2*(_a))*(_a))*(_a))*(_a))*(_a))*(_a))*(_a))*(_a))*(_a))*(_a))*(_a))*(_a))*(_a))  : sqrt(1.+(_a)) ) 

/* Produces a periodic function over x0 to xn */
#define PERIODIC_SHIFT(x,x0,xn) ( ((x) < (x0)) ? (1+((int) (fabs((x)-(x0))/((xn)-(x0))))) : (-((int) (fabs((x)-(x0))/((xn)-(x0))))) )
#define PERIODIC(x,x0,xn) ( (x) + ((xn)-(x0))*(PERIODIC_SHIFT((x),(x0),(xn))) )

#define ALLOC_ARRAY(_array,_npts)    {\
    if( ( (_array) = calloc( (_npts) , sizeof(*(_array)))) == NULL ) {	\
     fprintf(stdout,"%s(): Cannot allocate %s on line %d  of %s \n",__func__,#_array,__LINE__,__FILE__);  exit(0); \
   }\
 }

#define ALLOC_2D_ARRAY(_array,_n1,_n2) { \
 int _i1;\
 ALLOC_ARRAY(_array,_n1);\
 for(_i1=0; _i1<_n1; _i1++) { ALLOC_ARRAY(_array[_i1],_n2); }\
 }

#define ALLOC_3D_ARRAY(_array,_n1,_n2,_n3) { \
 int _i1, _i2;\
 ALLOC_2D_ARRAY(_array,_n1,_n2);\
 for(_i1=0; _i1<_n1; _i1++) for(_i2=0; _i2<_n2; _i2++) {  ALLOC_ARRAY(_array[_i1][_i2],_n3); }\
 }

#define ALLOC_4D_ARRAY(_array,_n1,_n2,_n3,_n4) { \
 int _i1, _i2, _i3;\
 ALLOC_3D_ARRAY(_array,_n1,_n2,_n3);\
 for(_i1=0; _i1<_n1; _i1++) for(_i2=0; _i2<_n2; _i2++) for(_i3=0; _i3<_n3; _i3++) {  ALLOC_ARRAY(_array[_i1][_i2][_i3],_n4); } \
 }

#define ALLOC_5D_ARRAY(_array,_n1,_n2,_n3,_n4,_n5) {     \
 int _i1, _i2, _i3,_i4;\
 ALLOC_4D_ARRAY(_array,_n1,_n2,_n3,_n4);\
 for(_i1=0; _i1<_n1; _i1++) for(_i2=0; _i2<_n2; _i2++) for(_i3=0; _i3<_n3; _i3++) for(_i4=0; _i4<_n4; _i4++) {  ALLOC_ARRAY(_array[_i1][_i2][_i3][_i4],_n5); }\
 }

#define DEALLOC_ARRAY(_array,_npts)    {FREE((_array)); }

#define DEALLOC_2D_ARRAY(_array,_n1,_n2) { \
 int _i1;\
 for(_i1=0; _i1<_n1; _i1++) { DEALLOC_ARRAY(_array[_i1],_n2); }\
 DEALLOC_ARRAY(_array,_n1);\
 }

#define DEALLOC_3D_ARRAY(_array,_n1,_n2,_n3) { \
 int _i1, _i2;\
 for(_i1=0; _i1<_n1; _i1++) for(_i2=0; _i2<_n2; _i2++) {  DEALLOC_ARRAY(_array[_i1][_i2],_n3); }\
 DEALLOC_2D_ARRAY(_array,_n1,_n2);\
 }

#define DEALLOC_4D_ARRAY(_array,_n1,_n2,_n3,_n4) { \
 int _i1, _i2, _i3;\
 for(_i1=0; _i1<_n1; _i1++) for(_i2=0; _i2<_n2; _i2++) for(_i3=0; _i3<_n3; _i3++) {  DEALLOC_ARRAY(_array[_i1][_i2][_i3],_n4); } \
 DEALLOC_3D_ARRAY(_array,_n1,_n2,_n3);\
 }

#define DEALLOC_5D_ARRAY(_array,_n1,_n2,_n3,_n4,_n5) {     \
 int _i1, _i2, _i3,_i4;\
 for(_i1=0; _i1<_n1; _i1++) for(_i2=0; _i2<_n2; _i2++) for(_i3=0; _i3<_n3; _i3++) for(_i4=0; _i4<_n4; _i4++) {  DEALLOC_ARRAY(_array[_i1][_i2][_i3][_i4],_n5); }\
 DEALLOC_4D_ARRAY(_array,_n1,_n2,_n3,_n4);\
 }

/* Safe version of DEALLOC_?D_ARRAY()  that make sure that arrays are allocated before looping through them : */
#define DEALLOC_2D_ARRAY_SAFE(_array,_n1,_n2) { \
 int _i1;\
 if( _array != NULL) { for(_i1=0; _i1<_n1; _i1++) { DEALLOC_ARRAY(_array[_i1],_n2); }} \
 DEALLOC_ARRAY(_array,_n1);\
 }

#define DEALLOC_3D_ARRAY_SAFE(_array,_n1,_n2,_n3) { \
 int _i1, _i2;\
 if( _array != NULL) { for(_i1=0; _i1<_n1; _i1++) if( _array[_i1] != NULL) { for(_i2=0; _i2<_n2; _i2++) {  DEALLOC_ARRAY(_array[_i1][_i2],_n3); }}} \
 DEALLOC_2D_ARRAY(_array,_n1,_n2);\
 }

#define DEALLOC_4D_ARRAY_SAFE(_array,_n1,_n2,_n3,_n4) { \
 int _i1, _i2, _i3;\
 if( _array != NULL) { for(_i1=0; _i1<_n1; _i1++) if( _array[_i1] != NULL) { for(_i2=0; _i2<_n2; _i2++) if( _array[_i1][_i2] != NULL) { for(_i3=0; _i3<_n3; _i3++) {  DEALLOC_ARRAY(_array[_i1][_i2][_i3],_n4); }}}} \
 DEALLOC_3D_ARRAY(_array,_n1,_n2,_n3);\
 }

#define DEALLOC_5D_ARRAY_SAFE(_array,_n1,_n2,_n3,_n4,_n5) {     \
 int _i1, _i2, _i3,_i4;\
 for(_i1=0; _i1<_n1; _i1++) for(_i2=0; _i2<_n2; _i2++) for(_i3=0; _i3<_n3; _i3++) for(_i4=0; _i4<_n4; _i4++) {  DEALLOC_ARRAY(_array[_i1][_i2][_i3][_i4],_n5); }\
 if( _array != NULL) { for(_i1=0; _i1<_n1; _i1++) if( _array[_i1] != NULL) { for(_i2=0; _i2<_n2; _i2++) if( _array[_i1][_i2] != NULL) { for(_i3=0; _i3<_n3; _i3++) if( _array[_i1][_i2][_i3] != NULL) {  for(_i4=0; _i4<_n4; _i4++) {  DEALLOC_ARRAY(_array[_i1][_i2][_i3][_i4],_n5); }}}}} \
 DEALLOC_4D_ARRAY(_array,_n1,_n2,_n3,_n4);\
 }


/* CONN_ID returns the connection's index referring to a cell at (i,j,k)=(it,jt,kt)  */
#if(METRIC_DIM==3) 
 //#define CONN_ID(it,jt,kt)  (-NG*(1+N3*(1+N2)) + ((kt) + N3*((jt) + N2*(it))))     /*  = ((kt-NG) + N3*((jt-NG) + N2*((it-NG))));  */
#define CONN_ID(it,jt,kt)  ((kt-NG+NDC) + (N3_CONN_glob)*((jt-NG+NDC) + (N2_CONN_glob)*((it-NG+NDC)))) 
#elif(METRIC_DIM==2) 
//#define CONN_ID(it,jt,kt)  (-NG*(1+N2) + (jt) + N2*(it))       /*   = (jt-NG) + N2*((it-NG)); */
#define CONN_ID(it,jt,kt) (((jt-NG+NDC) + (N2_CONN_glob)*((it-NG+NDC))))     /*   = (jt-NG) + N2*((it-NG)); */
#elif(METRIC_DIM==1) 
 //#define CONN_ID(it,jt,kt)  ((it)-NG)  
#define CONN_ID(it,jt,kt) ((((it-NG+NDC))))
#else
#define CONN_ID(it,jt,kt)  (0)
#endif

/* GEOM_ID returns geom_arr's index referring to a cell at (i,j,k,pos)=(it,jt,kt,pt)  */
#if(METRIC_DIM==3) 
#define GEOM_ID(it,jt,kt,pt)  ( (kt) + N3TOT*((jt) + N2TOT*((it) + N1TOT*(pt) )) )    
#elif(METRIC_DIM==2) 
#define GEOM_ID(it,jt,kt,pt)  ( (jt) + N2TOT*((it) + N1TOT*(pt) ) )    
#elif(METRIC_DIM==1) 
#define GEOM_ID(it,jt,kt,pt)  ( (it) + N1TOT*(pt) )    
#else
#define GEOM_ID(it,jt,kt,pt)  (0)
#endif

/* COORD_ID returns geom_arr's index referring to a cell at (i,j,k,pos)=(it,jt,kt,pt)  */
#define COORD_ID(it,jt,kt,pt)  ( (kt) + N3TOT*((jt) + N2TOT*((it) + N1TOT*(pt) )) )    

/* Macro to replace get_geometry() routine  that returns the pointer to the geometry structure for 
   a given cell (ii_g,jj_g,kk_g)  at position=dd_g  in the cell and at time level nn_g             e*/
#define get_geometry(ii_g,jj_g,kk_g,dd_g,nn_g,geom_struct_g)  { icurr=(ii_g); jcurr=(jj_g); kcurr=(kk_g); pcurr=(dd_g);  geom_struct_g = &(geom_arr[(nn_g)][GEOM_ID((ii_g),(jj_g),(kk_g),(dd_g))]) ;  }

/* Macro to replace x_of_xp and dx_dxp_calc and coord()  routines that returns the pointer to the 
   coordinate structure for a given cell (ii_g,jj_g,kk_g)  at position=dd_g  in the cell and at 
   time level nn_g  */
#if( DYNAMIC_COORDINATES )
# define get_coord(ii_g,jj_g,kk_g,dd_g,nn_g,coord_struct_g)  { icurr=(ii_g); jcurr=(jj_g); kcurr=(kk_g); pcurr=(dd_g);  coord_struct_g = &(coord_arr[(nn_g)][COORD_ID((ii_g),(jj_g),(kk_g),(dd_g))]) ;  }
#else 
# define get_coord(ii_g,jj_g,kk_g,dd_g,nn_g,coord_struct_g)  { icurr=(ii_g); jcurr=(jj_g); kcurr=(kk_g); pcurr=(dd_g);  coord_struct_g = &(coord_arr[0][COORD_ID((ii_g),(jj_g),(kk_g),(dd_g))]) ;  }
#endif


/*************************************************************************************/
/** TIMING MACROS and VARIABLES                                   ********************/
/*************************************************************************************/
#if( MAKE_TIMERS ) 
#include <sys/time.h>

#define TIMER_TYPE_TOTAL    (0)   /* the full update process including everything done in one timestep */
#define TIMER_TYPE_ADVANCE  (1)   /* timer for calls to advance() */
#define TIMER_TYPE_METRIC   (2)   /* timer for calls to routines that update the spacetime metric or connection */
#define TIMER_TYPE_BOUNDS   (3)   /* timer for routines involving bounds() and its associated MPI time  */
#define TIMER_TYPE_DIAG     (4)   /* timer for all calls to diag() during an update step (i.e. not beginning and end calls) */
#define TIMER_TYPE_MPI      (5)   /* timer for all MPI_Wait(), i.e. how long we are waiting on our neighbors */
#define N_TIMER_TYPES       (6)

#define TIMER_SECONDS      (0)
#define TIMER_MICROSECONDS (1)

#define BEG_TIMING(_itimer)  gettimeofday(&(beg_times[(_itimer)]), NULL); 

#define END_TIMING(_itimer)  { gettimeofday(&end_time, NULL);\
    elapsed_times[(_itimer)][TIMER_SECONDS     ] += end_time.tv_sec  - beg_times[(_itimer)].tv_sec ;\
    elapsed_times[(_itimer)][TIMER_MICROSECONDS] += end_time.tv_usec - beg_times[(_itimer)].tv_usec;\
  }

# if(USEMPI)
#define MY_MPI_Wait(_arg1,_arg2) { \
  BEG_TIMING(TIMER_TYPE_MPI);\
  MPI_Wait((_arg1),(_arg2));\
  END_TIMING(TIMER_TYPE_MPI);\
  }
# endif

#else
#define BEG_TIMING(_itimer)  {;}
#define END_TIMING(_itimer)  {;}

# if(USEMPI)
#define MY_MPI_Wait(_arg1,_arg2)  MPI_Wait((_arg1),(_arg2));
# endif

#endif

#if( GATHER_IO ) 
#define myH5_write_gfunc  myH5_write_gfunc_gather
#define myH5_write_int_gfunc  myH5_write_int_gfunc_gather
#define myH5_write_scalar2  myH5_write_scalar2_gather
#define myH5_read_scalar2  myH5_read_scalar2_gather
#define myH5_read_gfunc myH5_read_gfunc_gather
#define myH5_Fcreate myH5_Fcreate_gather
#define restart_check_hdf5 restart_check_hdf5_gather
#else 
#define myH5_write_gfunc  myH5_write_gfunc_orig
#define myH5_write_int_gfunc  myH5_write_int_gfunc_orig
#define myH5_write_scalar2  myH5_write_scalar2_orig
#define myH5_read_scalar2  myH5_read_scalar2_orig
#define myH5_read_gfunc myH5_read_gfunc_orig
#define myH5_Fcreate myH5_Fcreate_orig
#define restart_check_hdf5 restart_check_hdf5_orig
#endif


